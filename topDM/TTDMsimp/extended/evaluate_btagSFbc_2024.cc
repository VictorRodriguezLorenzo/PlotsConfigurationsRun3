#ifndef BTAGSF
#define BTAGSF
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "TVector2.h"
#include "Math/Vector4Dfwd.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "TString.h"

#include <iostream>
#include "ROOT/RVec.hxx"
#include "TH2D.h"
#include "TFile.h"
#include "correction.h"

using namespace ROOT;
using namespace ROOT::VecOps;
using correction::CorrectionSet;

class btagSFbc {

  public:
    btagSFbc(TString eff_map, const string year, TString algo_extension);
    ~btagSFbc();

    TH2F* h_bjet_eff;
    TH2F* h_cjet_eff;

    std::string year_;
    std::unique_ptr<CorrectionSet> cset;
    std::unique_ptr<CorrectionSet> cset_light;

    RVecF operator()(
                    RVecF         CleanJet_pt,
                    RVecF         CleanJet_eta,
                    RVecI         CleanJet_jetIdx,
                    unsigned int  nCleanJet,
                    RVecI         Jet_hadronFlavour,
                    RVecF         Jet_btag,
                    const string  WP,
                    RVec<string> systematic
                     ){

      RVecF results(systematic.size(), 1.0);

      std::string correctionKey;
      std::string wpKey;

      std::string correctionKey_alt;

      if (year_ == "2024_Summer24"){
        correctionKey = "UParTAK4_kinfit";
        wpKey = "UParTAK4_wp_values";

        correctionKey_alt = "robustParticleTransformer_mujets";
      }else{
        correctionKey   = "robustParticleTransformer_mujets";
        wpKey = "robustParticleTransformer_wp_values";
      }

      auto cset_deepJet_mujets_alt  = cset->at(correctionKey);
      if (year_ == "2024_Summer24")
        cset_deepJet_mujets_alt  = cset_light->at(correctionKey_alt);
      auto cset_deepJet_mujets  = cset->at(correctionKey);
      auto cset_deepJet_wps     = cset->at(wpKey);

      for (long unsigned int i=0; i<systematic.size(); i++){

        float btag_sf    = 1.;

        for (unsigned iJ{0}; iJ != nCleanJet; ++iJ) {

          double pt = ROOT::VecOps::Min(ROOT::RVecD{CleanJet_pt[iJ], 9999.9});

          if (CleanJet_pt[iJ] <= 30. || abs(CleanJet_eta[iJ]) >= 2.5) continue;
          if (Jet_btag[CleanJet_jetIdx[iJ]] > cset_deepJet_wps->evaluate({WP})) {
            if (Jet_hadronFlavour[CleanJet_jetIdx[iJ]] == 5) {
              btag_sf *= cset_deepJet_mujets->evaluate({systematic[i], WP, 5, abs(CleanJet_eta[iJ]), pt});
            }
            else if (Jet_hadronFlavour[CleanJet_jetIdx[iJ]] == 4) {

              if (year_ == "2024_Summer24"){
                //btag_sf *= cset_deepJet_mujets_alt->evaluate({systematic[i], WP, 4, abs(CleanJet_eta[iJ]), CleanJet_pt[iJ]});
                btag_sf *= cset_deepJet_mujets_alt->evaluate({"central", WP, 4, abs(CleanJet_eta[iJ]), pt});
              }else{
                btag_sf *= cset_deepJet_mujets->evaluate({systematic[i], WP, 4, abs(CleanJet_eta[iJ]), pt});
              }
            }
            else {
              continue;
            }
          }
          else {
            float btag_eff = getEff(pt, CleanJet_eta[iJ], Jet_hadronFlavour[CleanJet_jetIdx[iJ]]);
            if (btag_eff == 1.) {
              continue;
            }else {
              if (Jet_hadronFlavour[CleanJet_jetIdx[iJ]] == 5) {
                btag_sf *= (1-btag_eff*cset_deepJet_mujets->evaluate({systematic[i], WP, 5, abs(CleanJet_eta[iJ]), pt}))/(1-btag_eff);
              }
              else if (Jet_hadronFlavour[CleanJet_jetIdx[iJ]] == 4) {
                if (year_ == "2024_Summer24"){
                  //btag_sf *= (1-btag_eff*cset_deepJet_mujets_alt->evaluate({systematic[i], WP, 4, abs(CleanJet_eta[iJ]), CleanJet_pt[iJ]}))/(1-btag_eff);
                  btag_sf *= (1-btag_eff*cset_deepJet_mujets_alt->evaluate({"central", WP, 4, abs(CleanJet_eta[iJ]), pt}))/(1-btag_eff);
                }else{
                  btag_sf *= (1-btag_eff*cset_deepJet_mujets->evaluate({systematic[i], WP, 4, abs(CleanJet_eta[iJ]), pt}))/(1-btag_eff);
                }
              }
              else {
                continue;
              }
            }
          }
        }
        results[i] = btag_sf;
      }
      return results;
    }

  private:
    float getEff(float pt, float eta, int flavour);

};

btagSFbc::btagSFbc(TString eff_map, const string year, TString algo_extension = "") {
  std::string home;
  year_ = year;
  if (year == "2024_Summer24"){
    home = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/";
    cset = CorrectionSet::from_file(home + "/btagging_preliminary.json.gz");

    cset_light = CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/btagging.json.gz");
  }

  // deepJet -> bjet_eff
  // PNet -> bjet_pnet_eff
  // RobPartT -> bjet_parT_eff

  TFile *reff = TFile::Open(eff_map, "READ");
  h_bjet_eff  = (TH2F*)reff->Get("bjet"+algo_extension+"_eff")->Clone();
  h_cjet_eff  = (TH2F*)reff->Get("cjet"+algo_extension+"_eff")->Clone();
  h_bjet_eff->SetDirectory(0);
  h_cjet_eff->SetDirectory(0);
  reff->Close();
}

float btagSFbc::getEff(float pt, float eta, int flavour) {
  int xbin, ybin;
  float eff;
  if (flavour == 5) {
    xbin = h_bjet_eff->GetXaxis()->FindBin(pt);
    ybin = h_bjet_eff->GetYaxis()->FindBin(eta);
    eff  = h_bjet_eff->GetBinContent(xbin, ybin);
  }
  else if (flavour == 4) {
    xbin = h_cjet_eff->GetXaxis()->FindBin(pt);
    ybin = h_cjet_eff->GetYaxis()->FindBin(eta);
    eff  = h_cjet_eff->GetBinContent(xbin, ybin);
  }
  else {
    eff   = 1.;
  }
  return eff;
}

btagSFbc::~btagSFbc(){
  std::cout << "Cleaning up memory" << std::endl;
  h_bjet_eff->Delete();
  h_cjet_eff->Delete();
}

#endif
