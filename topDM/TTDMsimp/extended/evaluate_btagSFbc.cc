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
  
    std::unique_ptr<CorrectionSet> cset;
 
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
      
      auto cset_deepJet_mujets  = cset->at("robustParticleTransformer_mujets");
      auto cset_deepJet_wps     = cset->at("robustParticleTransformer_wp_values");

      for (long unsigned int i=0; i<systematic.size(); i++){
	
	float btag_sf    = 1.;
	
	for (unsigned iJ{0}; iJ != nCleanJet; ++iJ) {
	  if (CleanJet_pt[iJ] <= 30. || abs(CleanJet_eta[iJ]) >= 2.5) continue;
	  if (Jet_btag[CleanJet_jetIdx[iJ]] > cset_deepJet_wps->evaluate({WP})) {
	    if (Jet_hadronFlavour[CleanJet_jetIdx[iJ]] == 5) {
	      btag_sf *= cset_deepJet_mujets->evaluate({systematic[i], WP, 5, abs(CleanJet_eta[iJ]), CleanJet_pt[iJ]});
	    }
	    else if (Jet_hadronFlavour[CleanJet_jetIdx[iJ]] == 4) {
	      btag_sf *= cset_deepJet_mujets->evaluate({systematic[i], WP, 4, abs(CleanJet_eta[iJ]), CleanJet_pt[iJ]});
	    }
	    else {
	      continue;
	    }
	  }
	  else {
	    float btag_eff = getEff(CleanJet_pt[iJ], CleanJet_eta[iJ], Jet_hadronFlavour[CleanJet_jetIdx[iJ]]);
	    if (btag_eff == 1.) {
	      continue;
	    }else {
	      if (Jet_hadronFlavour[CleanJet_jetIdx[iJ]] == 5) {
		btag_sf *= (1-btag_eff*cset_deepJet_mujets->evaluate({systematic[i], WP, 5, abs(CleanJet_eta[iJ]), CleanJet_pt[iJ]}))/(1-btag_eff);
	      }
	      else if (Jet_hadronFlavour[CleanJet_jetIdx[iJ]] == 4) {
		btag_sf *= (1-btag_eff*cset_deepJet_mujets->evaluate({systematic[i], WP, 4, abs(CleanJet_eta[iJ]), CleanJet_pt[iJ]}))/(1-btag_eff);
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
  std::string home = "/afs/cern.ch/user/v/victorr/private/mkShapesRDF/mkShapesRDF/processor/data/jsonpog-integration/POG/BTV/" + year;  
  cset = CorrectionSet::from_file(home + "/btagging.json.gz");

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
