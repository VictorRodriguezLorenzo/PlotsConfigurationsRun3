#ifndef BTAGSFLIGHT
#define BTAGSFLIGHT
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

class btagSFlight {

  public:
     btagSFlight(TString eff_map, const string year, TString algo_extension);
    ~btagSFlight();

    TH2F* h_ljet_eff;
  
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
      
      auto cset_deepJet_mujets  = cset->at("robustParticleTransformer_light");
      auto cset_deepJet_wps     = cset->at("robustParticleTransformer_wp_values");      
      
      for (long unsigned int i=0; i<systematic.size(); i++){
	
	float btag_sf    = 1.;
	for (unsigned iJ{0}; iJ != nCleanJet; ++iJ) {
	  
	  if (CleanJet_pt[iJ] <= 30. || abs(CleanJet_eta[iJ]) >= 2.5) continue;
	  if (Jet_btag[CleanJet_jetIdx[iJ]] > cset_deepJet_wps->evaluate({WP})) {
	    if (Jet_hadronFlavour[CleanJet_jetIdx[iJ]] == 0) {
	      btag_sf *= cset_deepJet_mujets->evaluate({systematic[i], WP, 0, abs(CleanJet_eta[iJ]), CleanJet_pt[iJ]});
	    }
	    else {
	      continue;
	    }
	  }
	  else {
	    float btag_eff = getEff(CleanJet_pt[iJ], CleanJet_eta[iJ], Jet_hadronFlavour[CleanJet_jetIdx[iJ]]);
	    if (btag_eff == 1.) {
	      continue;
	    }
	    else {
	      if (Jet_hadronFlavour[CleanJet_jetIdx[iJ]] == 0) {
		btag_sf *= (1-btag_eff*cset_deepJet_mujets->evaluate({systematic[i], WP, 0, abs(CleanJet_eta[iJ]), CleanJet_pt[iJ]}))/(1-btag_eff);
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

btagSFlight::btagSFlight(TString eff_map, const string year, TString algo_extension = "") {
  std::string home = "/afs/cern.ch/user/v/victorr/private/mkShapesRDF/mkShapesRDF/processor/data/jsonpog-integration/POG/BTV/" + year;
  cset = CorrectionSet::from_file(home + "/btagging.json.gz");

  // deepJet -> ljet_eff
  // PNet -> ljet_pnet_eff
  // RobPartT -> ljet_parT_eff

  TFile *reff = TFile::Open(eff_map, "READ");
  h_ljet_eff  = (TH2F*)reff->Get("ljet"+algo_extension+"_eff")->Clone();
  h_ljet_eff->SetDirectory(0);
  reff->Close();
}

float btagSFlight::getEff(float pt, float eta, int flavour) {
  int xbin, ybin;
  float eff;
  if (flavour == 0) {
    xbin = h_ljet_eff->GetXaxis()->FindBin(pt);
    ybin = h_ljet_eff->GetYaxis()->FindBin(eta);
    eff  = h_ljet_eff->GetBinContent(xbin, ybin);
  }
  else {
    eff   = 1.;
  }
  return eff;
}

btagSFlight::~btagSFlight(){
  std::cout << "Cleaning up memory" << std::endl;
  h_ljet_eff->Delete();
}
 
#endif
