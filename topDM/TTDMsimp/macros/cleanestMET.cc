#include <vector>
#include "TLorentzVector.h"
#include "correction.h"
#include "ROOT/RVec.hxx"

using namespace ROOT;
using namespace ROOT::VecOps;

RVecF cleanestMET(
		  RVecF CleanJet_pt,
		  RVecF CleanJet_eta,
		  RVecF CleanJet_phi,
		  RVecF Jet_mass,
		  RVecF CleanJet_jetIdx,
		  RVecF CleanJet_mask,
		  float MET_pt,
		  float MET_phi,
		  RVecF Lepton_pt,
                  RVecF Lepton_eta,
                  RVecF Lepton_phi){

  TLorentzVector MET;
  TLorentzVector fakeJet;

  MET.SetPtEtaPhiM(MET_pt, 0.0, MET_phi, 0.0);
  for (unsigned int i=0; i<CleanJet_pt.size(); i++){
    if (CleanJet_mask[i])
      continue;

    // Remove jet in the horns from MET
    fakeJet.SetPtEtaPhiM(CleanJet_pt[i], CleanJet_eta[i], CleanJet_phi[i], Jet_mass[CleanJet_jetIdx[i]]);
    MET = MET - fakeJet;    
  }

  TLorentzVector L1;
  TLorentzVector L2;

  L1.SetPtEtaPhiM(Lepton_pt[0], Lepton_eta[0], Lepton_phi[0], 0.0);
  L2.SetPtEtaPhiM(Lepton_pt[1], Lepton_eta[1], Lepton_phi[1], 0.0);

  float dphillmet = abs(DeltaPhi((L1+L2).Phi(),MET.Phi()));
  
  return {(float)MET.Pt(), (float)MET.Phi(), dphillmet};		    
}
