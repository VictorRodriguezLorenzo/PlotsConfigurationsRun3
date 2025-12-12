#include <vector>
#include "TLorentzVector.h"
#include "correction.h"
#include "ROOT/RVec.hxx"

using namespace ROOT;
using namespace ROOT::VecOps;

RVecF cleanestKinematics(
			 RVecF CleanJet_pt,
			 RVecF CleanJet_eta,
			 RVecF CleanJet_phi,
			 RVecF Jet_mass,
			 RVecF CleanJet_jetIdx,
			 float MET_pt,
			 float MET_phi){


  TLorentzVector J1;
  TLorentzVector J2;

  if (CleanJet_pt.size()<=1)
    return {-999.9,-999.9,-999.9,-999.9,-999.9};
  
  J1.SetPtEtaPhiM(CleanJet_pt[0], CleanJet_eta[0], CleanJet_phi[0], Jet_mass[CleanJet_jetIdx[0]]);
  J2.SetPtEtaPhiM(CleanJet_pt[1], CleanJet_eta[1], CleanJet_phi[1], Jet_mass[CleanJet_jetIdx[1]]);

  float mjj = (J1+J2).M();
  float dphijj = DeltaPhi(CleanJet_phi[0], CleanJet_phi[1]);;
  float drjj = DeltaR(CleanJet_eta[0], CleanJet_eta[1], CleanJet_phi[0], CleanJet_phi[1]);
  float detajj = abs(CleanJet_eta[0] - CleanJet_eta[1]);
  float dphijjmet = std::abs(TVector2::Phi_mpi_pi((J1 + J2).Phi() - MET_phi));

  return {mjj, dphijj, drjj, detajj, dphijjmet};
  
}
