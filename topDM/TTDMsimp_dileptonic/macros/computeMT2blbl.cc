#ifndef COMPUTE_MT2BLBL_CC
#define COMPUTE_MT2BLBL_CC

#include <algorithm>
#include <cmath>
#include <vector>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"

// Reuse existing Lester MT2 implementation
#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/computeMT2.cc"

using ROOT::VecOps::RVec;

static inline bool _isGoodFloat(double x) {
  return std::isfinite(x);
}

static inline double _lepMassFromPdgId(int pdgId) {
  const int apdg = std::abs(pdgId);
  if (apdg == 13) return 0.105658;   // muon mass
  if (apdg == 11) return 0.000511;   // electron mass
  return 0.0;
}

static inline float _mt2TwoVisibleSystems(
    const TLorentzVector& vis1,
    const TLorentzVector& vis2,
    float met_pt,
    float met_phi
) {
  if (!_isGoodFloat(vis1.M()) || !_isGoodFloat(vis1.Px()) || !_isGoodFloat(vis1.Py()) ||
      !_isGoodFloat(vis2.M()) || !_isGoodFloat(vis2.Px()) || !_isGoodFloat(vis2.Py()) ||
      !_isGoodFloat(met_pt) || !_isGoodFloat(met_phi)) {
    return -999.f;
  }

  const double pxMiss = met_pt * std::cos(met_phi);
  const double pyMiss = met_pt * std::sin(met_phi);

  const double mt2 = asymm_mt2_lester_bisect::get_mT2(
      std::fabs(vis1.M()), vis1.Px(), vis1.Py(),
      std::fabs(vis2.M()), vis2.Px(), vis2.Py(),
      pxMiss, pyMiss,
      0.0, 0.0,
      0.05   // precision in GeV; faster than 0 and more than enough for plotting/MVA
  );

  if (!_isGoodFloat(mt2) || mt2 < 0) return -999.f;

  return static_cast<float>(mt2);
}

static inline float _computeMT2blbl_fixedPair(
    float l1_pt, float l1_eta, float l1_phi, int l1_pdgId,
    float l2_pt, float l2_eta, float l2_phi, int l2_pdgId,
    float b1_pt, float b1_eta, float b1_phi, float b1_mass,
    float b2_pt, float b2_eta, float b2_phi, float b2_mass,
    float met_pt, float met_phi
) {
  if (l1_pt <= 0 || l2_pt <= 0 || b1_pt <= 0 || b2_pt <= 0) return -999.f;

  TLorentzVector l1, l2, b1, b2;

  l1.SetPtEtaPhiM(l1_pt, l1_eta, l1_phi, _lepMassFromPdgId(l1_pdgId));
  l2.SetPtEtaPhiM(l2_pt, l2_eta, l2_phi, _lepMassFromPdgId(l2_pdgId));

  b1.SetPtEtaPhiM(b1_pt, b1_eta, b1_phi, std::max(0.f, b1_mass));
  b2.SetPtEtaPhiM(b2_pt, b2_eta, b2_phi, std::max(0.f, b2_mass));

  // Pairing A: (l1+b1, l2+b2)
  TLorentzVector visA1 = l1 + b1;
  TLorentzVector visA2 = l2 + b2;
  const float mt2_A = _mt2TwoVisibleSystems(visA1, visA2, met_pt, met_phi);

  // Pairing B: (l1+b2, l2+b1)
  TLorentzVector visB1 = l1 + b2;
  TLorentzVector visB2 = l2 + b1;
  const float mt2_B = _mt2TwoVisibleSystems(visB1, visB2, met_pt, met_phi);

  if (mt2_A < 0 && mt2_B < 0) return -999.f;
  if (mt2_A < 0) return mt2_B;
  if (mt2_B < 0) return mt2_A;

  // Standard choice: take the smaller pairing
  return std::min(mt2_A, mt2_B);
}

float computeMT2blbl(
    int nCleanJet,
    const RVec<float>& CleanJet_pt,
    const RVec<float>& CleanJet_eta,
    const RVec<float>& CleanJet_phi,
    const RVec<float>& CleanJet_mass,
    const RVec<int>& CleanJet_jetIdx,

    int nLepton,
    const RVec<float>& Lepton_pt,
    const RVec<float>& Lepton_eta,
    const RVec<float>& Lepton_phi,
    const RVec<int>& Lepton_pdgId,

    float PuppiMET_pt,
    float PuppiMET_phi,

    const RVec<float>& Jet_btag,
    float btagWP
) {
  if (nLepton < 2) return -999.f;
  if ((int)Lepton_pt.size() < 2 || (int)Lepton_eta.size() < 2 ||
      (int)Lepton_phi.size() < 2 || (int)Lepton_pdgId.size() < 2) {
    return -999.f;
  }

  if (nCleanJet < 2) return -999.f;

  std::vector<int> bCleanIdx;

  const int nCJ = std::min({
      nCleanJet,
      (int)CleanJet_pt.size(),
      (int)CleanJet_eta.size(),
      (int)CleanJet_phi.size(),
      (int)CleanJet_mass.size(),
      (int)CleanJet_jetIdx.size()
  });

  for (int i = 0; i < nCJ; ++i) {
    const int jetIdx = CleanJet_jetIdx[i];

    if (jetIdx < 0) continue;
    if (jetIdx >= (int)Jet_btag.size()) continue;

    if (Jet_btag[jetIdx] > btagWP) {
      bCleanIdx.push_back(i);
    }
  }

  // Strict MT2blbl needs two b-tagged jets
  if (bCleanIdx.size() < 2) return -999.f;

  float best = 1e9;
  bool found = false;

  // If there are >2 b-jets, take the minimum over all b-jet pairs.
  // This protects the ttbar endpoint against wrong b choices.
  for (unsigned int a = 0; a < bCleanIdx.size(); ++a) {
    for (unsigned int b = a + 1; b < bCleanIdx.size(); ++b) {
      const int ib1 = bCleanIdx[a];
      const int ib2 = bCleanIdx[b];

      const float val = _computeMT2blbl_fixedPair(
          Lepton_pt[0], Lepton_eta[0], Lepton_phi[0], Lepton_pdgId[0],
          Lepton_pt[1], Lepton_eta[1], Lepton_phi[1], Lepton_pdgId[1],

          CleanJet_pt[ib1], CleanJet_eta[ib1], CleanJet_phi[ib1], CleanJet_mass[ib1],
          CleanJet_pt[ib2], CleanJet_eta[ib2], CleanJet_phi[ib2], CleanJet_mass[ib2],

          PuppiMET_pt,
          PuppiMET_phi
      );

      if (val > -998.f && val < best) {
        best = val;
        found = true;
      }
    }
  }

  if (!found) return -999.f;

  return best;
}

#endif
