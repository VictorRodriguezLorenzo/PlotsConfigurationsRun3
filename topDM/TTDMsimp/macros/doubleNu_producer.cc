#ifndef DOUBLENU_PRODUCER_CC
#define DOUBLENU_PRODUCER_CC

#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp/macros/doubleNu_producer.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "TLorentzVector.h"
#include "TMath.h"
#include "ROOT/RVec.hxx"

using namespace ROOT;
using namespace ROOT::VecOps;
using namespace nuana;

// ---- helper ----
RVecF emptyResult() {
    return RVecF{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, 0};
}

RVecF doubleNu_producer(
        int nJet,
        RVecF CleanJet_pt, RVecF CleanJet_eta, RVecF CleanJet_phi, RVecF CleanJet_mass, RVecI Jet_jetIdx,
        int nLep,
        RVecF Lep_pt, RVecF Lep_eta, RVecF Lep_phi, RVecI Lep_pdgId,
        float PuppiMET_pt, float PuppiMET_phi,
	RVecF Jet_btagger, float bAlgo_WP
        ){

        // -------------------------
        // 1) Select leptons
        // -------------------------
        if (nLep < 2) return emptyResult();

        auto leptonMass = [](int pdgId) {
            const int absId = std::abs(pdgId);
            if (absId == 11) return 0.000511f; // electron mass (GeV)
            if (absId == 13) return 0.105658f; // muon mass (GeV)
            return 0.0f;
        };

        TLorentzVector l1, l2;
        const float l1_mass = leptonMass(Lep_pdgId[0]);
        const float l2_mass = leptonMass(Lep_pdgId[1]);
        l1.SetPtEtaPhiM(Lep_pt[0], Lep_eta[0], Lep_phi[0], l1_mass);
        l2.SetPtEtaPhiM(Lep_pt[1], Lep_eta[1], Lep_phi[1], l2_mass);

        // -------------------------
        // 2) Select b-jets
        // -------------------------
        if (nJet < 2) return emptyResult();

        std::vector<int> bjet_indices;
        
        for (size_t i = 0; i < CleanJet_pt.size(); ++i) {
            int jetIdx = Jet_jetIdx[i];
        
            // Guard against corrupt indices
            if (jetIdx < 0 || jetIdx >= CleanJet_pt.size()) continue;
        
            float pt  = CleanJet_pt[i];
            float eta = CleanJet_eta[i];
            float btag = Jet_btagger[jetIdx];
        
            if (pt > 30.0 && std::abs(eta) < 2.5 && btag > bAlgo_WP)
                bjet_indices.push_back(i);
        }
        
        if (bjet_indices.size() < 2)
            return emptyResult();
        
        // Build TLorentzVectors for the first two b-jets
        int b1_idx = bjet_indices[0];
        int b2_idx = bjet_indices[1];
        
        TLorentzVector bj1, bj2;
        bj1.SetPtEtaPhiM(
            CleanJet_pt[b1_idx], CleanJet_eta[b1_idx], 
            CleanJet_phi[b1_idx], CleanJet_mass[b1_idx]
        );
        
        bj2.SetPtEtaPhiM(
            CleanJet_pt[b2_idx], CleanJet_eta[b2_idx], 
            CleanJet_phi[b2_idx], CleanJet_mass[b2_idx]
        );


        // -------------------------
        // 3) MET
        // -------------------------
        double met_x = PuppiMET_pt * std::cos(PuppiMET_phi);
        double met_y = PuppiMET_pt * std::sin(PuppiMET_phi);

        // -------------------------
        // 4) Solve
        // -------------------------
        nuana::doubleNeutrinoSolution solver(bj1, bj2, l1, l2, met_x, met_y);

        if (!solver.hasValidSolution())
            return emptyResult();

        size_t idx = 0;

        auto kin = nuana::computeEventKinematics(
            bj1, bj2, l1, l2, met_x, met_y, solver, idx
        );

        // -------------------------
        // 5) Fill result
        // -------------------------
        RVecF out(10);
        out[0] = solver.nu1_px(idx);
        out[1] = solver.nu1_py(idx);
        out[2] = solver.nu2_px(idx);
        out[3] = solver.nu2_py(idx);

        out[4] = kin.top1.Pt();
        out[5] = kin.top2.Pt();

        out[6] = kin.chel;
        out[7] = kin.dphi_ttbar;
        out[8] = kin.pdark;
        out[9] = kin.valid ? 1.0f : 0.0f;

        return out;
}

#endif
