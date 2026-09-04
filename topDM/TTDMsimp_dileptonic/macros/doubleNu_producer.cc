#ifndef DOUBLENU_PRODUCER_CC
#define DOUBLENU_PRODUCER_CC

#include "doubleNu_producer.h"
#include "topDM_spin_observables.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include "TLorentzVector.h"
#include "TMath.h"
#include "ROOT/RVec.hxx"

using namespace ROOT;
using namespace ROOT::VecOps;
using namespace nuana;

namespace {

constexpr int kDoubleNuOutputSize = 25;

RVecF emptyResult() {
    RVecF out(kDoubleNuOutputSize, NAN);
    out[9] = 0.f;   // valid flag
    return out;
}


} // namespace

RVecF doubleNu_producer(
        int nCleanJet,
        RVecF CleanJet_pt, RVecF CleanJet_eta, RVecF CleanJet_phi, RVecF CleanJet_mass, RVecI CleanJet_jetIdx,
        int nLep,
        RVecF Lep_pt, RVecF Lep_eta, RVecF Lep_phi, RVecI Lep_pdgId,
        float PuppiMET_pt, float PuppiMET_phi,
        RVecF Jet_btagger, float bAlgo_WP
        ){

        // -------------------------
        // 1) Select leptons
        // -------------------------
        const int nSafeLep = std::min({
            nLep,
            static_cast<int>(Lep_pt.size()),
            static_cast<int>(Lep_eta.size()),
            static_cast<int>(Lep_phi.size()),
            static_cast<int>(Lep_pdgId.size())
        });
        if (nSafeLep < 2) return emptyResult();

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
        const int nSafeCleanJet = std::min({
            nCleanJet,
            static_cast<int>(CleanJet_pt.size()),
            static_cast<int>(CleanJet_eta.size()),
            static_cast<int>(CleanJet_phi.size()),
            static_cast<int>(CleanJet_mass.size()),
            static_cast<int>(CleanJet_jetIdx.size())
        });
        if (nSafeCleanJet < 2) return emptyResult();

        std::vector<int> bjet_indices;
        for (int i = 0; i < nSafeCleanJet; ++i) {
            const int jetIdx = CleanJet_jetIdx[i];
            if (jetIdx < 0 || jetIdx >= static_cast<int>(Jet_btagger.size())) continue;

            const float pt = CleanJet_pt[i];
            const float eta = CleanJet_eta[i];
            const float btag = Jet_btagger[jetIdx];
            if (!std::isfinite(pt) || !std::isfinite(eta) || !std::isfinite(btag)) continue;

            if (pt > 30.0 && std::abs(eta) < 2.5 && btag > bAlgo_WP)
                bjet_indices.push_back(i);
        }

        if (bjet_indices.size() < 2) return emptyResult();

        // Keep the original b-jet choice (first two selected tagged jets).
        const int b1_idx = bjet_indices[0];
        const int b2_idx = bjet_indices[1];

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
        const double met_x = PuppiMET_pt * std::cos(PuppiMET_phi);
        const double met_y = PuppiMET_pt * std::sin(PuppiMET_phi);

        // -------------------------
        // 4) Solve
        // -------------------------
        nuana::doubleNeutrinoSolution solver(bj1, bj2, l1, l2, met_x, met_y);
        const size_t idx = 0;

        // Reconstruct the event using the pairing selected by the solver and
        // the full neutrino momenta (including pz).
        auto kin = nuana::computeEventKinematics(
            bj1, bj2, l1, l2, met_x, met_y, solver, idx
        );

        if (!kin.valid) return emptyResult();

        // -------------------------
        // 5) Fill reconstructed quantities (indices 0-9)
        // -------------------------
        RVecF out(kDoubleNuOutputSize, NAN);
        out[0] = solver.nu1_px(idx);
        out[1] = solver.nu1_py(idx);
        out[2] = solver.nu2_px(idx);
        out[3] = solver.nu2_py(idx);
        out[4] = kin.top1.Pt();
        out[5] = kin.top2.Pt();
        out[6] = kin.chel;
        out[7] = kin.dphi_ttbar;
        out[8] = kin.pdark;
        out[9] = 1.f;

        // -------------------------
        // 6) Spin observables from the SAME corrected top reconstruction
        // -------------------------
        const bool swapped = solver.selectedPairingSwapped();
        const TLorentzVector& pairedL1 = swapped ? l2 : l1;
        const TLorentzVector& pairedL2 = swapped ? l1 : l2;
        const int pairedPdg1 = swapped ? Lep_pdgId[1] : Lep_pdgId[0];
        const int pairedPdg2 = swapped ? Lep_pdgId[0] : Lep_pdgId[1];

        // PDG convention: l+ has negative pdgId, l- has positive pdgId.
        const TLorentzVector* top = nullptr;
        const TLorentzVector* antitop = nullptr;
        const TLorentzVector* lplus = nullptr;
        const TLorentzVector* lminus = nullptr;

        if (pairedPdg1 < 0 && pairedPdg2 > 0) {
            top = &kin.top1;
            antitop = &kin.top2;
            lplus = &pairedL1;
            lminus = &pairedL2;
        } else if (pairedPdg2 < 0 && pairedPdg1 > 0) {
            top = &kin.top2;
            antitop = &kin.top1;
            lplus = &pairedL2;
            lminus = &pairedL1;
        } else {
            return out;
        }

        const RVecF spin = topdm::spinObservables(*top, *antitop, *lplus, *lminus);
        out[10] = spin[topdm::kCosPhi];
        out[11] = spin[topdm::kDevt];
        out[12] = spin[topdm::kXkk];
        out[13] = spin[topdm::kXrr];
        out[14] = spin[topdm::kXnn];
        out[15] = spin[topdm::kXkr];
        out[16] = spin[topdm::kXkn];
        out[17] = spin[topdm::kXrk];
        out[18] = spin[topdm::kXrn];
        out[19] = spin[topdm::kXnk];
        out[20] = spin[topdm::kXnr];
        out[21] = spin[topdm::kMttReco];
        out[22] = spin[topdm::kPtttReco];
        out[23] = spin[topdm::kYttReco];
        out[24] = spin[topdm::kAbsDeltaYtt];

        return out;
}
#endif
