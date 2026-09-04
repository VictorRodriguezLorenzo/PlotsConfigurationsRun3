#ifndef TOPDM_RESTFRAME_VARS_CC
#define TOPDM_RESTFRAME_VARS_CC

#include <cmath>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "ROOT/RVec.hxx"

using namespace ROOT::VecOps;

static inline bool is_good_float_rf(float x) {
    return std::isfinite(x);
}

static inline bool is_good_positive_rf(float x) {
    return std::isfinite(x) && x > 0.f;
}

static inline float lep_mass_rf(int pdgId) {
    const int apdg = std::abs(pdgId);
    if (apdg == 11) return 0.000511f;
    if (apdg == 13) return 0.105658f;
    return 0.0f;
}

static inline float dphi_abs_rf(float phi1, float phi2) {
    if (!std::isfinite(phi1) || !std::isfinite(phi2)) return -999.f;

    float dphi = std::fabs(phi1 - phi2);
    while (dphi > M_PI) dphi = std::fabs(dphi - 2.0f * M_PI);
    return dphi;
}

static inline bool make_p4_rf(
    TLorentzVector& p4,
    float pt,
    float eta,
    float phi,
    float mass
) {
    if (!std::isfinite(pt)   || pt <= 0.f) return false;
    if (!std::isfinite(eta)) return false;
    if (!std::isfinite(phi)) return false;
    if (!std::isfinite(mass)) mass = 0.f;
    if (mass < 0.f) mass = 0.f;

    p4.SetPtEtaPhiM(pt, eta, phi, mass);

    if (!std::isfinite(p4.E()))  return false;
    if (!std::isfinite(p4.Px())) return false;
    if (!std::isfinite(p4.Py())) return false;
    if (!std::isfinite(p4.Pz())) return false;

    return true;
}

static inline TLorentzVector boost_to_rest_frame_rf(
    const TLorentzVector& obj,
    const TLorentzVector& frame
) {
    TLorentzVector out = obj;

    if (!std::isfinite(frame.E()) || frame.E() <= 0.) return out;
    if (!std::isfinite(frame.M()) || frame.M() <= 1e-6) return out;

    TVector3 beta = -frame.BoostVector();

    if (!std::isfinite(beta.X())) return out;
    if (!std::isfinite(beta.Y())) return out;
    if (!std::isfinite(beta.Z())) return out;
    if (beta.Mag2() >= 1.0) return out;

    out.Boost(beta);

    if (!std::isfinite(out.E()))  return obj;
    if (!std::isfinite(out.Px())) return obj;
    if (!std::isfinite(out.Py())) return obj;
    if (!std::isfinite(out.Pz())) return obj;

    return out;
}

static inline float angle_safe_rf(
    const TLorentzVector& a,
    const TLorentzVector& b
) {
    if (a.Vect().Mag() <= 1e-6) return -999.f;
    if (b.Vect().Mag() <= 1e-6) return -999.f;

    const float angle = a.Vect().Angle(b.Vect());
    if (!std::isfinite(angle)) return -999.f;

    return angle;
}

static inline float cos_to_lab_axis_rf(
    const TLorentzVector& boosted_obj,
    const TLorentzVector& frame_lab
) {
    TVector3 axis = frame_lab.Vect();

    if (axis.Mag() <= 1e-6) return -999.f;
    if (boosted_obj.Vect().Mag() <= 1e-6) return -999.f;

    const float angle = boosted_obj.Vect().Angle(axis);
    if (!std::isfinite(angle)) return -999.f;

    const float c = std::cos(angle);
    if (!std::isfinite(c)) return -999.f;

    return c;
}

RVecF topDM_restFrame_vars(
    int nCleanJet,
    RVecF CleanJet_pt,
    RVecF CleanJet_eta,
    RVecF CleanJet_phi,
    RVecF CleanJet_mass,
    RVecI CleanJet_jetIdx,

    int nLepton,
    RVecF Lepton_pt,
    RVecF Lepton_eta,
    RVecF Lepton_phi,
    RVecI Lepton_pdgId,

    float PuppiMET_pt,
    float PuppiMET_phi,

    RVecF Jet_btagger,
    float bAlgo_WP
) {
    /*
      Output map:
       0 angle_ll_llbb_rf
       1 dphi_ll_llbb_rf
       2 cos_l1_llbb_rf
       3 cos_l2_llbb_rf

       4 angle_ll_llmet_rf
       5 dphi_ll_llmet_rf
       6 cos_l1_llmet_rf
       7 cos_l2_llmet_rf
    */

    RVecF out(8, -999.f);

    // Fully guard lepton sizes, not only nLepton.
    const int nLep = std::min({
        nLepton,
        (int)Lepton_pt.size(),
        (int)Lepton_eta.size(),
        (int)Lepton_phi.size(),
        (int)Lepton_pdgId.size()
    });

    if (nLep < 2) return out;

    TLorentzVector l1, l2;

    const bool good_l1 = make_p4_rf(
        l1,
        Lepton_pt[0],
        Lepton_eta[0],
        Lepton_phi[0],
        lep_mass_rf(Lepton_pdgId[0])
    );

    const bool good_l2 = make_p4_rf(
        l2,
        Lepton_pt[1],
        Lepton_eta[1],
        Lepton_phi[1],
        lep_mass_rf(Lepton_pdgId[1])
    );

    if (!good_l1 || !good_l2) return out;

    TLorentzVector met;

    if (
        std::isfinite(PuppiMET_pt) &&
        PuppiMET_pt >= 0.f &&
        std::isfinite(PuppiMET_phi)
    ) {
        met.SetPtEtaPhiM(PuppiMET_pt, 0.0, PuppiMET_phi, 0.0);
    } else {
        return out;
    }

    TLorentzVector llmet = l1 + l2 + met;

    // Approximate ll+MET rest-frame variables
    if (
        std::isfinite(llmet.M()) &&
        std::isfinite(llmet.E()) &&
        llmet.M() > 1e-6 &&
        llmet.E() > 0.
    ) {
        TLorentzVector l1_llmet = boost_to_rest_frame_rf(l1, llmet);
        TLorentzVector l2_llmet = boost_to_rest_frame_rf(l2, llmet);

        out[4] = angle_safe_rf(l1_llmet, l2_llmet);
        out[5] = dphi_abs_rf(l1_llmet.Phi(), l2_llmet.Phi());
        out[6] = cos_to_lab_axis_rf(l1_llmet, llmet);
        out[7] = cos_to_lab_axis_rf(l2_llmet, llmet);
    }

    // Select b jets safely.
    std::vector<int> bidx;

    const int nCJ = std::min({
        nCleanJet,
        (int)CleanJet_pt.size(),
        (int)CleanJet_eta.size(),
        (int)CleanJet_phi.size(),
        (int)CleanJet_mass.size(),
        (int)CleanJet_jetIdx.size()
    });

    for (int i = 0; i < nCJ; ++i) {
        if (!std::isfinite(CleanJet_pt[i])) continue;
        if (!std::isfinite(CleanJet_eta[i])) continue;
        if (!std::isfinite(CleanJet_phi[i])) continue;

        if (CleanJet_pt[i] < 30.f) continue;
        if (std::abs(CleanJet_eta[i]) > 2.5f) continue;

        const int jetIdx = CleanJet_jetIdx[i];

        // This is the important protection against segfaults.
        if (jetIdx < 0) continue;
        if (jetIdx >= (int)Jet_btagger.size()) continue;
        if (!std::isfinite(Jet_btagger[jetIdx])) continue;

        if (Jet_btagger[jetIdx] > bAlgo_WP) {
            bidx.push_back(i);
        }
    }

    if (bidx.size() < 2) return out;

    // Do not trust CleanJet order after JES/JER variations.
    // Sort selected b jets by the corrected CleanJet_pt.
    std::sort(
        bidx.begin(),
        bidx.end(),
        [&](int a, int b) {
            return CleanJet_pt[a] > CleanJet_pt[b];
        }
    );

    TLorentzVector b1, b2;

    const bool good_b1 = make_p4_rf(
        b1,
        CleanJet_pt[bidx[0]],
        CleanJet_eta[bidx[0]],
        CleanJet_phi[bidx[0]],
        CleanJet_mass[bidx[0]]
    );

    const bool good_b2 = make_p4_rf(
        b2,
        CleanJet_pt[bidx[1]],
        CleanJet_eta[bidx[1]],
        CleanJet_phi[bidx[1]],
        CleanJet_mass[bidx[1]]
    );

    if (!good_b1 || !good_b2) return out;

    TLorentzVector llbb = l1 + l2 + b1 + b2;

    // Visible ttbar-like rest-frame variables
    if (
        std::isfinite(llbb.M()) &&
        std::isfinite(llbb.E()) &&
        llbb.M() > 1e-6 &&
        llbb.E() > 0.
    ) {
        TLorentzVector l1_llbb = boost_to_rest_frame_rf(l1, llbb);
        TLorentzVector l2_llbb = boost_to_rest_frame_rf(l2, llbb);

        out[0] = angle_safe_rf(l1_llbb, l2_llbb);
        out[1] = dphi_abs_rf(l1_llbb.Phi(), l2_llbb.Phi());
        out[2] = cos_to_lab_axis_rf(l1_llbb, llbb);
        out[3] = cos_to_lab_axis_rf(l2_llbb, llbb);
    }

    return out;
}

#endif
