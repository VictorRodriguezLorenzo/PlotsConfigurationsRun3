#ifndef TOPDM_DNN_VARS_CC
#define TOPDM_DNN_VARS_CC

#include <cmath>
#include <vector>
#include <algorithm>
#include <initializer_list>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"

#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/computeMT2.cc"

using namespace ROOT::VecOps;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline bool is_good_float_topdm(float x) {
    return std::isfinite(x);
}

static inline bool is_good_pt_topdm(float x) {
    return std::isfinite(x) && x > 0.f;
}

static inline float dphi_abs_topdm(float phi1, float phi2) {
    if (!std::isfinite(phi1) || !std::isfinite(phi2)) return -999.f;

    float dphi = std::fmod(phi1 - phi2, 2.0f * M_PI);

    if (dphi > M_PI) dphi -= 2.0f * M_PI;
    if (dphi <= -M_PI) dphi += 2.0f * M_PI;

    return std::fabs(dphi);
}

static inline float dr_topdm(float eta1, float phi1, float eta2, float phi2) {
    if (!std::isfinite(eta1) || !std::isfinite(eta2)) return -999.f;

    const float dphi = dphi_abs_topdm(phi1, phi2);
    if (dphi < -998.f) return -999.f;

    const float deta = eta1 - eta2;
    const float dr2 = deta * deta + dphi * dphi;

    if (!std::isfinite(dr2) || dr2 < 0.f) return -999.f;

    return std::sqrt(dr2);
}

static inline float lep_mass_topdm(int pdgId) {
    const int apdg = std::abs(pdgId);
    if (apdg == 11) return 0.000511f;
    if (apdg == 13) return 0.105658f;
    return 0.0f;
}

static inline bool make_p4_topdm(
    TLorentzVector& p4,
    float pt,
    float eta,
    float phi,
    float mass
) {
    if (!std::isfinite(pt) || pt <= 0.f) return false;
    if (!std::isfinite(eta)) return false;
    if (!std::isfinite(phi)) return false;
    if (std::fabs(eta) > 10.f) return false;

    if (!std::isfinite(mass)) mass = 0.f;
    if (mass < 0.f) mass = 0.f;

    p4.SetPtEtaPhiM(pt, eta, phi, mass);

    if (!std::isfinite(p4.E()))  return false;
    if (!std::isfinite(p4.Px())) return false;
    if (!std::isfinite(p4.Py())) return false;
    if (!std::isfinite(p4.Pz())) return false;

    return true;
}

static inline float mt2_from_visible_systems_topdm(
    const TLorentzVector& v1,
    const TLorentzVector& v2,
    float met_pt,
    float met_phi
) {
    if (!std::isfinite(v1.E()) || !std::isfinite(v2.E())) return -999.f;
    if (!std::isfinite(v1.M()) || !std::isfinite(v2.M())) return -999.f;
    if (!std::isfinite(v1.Px()) || !std::isfinite(v1.Py())) return -999.f;
    if (!std::isfinite(v2.Px()) || !std::isfinite(v2.Py())) return -999.f;

    if (!std::isfinite(met_pt) || met_pt < 0.f) return -999.f;
    if (!std::isfinite(met_phi)) return -999.f;

    const double pxMiss = met_pt * std::cos(met_phi);
    const double pyMiss = met_pt * std::sin(met_phi);

    if (!std::isfinite(pxMiss) || !std::isfinite(pyMiss)) return -999.f;

    const double mt2 = asymm_mt2_lester_bisect::get_mT2(
        std::fabs(v1.M()), v1.Px(), v1.Py(),
        std::fabs(v2.M()), v2.Px(), v2.Py(),
        pxMiss, pyMiss,
        0.0, 0.0,
        0.05
    );

    if (!std::isfinite(mt2) || mt2 < 0.) return -999.f;

    return static_cast<float>(mt2);
}

static inline float transverse_mass_visible_met_topdm(
    const TLorentzVector& vis,
    float met_pt,
    float met_phi
) {
    if (!std::isfinite(vis.E()))  return -999.f;
    if (!std::isfinite(vis.Pt())) return -999.f;
    if (!std::isfinite(vis.M()))  return -999.f;

    if (!std::isfinite(met_pt) || met_pt < 0.f) return -999.f;
    if (!std::isfinite(met_phi)) return -999.f;

    const double met_px = met_pt * std::cos(met_phi);
    const double met_py = met_pt * std::sin(met_phi);

    if (!std::isfinite(met_px) || !std::isfinite(met_py)) return -999.f;

    const double mass = std::fabs(vis.M());
    const double et_vis = std::sqrt(std::max(0.0, mass * mass + vis.Pt() * vis.Pt()));

    const double px_tot = vis.Px() + met_px;
    const double py_tot = vis.Py() + met_py;

    const double mt2 = std::pow(et_vis + met_pt, 2) - px_tot * px_tot - py_tot * py_tot;

    if (!std::isfinite(mt2) || mt2 < 0.) return -999.f;

    return static_cast<float>(std::sqrt(mt2));
}

RVecF topDM_DNN_vars(
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
       0  dphi_met_ll
       1  ht
       2  st
       3  met_over_sqrt_ht
       4  met_over_st
       5  dphi_min_j_met

       6  pt_llb
       7  dphi_met_llb

       8  pt_llbb
       9  dphi_met_llbb
      10  m_llbb
      11  mT_llbb_met
      12  met_over_pt_llbb

      13  mt2_bell_l
      14  mbl_min
      15  mbl_max
      16  mtbl

      17  mbb
      18  drbb
      19  ptbb

      20  max_nonleading_btag
      21  pt_b2

      22  nForwardJet
      23  leadingForwardJet_pt
      24  leadingForwardJet_absEta
      25  deta_forwardJet_b
      26  dphi_forwardJet_met
    */

    RVecF out(27, -999.f);

    // ------------------------------------------------------------
    // Fully guard lepton access
    // ------------------------------------------------------------

    int nLep = std::min({
        nLepton,
        (int)Lepton_pt.size(),
        (int)Lepton_eta.size(),
        (int)Lepton_phi.size(),
        (int)Lepton_pdgId.size()
    });

    if (nLep < 2) return out;

    TLorentzVector l1, l2;

    const bool good_l1 = make_p4_topdm(
        l1,
        Lepton_pt[0],
        Lepton_eta[0],
        Lepton_phi[0],
        lep_mass_topdm(Lepton_pdgId[0])
    );

    const bool good_l2 = make_p4_topdm(
        l2,
        Lepton_pt[1],
        Lepton_eta[1],
        Lepton_phi[1],
        lep_mass_topdm(Lepton_pdgId[1])
    );

    if (!good_l1 || !good_l2) return out;

    if (!std::isfinite(PuppiMET_pt) || PuppiMET_pt < 0.f) return out;
    if (!std::isfinite(PuppiMET_phi)) return out;

    TLorentzVector ll = l1 + l2;

    if (std::isfinite(ll.Phi())) {
        out[0] = dphi_abs_topdm(PuppiMET_phi, ll.Phi());
    }

    // ------------------------------------------------------------
    // Fully guard CleanJet access
    // ------------------------------------------------------------

    int nCJ = std::min({
        nCleanJet,
        (int)CleanJet_pt.size(),
        (int)CleanJet_eta.size(),
        (int)CleanJet_phi.size(),
        (int)CleanJet_mass.size(),
        (int)CleanJet_jetIdx.size()
    });

    if (nCJ < 0) nCJ = 0;

    std::vector<int> bidx;

    float ht = 0.f;
    float dphiMinJMet = 999.f;

    int leadFwdIdx = -1;
    float leadFwdPt = -1.f;
    int nFwd = 0;

    for (int i = 0; i < nCJ; ++i) {
        const float pt  = CleanJet_pt[i];
        const float eta = CleanJet_eta[i];
        const float phi = CleanJet_phi[i];

        if (!std::isfinite(pt)) continue;
        if (!std::isfinite(eta)) continue;
        if (!std::isfinite(phi)) continue;
        if (std::fabs(eta) > 10.f) continue;

        if (pt > 30.f) {
            ht += pt;

            const float dphi_j_met = dphi_abs_topdm(phi, PuppiMET_phi);
            if (dphi_j_met > -998.f) {
                dphiMinJMet = std::min(dphiMinJMet, dphi_j_met);
            }
        }

        if (pt > 30.f && std::abs(eta) > 2.5f && std::abs(eta) < 4.7f) {
            nFwd++;

            if (pt > leadFwdPt) {
                leadFwdPt = pt;
                leadFwdIdx = i;
            }
        }

        const int jetIdx = CleanJet_jetIdx[i];

        // Important protection against segfaults after JES/JER shifts
        if (jetIdx < 0) continue;
        if (jetIdx >= (int)Jet_btagger.size()) continue;
        if (!std::isfinite(Jet_btagger[jetIdx])) continue;

        const float btag = Jet_btagger[jetIdx];

        if (pt > 30.f && std::abs(eta) < 2.5f && btag > bAlgo_WP) {
            bidx.push_back(i);
        }
    }

    // Do not trust CleanJet ordering after JES/JER variations.
    // Sort selected b jets by corrected CleanJet_pt.
    std::sort(
        bidx.begin(),
        bidx.end(),
        [&](int a, int b) {
            return CleanJet_pt[a] > CleanJet_pt[b];
        }
    );

    // ------------------------------------------------------------
    // Inclusive jet / MET variables
    // ------------------------------------------------------------

    out[1] = ht;
    out[2] = ht + Lepton_pt[0] + Lepton_pt[1] + PuppiMET_pt;
    out[3] = (ht > 0.f) ? PuppiMET_pt / std::sqrt(ht) : -999.f;
    out[4] = (out[2] > 0.f) ? PuppiMET_pt / out[2] : -999.f;
    out[5] = (dphiMinJMet < 998.f) ? dphiMinJMet : -999.f;

    out[22] = (float)nFwd;

    if (leadFwdIdx >= 0 && leadFwdIdx < nCJ) {
        out[23] = CleanJet_pt[leadFwdIdx];
        out[24] = std::abs(CleanJet_eta[leadFwdIdx]);
        out[26] = dphi_abs_topdm(CleanJet_phi[leadFwdIdx], PuppiMET_phi);
    }

    // ------------------------------------------------------------
    // One-b variables
    // ------------------------------------------------------------

    if (bidx.size() >= 1) {
        const int ib1 = bidx[0];

        if (ib1 >= 0 && ib1 < nCJ) {
            TLorentzVector b1;

            const bool good_b1 = make_p4_topdm(
                b1,
                CleanJet_pt[ib1],
                CleanJet_eta[ib1],
                CleanJet_phi[ib1],
                CleanJet_mass[ib1]
            );

            if (good_b1) {
                TLorentzVector llb = ll + b1;

                if (std::isfinite(llb.Pt())) {
                    out[6] = llb.Pt();
                }

                if (std::isfinite(llb.Phi())) {
                    out[7] = dphi_abs_topdm(PuppiMET_phi, llb.Phi());
                }

                if (leadFwdIdx >= 0 && leadFwdIdx < nCJ) {
                    out[25] = std::abs(CleanJet_eta[leadFwdIdx] - CleanJet_eta[ib1]);
                }
            }

            // mbl min/max over all selected b jets
            float mblMin = 999999.f;
            float mblMax = -999.f;

            for (unsigned int ib = 0; ib < bidx.size(); ++ib) {
                const int idx = bidx[ib];

                if (idx < 0 || idx >= nCJ) continue;

                TLorentzVector b;

                const bool good_b = make_p4_topdm(
                    b,
                    CleanJet_pt[idx],
                    CleanJet_eta[idx],
                    CleanJet_phi[idx],
                    CleanJet_mass[idx]
                );

                if (!good_b) continue;

                const float m1 = (b + l1).M();
                const float m2 = (b + l2).M();

                if (std::isfinite(m1)) {
                    mblMin = std::min(mblMin, m1);
                    mblMax = std::max(mblMax, m1);
                }

                if (std::isfinite(m2)) {
                    mblMin = std::min(mblMin, m2);
                    mblMax = std::max(mblMax, m2);
                }
            }

            if (mblMin < 999998.f) out[14] = mblMin;
            if (mblMax > -998.f)   out[15] = mblMax;

            // MT2^{b l, l}; useful for 1b / tW-like topology
            float bestBellL = 999999.f;

            for (unsigned int ib = 0; ib < bidx.size(); ++ib) {
                const int idx = bidx[ib];

                if (idx < 0 || idx >= nCJ) continue;

                TLorentzVector b;

                const bool good_b = make_p4_topdm(
                    b,
                    CleanJet_pt[idx],
                    CleanJet_eta[idx],
                    CleanJet_phi[idx],
                    CleanJet_mass[idx]
                );

                if (!good_b) continue;

                const float v1 = mt2_from_visible_systems_topdm(
                    b + l1,
                    l2,
                    PuppiMET_pt,
                    PuppiMET_phi
                );

                const float v2 = mt2_from_visible_systems_topdm(
                    b + l2,
                    l1,
                    PuppiMET_pt,
                    PuppiMET_phi
                );

                if (v1 > -998.f) bestBellL = std::min(bestBellL, v1);
                if (v2 > -998.f) bestBellL = std::min(bestBellL, v2);
            }

            if (bestBellL < 999998.f) out[13] = bestBellL;

            // max b-tag score among central jets except the leading selected b
            float maxNonLeadingBtag = -999.f;

            for (int i = 0; i < nCJ; ++i) {
                if (i == ib1) continue;

                if (!std::isfinite(CleanJet_pt[i])) continue;
                if (!std::isfinite(CleanJet_eta[i])) continue;

                if (CleanJet_pt[i] < 20.f) continue;
                if (std::abs(CleanJet_eta[i]) > 2.5f) continue;

                const int jetIdx = CleanJet_jetIdx[i];

                if (jetIdx < 0) continue;
                if (jetIdx >= (int)Jet_btagger.size()) continue;
                if (!std::isfinite(Jet_btagger[jetIdx])) continue;

                maxNonLeadingBtag = std::max(maxNonLeadingBtag, Jet_btagger[jetIdx]);
            }

            out[20] = maxNonLeadingBtag;
        }
    }

    // ------------------------------------------------------------
    // Two-b variables
    // ------------------------------------------------------------

    if (bidx.size() >= 2) {
        const int ib1 = bidx[0];
        const int ib2 = bidx[1];

        if (ib1 >= 0 && ib1 < nCJ && ib2 >= 0 && ib2 < nCJ) {
            out[21] = CleanJet_pt[ib2];

            TLorentzVector b1, b2;

            const bool good_b1 = make_p4_topdm(
                b1,
                CleanJet_pt[ib1],
                CleanJet_eta[ib1],
                CleanJet_phi[ib1],
                CleanJet_mass[ib1]
            );

            const bool good_b2 = make_p4_topdm(
                b2,
                CleanJet_pt[ib2],
                CleanJet_eta[ib2],
                CleanJet_phi[ib2],
                CleanJet_mass[ib2]
            );

            if (good_b1 && good_b2) {
                TLorentzVector bb = b1 + b2;
                TLorentzVector llbb = ll + bb;

                if (std::isfinite(llbb.Pt())) {
                    out[8] = llbb.Pt();
                }

                if (std::isfinite(llbb.Phi())) {
                    out[9] = dphi_abs_topdm(PuppiMET_phi, llbb.Phi());
                }

                if (std::isfinite(llbb.M())) {
                    out[10] = llbb.M();
                }

                out[11] = transverse_mass_visible_met_topdm(
                    llbb,
                    PuppiMET_pt,
                    PuppiMET_phi
                );

                out[12] = (std::isfinite(llbb.Pt()) && llbb.Pt() > 0.f)
                    ? PuppiMET_pt / llbb.Pt()
                    : -999.f;

                if (std::isfinite(bb.M())) {
                    out[17] = bb.M();
                }

                out[18] = dr_topdm(
                    CleanJet_eta[ib1],
                    CleanJet_phi[ib1],
                    CleanJet_eta[ib2],
                    CleanJet_phi[ib2]
                );

                if (std::isfinite(bb.Pt())) {
                    out[19] = bb.Pt();
                }
            }

            // Exact m_bl minimax, generalized over all b-jet pairs
            float bestMtbl = 999999.f;

            for (unsigned int a = 0; a < bidx.size(); ++a) {
                for (unsigned int b = a + 1; b < bidx.size(); ++b) {
                    const int ia = bidx[a];
                    const int ib = bidx[b];

                    if (ia < 0 || ia >= nCJ) continue;
                    if (ib < 0 || ib >= nCJ) continue;

                    TLorentzVector ba, bb_;

                    const bool good_ba = make_p4_topdm(
                        ba,
                        CleanJet_pt[ia],
                        CleanJet_eta[ia],
                        CleanJet_phi[ia],
                        CleanJet_mass[ia]
                    );

                    const bool good_bb = make_p4_topdm(
                        bb_,
                        CleanJet_pt[ib],
                        CleanJet_eta[ib],
                        CleanJet_phi[ib],
                        CleanJet_mass[ib]
                    );

                    if (!good_ba || !good_bb) continue;

                    const float m_a_l1 = (ba  + l1).M();
                    const float m_b_l2 = (bb_ + l2).M();
                    const float m_a_l2 = (ba  + l2).M();
                    const float m_b_l1 = (bb_ + l1).M();

                    if (
                        !std::isfinite(m_a_l1) ||
                        !std::isfinite(m_b_l2) ||
                        !std::isfinite(m_a_l2) ||
                        !std::isfinite(m_b_l1)
                    ) {
                        continue;
                    }

                    const float pairingA = std::max(m_a_l1, m_b_l2);
                    const float pairingB = std::max(m_a_l2, m_b_l1);

                    bestMtbl = std::min(bestMtbl, std::min(pairingA, pairingB));
                }
            }

            if (bestMtbl < 999998.f) out[16] = bestMtbl;
        }
    }

    return out;
}

#endif
