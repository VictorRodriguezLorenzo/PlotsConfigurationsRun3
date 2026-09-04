#ifndef TTZ_LEPTONS_CC
#define TTZ_LEPTONS_CC

#include <cmath>
#include <iostream>
#include <vector>
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"

using namespace ROOT::VecOps;

static const float mEle = 0.000511;
static const float mMu  = 0.105658;
static const float mZ   = 91.1876;

// Debug flag (off by default)
static bool DEBUG_Z = false;

// Function to toggle debugging
extern "C" void setZDebug(bool x) {
    DEBUG_Z = x;
}

// --------------------------------------------------------
// Nominal lepton masses
// --------------------------------------------------------
float getMass(int pdgid) {
    int apdg = std::abs(pdgid);
    if (apdg == 11) return mEle;
    if (apdg == 13) return mMu;
    return 0.0;
}

// --------------------------------------------------------
// Build TLorentzVector with nominal mass
// --------------------------------------------------------
TLorentzVector makeLV(float pt, float eta, float phi, int pdgId) {
    TLorentzVector lv;
    float mass = getMass(pdgId);
    lv.SetPtEtaPhiM(pt, eta, phi, mass);
    return lv;
}

// --------------------------------------------------------
// Struct containing all Z information
// --------------------------------------------------------
struct ZInfo {
    int z1, z2, other;
    float mll;
};

// --------------------------------------------------------
// Compute best OSSF Z pair + remaining lepton
// --------------------------------------------------------
ZInfo computeZInfo(int nLep, RVecF pt, RVecF eta, RVecF phi, RVecI pdgId)
{
    ZInfo out = {-1, -1, -1, -1.0};

    if (DEBUG_Z) {
        std::cout << "\n========== computeZInfo DEBUG ==========\n";
        std::cout << "nLep = " << nLep << "\n";
        for (int i = 0; i < nLep; i++) {
            std::cout << "  Lep " << i
                      << ": pt=" << pt[i]
                      << " eta=" << eta[i]
                      << " phi=" << phi[i]
                      << " pdgId=" << pdgId[i]
                      << " mass=" << getMass(pdgId[i])
                      << "\n";
        }
    }

    if (nLep < 3) {
        if (DEBUG_Z) std::cout << "Not enough leptons.\n";
        return out;
    }

    float bestDiff = 1e9;
    int best_i = -1, best_j = -1;
    float best_mll = -1.0;

    // Loop over OSSF pairs
    for (int i = 0; i < nLep; i++) {
        for (int j = i + 1; j < nLep; j++) {

            if (pdgId[i] != -pdgId[j]) continue;  // require OSSF

            TLorentzVector li = makeLV(pt[i], eta[i], phi[i], pdgId[i]);
            TLorentzVector lj = makeLV(pt[j], eta[j], phi[j], pdgId[j]);

            float mll = (li + lj).M();
            float diff = std::abs(mll - mZ);

            if (DEBUG_Z) {
                std::cout << "  Pair (" << i << "," << j << ") is OSSF: "
                          << "mll=" << mll
                          << " |mll - mZ|=" << diff << "\n";
            }

            if (diff < bestDiff) {
                bestDiff = diff;
                best_i = i;
                best_j = j;
                best_mll = mll;
            }
        }
    }

    if (best_i < 0) {
        if (DEBUG_Z) std::cout << "No OSSF pair found.\n";
        return out;
    }

    out.z1 = best_i;
    out.z2 = best_j;
    out.mll = best_mll;

    if (DEBUG_Z) {
        std::cout << "Best Z pair: (" << best_i << "," << best_j
                  << ")  mll=" << best_mll << "\n";
    }

    // Remaining lepton
    for (int k = 0; k < nLep; k++) {
        if (k != best_i && k != best_j) {
            out.other = k;
            break;
        }
    }

    if (DEBUG_Z) {
        std::cout << "Remaining lepton (other) index: " << out.other << "\n";
        std::cout << "==========================================\n";
    }

    return out;
}

extern "C" {

// --------------------------------------------------------
// Aliases-accessible functions
// --------------------------------------------------------
int getZLep1Index(int nLep, RVecF pt, RVecF eta, RVecF phi, RVecI pdgId) {
    return computeZInfo(nLep, pt, eta, phi, pdgId).z1;
}

int getZLep2Index(int nLep, RVecF pt, RVecF eta, RVecF phi, RVecI pdgId) {
    return computeZInfo(nLep, pt, eta, phi, pdgId).z2;
}

int getOtherLepIndex(int nLep, RVecF pt, RVecF eta, RVecF phi, RVecI pdgId) {
    return computeZInfo(nLep, pt, eta, phi, pdgId).other;
}

float getZLepMll(int nLep, RVecF pt, RVecF eta, RVecF phi, RVecI pdgId) {
    return computeZInfo(nLep, pt, eta, phi, pdgId).mll;
}

} // extern "C"

#endif
