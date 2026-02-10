#ifndef TOPNESS_PRODUCER_MINUIT_CC
#define TOPNESS_PRODUCER_MINUIT_CC

#include <cmath>
#include <vector>
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace ROOT;
using namespace ROOT::VecOps;

// --------------------
// Topness functor
// --------------------
struct TopnessFunctor {

    TLorentzVector b1, b2, lep;
    double nu_px, nu_py;

    double operator()(const double *x) {
        const double nu_pz = x[0];

        const double mW  = 80.379;
        const double mt  = 172.5;
        const double aW  = 5.0;
        const double at  = 15.0;
        const double aCM = 1000.0;

        TLorentzVector nu;
        double nu_E = std::sqrt(nu_px*nu_px + nu_py*nu_py + nu_pz*nu_pz);
        nu.SetPxPyPzE(nu_px, nu_py, nu_pz, nu_E);

        TLorentzVector W    = lep + nu;
        TLorentzVector top1 = b1 + W;
        TLorentzVector top2 = b2 + W;
        TLorentzVector sum  = lep + nu + b1 + b2;

        double S =
            std::pow(mW*mW - W.M2(), 2) / std::pow(aW, 4) +
            std::pow(mt*mt - top1.M2(), 2) / std::pow(at, 4) +
            std::pow(mt*mt - top2.M2(), 2) / std::pow(at, 4) +
            std::pow(4*mt*mt - sum.M2(), 2) / std::pow(aCM, 4);

        return S;
    }
};

// --------------------
// Producer
// --------------------
RVecF topness_producer_minuit(
    int nCleanJet,
    RVecF CleanJet_pt, RVecF CleanJet_eta, RVecF CleanJet_phi, RVecF CleanJet_mass, RVecI CleanJet_jetIdx,
    int nLep,
    RVecF Lep_pt, RVecF Lep_eta, RVecF Lep_phi, RVecI Lep_pdgId,
    float PuppiMET_pt, float PuppiMET_phi,
    RVecF Jet_btagger, float bAlgo_WP
) {

    if (nLep < 1 || nCleanJet < 2)
        return RVecF{NAN};

    // --------------------
    // Lepton
    // --------------------
    double lepMass = (std::abs(Lep_pdgId[0]) == 13 ? 0.105658 : 0.000511);
    TLorentzVector lep;
    lep.SetPtEtaPhiM(Lep_pt[0], Lep_eta[0], Lep_phi[0], lepMass);

    // --------------------
    // b-jets
    // --------------------
    std::vector<int> bjets;
    for (size_t i = 0; i < CleanJet_pt.size(); ++i) {
        int idx = CleanJet_jetIdx[i];
        if (idx < 0 || idx >= (int)Jet_btagger.size()) continue;
        if (CleanJet_pt[i] > 20.0 &&
            std::abs(CleanJet_eta[i]) < 2.5 &&
            Jet_btagger[idx] > bAlgo_WP)
            bjets.push_back(i);
    }

    if (bjets.size() < 2)
        return RVecF{NAN};

    TLorentzVector b1, b2;
    b1.SetPtEtaPhiM(
        CleanJet_pt[bjets[0]], CleanJet_eta[bjets[0]],
        CleanJet_phi[bjets[0]], CleanJet_mass[bjets[0]]
    );
    b2.SetPtEtaPhiM(
        CleanJet_pt[bjets[1]], CleanJet_eta[bjets[1]],
        CleanJet_phi[bjets[1]], CleanJet_mass[bjets[1]]
    );

    // --------------------
    // MET
    // --------------------
    double nu_px = PuppiMET_pt * std::cos(PuppiMET_phi);
    double nu_py = PuppiMET_pt * std::sin(PuppiMET_phi);

    // --------------------
    // Build functor
    // --------------------
    TopnessFunctor f;
    f.b1 = b1;
    f.b2 = b2;
    f.lep = lep;
    f.nu_px = nu_px;
    f.nu_py = nu_py;

    ROOT::Math::Functor functor(f, 1);

    // --------------------
    // Minuit
    // --------------------
    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad")
    );

    min->SetFunction(functor);
    min->SetLimitedVariable(0, "nu_pz", 0.0, 10.0, -1000.0, 1000.0);
    min->SetMaxFunctionCalls(10000);
    min->SetTolerance(1e-3);
    min->Minimize();

    double minS = min->MinValue();

    if (!std::isfinite(minS) || minS <= 0)
        return RVecF{NAN};

    return RVecF{static_cast<float>(std::log(minS))};
}

#endif

