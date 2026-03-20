#ifndef TOPNESS_PRODUCER_MINUIT_CC
#define TOPNESS_PRODUCER_MINUIT_CC

#include <cmath>
#include <vector>
#include <memory>
#include <limits>
#include <algorithm>

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

        const double pWx   = x[0];
        const double pWy   = x[1];
        const double pWz   = x[2];
        const double nu_pz = x[3];

        const double mW  = 80.379;
        const double mt  = 172.5;

        const double aW  = 5.0;
        const double at  = 15.0;
        const double aCM = 1000.0;

        // --------------------
        // neutrino
        // --------------------
        TLorentzVector nu;
        double nu_E = std::sqrt(nu_px*nu_px + nu_py*nu_py + nu_pz*nu_pz);
        nu.SetPxPyPzE(nu_px, nu_py, nu_pz, nu_E);

        // --------------------
        // lost W boson
        // --------------------
        TLorentzVector W;
        double EW = std::sqrt(pWx*pWx + pWy*pWy + pWz*pWz + mW*mW);
        W.SetPxPyPzE(pWx, pWy, pWz, EW);

        // --------------------
        // top hypotheses
        // --------------------
        TLorentzVector top1 = b1 + lep + nu;
        TLorentzVector top2 = b2 + W;

        // --------------------
        // CM system (5 particles)
        // --------------------
        TLorentzVector sum = b1 + b2 + lep + nu + W;

        // --------------------
        // S function
        // --------------------
        double S =
            std::pow(mW*mW - W.M2(), 2) / std::pow(aW,4) +
            std::pow(mt*mt - top1.M2(), 2) / std::pow(at,4) +
            std::pow(mt*mt - top2.M2(), 2) / std::pow(at,4) +
            std::pow(4*mt*mt - sum.M2(), 2) / std::pow(aCM,4);

        return S;
    }
};

// --------------------
// Producer
// --------------------
double topness_producer_minuit(

    int nCleanJet,
    RVecF CleanJet_pt,
    RVecF CleanJet_eta,
    RVecF CleanJet_phi,
    RVecF CleanJet_mass,

    int nLep,
    RVecF Lep_pt,
    RVecF Lep_eta,
    RVecF Lep_phi,
    RVecI Lep_pdgId,

    float PuppiMET_pt,
    float PuppiMET_phi,

    RVecF Jet_btagger,
    float bAlgo_WP
) {

    if (nLep < 1 || nCleanJet < 2)
        return NAN;

    // --------------------
    // Lepton
    // --------------------
    double lepMass = (std::abs(Lep_pdgId[0]) == 13 ? 0.105658 : 0.000511);

    TLorentzVector lep;
    lep.SetPtEtaPhiM(
        Lep_pt[0],
        Lep_eta[0],
        Lep_phi[0],
        lepMass
    );

    // --------------------
    // select b jets
    // --------------------
    std::vector<int> bjets;
    for (size_t i = 0; i < CleanJet_pt.size(); ++i) {
        if (i < 0 || i >= (int)Jet_btagger.size()) continue;

        if (CleanJet_pt[i] > 30.0 &&
            std::abs(CleanJet_eta[i]) < 2.5 &&
            Jet_btagger[i] > bAlgo_WP)
        {
            bjets.push_back(i);
        }
    }

    // --------------------
    // MET equals neutrino px,py
    // --------------------
    double nu_px = PuppiMET_pt * std::cos(PuppiMET_phi);
    double nu_py = PuppiMET_pt * std::sin(PuppiMET_phi);

    // --------------------
    // build jet pairs
    // --------------------
    std::vector<std::pair<int,int>> jetPairs;

    if (bjets.size() >= 2) {
        jetPairs.push_back({bjets[0], bjets[1]});
    } else if (bjets.size() == 1) {
        // only one b-jet: pair it with the two hardest untagged jets
        std::vector<int> untagged;
        for (size_t i = 0; i < CleanJet_pt.size(); ++i) {
            if ((int)i == bjets[0]) continue;
            if (CleanJet_pt[i] > 30.0 && std::abs(CleanJet_eta[i]) < 2.5)
                untagged.push_back(i);
        }

        // sort untagged jets by pt descending
        std::sort(untagged.begin(), untagged.end(), [&](int a,int b){
            return CleanJet_pt[a] > CleanJet_pt[b];
        });

        for (int i = 0; i < std::min(2,(int)untagged.size()); ++i) {
            jetPairs.push_back({bjets[0], untagged[i]});
        }

    } else {
        // no b-jets: cannot proceed
        return NAN;
    }

    double bestS = std::numeric_limits<double>::infinity();

    // --------------------
    // loop over jetPairs and permutations
    // --------------------
    for (auto &[i1,i2] : jetPairs) {

        TLorentzVector jet1, jet2;
        jet1.SetPtEtaPhiM(CleanJet_pt[i1], CleanJet_eta[i1], CleanJet_phi[i1], CleanJet_mass[i1]);
        jet2.SetPtEtaPhiM(CleanJet_pt[i2], CleanJet_eta[i2], CleanJet_phi[i2], CleanJet_mass[i2]);

        for (int perm = 0; perm < 2; perm++) {

            TLorentzVector b1 = (perm == 0 ? jet1 : jet2);
            TLorentzVector b2 = (perm == 0 ? jet2 : jet1);

            TopnessFunctor f;
            f.b1 = b1;
            f.b2 = b2;
            f.lep = lep;
            f.nu_px = nu_px;
            f.nu_py = nu_py;

            ROOT::Math::Functor functor(f,4);

            std::unique_ptr<ROOT::Math::Minimizer> min(
                ROOT::Math::Factory::CreateMinimizer("Minuit2","Simplex") // Nelder-Mead
            );
            
            min->SetFunction(functor);
            
            // parameters: pWx, pWy, pWz, nu_pz
            min->SetVariable(0,"pWx",0.0,10.0);
            min->SetVariable(1,"pWy",0.0,10.0);
            min->SetVariable(2,"pWz",0.0,10.0);
            min->SetVariable(3,"nu_pz",0.0,10.0);
            
            min->SetMaxFunctionCalls(10);
            min->SetMaxIterations(10);
            
            min->Minimize();

            double S = min->MinValue();

            if (std::isfinite(S) && S < bestS)
                bestS = S;
        }
    }

    if (!std::isfinite(bestS))
        return NAN;

    float topness = std::log(bestS);

    return topness;
}

#endif
