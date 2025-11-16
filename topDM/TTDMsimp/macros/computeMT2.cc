#ifndef COMPUTE_MT2_CC
#define COMPUTE_MT2_CC

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TMinuit.h"
#include "TLorentzVector.h"
#include "TVector2.h"

struct MT2Input {
    double px1, py1;
    double px2, py2;
    double metx, mety;
};

static MT2Input gMT2;

void functionMT2(int&, double*, double& f, double par[], int) {

    double met1  = par[0];
    double phi1  = par[1];

    double metx1 = met1 * std::cos(phi1);
    double mety1 = met1 * std::sin(phi1);

    double metx2 = gMT2.metx - metx1;
    double mety2 = gMT2.mety - mety1;

    double met2  = std::sqrt(metx2 * metx2 + mety2 * mety2);

    double p1 = std::sqrt(gMT2.px1 * gMT2.px1 + gMT2.py1 * gMT2.py1);
    double p2 = std::sqrt(gMT2.px2 * gMT2.px2 + gMT2.py2 * gMT2.py2);

    // Δφ
    double dphi1 = TVector2::Phi_mpi_pi(std::atan2(gMT2.py1, gMT2.px1) - phi1);
    double phi2  = std::atan2(mety2, metx2); 
    double dphi2 = TVector2::Phi_mpi_pi(std::atan2(gMT2.py2, gMT2.px2) - phi2);

    // MT formulas
    double mt1 = 2.0 * p1 * met1 * (1.0 - std::cos(dphi1));
    double mt2 = 2.0 * p2 * met2 * (1.0 - std::cos(dphi2));

    double maxmt2 = mt1 > mt2 ? mt1 : mt2;

    f = std::sqrt(maxmt2);  
}

// ---------------------------------------------------------
// Compute MT2 from four vectors + MET
// ---------------------------------------------------------
float computeMT2(
    float l1_pt, float l1_eta, float l1_phi,
    float l2_pt, float l2_eta, float l2_phi,
    float met_pt, float met_phi
) {

    TLorentzVector L1, L2;
    L1.SetPtEtaPhiM(l1_pt, l1_eta, l1_phi, 0.);
    L2.SetPtEtaPhiM(l2_pt, l2_eta, l2_phi, 0.);

    TVector2 MET;
    MET.SetMagPhi(met_pt, met_phi);

    // Fill the structure for Minuit
    gMT2.px1  = L1.Px();
    gMT2.py1  = L1.Py();
    gMT2.px2  = L2.Px();
    gMT2.py2  = L2.Py();
    gMT2.metx = MET.X();
    gMT2.mety = MET.Y();

    // Setup Minuit
    TMinuit minuit(2);
    minuit.SetFCN(functionMT2);
    minuit.SetPrintLevel(-1);

    const double met = MET.Mod();
    double par[2]      = {met * 0.5, 0.0};
    double step[2]     = {met / 100., 0.01};
    double minVal[2]   = {0.0, -M_PI};
    double maxVal[2]   = {2.0 * met, M_PI};
    const char *names[2] = {"met1", "phi1"};

    for (int i = 0; i < 2; i++) {
        minuit.DefineParameter(i, names[i], par[i], step[i], minVal[i], maxVal[i]);
    }

    minuit.Migrad();

    // Minimized value
    double mt2 = minuit.fAmin;

    return mt2;
}

#endif

