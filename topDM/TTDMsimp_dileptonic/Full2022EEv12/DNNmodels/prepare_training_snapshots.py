import ROOT
import os
import sys

ROOT.EnableImplicitMT()

# ============================================================
# INPUT ARGUMENT
# ============================================================
# Condor passes the file path
input_file = sys.argv[1]

outdir = "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022EEv12/DNNmodels/files_for_training"
os.makedirs(outdir, exist_ok=True)

outfile = os.path.join(
    outdir,
    os.path.basename(input_file).replace(".root", "_snapshot.root")
)

# ============================================================
# VARIABLE LIST
# ============================================================
var = [
'lep_pt1',
'lep_pt2',
'lep_eta1',
'lep_eta2',

'mll',
'ptll',
'drll',
'detall',
'dphill',
'yll',

'PuppiMET_pt',
'PuppiMET_phi',
'dphilmet',
'dphilmet1',
'dphilmet2',
'dphillmet',

'mtw1',
'mtw2',
'mth',
'mTi',
'mR',
'mT2',
'mTe',

'recoil',
'upara',
'uperp',
'pTWW',

'mcoll',
'mcollWW',
'choiMass',

'nbjet_jet_ratio',
'njet',
'ht',
'vht_pt',
'dphijet1met',
'dphijet2met',
'dphijjmet',

'chel',
'pdark',
'dphi_ttbar',
'dphi_met_llb'
]

#var = [
#  'dphill',
#  'PuppiMET_pt',
#  'mT2',
#  'pdark',
#  'chel',
#  'dphi_ttbar',
#  'dphi_met_llb',
#]

# ============================================================
# DECLARATIONS
# ============================================================

ROOT.gROOT.ProcessLine(
"""
template
<typename container>
float Alt(container c, int index, float alt){
    if (index < c.size()) return c[index];
    else return alt;
}
"""
)

ROOT.gInterpreter.Declare(
"""
#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/computeMT2.cc"
#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/jet_horns.cc"
#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/doubleNu_producer.h"

#include <cmath>
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"

using namespace ROOT;
using namespace ROOT::VecOps;
using namespace nuana;

// ---- helper ----
RVecF emptyResult() {
    return RVecF{-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, 0};
}

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
        if (nCleanJet < 2) return emptyResult();

        std::vector<int> bjet_indices;

        for (size_t i = 0; i < CleanJet_pt.size(); ++i) {
            int jetIdx = CleanJet_jetIdx[i];

            // Guard against corrupt indices
            if (jetIdx < 0 || jetIdx >= (int)Jet_btagger.size()) continue;

            float pt  = CleanJet_pt[i];
            float eta = CleanJet_eta[i];
            float btag = Jet_btagger[jetIdx];

            if (pt > 20.0 && std::abs(eta) < 2.5 && btag > bAlgo_WP)
                bjet_indices.push_back(i);
        }

        if (bjet_indices.size() < 2) return emptyResult();

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
        out[9] = kin.valid ? 1 : 0;
//      std::cout << "out[9] = " << out[9] << std::endl;
        if (kin.valid) {
            return out;
        } else {
            return RVecF{-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, 0};
        }
}


"""

)

# ============================================================
# BUILD DATAFRAME FROM SINGLE FILE
# ============================================================

df = ROOT.RDataFrame("Events", input_file)

# Algo / WP / WP cut
btagging_WPs = {
    "DeepFlavB" : {
        "loose"    : "0.0614",
        "medium"   : "0.3196",
        "tight"    : "0.73",
        "xtight"   : "0.8184",
        "xxtight"  : "0.9542",
    },
    "RobustParTAK4B" : {
        "loose"    : "0.0897",
        "medium"   : "0.451",
        "tight"    : "0.8604",
        "xtight"   : "0.9234",
        "xxtight"  : "0.9893",
    },
    "PNetB" : {
        "loose"    : "0.0499",
        "medium"   : "0.2605",
        "tight"    : "0.6915",
        "xtight"   : "0.8033",
        "xxtight"  : "0.9664",
    }
}

# Algo / SF name
btagging_SFs = {
    "DeepFlavB"      : "deepjet",
    "RobustParTAK4B" : "RobustParT",
    "PNetB"          : "partNet",
}

# Algorithm and WP selection
bAlgo = 'RobustParTAK4B' # ['DeepFlavB','RobustParTAK4B','PNetB'] 
bWP    = 'medium'     # ['loose','medium','tight','xtight','xxtight']

df = df.Define(
"doubleNu",
"doubleNu_producer(nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, PuppiMET_pt, PuppiMET_phi, Jet_btag{0}, {1})".format(bAlgo, btagging_WPs[bAlgo][bWP])
)

df = df.Define("lep_eta1", "Lepton_eta[0]")
df = df.Define("lep_eta2", "Lepton_eta[1]")
df = df.Define("lep_pt1", "Lepton_pt[0]")
df = df.Define("lep_pt2", "Lepton_pt[1]")
df = df.Define("nbjet_jet_ratio", 'Sum(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{}, CleanJet_jetIdx) > {})/njet'.format(bAlgo, btagging_WPs[bAlgo][bWP]))
df = df.Define("mT2",'computeMT2(Lepton_pt[0],Lepton_eta[0],Lepton_phi[0],Lepton_pt[1],Lepton_eta[1],Lepton_phi[1],PuppiMET_pt,PuppiMET_phi)')
df = df.Define("pdark","doubleNu[8]")
df = df.Define("chel","doubleNu[6]")
df = df.Define("dphi_ttbar","doubleNu[7]")
df = df.Define("tt_reco","doubleNu[9]")
df = df.Define("bjet_idx",'(nCleanJet > 0 && Jet_btag{algo}[CleanJet_jetIdx[0]] > {wp}) ? 0 : (nCleanJet > 1 && Jet_btag{algo}[CleanJet_jetIdx[1]] > {wp}) ? 1 : (nCleanJet > 2 && Jet_btag{algo}[CleanJet_jetIdx[2]] > {wp}) ? 2 : -1'.format(algo=bAlgo, wp=btagging_WPs[bAlgo][bWP]))
df = df.Define("dphi_met_llb","abs(DeltaPhi(PuppiMET_phi, atan2(Lepton_pt[0]*sin(Lepton_phi[0]) + Lepton_pt[1]*sin(Lepton_phi[1]) + CleanJet_pt[bjet_idx]*sin(CleanJet_phi[bjet_idx]), Lepton_pt[0]*cos(Lepton_phi[0]) + Lepton_pt[1]*cos(Lepton_phi[1]) + CleanJet_pt[bjet_idx]*cos(CleanJet_phi[bjet_idx]))))")
df = df.Define("noJetInHorn","Jet_inHorns(CleanJet_pt, CleanJet_eta)")
df = df.Define("bReq",'Sum(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{}, CleanJet_jetIdx) > {}) >= 1'.format(bAlgo, btagging_WPs[bAlgo][bWP]))

# ============================================================
# FILTER
# ============================================================

df = df.Filter(
"((abs(Lepton_pdgId[0]) == 11 || abs(Lepton_pdgId[0]) == 13) && "
"(abs(Lepton_pdgId[1]) == 11 || abs(Lepton_pdgId[1]) == 13)) && "
"Lepton_pt[0]>25 && Lepton_pt[1]>20 && Alt(Lepton_pt,2, 0)<10 && "
"abs(Lepton_eta[0]) < 2.4 && abs(Lepton_eta[1]) < 2.4 && "
"mll > 20 && noJetInHorn && bReq && tt_reco && "
"mT2 > 80 && (!(abs(Lepton_pdgId[0]) == abs(Lepton_pdgId[1])) || abs(91.1876 - mll) > 15)"
)

# ============================================================
# SNAPSHOT OUTPUT
# ============================================================

df.Snapshot("Events", outfile, var)

print("Wrote:", outfile)
