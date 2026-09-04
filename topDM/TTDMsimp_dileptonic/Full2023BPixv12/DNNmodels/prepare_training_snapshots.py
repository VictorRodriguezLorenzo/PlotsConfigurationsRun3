import ROOT
import os
import sys
import re

ROOT.EnableImplicitMT()

# ============================================================
# INPUT ARGUMENT
# ============================================================

PARAMETRIC = True

input_file = sys.argv[1]

outdir = "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2023BPixv12/DNNmodels/files_for_training"
os.makedirs(outdir, exist_ok=True)

outfile = os.path.join(
    outdir,
    os.path.basename(input_file).replace(".root", "_snapshot.root")
)

# ============================================================
# VARIABLE LIST
# ============================================================

OLD_COMMENTED_VARS = [
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
    'nbjets',
    'ht',
    'vht_pt',
    'dphijet1met',
    'dphijet2met',
    'dphijjmet',

    'chel',
    'pdark',
    'dphi_ttbar',
    'dphi_met_llb',
]

TOPDM_EXTRA_VARS = [
    'mt2blbl',

    'dphi_met_ll',
    'st',
    'met_over_sqrt_ht',
    'met_over_st',
    'dphi_min_j_met',

    'pt_llb',
    'dphi_met_llb_safe',

    'pt_llbb',
    'dphi_met_llbb',
    'm_llbb',
    'mT_llbb_met',
    'met_over_pt_llbb',

    'mt2_bell_l',
    'mbl_min',
    'mbl_max',
    'mtbl',

    'mbb',
    'drbb',
    'ptbb',

    'max_nonleading_btag',
    'pt_b2',

    'nForwardJet',
    'leadingForwardJet_pt',
    'leadingForwardJet_absEta',
    'deta_forwardJet_b',
    'dphi_forwardJet_met',

    'top1_pt_reco',
    'top2_pt_reco',
    'tt_reco',
    'pdark_over_met',
]

RESTFRAME_VARS = [
    'angle_ll_llbb_rf',
    'dphi_ll_llbb_rf',
    'cos_l1_llbb_rf',
    'cos_l2_llbb_rf',

    'angle_ll_llmet_rf',
    'dphi_ll_llmet_rf',
    'cos_l1_llmet_rf',
    'cos_l2_llmet_rf',
]

# Unique list, preserving order
var = list(dict.fromkeys(OLD_COMMENTED_VARS + TOPDM_EXTRA_VARS + RESTFRAME_VARS))

if PARAMETRIC:
    var.append("mPhi")

# ============================================================
# DECLARATIONS
# ============================================================

ROOT.gROOT.ProcessLine(
"""
template <typename container>
float Alt(container c, int index, float alt) {
    if (index >= 0 && index < (int)c.size()) return c[index];
    else return alt;
}
"""
)

ROOT.gInterpreter.Declare(
"""
#include <cmath>
#include <algorithm>
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"

using namespace ROOT;
using namespace ROOT::VecOps;

#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/computeMT2.cc"
#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/computeMT2blbl.cc"
#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/topDM_vars.cc"
#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/topDM_restFrame_vars.cc"

#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/jet_horns.cc"
#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/doubleNu_producer.h"

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

        if (nLep < 2) return emptyResult();

        auto leptonMass = [](int pdgId) {
            const int absId = std::abs(pdgId);
            if (absId == 11) return 0.000511f;
            if (absId == 13) return 0.105658f;
            return 0.0f;
        };

        TLorentzVector l1, l2;
        l1.SetPtEtaPhiM(Lep_pt[0], Lep_eta[0], Lep_phi[0], leptonMass(Lep_pdgId[0]));
        l2.SetPtEtaPhiM(Lep_pt[1], Lep_eta[1], Lep_phi[1], leptonMass(Lep_pdgId[1]));

        if (nCleanJet < 2) return emptyResult();

        std::vector<int> bjet_indices;

        for (size_t i = 0; i < CleanJet_pt.size(); ++i) {
            int jetIdx = CleanJet_jetIdx[i];

            if (jetIdx < 0 || jetIdx >= (int)Jet_btagger.size()) continue;

            float pt   = CleanJet_pt[i];
            float eta  = CleanJet_eta[i];
            float btag = Jet_btagger[jetIdx];

            if (pt > 30.0 && std::abs(eta) < 2.5 && btag > bAlgo_WP)
                bjet_indices.push_back(i);
        }

        if (bjet_indices.size() < 2) return emptyResult();

        int b1_idx = bjet_indices[0];
        int b2_idx = bjet_indices[1];

        TLorentzVector bj1, bj2;
        bj1.SetPtEtaPhiM(CleanJet_pt[b1_idx], CleanJet_eta[b1_idx], CleanJet_phi[b1_idx], CleanJet_mass[b1_idx]);
        bj2.SetPtEtaPhiM(CleanJet_pt[b2_idx], CleanJet_eta[b2_idx], CleanJet_phi[b2_idx], CleanJet_mass[b2_idx]);

        double met_x = PuppiMET_pt * std::cos(PuppiMET_phi);
        double met_y = PuppiMET_pt * std::sin(PuppiMET_phi);

        nuana::doubleNeutrinoSolution solver(bj1, bj2, l1, l2, met_x, met_y);

        size_t idx = 0;

        auto kin = nuana::computeEventKinematics(
            bj1, bj2, l1, l2, met_x, met_y, solver, idx
        );

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

mphi_match = re.search(r"m[Pp]hi[_-]?(\d+)", os.path.basename(input_file))
mphi_value = float(mphi_match.group(1)) if mphi_match else -999.0

# ============================================================
# B-TAGGING CONFIG
# ============================================================

btagging_WPs = {
    "DeepFlavB" : {
        "loose"    : "0.0583",
        "medium"   : "0.3086",
        "tight"    : "0.7183",
        "xtight"   : "0.8111",
        "xxtight"  : "0.9512",
    },
    "RobustParTAK4B" : {
        "loose"    : "0.0849",
        "medium"   : "0.4319",
        "tight"    : "0.8482",
        "xtight"   : "0.9151",
        "xxtight"  : "0.9874",
    },
    "PNetB" : {
        "loose"    : "0.047",
        "medium"   : "0.245",
        "tight"    : "0.6734",
        "xtight"   : "0.7862",
        "xxtight"  : "0.961",
    }
}

# Algo / SF name
btagging_SFs = {
    "DeepFlavB"      : "deepjet",
    "RobustParTAK4B" : "RobustParT",
    "PNetB"          : "partNet",
}

bAlgo = 'PNetB'
bWP   = 'medium'

# ============================================================
# SMALL HELPERS FOR DEFINING COLUMNS
# ============================================================

def _columns(df_):
    return set(str(c) for c in df_.GetColumnNames())

def has_column(df_, name):
    return name in _columns(df_)

def define_if_missing(df_, name, expr):
    if has_column(df_, name):
        return df_
    return df_.Define(name, expr)

# ============================================================
# DEFINE VARIABLES
# ============================================================

df = df.Define(
    "doubleNu",
    "doubleNu_producer(nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, "
    "nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, "
    "PuppiMET_pt, PuppiMET_phi, Jet_btag{0}, {1})".format(
        bAlgo,
        btagging_WPs[bAlgo][bWP]
    )
)

# Basic lepton variables
df = define_if_missing(df, "lep_eta1", "Alt(Lepton_eta, 0, -999.f)")
df = define_if_missing(df, "lep_eta2", "Alt(Lepton_eta, 1, -999.f)")
df = define_if_missing(df, "lep_pt1",  "Alt(Lepton_pt,  0, -999.f)")
df = define_if_missing(df, "lep_pt2",  "Alt(Lepton_pt,  1, -999.f)")

df = define_if_missing(
    df,
    "mll",
    "sqrt(std::max(0.0, 2.0 * Lepton_pt[0] * Lepton_pt[1] * "
    "(cosh(Lepton_eta[0] - Lepton_eta[1]) - cos(DeltaPhi(Lepton_phi[0], Lepton_phi[1])))))"
)

df = define_if_missing(df, "dphill", "abs(DeltaPhi(Lepton_phi[0], Lepton_phi[1]))")
df = define_if_missing(df, "detall", "Lepton_eta[0] - Lepton_eta[1]")
df = define_if_missing(
    df,
    "drll",
    "sqrt(pow(Lepton_eta[0]-Lepton_eta[1], 2) + pow(DeltaPhi(Lepton_phi[0], Lepton_phi[1]), 2))"
)

df = define_if_missing(
    df,
    "ptll",
    "sqrt(pow(Lepton_pt[0]*cos(Lepton_phi[0]) + Lepton_pt[1]*cos(Lepton_phi[1]), 2) + "
    "pow(Lepton_pt[0]*sin(Lepton_phi[0]) + Lepton_pt[1]*sin(Lepton_phi[1]), 2))"
)

yll_expr = (
    "((Lepton_pt[0]*cosh(Lepton_eta[0]) + Lepton_pt[1]*cosh(Lepton_eta[1]) - "
    "  Lepton_pt[0]*sinh(Lepton_eta[0]) - Lepton_pt[1]*sinh(Lepton_eta[1])) > 0 && "
    " (Lepton_pt[0]*cosh(Lepton_eta[0]) + Lepton_pt[1]*cosh(Lepton_eta[1]) + "
    "  Lepton_pt[0]*sinh(Lepton_eta[0]) + Lepton_pt[1]*sinh(Lepton_eta[1])) > 0) ? "
    "0.5 * log( "
    "(Lepton_pt[0]*cosh(Lepton_eta[0]) + Lepton_pt[1]*cosh(Lepton_eta[1]) + "
    " Lepton_pt[0]*sinh(Lepton_eta[0]) + Lepton_pt[1]*sinh(Lepton_eta[1])) / "
    "(Lepton_pt[0]*cosh(Lepton_eta[0]) + Lepton_pt[1]*cosh(Lepton_eta[1]) - "
    " Lepton_pt[0]*sinh(Lepton_eta[0]) - Lepton_pt[1]*sinh(Lepton_eta[1])) ) : -999.f"
)
df = define_if_missing(df, "yll", yll_expr)

# Jets, b jets, and ratios
df = define_if_missing(df, "njets", "Sum(CleanJet_pt > 20.)")
df = define_if_missing(df, "njet", "njets")

df = define_if_missing(
    df,
    "nbjets",
    "Sum(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && "
    "Take(Jet_btag{0}, CleanJet_jetIdx) > {1})".format(
        bAlgo,
        btagging_WPs[bAlgo][bWP]
    )
)

df = define_if_missing(
    df,
    "nbjet_jet_ratio",
    "(njets > 0) ? float(nbjets)/float(njets) : 0.f"
)

df = define_if_missing(
    df,
    "bReq",
    "Sum(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && "
    "Take(Jet_btag{0}, CleanJet_jetIdx) > {1}) >= 1".format(
        bAlgo,
        btagging_WPs[bAlgo][bWP]
    )
)

# MET and angular variables
df = define_if_missing(df, "dphilmet1", "abs(DeltaPhi(Lepton_phi[0], PuppiMET_phi))")
df = define_if_missing(df, "dphilmet2", "abs(DeltaPhi(Lepton_phi[1], PuppiMET_phi))")
df = define_if_missing(df, "dphilmet", "std::min(dphilmet1, dphilmet2)")

df = define_if_missing(
    df,
    "dphillmet",
    "abs(DeltaPhi(PuppiMET_phi, atan2(Lepton_pt[0]*sin(Lepton_phi[0]) + Lepton_pt[1]*sin(Lepton_phi[1]), "
    "Lepton_pt[0]*cos(Lepton_phi[0]) + Lepton_pt[1]*cos(Lepton_phi[1]))))"
)

df = define_if_missing(
    df,
    "dphijet1met",
    "(nCleanJet > 0) ? abs(DeltaPhi(CleanJet_phi[0], PuppiMET_phi)) : -999.f"
)

df = define_if_missing(
    df,
    "dphijet2met",
    "(nCleanJet > 1) ? abs(DeltaPhi(CleanJet_phi[1], PuppiMET_phi)) : -999.f"
)

df = define_if_missing(
    df,
    "dphijjmet",
    "(nCleanJet > 1) ? abs(DeltaPhi(PuppiMET_phi, atan2(CleanJet_pt[0]*sin(CleanJet_phi[0]) + CleanJet_pt[1]*sin(CleanJet_phi[1]), "
    "CleanJet_pt[0]*cos(CleanJet_phi[0]) + CleanJet_pt[1]*cos(CleanJet_phi[1])))) : -999.f"
)

# Transverse masses and old HWW-style variables
df = define_if_missing(
    df,
    "mtw1",
    "sqrt(2.0 * Lepton_pt[0] * PuppiMET_pt * (1.0 - cos(DeltaPhi(Lepton_phi[0], PuppiMET_phi))))"
)

df = define_if_missing(
    df,
    "mtw2",
    "sqrt(2.0 * Lepton_pt[1] * PuppiMET_pt * (1.0 - cos(DeltaPhi(Lepton_phi[1], PuppiMET_phi))))"
)

df = define_if_missing(
    df,
    "mth",
    "sqrt(2.0 * ptll * PuppiMET_pt * (1.0 - cos(dphillmet)))"
)

df = define_if_missing(df, "mTi", "-999.f")
df = define_if_missing(df, "mR", "-999.f")
df = define_if_missing(df, "mTe", "-999.f")

# MT2 and ttbar reconstruction variables
df = define_if_missing(
    df,
    "mT2",
    "computeMT2(Lepton_pt[0], Lepton_eta[0], Lepton_phi[0], "
    "Lepton_pt[1], Lepton_eta[1], Lepton_phi[1], "
    "PuppiMET_pt, PuppiMET_phi)"
)

df = define_if_missing(
    df,
    "mt2blbl",
    "computeMT2blbl(nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, "
    "nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, "
    "PuppiMET_pt, PuppiMET_phi, Jet_btag{0}, {1})".format(
        bAlgo,
        btagging_WPs[bAlgo][bWP]
    )
)

df = define_if_missing(df, "pdark", "doubleNu[8]")
df = define_if_missing(df, "chel", "doubleNu[6]")
df = define_if_missing(df, "dphi_ttbar", "doubleNu[7]")
df = define_if_missing(df, "tt_reco", "doubleNu[9]")
df = define_if_missing(df, "top1_pt_reco", "doubleNu[4]")
df = define_if_missing(df, "top2_pt_reco", "doubleNu[5]")

df = define_if_missing(
    df,
    "pdark_over_met",
    "(PuppiMET_pt > 0 && pdark > -998) ? pdark/PuppiMET_pt : -999.f"
)

# Old dphi_met_llb, kept because you asked for ALL variables
df = define_if_missing(
    df,
    "bjet_idx",
    "(nCleanJet > 0 && CleanJet_jetIdx[0] >= 0 && CleanJet_jetIdx[0] < (int)Jet_btag{algo}.size() && Jet_btag{algo}[CleanJet_jetIdx[0]] > {wp}) ? 0 : "
    "(nCleanJet > 1 && CleanJet_jetIdx[1] >= 0 && CleanJet_jetIdx[1] < (int)Jet_btag{algo}.size() && Jet_btag{algo}[CleanJet_jetIdx[1]] > {wp}) ? 1 : "
    "(nCleanJet > 2 && CleanJet_jetIdx[2] >= 0 && CleanJet_jetIdx[2] < (int)Jet_btag{algo}.size() && Jet_btag{algo}[CleanJet_jetIdx[2]] > {wp}) ? 2 : -1".format(
        algo=bAlgo,
        wp=btagging_WPs[bAlgo][bWP]
    )
)

df = define_if_missing(
    df,
    "dphi_met_llb",
    "(bjet_idx >= 0) ? abs(DeltaPhi(PuppiMET_phi, atan2("
    "Lepton_pt[0]*sin(Lepton_phi[0]) + Lepton_pt[1]*sin(Lepton_phi[1]) + CleanJet_pt[bjet_idx]*sin(CleanJet_phi[bjet_idx]), "
    "Lepton_pt[0]*cos(Lepton_phi[0]) + Lepton_pt[1]*cos(Lepton_phi[1]) + CleanJet_pt[bjet_idx]*cos(CleanJet_phi[bjet_idx])))) : -999.f"
)

# New topDM variables
df = df.Define(
    "topDMVars",
    "topDM_DNN_vars(nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, "
    "nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, "
    "PuppiMET_pt, PuppiMET_phi, Jet_btag{0}, {1})".format(
        bAlgo,
        btagging_WPs[bAlgo][bWP]
    )
)

_topdm_defines = {
    'dphi_met_ll'              : 0,
    'ht'                       : 1,
    'st'                       : 2,
    'met_over_sqrt_ht'         : 3,
    'met_over_st'              : 4,
    'dphi_min_j_met'           : 5,

    'pt_llb'                   : 6,
    'dphi_met_llb_safe'        : 7,

    'pt_llbb'                  : 8,
    'dphi_met_llbb'            : 9,
    'm_llbb'                   : 10,
    'mT_llbb_met'              : 11,
    'met_over_pt_llbb'         : 12,

    'mt2_bell_l'               : 13,
    'mbl_min'                  : 14,
    'mbl_max'                  : 15,
    'mtbl'                     : 16,

    'mbb'                      : 17,
    'drbb'                     : 18,
    'ptbb'                     : 19,

    'max_nonleading_btag'      : 20,
    'pt_b2'                    : 21,

    'nForwardJet'              : 22,
    'leadingForwardJet_pt'     : 23,
    'leadingForwardJet_absEta' : 24,
    'deta_forwardJet_b'        : 25,
    'dphi_forwardJet_met'      : 26,
}

for _name, _idx in _topdm_defines.items():
    df = define_if_missing(df, _name, "topDMVars[{}]".format(_idx))

# Rest-frame variables
df = df.Define(
    "topDMRestFrameVars",
    "topDM_restFrame_vars(nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, "
    "nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, "
    "PuppiMET_pt, PuppiMET_phi, Jet_btag{0}, {1})".format(
        bAlgo,
        btagging_WPs[bAlgo][bWP]
    )
)

_restframe_defines = {
    'angle_ll_llbb_rf' : 0,
    'dphi_ll_llbb_rf'  : 1,
    'cos_l1_llbb_rf'   : 2,
    'cos_l2_llbb_rf'   : 3,

    'angle_ll_llmet_rf': 4,
    'dphi_ll_llmet_rf' : 5,
    'cos_l1_llmet_rf'  : 6,
    'cos_l2_llmet_rf'  : 7,
}

for _name, _idx in _restframe_defines.items():
    df = define_if_missing(df, _name, "topDMRestFrameVars[{}]".format(_idx))

# Recoil / HWW-style variables.
# If they already exist in NanoLatino, they are kept.
# If not, these are simple fallbacks so the snapshot does not crash.
df = define_if_missing(
    df,
    "recoil",
    "sqrt(pow(ptll*cos(atan2(Lepton_pt[0]*sin(Lepton_phi[0]) + Lepton_pt[1]*sin(Lepton_phi[1]), "
    "Lepton_pt[0]*cos(Lepton_phi[0]) + Lepton_pt[1]*cos(Lepton_phi[1]))) + PuppiMET_pt*cos(PuppiMET_phi), 2) + "
    "pow(ptll*sin(atan2(Lepton_pt[0]*sin(Lepton_phi[0]) + Lepton_pt[1]*sin(Lepton_phi[1]), "
    "Lepton_pt[0]*cos(Lepton_phi[0]) + Lepton_pt[1]*cos(Lepton_phi[1]))) + PuppiMET_pt*sin(PuppiMET_phi), 2))"
)

df = define_if_missing(df, "pTWW", "recoil")

df = define_if_missing(df, "upara", "-PuppiMET_pt * cos(dphillmet)")
df = define_if_missing(df, "uperp", "PuppiMET_pt * sin(dphillmet)")

df = define_if_missing(
    df,
    "vht_pt",
    "sqrt(pow(Sum(CleanJet_pt * cos(CleanJet_phi)) + Lepton_pt[0]*cos(Lepton_phi[0]) + Lepton_pt[1]*cos(Lepton_phi[1]) + PuppiMET_pt*cos(PuppiMET_phi), 2) + "
    "pow(Sum(CleanJet_pt * sin(CleanJet_phi)) + Lepton_pt[0]*sin(Lepton_phi[0]) + Lepton_pt[1]*sin(Lepton_phi[1]) + PuppiMET_pt*sin(PuppiMET_phi), 2))"
)

for _missing_var in [
    "mcoll",
    "mcollWW",
    "choiMass",
]:
    df = define_if_missing(df, _missing_var, "-999.f")

# Event cleaning
df = define_if_missing(
    df,
    "noJetInHorn_pT20",
    "Jet_inHorns(CleanJet_pt, CleanJet_eta, true)"
)

if PARAMETRIC:
    df = df.Define("mPhi", str(mphi_value))

# ============================================================
# FILTER
# ============================================================

df = df.Filter(
    "((abs(Lepton_pdgId[0]) == 11 || abs(Lepton_pdgId[0]) == 13) && "
    "(abs(Lepton_pdgId[1]) == 11 || abs(Lepton_pdgId[1]) == 13)) && "
    "Lepton_pt[0] > 25 && Lepton_pt[1] > 20 && Alt(Lepton_pt, 2, 0) < 10 && "
    "((abs(Lepton_pdgId[0]) == 11 && abs(Lepton_eta[0]) < 2.5) || "
    " (abs(Lepton_pdgId[0]) == 13 && abs(Lepton_eta[0]) < 2.4)) && "
    "((abs(Lepton_pdgId[1]) == 11 && abs(Lepton_eta[1]) < 2.5) || "
    " (abs(Lepton_pdgId[1]) == 13 && abs(Lepton_eta[1]) < 2.4)) && "
    "mll > 20 && noJetInHorn_pT20 && bReq && "
    "mT2 > 80 && "
    "(!(abs(Lepton_pdgId[0]) == abs(Lepton_pdgId[1])) || abs(91.1876 - mll) > 15)"
)

# ============================================================
# SNAPSHOT OUTPUT
# ============================================================

print("Writing variables:")
for v in var:
    print("  ", v)

df.Snapshot("Events", outfile, var)

print("Wrote:", outfile)
