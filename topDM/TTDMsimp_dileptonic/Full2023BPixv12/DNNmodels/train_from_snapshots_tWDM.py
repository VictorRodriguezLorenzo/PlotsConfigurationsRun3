import ROOT
import glob
import pandas as pd
import numpy as np
import os
import re
import matplotlib.pyplot as plt
import joblib

import tensorflow as tf

from tensorflow.keras import Input
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization
from tensorflow.keras import callbacks, regularizers
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.metrics import roc_auc_score

# ============================================================
# CONTROL PARAMETERS
# ============================================================

# Run both mediator hypotheses automatically when the script is launched
# normally. Each child process executes the same training pipeline with a
# different signal sample and keeps all outputs separated by hypothesis.
SIGNAL_TYPES = ("ps", "s")
SIGNAL_TYPE = os.environ.get("TWDM_SIGNAL_TYPE")

if SIGNAL_TYPE is None:
    import subprocess
    import sys

    for signal_type in SIGNAL_TYPES:
        print("\n" + "=" * 80)
        print(f"Starting tWDM training for signal type: {signal_type}")
        print("=" * 80)

        env = os.environ.copy()
        env["TWDM_SIGNAL_TYPE"] = signal_type
        subprocess.run([sys.executable, os.path.abspath(__file__)], env=env, check=True)

    print("\nFinished tWDM trainings for: ps, s")
    sys.exit(0)

if SIGNAL_TYPE not in SIGNAL_TYPES:
    raise ValueError(
        f"Unknown TWDM_SIGNAL_TYPE={SIGNAL_TYPE!r}. "
        f"Expected one of {SIGNAL_TYPES}."
    )

ANALYSIS_NAME = f"model_DNN_tWDM_{SIGNAL_TYPE}"
save_model = True

loaded_model = False

# Keep the trained .keras files next to this training script, with the exact
# names used by the scalar/pseudoscalar inference modules.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_NAME = os.path.join(SCRIPT_DIR, f"./Models/model_DNN_tWDM_{SIGNAL_TYPE}.keras")

PARAMETRIC = True

SNAPSHOT_DIR = "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2023BPixv12/DNNmodels/files_for_training"

# ============================================================
# VARIABLES
# ============================================================

#var = [
#    'vht_pt',
#    'mtw2',
#    'mtw1',
#    'mt2blbl',
#    'pt_llbb',
#    'ptbb',
#    'lep_pt1',
#    'mth',
#    'dphi_min_j_met',
#    'PuppiMET_pt',
#    'pdark',
#    'm_llbb',
#    'dphilmet',
#    'pt_b2',
#    'ptll',
#    'chel',
#    'met_over_st',
#    'pdark_over_met',
#    'angle_ll_llbb_rf',
#    'dphi_ll_llmet_rf',
#    'cos_l1_llmet_rf',
#]

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

var = list(dict.fromkeys(OLD_COMMENTED_VARS + TOPDM_EXTRA_VARS + RESTFRAME_VARS))

if PARAMETRIC:
    var.append("mPhi")

# ============================================================
# FILE DISCOVERY
# ============================================================
# ============================================================
# TRAINING / PARAMETRIC CONTROLS
# ============================================================

# Main training recommended here:
#   tWDM vs (ttbar + single top + visible ttV)
# with TTNuNu excluded and reserved for the dedicated classifier.
#
# Set to "ttnunu" to train the second classifier with exactly the same
# balancing, preprocessing, validation and plotting machinery.
TRAINING_TARGET = "main"          # "main" or "ttnunu"

# Keep "vanilla" first for a clean comparison with the current model.
# Set to "affine" later to test repeated affine conditioning with the
# same dataset and the same balanced training procedure.
MODEL_ARCHITECTURE = "affine"   # "vanilla" or "affine"

SEED = 12345
TEST_SIZE = 0.20
BATCH_SIZE = 1024
MAX_EPOCHS = 300
PREDICT_BATCH_SIZE = 4096

# ------------------------------------------------------------------
# BALANCING STRATEGY
# ------------------------------------------------------------------
# The pNN reference recommends:
#   * 50/50 signal-background class balance;
#   * equal representation of all signal mass hypotheses in each batch;
#   * explicit representation of all background processes;
#   * background mPhi sampled from the same discrete set as the signal.
#
# A strict paper-style full balance would give every background process the
# same fraction. For this analysis that would strongly over-emphasize the
# rare backgrounds and can reduce the ttDM-vs-TTTo2L2Nu performance.
# Therefore the default is an analysis-aware adaptation: all groups remain
# present in every batch, TTTo2L2Nu remains the dominant adversary, and the
# total 15% ttV share is split equally among other ttZ, ttW and ttH so that
# raw MC statistics do not make one visible-ttV process dominate the others.
#
# Available modes:
#   "analysis_weighted" : recommended default for the main DNN
#   "paper_full"        : strict equal balance over background groups
BACKGROUND_BALANCE_MODE = "analysis_weighted"

MAIN_BACKGROUND_BATCH_FRACTIONS = {
    "tt2l":        0.50,
    "other_tt":    0.08,
    "single_top":  0.12,
    "other_ttZ":   0.05,
    "ttW":         0.05,
    "ttH":         0.05,
    "TTNuNu":      0.15,
}

# Backgrounds receive one of the real signal hypotheses, sampled again
# whenever they enter a training mini-batch (identical-sampled strategy).
DYNAMIC_BACKGROUND_MPHI = True

# Define one generator epoch from the signal exposure, not from the size of
# the huge raw background pool. With 50% signal batches, this makes the total
# number of signal draws per epoch approximately equal to the number of
# available training signal events.
STEPS_PER_EPOCH_MODE = "signal_once"  # "signal_once" or "all_events"

# DNN-tail diagnostics. These are raw-event efficiencies unless a nominal
# event weight is added to the snapshots and handled explicitly.
TAIL_THRESHOLDS = [0.80, 0.90, 0.95]

# ------------------------------------------------------------------
# SUMMARY DIAGNOSTIC PLOTS
# ------------------------------------------------------------------
# These plots turn the numerical diagnostics into the physics summary needed
# to decide whether the new training preserves ttDM-vs-tt(2l) separation,
# learns visible ttV, and whether a dedicated TTNuNu classifier adds genuinely
# new information.
MAKE_SUMMARY_DIAGNOSTIC_PLOTS = True
SUMMARY_PLOT_DPI = 350
SUMMARY_MASSES = [10.0, 400.0, 800.0, 1000.0]

# Reference AUC values from the earlier dedicated ttDM-vs-TTNuNu training.
# They are used only for an optional comparison plot and do not affect training.
# Set this dictionary to {} to disable the comparison.
DEDICATED_TTNUNU_AUC_REFERENCE = {
    10.0: 0.506,
    50.0: 0.545,
    100.0: 0.564,
    150.0: 0.572,
    200.0: 0.613,
    250.0: 0.618,
    300.0: 0.615,
    350.0: 0.693,
    400.0: 0.669,
    500.0: 0.649,
    600.0: 0.751,
    700.0: 0.722,
    800.0: 0.779,
    1000.0: 0.773,
}

# Optional pruning kept from the current code.
DROP_CONSTANT_FEATURES = True
DROP_HIGH_INVALID_FEATURES = False
INVALID_THRESHOLD = 0.95

np.random.seed(SEED)
tf.keras.utils.set_random_seed(SEED)

# ============================================================
# FILE DISCOVERY
# ============================================================

files = sorted(glob.glob(os.path.join(SNAPSHOT_DIR, "*.root")))

if len(files) == 0:
    raise RuntimeError(f"No ROOT snapshots found in: {SNAPSHOT_DIR}")

# Keep sample matching explicit and physics-driven. The training groups are
# coarser than the evaluation groups on purpose: mini-batch balancing should
# not depend on how a physics process happens to be split into datasets.
TT2L_PATTERNS = [
    "TTTo2L2Nu",
]

OTHER_TT_PATTERNS = [
    "TTToSemiLeptonic",
]

SINGLE_TOP_PATTERNS = [
    "ST_t-channel_top",
    "ST_t-channel_antitop",
    "ST_s-channel_plus",
    "ST_s-channel_minus",
    "TWminusto2L2Nu",
    "TbarWplusto2L2Nu",
    "ST_tW_top",
    "ST_tW_antitop",
]

OTHER_TTZ_PATTERNS = [
    "TTLL_MLL-4to50",
    "TTLL_MLL-50",
    "TTZ-ZtoQQ",
]

TTW_PATTERNS = [
    "TTLNu",
    "TTW",
]

TTH_PATTERNS = [
    "ttH",
    "TTH",
]

TTNUNU_PATTERNS = [
    "TTNuNu",
]


def matches_any(path, patterns):
    name = os.path.basename(path).lower()
    return any(pattern.lower() in name for pattern in patterns)


def extract_sample_tag(path):
    basename = os.path.basename(path)
    match = re.search(r"nanoLatino_(.*)__part\d+_snapshot", basename)
    if match:
        return match.group(1)
    return basename.replace("nanoLatino_", "").replace("_snapshot.root", "")


def extract_mphi(sample_tag):
    match = re.search(r"m[Pp]hi[_-]?(\d+)", sample_tag)
    return float(match.group(1)) if match else None


def evaluation_group_from_file(path):
    if matches_any(path, TTNUNU_PATTERNS):
        return "TTNuNu"
    if matches_any(path, OTHER_TTZ_PATTERNS):
        return "other_ttZ"
    if matches_any(path, TTW_PATTERNS):
        return "ttW"
    if matches_any(path, TTH_PATTERNS):
        return "ttH"
    if matches_any(path, TT2L_PATTERNS):
        return "tt2l"
    if matches_any(path, OTHER_TT_PATTERNS):
        return "other_tt"
    if matches_any(path, SINGLE_TOP_PATTERNS):
        return "single_top"
    return None


SIGNAL_PATTERN = f"TWto2LDMsimpSpin0_{SIGNAL_TYPE}_"
files_sig = [f for f in files if SIGNAL_PATTERN in os.path.basename(f)]

files_ttnunu = [
    f for f in files
    if matches_any(f, TTNUNU_PATTERNS) and f not in files_sig
]
files_tt2l = [
    f for f in files
    if matches_any(f, TT2L_PATTERNS)
    and f not in files_sig
    and f not in files_ttnunu
]
files_other_tt = [
    f for f in files
    if matches_any(f, OTHER_TT_PATTERNS)
    and f not in files_sig
    and f not in files_ttnunu
]
files_single_top = [
    f for f in files
    if matches_any(f, SINGLE_TOP_PATTERNS)
    and f not in files_sig
    and f not in files_ttnunu
]
files_other_ttz = [
    f for f in files
    if matches_any(f, OTHER_TTZ_PATTERNS)
    and f not in files_sig
    and f not in files_ttnunu
]
files_ttw = [
    f for f in files
    if matches_any(f, TTW_PATTERNS)
    and f not in files_sig
    and f not in files_ttnunu
]
files_tth = [
    f for f in files
    if matches_any(f, TTH_PATTERNS)
    and f not in files_sig
    and f not in files_ttnunu
]

# A visible ttV training group: all ttV except TTNuNu.
files_visible_ttv = sorted(set(files_other_ttz + files_ttw + files_tth))

if TRAINING_TARGET == "main":
    training_background_files = {
        "tt2l": files_tt2l,
        "other_tt": files_other_tt,
        "single_top": files_single_top,
        "other_ttZ": files_other_ttz,
        "ttW": files_ttw,
        "ttH": files_tth,
    }

    if BACKGROUND_BALANCE_MODE == "analysis_weighted":
        ACTIVE_BACKGROUND_BATCH_FRACTIONS = MAIN_BACKGROUND_BATCH_FRACTIONS
    elif BACKGROUND_BALANCE_MODE == "paper_full":
        ACTIVE_BACKGROUND_BATCH_FRACTIONS = None
    else:
        raise ValueError(
            f"Unknown BACKGROUND_BALANCE_MODE={BACKGROUND_BALANCE_MODE!r}. "
            "Use 'analysis_weighted' or 'paper_full'."
        )
elif TRAINING_TARGET == "ttnunu":
    training_background_files = {
        "TTNuNu": files_ttnunu,
    }
    ACTIVE_BACKGROUND_BATCH_FRACTIONS = None
else:
    raise ValueError(
        f"Unknown TRAINING_TARGET={TRAINING_TARGET!r}. "
        "Use 'main' or 'ttnunu'."
    )

# Remove empty groups explicitly. This makes the error message below useful
# while still allowing the script to run before all optional samples exist.
training_background_files = {
    group: group_files
    for group, group_files in training_background_files.items()
    if len(group_files) > 0
}

# Keep the configured analysis-aware fractions consistent if an optional
# process group is absent from the snapshots.
if ACTIVE_BACKGROUND_BATCH_FRACTIONS is not None:
    ACTIVE_BACKGROUND_BATCH_FRACTIONS = {
        group: fraction
        for group, fraction in ACTIVE_BACKGROUND_BATCH_FRACTIONS.items()
        if group in training_background_files
    }

print("")
print("Input snapshot summary")
print(f"  Signal files:          {len(files_sig)}")
print(f"  tt(2l) files:          {len(files_tt2l)}")
print(f"  other tt files:        {len(files_other_tt)}")
print(f"  single-top files:      {len(files_single_top)}")
print(f"  other ttZ files:       {len(files_other_ttz)}")
print(f"  ttW files:             {len(files_ttw)}")
print(f"  ttH files:             {len(files_tth)}")
print(f"  visible ttV files:     {len(files_visible_ttv)}")
print(f"  TTNuNu files:          {len(files_ttnunu)}")
print(f"  training target:       {TRAINING_TARGET}")
print(f"  model architecture:    {MODEL_ARCHITECTURE}")
print(f"  background balance:    {BACKGROUND_BALANCE_MODE if TRAINING_TARGET == 'main' else 'single-background'}")

if len(files_sig) == 0:
    raise RuntimeError("No signal files found.")
if len(training_background_files) == 0:
    raise RuntimeError(
        f"No background files found for TRAINING_TARGET={TRAINING_TARGET!r}."
    )

print("\nBackground groups used for training:")
for group, group_files in training_background_files.items():
    print(f"  {group:15s}: {len(group_files)} files")

print("\nTTNuNu is excluded from the main DNN training:", TRAINING_TARGET == "main")

# ============================================================
# GROUP SIGNAL FILES BY MASS SAMPLE
# ============================================================

signal_file_groups = {}

for f in files_sig:
    tag = extract_sample_tag(f)
    signal_file_groups.setdefault(tag, []).append(f)

signal_hypotheses = sorted({
    extract_mphi(tag)
    for tag in signal_file_groups
    if extract_mphi(tag) is not None
})

if not signal_hypotheses:
    raise RuntimeError("No signal mPhi hypotheses found.")

print("\nSignal hypotheses:")
print("  ", signal_hypotheses)
print("Signal sample groups:")
for tag, group_files in sorted(signal_file_groups.items()):
    print(f"  {tag}: {len(group_files)} files")

# ============================================================
# ROOT -> PANDAS HELPERS
# ============================================================

features_to_read = [v for v in var if v != "mPhi"]

def dataframe_from_files(file_list, columns, description):
    if len(file_list) == 0:
        return pd.DataFrame(columns=columns)

    chain = ROOT.TChain("Events")
    for f in file_list:
        chain.Add(f)

    df = ROOT.RDataFrame(chain)

    # Require 1 b jet
    df = df.Filter("nbjets == 1")

    n_entries = df.Count().GetValue()
    print(
        f"Reading {description}: {len(file_list)} files, "
        f"{n_entries} entries after nbjets == 1"
    )

    if n_entries <= 0:
        return pd.DataFrame(columns=columns)

    return pd.DataFrame(df.AsNumpy(columns))

# ============================================================
# BUILD SIGNAL DATAFRAME
# ============================================================

signal_frames = []

for sample_tag, file_list in sorted(signal_file_groups.items()):
    mphi = extract_mphi(sample_tag)
    if mphi is None:
        raise RuntimeError(f"Could not infer mPhi from signal tag: {sample_tag}")

    frame = dataframe_from_files(
        file_list,
        features_to_read,
        f"signal {sample_tag}",
    )

    if len(frame) == 0:
        continue

    frame["isSignal"] = 1
    frame["isBkg"] = 0
    frame["mPhi_true"] = float(mphi)
    frame["bkg_group"] = "signal"
    frame["eval_group"] = "signal"
    frame["sample_tag"] = sample_tag
    signal_frames.append(frame)

if len(signal_frames) == 0:
    raise RuntimeError("No signal events were read from the snapshots.")

Sig = pd.concat(signal_frames, ignore_index=True)

# ============================================================
# BUILD TRAINING BACKGROUND DATAFRAME
# ============================================================

background_frames = []

for bkg_group, group_files in training_background_files.items():
    # Read each evaluation subgroup separately so we retain fine labels for
    # tail-composition studies, while training still balances the coarser group.
    files_by_eval_group = {}
    for f in group_files:
        eval_group = evaluation_group_from_file(f)
        if eval_group is None:
            eval_group = bkg_group
        files_by_eval_group.setdefault(eval_group, []).append(f)

    for eval_group, eval_files in sorted(files_by_eval_group.items()):
        frame = dataframe_from_files(
            eval_files,
            features_to_read,
            f"background {bkg_group}/{eval_group}",
        )

        if len(frame) == 0:
            continue

        frame["isSignal"] = 0
        frame["isBkg"] = 1
        frame["mPhi_true"] = np.nan
        frame["bkg_group"] = bkg_group
        frame["eval_group"] = eval_group
        frame["sample_tag"] = eval_group
        background_frames.append(frame)

if len(background_frames) == 0:
    raise RuntimeError("No background events were read from the snapshots.")

Bkg = pd.concat(background_frames, ignore_index=True)

# ============================================================
# READ TTNuNu FOR INDEPENDENT EVALUATION
# ============================================================

TTNuNu = dataframe_from_files(
    files_ttnunu,
    features_to_read,
    "TTNuNu evaluation sample",
)

if len(TTNuNu) > 0:
    TTNuNu["isSignal"] = 0
    TTNuNu["isBkg"] = 1
    TTNuNu["mPhi_true"] = np.nan
    TTNuNu["bkg_group"] = "TTNuNu"
    TTNuNu["eval_group"] = "TTNuNu"
    TTNuNu["sample_tag"] = "TTNuNu"

# ============================================================
# CLEAN DATAFRAMES / OPTIONAL FEATURE PRUNING
# ============================================================


def clean_frame(frame, name):
    if len(frame) == 0:
        return frame.copy()

    frame = frame.copy()
    frame[features_to_read] = frame[features_to_read].replace(
        [np.inf, -np.inf],
        np.nan,
    )

    print(f"\nNaNs before cleaning: {name}")
    print(frame[features_to_read].isna().sum())

    frame.dropna(subset=features_to_read, inplace=True)
    frame.reset_index(drop=True, inplace=True)

    print(f"Entries after cleaning {name}: {len(frame)}")
    return frame


Sig = clean_frame(Sig, "signal")
Bkg = clean_frame(Bkg, "training background")
TTNuNu = clean_frame(TTNuNu, "TTNuNu")

df_all = pd.concat([Sig, Bkg], ignore_index=True)

if DROP_CONSTANT_FEATURES:
    constant_vars = [
        v for v in features_to_read
        if df_all[v].nunique(dropna=True) <= 1
    ]

    if constant_vars:
        print("\nDropping constant variables:")
        for v in constant_vars:
            print("  ", v)

        features_to_read = [v for v in features_to_read if v not in constant_vars]

if DROP_HIGH_INVALID_FEATURES:
    high_invalid_vars = []

    for v in features_to_read:
        invalid_fraction = ((df_all[v] <= -998) | df_all[v].isna()).mean()
        if invalid_fraction > INVALID_THRESHOLD:
            high_invalid_vars.append(v)

    if high_invalid_vars:
        print("\nDropping high-invalid-fraction variables:")
        for v in high_invalid_vars:
            print("  ", v)

        features_to_read = [v for v in features_to_read if v not in high_invalid_vars]

# Rebuild the final model input list after pruning.
var = list(features_to_read)
if PARAMETRIC:
    var.append("mPhi")

print("\nFinal training variables:")
for v in var:
    print("  ", v)

print("Number of DNN inputs:", len(var))
print("Final signal entries:", len(Sig))
print("Final training-background entries:", len(Bkg))
print("Signal fraction in the unbalanced full dataset:", df_all["isSignal"].mean())

print("\nSignal entries per mass:")
print(Sig.groupby("mPhi_true").size())

print("\nTraining-background entries per balance group:")
print(Bkg.groupby("bkg_group").size())

print("\nTraining-background entries per evaluation group:")
print(Bkg.groupby("eval_group").size())

os.makedirs("./Models", exist_ok=True)
os.makedirs("./Plots", exist_ok=True)

with open(f"./Models/features_{ANALYSIS_NAME}.txt", "w") as f:
    for feature in var:
        f.write(feature + "\n")

joblib.dump(var, f"./Models/features_{ANALYSIS_NAME}.pkl")
joblib.dump(signal_hypotheses, f"./Models/mPhi_hypotheses_{ANALYSIS_NAME}.pkl")

# ============================================================
# TRAIN / TEST SPLIT
# ============================================================

# Preserve every signal mass and every background training group in both
# partitions. This is more informative than stratifying only on isSignal.
df_all = df_all.copy()
df_all["split_stratum"] = "bkg_" + df_all["bkg_group"].astype(str)
signal_split_mask = df_all["isSignal"] == 1
df_all.loc[signal_split_mask, "split_stratum"] = (
    "sig_"
    + df_all.loc[signal_split_mask, "mPhi_true"].astype(int).astype(str)
)

stratum_counts = df_all["split_stratum"].value_counts()
rare_strata = stratum_counts[stratum_counts < 2]

if len(rare_strata) > 0:
    print("\nWARNING: Some strata have fewer than two events; falling back to class-only stratification:")
    print(rare_strata)
    stratify_labels = df_all["isSignal"]
else:
    stratify_labels = df_all["split_stratum"]

train_df, test_df = train_test_split(
    df_all,
    test_size=TEST_SIZE,
    random_state=SEED,
    stratify=stratify_labels,
)

train_df = train_df.reset_index(drop=True)
test_df = test_df.reset_index(drop=True)

print("\nTrain/test split")
print(f"  Train entries: {len(train_df)}")
print(f"  Test entries:  {len(test_df)}")
print("\nTrain composition:")
print(train_df.groupby(["isSignal", "bkg_group"]).size())
print("\nTest composition:")
print(test_df.groupby(["isSignal", "bkg_group"]).size())

# ============================================================
# PARAMETRIC INPUT BUILDING / PREPROCESSING
# ============================================================


def build_model_input_frame(
    frame,
    rng=None,
    background_mass=None,
):
    """
    Build the actual DNN inputs.

    Signal events always receive their generated mPhi.
    Background events receive either:
      - a fixed tested hypothesis (background_mass), or
      - a fresh random value sampled from the discrete signal hypotheses.
    """
    X = frame[features_to_read].copy()

    if not PARAMETRIC:
        return X

    masses = frame["mPhi_true"].to_numpy(dtype=float, copy=True)
    bkg_mask = frame["isSignal"].to_numpy(dtype=int) == 0

    if np.any(bkg_mask):
        if background_mass is not None:
            masses[bkg_mask] = float(background_mass)
        else:
            if rng is None:
                rng = np.random.default_rng()
            masses[bkg_mask] = rng.choice(
                signal_hypotheses,
                size=np.count_nonzero(bkg_mask),
                replace=True,
            )

    if np.isnan(masses).any():
        raise RuntimeError("NaN mPhi values remain after parametric input assignment.")

    X["mPhi"] = masses
    return X[var]


# Fit one scaler object, preserving compatibility with the current inference
# pattern (single feature vector including mPhi). The background masses used
# only for fitting the scaler follow the same discrete signal-hypothesis law.
preprocessing_rng = np.random.default_rng(SEED + 100)
X_train_for_scaler = build_model_input_frame(
    train_df,
    rng=preprocessing_rng,
)

scaler = StandardScaler()
scaler.fit(X_train_for_scaler[var])

joblib.dump(scaler, f"./Models/scaler_{ANALYSIS_NAME}.pkl")

# Fixed, unbalanced validation sample. Background masses are sampled from the
# same discrete hypothesis set, but only once for reproducible val_loss/val_auc.
validation_rng = np.random.default_rng(SEED + 200)
X_test_eval = build_model_input_frame(
    test_df,
    rng=validation_rng,
)
X_test_scaled = scaler.transform(X_test_eval[var])
Y_test = test_df[["isSignal"]].copy()

# A similarly fixed training view is only used for train-vs-test diagnostic
# plots and permutation importance. It is NOT used by the mini-batch trainer.
train_eval_rng = np.random.default_rng(SEED + 300)
X_train_eval = build_model_input_frame(
    train_df,
    rng=train_eval_rng,
)
X_train_scaled = scaler.transform(X_train_eval[var])
Y_train = train_df[["isSignal"]].copy()

# ============================================================
# CLASS- AND MASS-BALANCED PARAMETRIC BATCH GENERATOR
# ============================================================


def integer_partition(total, n_parts):
    if n_parts <= 0:
        raise ValueError("n_parts must be positive.")
    base = total // n_parts
    remainder = total % n_parts
    return [base + (i < remainder) for i in range(n_parts)]


def normalized_group_fractions(groups, requested_fractions=None):
    groups = list(groups)

    if requested_fractions is None:
        return {group: 1.0 / len(groups) for group in groups}

    unknown = set(requested_fractions) - set(groups)
    if unknown:
        raise ValueError(f"Unknown background groups in fractions: {sorted(unknown)}")

    fractions = {
        group: float(requested_fractions.get(group, 0.0))
        for group in groups
    }

    total = sum(fractions.values())
    if total <= 0:
        raise ValueError("BACKGROUND_BATCH_FRACTIONS must sum to a positive value.")

    return {group: value / total for group, value in fractions.items()}


def counts_from_fractions(total, fractions):
    groups = list(fractions)
    raw = np.array([fractions[g] * total for g in groups], dtype=float)
    counts = np.floor(raw).astype(int)

    remaining = total - counts.sum()
    if remaining > 0:
        order = np.argsort(raw - counts)[::-1]
        for i in order[:remaining]:
            counts[i] += 1

    return {group: int(count) for group, count in zip(groups, counts)}


def full_balanced_batch_generator(
    frame,
    scaler,
    batch_size=1024,
    background_fractions=None,
    seed=12345,
):
    """
    Full balance:
      1) 50% signal / 50% background;
      2) signal half balanced over all mPhi hypotheses;
      3) every relevant background group is explicitly represented;
         group fractions are either strict equal balance (paper_full) or the
         analysis-aware fractions configured above;
      4) every selected background event receives a freshly sampled mPhi.
    """
    rng = np.random.default_rng(seed)

    signal_half = batch_size // 2
    background_half = batch_size - signal_half

    signal_indices_by_mass = {}
    for mass in signal_hypotheses:
        idx = frame.index[
            (frame["isSignal"] == 1)
            & np.isclose(frame["mPhi_true"], mass)
        ].to_numpy()

        if len(idx) == 0:
            raise RuntimeError(f"No training signal events for mPhi={mass}.")

        signal_indices_by_mass[mass] = idx

    background_groups = sorted(
        frame.loc[frame["isSignal"] == 0, "bkg_group"].unique()
    )

    if len(background_groups) == 0:
        raise RuntimeError("No background groups in the training dataframe.")

    background_indices_by_group = {}
    for group in background_groups:
        idx = frame.index[
            (frame["isSignal"] == 0)
            & (frame["bkg_group"] == group)
        ].to_numpy()

        if len(idx) == 0:
            raise RuntimeError(f"No training background events for group={group}.")

        background_indices_by_group[group] = idx

    sig_counts = dict(zip(
        signal_hypotheses,
        integer_partition(signal_half, len(signal_hypotheses)),
    ))

    fractions = normalized_group_fractions(
        background_groups,
        requested_fractions=background_fractions,
    )
    bkg_counts = counts_from_fractions(background_half, fractions)

    print("\nBalanced mini-batch composition")
    print(f"  Batch size:       {batch_size}")
    print(f"  Signal half:      {signal_half}")
    print(f"  Background half:  {background_half}")
    print("  Signal events per mass:")
    for mass, count in sig_counts.items():
        print(f"    mPhi={mass:g}: {count}")
    print("  Background events per group:")
    for group, count in bkg_counts.items():
        print(f"    {group}: {count}")

    while True:
        selected_indices = []

        # Signal balance: every mass is present in every mini-batch.
        for mass, count in sig_counts.items():
            pool = signal_indices_by_mass[mass]
            selected_indices.extend(
                rng.choice(pool, size=count, replace=True).tolist()
            )

        # Background balance: every process group is present in every batch.
        for group, count in bkg_counts.items():
            pool = background_indices_by_group[group]
            selected_indices.extend(
                rng.choice(pool, size=count, replace=True).tolist()
            )

        selected_indices = np.asarray(selected_indices, dtype=int)
        batch_frame = frame.loc[selected_indices].copy()

        # Dynamic identical-distribution assignment for background mPhi.
        X_batch = build_model_input_frame(
            batch_frame,
            rng=rng if DYNAMIC_BACKGROUND_MPHI else np.random.default_rng(SEED),
        )
        y_batch = batch_frame["isSignal"].to_numpy(dtype=np.float32)

        X_batch_scaled = scaler.transform(X_batch[var]).astype(np.float32)

        permutation = rng.permutation(len(y_batch))
        yield X_batch_scaled[permutation], y_batch[permutation]


train_generator = full_balanced_batch_generator(
    train_df,
    scaler,
    batch_size=BATCH_SIZE,
    background_fractions=ACTIVE_BACKGROUND_BATCH_FRACTIONS,
    seed=SEED + 400,
)

if STEPS_PER_EPOCH_MODE == "signal_once":
    n_train_signal = int((train_df["isSignal"] == 1).sum())
    signal_per_batch = BATCH_SIZE // 2
    steps_per_epoch = max(1, int(np.ceil(n_train_signal / signal_per_batch)))
elif STEPS_PER_EPOCH_MODE == "all_events":
    steps_per_epoch = max(1, len(train_df) // BATCH_SIZE)
else:
    raise ValueError(
        f"Unknown STEPS_PER_EPOCH_MODE={STEPS_PER_EPOCH_MODE!r}. "
        "Use 'signal_once' or 'all_events'."
    )


def print_epoch_sampling_diagnostics(frame, steps, batch_size, background_fractions):
    signal_half = batch_size // 2
    background_half = batch_size - signal_half

    signal_counts = dict(zip(
        signal_hypotheses,
        integer_partition(signal_half, len(signal_hypotheses)),
    ))

    background_groups = sorted(
        frame.loc[frame["isSignal"] == 0, "bkg_group"].unique()
    )
    fractions = normalized_group_fractions(
        background_groups,
        requested_fractions=background_fractions,
    )
    background_counts = counts_from_fractions(background_half, fractions)

    print("\nEpoch sampling diagnostics")
    print(f"  steps_per_epoch: {steps}")
    print(f"  mode:            {STEPS_PER_EPOCH_MODE}")
    print("  Signal draws per epoch:")

    for mass, per_batch in signal_counts.items():
        pool = int(np.sum(
            (frame["isSignal"] == 1)
            & np.isclose(frame["mPhi_true"], mass)
        ))
        draws = per_batch * steps
        reuse = draws / pool if pool > 0 else np.nan
        print(
            f"    mPhi={mass:g}: draws={draws:6d}, "
            f"pool={pool:6d}, draws/pool={reuse:6.2f}"
        )

    print("  Background draws per epoch:")
    for group, per_batch in background_counts.items():
        pool = int(np.sum(
            (frame["isSignal"] == 0)
            & (frame["bkg_group"] == group)
        ))
        draws = per_batch * steps
        reuse = draws / pool if pool > 0 else np.nan
        print(
            f"    {group:12s}: draws={draws:6d}, "
            f"pool={pool:7d}, draws/pool={reuse:6.2f}"
        )


print_epoch_sampling_diagnostics(
    train_df,
    steps_per_epoch,
    BATCH_SIZE,
    ACTIVE_BACKGROUND_BATCH_FRACTIONS,
)

# ============================================================
# PER-MASS VALIDATION CALLBACK
# ============================================================


def evaluate_mass_hypothesis(
    model,
    frame,
    mass,
    scaler,
    batch_size=PREDICT_BATCH_SIZE,
):
    signal_mass = frame[
        (frame["isSignal"] == 1)
        & np.isclose(frame["mPhi_true"], mass)
    ]
    backgrounds = frame[frame["isSignal"] == 0]

    evaluation_frame = pd.concat(
        [signal_mass, backgrounds],
        ignore_index=True,
    )

    if len(signal_mass) == 0 or len(backgrounds) == 0:
        return None

    X_mass = build_model_input_frame(
        evaluation_frame,
        background_mass=mass,
    )
    X_mass_scaled = scaler.transform(X_mass[var])
    y_mass = evaluation_frame["isSignal"].to_numpy(dtype=int)

    scores = model.predict(
        X_mass_scaled,
        batch_size=batch_size,
        verbose=0,
    ).reshape(-1)

    auc_mass = roc_auc_score(y_mass, scores)
    fpr_mass, tpr_mass, thresholds_mass = metrics.roc_curve(y_mass, scores)

    # Significance-ratio diagnostic from the pNN reference. For evaluation
    # only, signal and background totals are normalized to the same value;
    # this has no physical luminosity or cross-section meaning. With equal
    # total class weights the ratio simplifies to tpr/sqrt(tpr + fpr).
    denominator = np.sqrt(tpr_mass + fpr_mass)
    sigma_ratio_curve = np.divide(
        tpr_mass,
        denominator,
        out=np.zeros_like(tpr_mass, dtype=float),
        where=denominator > 0,
    )
    best_sigma_index = int(np.argmax(sigma_ratio_curve))
    sigma_ratio = float(sigma_ratio_curve[best_sigma_index])
    best_cut = float(thresholds_mass[best_sigma_index])

    return {
        "mass": float(mass),
        "auc": float(auc_mass),
        "sigma_ratio": sigma_ratio,
        "best_cut": best_cut,
        "fpr": fpr_mass,
        "tpr": tpr_mass,
        "thresholds": thresholds_mass,
        "scores": scores,
        "labels": y_mass,
    }


class PerMassValidationCallback(callbacks.Callback):
    def __init__(
        self,
        validation_frame,
        scaler,
        masses,
        batch_size=PREDICT_BATCH_SIZE,
    ):
        super().__init__()
        self.validation_frame = validation_frame
        self.scaler = scaler
        self.masses = list(masses)
        self.batch_size = batch_size
        self.rows = []

    def on_epoch_end(self, epoch, logs=None):
        logs = logs if logs is not None else {}

        aucs = []
        epoch_row = {"epoch": epoch + 1}

        for mass in self.masses:
            result = evaluate_mass_hypothesis(
                self.model,
                self.validation_frame,
                mass,
                self.scaler,
                batch_size=self.batch_size,
            )

            if result is None:
                continue

            auc_mass = result["auc"]
            aucs.append(auc_mass)
            epoch_row[f"auc_mPhi_{mass:g}"] = auc_mass
            epoch_row[f"sigmaRatio_mPhi_{mass:g}"] = result["sigma_ratio"]
            epoch_row[f"bestCut_mPhi_{mass:g}"] = result["best_cut"]

        if len(aucs) == 0:
            raise RuntimeError("Could not compute any per-mass validation AUC.")

        mean_auc = float(np.mean(aucs))
        min_auc = float(np.min(aucs))

        logs["val_mean_mass_auc"] = mean_auc
        logs["val_min_mass_auc"] = min_auc

        epoch_row["val_mean_mass_auc"] = mean_auc
        epoch_row["val_min_mass_auc"] = min_auc
        self.rows.append(epoch_row)

        print(
            f" - val_mean_mass_auc: {mean_auc:.4f}"
            f" - val_min_mass_auc: {min_auc:.4f}"
        )


per_mass_validation = PerMassValidationCallback(
    validation_frame=test_df,
    scaler=scaler,
    masses=signal_hypotheses,
)

# ============================================================
# MODEL DEFINITIONS
# ============================================================

# The vanilla network architecture is intentionally kept identical to the
# current code so the first comparison isolates the effect of the improved
# pNN training procedure. The reference paper explicitly notes that layer
# widths and related architecture hyperparameters are dataset dependent.


class AffineConditioning(tf.keras.layers.Layer):
    def __init__(self, units, l2_weight=1e-5, l2_bias=1e-6, **kwargs):
        super().__init__(**kwargs)
        self.units = units
        self.l2_weight = l2_weight
        self.l2_bias = l2_bias

        self.scale = Dense(
            units,
            activation=None,
            kernel_regularizer=regularizers.l2(l2_weight),
            bias_regularizer=regularizers.l2(l2_bias),
        )
        self.bias = Dense(
            units,
            activation=None,
            kernel_regularizer=regularizers.l2(l2_weight),
            bias_regularizer=regularizers.l2(l2_bias),
        )

    def call(self, inputs):
        hidden, mass = inputs
        return hidden * self.scale(mass) + self.bias(mass)

    def get_config(self):
        config = super().get_config()
        config.update({
            "units": self.units,
            "l2_weight": self.l2_weight,
            "l2_bias": self.l2_bias,
        })
        return config


def build_vanilla_model(n_inputs):
    model = tf.keras.Sequential([
        Input(shape=(n_inputs,)),

        Dense(128, activation="relu", kernel_regularizer=regularizers.l2(1e-4)),
        BatchNormalization(),
        Dropout(0.3),

        Dense(64, activation="relu", kernel_regularizer=regularizers.l2(1e-4)),
        BatchNormalization(),
        Dropout(0.3),

        Dense(32, activation="relu", kernel_regularizer=regularizers.l2(1e-4)),
        BatchNormalization(),
        Dropout(0.2),

        Dense(16, activation="relu", kernel_regularizer=regularizers.l2(1e-4)),
        Dense(1, activation="sigmoid"),
    ])

    return model


def build_affine_model(n_inputs):
    if not PARAMETRIC or "mPhi" not in var:
        raise RuntimeError("Affine model requires PARAMETRIC=True and mPhi in the inputs.")

    inputs = Input(shape=(n_inputs,), name="inputs")

    # mPhi is the final feature by construction. Keep the same one-vector input
    # interface as the current inference code and split it only inside the model.
    physics_inputs = inputs[:, :-1]
    mass_input = inputs[:, -1:]

    h = Dense(128, activation="relu", kernel_regularizer=regularizers.l2(1e-4))(physics_inputs)
    h = AffineConditioning(128, name="affine_128")([h, mass_input])
    h = Dropout(0.3)(h)

    h = Dense(64, activation="relu", kernel_regularizer=regularizers.l2(1e-4))(h)
    h = AffineConditioning(64, name="affine_64")([h, mass_input])
    h = Dropout(0.3)(h)

    h = Dense(32, activation="relu", kernel_regularizer=regularizers.l2(1e-4))(h)
    h = AffineConditioning(32, name="affine_32")([h, mass_input])
    h = Dropout(0.2)(h)

    h = Dense(16, activation="relu", kernel_regularizer=regularizers.l2(1e-4))(h)
    h = AffineConditioning(16, name="affine_16")([h, mass_input])

    output = Dense(1, activation="sigmoid", name="output")(h)
    return tf.keras.Model(inputs=inputs, outputs=output)


if MODEL_ARCHITECTURE == "vanilla":
    model = build_vanilla_model(len(var))
elif MODEL_ARCHITECTURE == "affine":
    model = build_affine_model(len(var))
else:
    raise ValueError(
        f"Unknown MODEL_ARCHITECTURE={MODEL_ARCHITECTURE!r}. "
        "Use 'vanilla' or 'affine'."
    )

model.compile(
    loss="binary_crossentropy",
    optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
    metrics=[tf.keras.metrics.AUC(name="auc")],
)

model.summary()

# The custom per-mass callback must come before EarlyStopping and
# ReduceLROnPlateau so the monitored metric is inserted into logs first.
early = callbacks.EarlyStopping(
    monitor="val_mean_mass_auc",
    mode="max",
    patience=15,
    restore_best_weights=True,
    verbose=1,
)

reduce_lr = callbacks.ReduceLROnPlateau(
    monitor="val_mean_mass_auc",
    mode="max",
    factor=0.5,
    patience=8,
    min_lr=1e-6,
    verbose=1,
)

# ============================================================
# TRAIN
# ============================================================

if loaded_model:
    print(f"\nLoading model: {MODEL_NAME}")
    model = tf.keras.models.load_model(
        MODEL_NAME,
        custom_objects={"AffineConditioning": AffineConditioning},
        compile=False,
    )
    model.compile(
        loss="binary_crossentropy",
        optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
        metrics=[tf.keras.metrics.AUC(name="auc")],
    )
    training = None
else:
    training = model.fit(
        train_generator,
        steps_per_epoch=steps_per_epoch//5,
        epochs=MAX_EPOCHS,
        # Original, unbalanced validation distribution. The background mass
        # hypothesis is sampled from the same discrete signal set.
        validation_data=(X_test_scaled, Y_test["isSignal"].to_numpy()),
        callbacks=[
            per_mass_validation,
            early,
            reduce_lr,
        ],
        verbose=2,
    )

# ============================================================
# SAVE TRAINING HISTORY / CONFIGURATION
# ============================================================

if training is not None:
    history_df = pd.DataFrame(training.history)
    history_df.to_csv(
        f"./Plots/TrainingHistory_{ANALYSIS_NAME}.csv",
        index=False,
    )

    loss = training.history["loss"]
    val_loss = training.history.get("val_loss", [])
    epochs = range(1, len(loss) + 1)

    plt.figure()
    plt.plot(epochs, loss)
    if len(val_loss) == len(loss):
        plt.plot(epochs, val_loss)
        plt.legend(["Train Loss", "Validation Loss"])
    else:
        plt.legend(["Train Loss"])
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("Training and Validation Loss")
    plt.grid(True)
    plt.yscale("log")
    plt.savefig(
        f"./Plots/loss_{ANALYSIS_NAME}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

if len(per_mass_validation.rows) > 0:
    per_mass_history_df = pd.DataFrame(per_mass_validation.rows)
    per_mass_history_df.to_csv(
        f"./Plots/PerMassValidationHistory_{ANALYSIS_NAME}.csv",
        index=False,
    )

    plt.figure(figsize=(8, 6))
    plt.plot(
        per_mass_history_df["epoch"],
        per_mass_history_df["val_mean_mass_auc"],
        label="Mean per-mass AUC",
    )
    plt.plot(
        per_mass_history_df["epoch"],
        per_mass_history_df["val_min_mass_auc"],
        label="Minimum per-mass AUC",
    )
    plt.xlabel("Epoch")
    plt.ylabel("AUC")
    plt.legend()
    plt.grid(True)
    plt.savefig(
        f"./Plots/PerMassValidation_{ANALYSIS_NAME}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

training_config = {
    "analysis_name": ANALYSIS_NAME,
    "signal_type": SIGNAL_TYPE,
    "signal_pattern": SIGNAL_PATTERN,
    "training_target": TRAINING_TARGET,
    "model_architecture": MODEL_ARCHITECTURE,
    "parametric": PARAMETRIC,
    "signal_hypotheses": signal_hypotheses,
    "batch_size": BATCH_SIZE,
    "steps_per_epoch": steps_per_epoch,
    "steps_per_epoch_mode": STEPS_PER_EPOCH_MODE,
    "test_size": TEST_SIZE,
    "dynamic_background_mPhi": DYNAMIC_BACKGROUND_MPHI,
    "background_balance_mode": BACKGROUND_BALANCE_MODE,
    "background_batch_fractions": ACTIVE_BACKGROUND_BATCH_FRACTIONS,
    "training_background_groups": sorted(training_background_files),
    "features": var,
}

import json
with open(f"./Models/config_{ANALYSIS_NAME}.json", "w") as f:
    json.dump(training_config, f, indent=2)

# ============================================================
# SAVE MODEL / SCALER
# ============================================================

if save_model and not loaded_model:
    print("\nSaving model and preprocessing objects")
    print(f"  Model: {MODEL_NAME}")
    model.save(MODEL_NAME)

# ============================================================
# PREDICTIONS ON THE FIXED UNBALANCED TRAIN/TEST VIEWS
# ============================================================

y_pred = model.predict(
    X_test_scaled,
    batch_size=PREDICT_BATCH_SIZE,
    verbose=0,
)
y_pred_t = model.predict(
    X_train_scaled,
    batch_size=PREDICT_BATCH_SIZE,
    verbose=0,
)

y_pred_L = y_pred
y_pred_t_L = y_pred_t

discriminant = np.squeeze(np.asarray(y_pred_L))
true_labels = np.squeeze(np.asarray(Y_test["isSignal"]))

discriminant0 = discriminant[np.array(true_labels == 0)]
discriminant1 = discriminant[np.array(true_labels == 1)]

discriminant_t = np.squeeze(np.asarray(y_pred_t_L))
true_labels_t = np.squeeze(np.asarray(Y_train["isSignal"]))

discriminant0_t = discriminant_t[np.array(true_labels_t == 0)]
discriminant1_t = discriminant_t[np.array(true_labels_t == 1)]

binning = np.linspace(0, 1, 51)

# ============================================================
# DNN OUTPUT PLOTS
# ============================================================

print("\nPlotting discriminant distributions")

plt.clf()
plt.figure(figsize=(6, 6))
plt.subplot(111)

pdf0, bins0, _ = plt.hist(
    discriminant0,
    bins=binning,
    alpha=0.0,
    histtype="stepfilled",
    linewidth=1,
    edgecolor="r",
    density=True,
)

pdf1, bins1, _ = plt.hist(
    discriminant1,
    bins=binning,
    alpha=0.0,
    histtype="stepfilled",
    linewidth=1,
    edgecolor="b",
    density=True,
)

plt.hist(
    discriminant0_t,
    bins=binning,
    alpha=0.3,
    histtype="stepfilled",
    linewidth=2,
    edgecolor="r",
    label="Backgrounds (train)",
    density=True,
)

plt.hist(
    discriminant1_t,
    bins=binning,
    alpha=0.3,
    histtype="stepfilled",
    linewidth=2,
    edgecolor="b",
    label="ttDM signal (train)",
    density=True,
)

plt.scatter(
    bins0[:-1] + 0.5 * (bins0[1:] - bins0[:-1]),
    pdf0,
    marker=".",
    s=30,
    alpha=0.8,
    label="Backgrounds",
)

plt.scatter(
    bins1[:-1] + 0.5 * (bins1[1:] - bins1[:-1]),
    pdf1,
    marker=".",
    s=30,
    alpha=0.8,
    label="ttDM signal",
)

plt.legend(loc="upper center")
plt.ylabel("Density", fontsize=12)
plt.xlabel("DNN discriminant", fontsize=12)
plt.savefig(
    f"./Plots/Discriminant_distribution_{ANALYSIS_NAME}.png",
    dpi=600,
)
plt.close()

plt.clf()
plt.figure(figsize=(6, 6))
plt.subplot(111)

pdf0, bins0, _ = plt.hist(
    discriminant0,
    bins=binning,
    alpha=0.0,
    histtype="stepfilled",
    linewidth=1,
    edgecolor="r",
    density=True,
)

pdf1, bins1, _ = plt.hist(
    discriminant1,
    bins=binning,
    alpha=0.0,
    histtype="stepfilled",
    linewidth=1,
    edgecolor="b",
    density=True,
)

plt.hist(
    discriminant0_t,
    bins=binning,
    alpha=0.3,
    histtype="stepfilled",
    linewidth=2,
    edgecolor="r",
    label="Backgrounds (train)",
    density=True,
)

plt.hist(
    discriminant1_t,
    bins=binning,
    alpha=0.3,
    histtype="stepfilled",
    linewidth=2,
    edgecolor="b",
    label="ttDM signal (train)",
    density=True,
)

plt.scatter(
    bins0[:-1] + 0.5 * (bins0[1:] - bins0[:-1]),
    pdf0,
    marker=".",
    s=30,
    alpha=0.8,
    label="Backgrounds",
)

plt.scatter(
    bins1[:-1] + 0.5 * (bins1[1:] - bins1[:-1]),
    pdf1,
    marker=".",
    s=30,
    alpha=0.8,
    label="ttDM signal",
)

plt.legend(loc="upper center")
plt.yscale("log")
plt.ylabel("Density", fontsize=12)
plt.xlabel("DNN discriminant", fontsize=12)
plt.savefig(
    f"./Plots/Log_Discriminant_distribution_{ANALYSIS_NAME}.png",
    dpi=600,
)
plt.close()

# ============================================================
# GLOBAL ROC ON FIXED UNBALANCED TEST VIEW
# ============================================================

print("\nPlotting global ROC")

fpr, tpr, thresholds = metrics.roc_curve(Y_test["isSignal"], y_pred_L)
auc = metrics.auc(fpr, tpr)

plt.clf()
plt.figure(figsize=(6, 6))
plt.subplot(111)
plt.plot(fpr, tpr, label="ROC curve")
plt.plot([0, 1], [0, 1], lw=2, linestyle="--", label="Random guess")
plt.legend(loc="lower right")
plt.xlabel("False Positive rate", fontsize=12)
plt.ylabel("True Positive rate", fontsize=12)
plt.text(0.55, 0.3, f"AUC = {auc:.3f}", fontsize=12)
plt.axvline(x=0, linestyle="--", linewidth=0.5)
plt.axhline(y=1, linestyle="--", linewidth=0.5)
plt.savefig(f"./Plots/ROC_{ANALYSIS_NAME}.png", dpi=600)
plt.close()

print("The global AUC of the model is:", auc)

# ============================================================
# CORRECT PER-MASS ROCS
# ============================================================

# For each tested mPhi, use only signal generated at that mass and ALL held-out
# backgrounds, with every background event evaluated at that same mPhi.
print("\nPlotting individual ROCs for different mPhi hypotheses")

roc_curves_by_mass = []
per_mass_auc_rows = []

for mass in signal_hypotheses:
    result = evaluate_mass_hypothesis(
        model,
        test_df,
        mass,
        scaler,
    )

    if result is None:
        print(f"Skipping mPhi={mass}: could not build both classes")
        continue

    roc_curves_by_mass.append(
        (mass, result["fpr"], result["tpr"], result["auc"])
    )
    per_mass_auc_rows.append({
        "mPhi": mass,
        "auc": result["auc"],
        "sigma_ratio": result["sigma_ratio"],
        "best_cut": result["best_cut"],
    })

if roc_curves_by_mass:
    plt.clf()
    plt.figure(figsize=(8, 7))
    plt.subplot(111)

    color_map = plt.get_cmap("tab20", len(roc_curves_by_mass))

    print("AUC per mPhi hypothesis:")
    for i, (mass, fpr_mass, tpr_mass, auc_mass) in enumerate(roc_curves_by_mass):
        diagnostic = next(
            row for row in per_mass_auc_rows
            if np.isclose(row["mPhi"], mass)
        )
        print(
            f"  mPhi={mass:g}: AUC={auc_mass:.4f}, "
            f"sigma_ratio={diagnostic['sigma_ratio']:.4f}, "
            f"best_cut={diagnostic['best_cut']:.4f}"
        )
        plt.plot(
            fpr_mass,
            tpr_mass,
            color=color_map(i),
            label=rf"$m_\Phi={mass:g}$ (AUC={auc_mass:.3f})",
        )

    plt.plot([0, 1], [0, 1], lw=2, linestyle="--", label="Random guess")
    plt.xlabel("False Positive rate", fontsize=12)
    plt.ylabel("True Positive rate", fontsize=12)
    plt.legend(loc="lower right", fontsize=9)
    plt.axvline(x=0, linestyle="--", linewidth=0.5)
    plt.axhline(y=1, linestyle="--", linewidth=0.5)
    plt.savefig(f"./Plots/ROC_{ANALYSIS_NAME}_mPhi.png", dpi=600)
    plt.close()

    pd.DataFrame(per_mass_auc_rows).to_csv(
        f"./Plots/Metrics_perMass_{ANALYSIS_NAME}.csv",
        index=False,
    )
else:
    print("No per-mPhi ROC could be produced.")

# ============================================================
# PER-MASS AUC AGAINST EACH BACKGROUND GROUP
# ============================================================

# This is the key diagnostic for checking that adding rare backgrounds does
# not degrade the original ttDM-vs-TTTo2L2Nu separation.
group_auc_frames = [test_df[test_df["isSignal"] == 0].copy()]

if len(TTNuNu) > 0 and TRAINING_TARGET == "main":
    if len(TTNuNu) > 1:
        _, ttnunu_group_test = train_test_split(
            TTNuNu,
            test_size=TEST_SIZE,
            random_state=SEED,
        )
    else:
        ttnunu_group_test = TTNuNu.copy()
    group_auc_frames.append(ttnunu_group_test)

background_group_auc_df = pd.concat(group_auc_frames, ignore_index=True)
background_group_auc_rows = []

print("\nComputing per-mass AUC against each background group")

for mass in signal_hypotheses:
    signal_mass = test_df[
        (test_df["isSignal"] == 1)
        & np.isclose(test_df["mPhi_true"], mass)
    ].copy()

    if len(signal_mass) == 0:
        continue

    for group, group_frame in background_group_auc_df.groupby("eval_group"):
        if len(group_frame) == 0:
            continue

        evaluation_frame = pd.concat(
            [signal_mass, group_frame],
            ignore_index=True,
        )

        X_group_auc = build_model_input_frame(
            evaluation_frame,
            background_mass=mass,
        )
        X_group_auc_scaled = scaler.transform(X_group_auc[var])
        y_group_auc = evaluation_frame["isSignal"].to_numpy(dtype=int)

        scores_group_auc = model.predict(
            X_group_auc_scaled,
            batch_size=PREDICT_BATCH_SIZE,
            verbose=0,
        ).reshape(-1)

        auc_group = roc_auc_score(y_group_auc, scores_group_auc)
        background_group_auc_rows.append({
            "mPhi": mass,
            "background_group": group,
            "auc": float(auc_group),
            "n_signal": len(signal_mass),
            "n_background": len(group_frame),
        })

if background_group_auc_rows:
    background_group_auc_table = pd.DataFrame(background_group_auc_rows)
    background_group_auc_table.to_csv(
        f"./Plots/AUC_perMass_perBackground_{ANALYSIS_NAME}.csv",
        index=False,
    )

    print("\nPer-mass AUC against each background group:")
    print(
        background_group_auc_table.pivot(
            index="mPhi",
            columns="background_group",
            values="auc",
        ).to_string()
    )

# ============================================================
# BACKGROUND COMPOSITION IN THE DNN TAIL
# ============================================================

# Use held-out events for every background that participated in training.
# TTNuNu was excluded from the main DNN, so take an independent reproducible
# test-sized subset for a statistically comparable diagnostic.
tail_frames = [test_df[test_df["isSignal"] == 0].copy()]

if len(TTNuNu) > 0 and TRAINING_TARGET == "main":
    if len(TTNuNu) > 1:
        _, ttnunu_test = train_test_split(
            TTNuNu,
            test_size=TEST_SIZE,
            random_state=SEED,
        )
    else:
        ttnunu_test = TTNuNu.copy()
    tail_frames.append(ttnunu_test)

background_tail_df = pd.concat(tail_frames, ignore_index=True)

tail_rows = []

print("\nStudying background composition in the DNN tail")

for mass in signal_hypotheses:
    for group, group_frame in background_tail_df.groupby("eval_group"):
        if len(group_frame) == 0:
            continue

        X_group = build_model_input_frame(
            group_frame,
            background_mass=mass,
        )
        X_group_scaled = scaler.transform(X_group[var])

        scores = model.predict(
            X_group_scaled,
            batch_size=PREDICT_BATCH_SIZE,
            verbose=0,
        ).reshape(-1)

        row = {
            "mPhi": mass,
            "group": group,
            "n_total": len(scores),
            "mean_score": float(np.mean(scores)),
        }

        for threshold in TAIL_THRESHOLDS:
            suffix = str(threshold).replace(".", "p")
            selected = scores > threshold
            row[f"n_gt_{suffix}"] = int(np.sum(selected))
            row[f"eff_gt_{suffix}"] = float(np.mean(selected))

        tail_rows.append(row)

if len(tail_rows) > 0:
    tail_table = pd.DataFrame(tail_rows)
    tail_table.to_csv(
        f"./Plots/BackgroundTail_{ANALYSIS_NAME}.csv",
        index=False,
    )

    for mass in signal_hypotheses:
        mass_table = tail_table[np.isclose(tail_table["mPhi"], mass)].copy()
        if len(mass_table) == 0:
            continue

        tightest_suffix = str(max(TAIL_THRESHOLDS)).replace(".", "p")
        mass_table = mass_table.sort_values(
            f"eff_gt_{tightest_suffix}",
            ascending=False,
        )

        print(f"\nBackground tail efficiencies for mPhi={mass:g}")
        print(mass_table.to_string(index=False))

# ============================================================
# SUMMARY DIAGNOSTIC PLOTS
# ============================================================


def save_summary_plot(filename):
    """Save the current figure in the common Plots directory."""
    plt.tight_layout()
    plt.savefig(
        f"./Plots/{filename}_{ANALYSIS_NAME}.png",
        dpi=SUMMARY_PLOT_DPI,
        bbox_inches="tight",
    )
    plt.close()


def ordered_existing_groups(groups):
    preferred_order = [
        "tt2l",
        "other_tt",
        "single_top",
        "other_ttZ",
        "ttW",
        "ttH",
        "TTNuNu",
    ]
    groups = list(groups)
    ordered = [group for group in preferred_order if group in groups]
    ordered.extend(group for group in groups if group not in ordered)
    return ordered


def nearest_available_mass(requested_mass, available_masses):
    available = np.asarray(sorted(available_masses), dtype=float)
    if len(available) == 0:
        return None
    return float(available[np.argmin(np.abs(available - requested_mass))])


if MAKE_SUMMARY_DIAGNOSTIC_PLOTS:
    print("\nCreating summary diagnostic plots")

    # --------------------------------------------------------
    # 1. Global per-mass performance
    # --------------------------------------------------------
    metrics_per_mass_df = pd.DataFrame(per_mass_auc_rows).sort_values("mPhi")

    if len(metrics_per_mass_df) > 0:
        metrics_per_mass_df.to_csv(
            f"./Plots/DiagnosticSummary_perMass_{ANALYSIS_NAME}.csv",
            index=False,
        )

        plt.figure(figsize=(8, 6))
        plt.plot(
            metrics_per_mass_df["mPhi"],
            metrics_per_mass_df["auc"],
            marker="o",
            linewidth=2,
            label="ROC AUC",
        )
        plt.plot(
            metrics_per_mass_df["mPhi"],
            metrics_per_mass_df["sigma_ratio"],
            marker="s",
            linewidth=2,
            label=r"$\sigma_{\rm ratio}$",
        )
        plt.xlabel(r"$m_\Phi$ [GeV]")
        plt.ylabel("Performance metric")
        plt.ylim(0.0, 1.02)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.title("Main DNN performance versus signal mass")
        save_summary_plot("PerformanceMetrics_vs_mPhi")

        plt.figure(figsize=(8, 6))
        plt.plot(
            metrics_per_mass_df["mPhi"],
            metrics_per_mass_df["best_cut"],
            marker="o",
            linewidth=2,
        )
        plt.xlabel(r"$m_\Phi$ [GeV]")
        plt.ylabel("Best DNN threshold")
        plt.ylim(0.0, 1.0)
        plt.grid(True, alpha=0.3)
        plt.title(r"Threshold maximizing $\sigma_{\rm ratio}$")
        save_summary_plot("BestCut_vs_mPhi")

    # --------------------------------------------------------
    # 2. Per-background AUC: full picture and heatmap
    # --------------------------------------------------------
    if "background_group_auc_table" in globals() and len(background_group_auc_table) > 0:
        auc_pivot = background_group_auc_table.pivot(
            index="mPhi",
            columns="background_group",
            values="auc",
        ).sort_index()

        auc_groups = ordered_existing_groups(auc_pivot.columns)
        auc_pivot = auc_pivot[auc_groups]

        plt.figure(figsize=(10, 7))
        for group in auc_groups:
            plt.plot(
                auc_pivot.index,
                auc_pivot[group],
                marker="o",
                linewidth=2,
                label=group,
            )
        plt.xlabel(r"$m_\Phi$ [GeV]")
        plt.ylabel("ROC AUC")
        plt.ylim(0.45, 1.01)
        plt.grid(True, alpha=0.3)
        plt.legend(ncol=2)
        plt.title("ttDM discrimination against each background")
        save_summary_plot("AUC_perBackground_vs_mPhi")

        # Annotated heatmap: one compact view of the complete table.
        heatmap_values = auc_pivot.T.to_numpy(dtype=float)
        plt.figure(figsize=(13, 6))
        image = plt.imshow(
            heatmap_values,
            aspect="auto",
            origin="lower",
            vmin=0.45,
            vmax=1.0,
        )
        plt.colorbar(image, label="ROC AUC")
        plt.xticks(
            np.arange(len(auc_pivot.index)),
            [f"{mass:g}" for mass in auc_pivot.index],
            rotation=45,
            ha="right",
        )
        plt.yticks(np.arange(len(auc_groups)), auc_groups)
        plt.xlabel(r"$m_\Phi$ [GeV]")
        plt.ylabel("Background")
        plt.title("Per-mass discrimination matrix")

        for row_idx in range(heatmap_values.shape[0]):
            for col_idx in range(heatmap_values.shape[1]):
                value = heatmap_values[row_idx, col_idx]
                if np.isfinite(value):
                    plt.text(
                        col_idx,
                        row_idx,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        fontsize=8,
                    )
        save_summary_plot("AUC_perBackground_heatmap")

        # ----------------------------------------------------
        # 3. Dedicated protection plot for tt(2l)
        # ----------------------------------------------------
        if "tt2l" in auc_pivot.columns:
            plt.figure(figsize=(8, 6))
            plt.plot(
                auc_pivot.index,
                auc_pivot["tt2l"],
                marker="o",
                linewidth=2,
            )
            plt.xlabel(r"$m_\Phi$ [GeV]")
            plt.ylabel("ROC AUC")
            plt.ylim(0.5, 1.01)
            plt.grid(True, alpha=0.3)
            plt.title(r"Main task protection: ttDM vs $t\bar{t}(2\ell)$")
            save_summary_plot("AUC_ttDM_vs_tt2l")

        # ----------------------------------------------------
        # 4. Difficulty ranking at representative masses
        # ----------------------------------------------------
        available_auc_masses = list(auc_pivot.index)
        used_masses = set()
        for requested_mass in SUMMARY_MASSES:
            mass = nearest_available_mass(requested_mass, available_auc_masses)
            if mass is None or mass in used_masses:
                continue
            used_masses.add(mass)

            ranking = auc_pivot.loc[mass].dropna().sort_values(ascending=True)
            plt.figure(figsize=(8, max(5, 0.55 * len(ranking))))
            positions = np.arange(len(ranking))
            plt.barh(positions, ranking.values)
            plt.yticks(positions, ranking.index)
            plt.xlabel("ROC AUC")
            plt.xlim(0.45, 1.0)
            plt.grid(True, axis="x", alpha=0.3)
            plt.title(rf"Background difficulty at $m_\Phi={mass:g}$ GeV")
            for pos, value in zip(positions, ranking.values):
                plt.text(value + 0.005, pos, f"{value:.3f}", va="center")
            save_summary_plot(f"BackgroundDifficulty_mPhi{mass:g}")

        # ----------------------------------------------------
        # 5. Main-DNN vs dedicated-TTNuNu AUC comparison
        # ----------------------------------------------------
        if "TTNuNu" in auc_pivot.columns and DEDICATED_TTNUNU_AUC_REFERENCE:
            comparison_rows = []
            for mass, main_auc in auc_pivot["TTNuNu"].dropna().items():
                mass_float = float(mass)
                if mass_float in DEDICATED_TTNUNU_AUC_REFERENCE:
                    dedicated_auc = DEDICATED_TTNUNU_AUC_REFERENCE[mass_float]
                    comparison_rows.append({
                        "mPhi": mass_float,
                        "main_auc": float(main_auc),
                        "dedicated_auc": float(dedicated_auc),
                        "dedicated_minus_main": float(dedicated_auc - main_auc),
                    })

            ttnunu_comparison_df = pd.DataFrame(comparison_rows).sort_values("mPhi")
            if len(ttnunu_comparison_df) > 0:
                ttnunu_comparison_df.to_csv(
                    f"./Plots/TTNuNu_MainVsDedicated_{ANALYSIS_NAME}.csv",
                    index=False,
                )

                plt.figure(figsize=(9, 6))
                plt.plot(
                    ttnunu_comparison_df["mPhi"],
                    ttnunu_comparison_df["main_auc"],
                    marker="o",
                    linewidth=2,
                    label="Main DNN vs TTNuNu",
                )
                plt.plot(
                    ttnunu_comparison_df["mPhi"],
                    ttnunu_comparison_df["dedicated_auc"],
                    marker="s",
                    linewidth=2,
                    label="Dedicated TTNuNu DNN",
                )
                plt.xlabel(r"$m_\Phi$ [GeV]")
                plt.ylabel("ROC AUC against TTNuNu")
                plt.ylim(0.45, 0.85)
                plt.grid(True, alpha=0.3)
                plt.legend()
                plt.title("Does the dedicated TTNuNu classifier add separation?")
                save_summary_plot("TTNuNu_MainVsDedicated_AUC")

                plt.figure(figsize=(9, 6))
                plt.plot(
                    ttnunu_comparison_df["mPhi"],
                    ttnunu_comparison_df["dedicated_minus_main"],
                    marker="o",
                    linewidth=2,
                )
                plt.axhline(0.0, linestyle="--", linewidth=1)
                plt.xlabel(r"$m_\Phi$ [GeV]")
                plt.ylabel("AUC(dedicated) - AUC(main)")
                plt.grid(True, alpha=0.3)
                plt.title("Additional TTNuNu discrimination from the second DNN")
                save_summary_plot("TTNuNu_DedicatedGain_vs_mPhi")

    # --------------------------------------------------------
    # 6. Background scores and high-DNN tails
    # --------------------------------------------------------
    if "tail_table" in globals() and len(tail_table) > 0:
        tail_groups = ordered_existing_groups(tail_table["group"].unique())

        plt.figure(figsize=(10, 7))
        for group in tail_groups:
            group_data = tail_table[tail_table["group"] == group].sort_values("mPhi")
            plt.plot(
                group_data["mPhi"],
                group_data["mean_score"],
                marker="o",
                linewidth=2,
                label=group,
            )
        plt.xlabel(r"$m_\Phi$ [GeV]")
        plt.ylabel("Mean main-DNN score")
        plt.ylim(0.0, 1.0)
        plt.grid(True, alpha=0.3)
        plt.legend(ncol=2)
        plt.title("How signal-like each background looks")
        save_summary_plot("MeanScore_perBackground_vs_mPhi")

        for threshold in TAIL_THRESHOLDS:
            suffix = str(threshold).replace(".", "p")
            eff_column = f"eff_gt_{suffix}"

            plt.figure(figsize=(10, 7))
            for group in tail_groups:
                group_data = tail_table[
                    tail_table["group"] == group
                ].sort_values("mPhi")
                plt.plot(
                    group_data["mPhi"],
                    group_data[eff_column],
                    marker="o",
                    linewidth=2,
                    label=group,
                )
            plt.xlabel(r"$m_\Phi$ [GeV]")
            plt.ylabel(rf"Background efficiency for DNN $>{threshold:.2f}$")
            plt.yscale("log")
            plt.grid(True, which="both", alpha=0.3)
            plt.legend(ncol=2)
            plt.title(rf"Background survival in the DNN $>{threshold:.2f}$ tail")
            save_summary_plot(f"TailEfficiency_DNNgt{suffix}_vs_mPhi")

        # Dedicated tt(2l) tail plot: directly answers whether the new
        # background balancing has damaged the dominant-background rejection.
        tt2l_tail = tail_table[tail_table["group"] == "tt2l"].sort_values("mPhi")
        if len(tt2l_tail) > 0:
            plt.figure(figsize=(9, 6))
            for threshold in TAIL_THRESHOLDS:
                suffix = str(threshold).replace(".", "p")
                plt.plot(
                    tt2l_tail["mPhi"],
                    tt2l_tail[f"eff_gt_{suffix}"],
                    marker="o",
                    linewidth=2,
                    label=rf"DNN $>{threshold:.2f}$",
                )
            plt.xlabel(r"$m_\Phi$ [GeV]")
            plt.ylabel(r"$t\bar{t}(2\ell)$ efficiency")
            plt.yscale("log")
            plt.grid(True, which="both", alpha=0.3)
            plt.legend()
            plt.title(r"Protection of the $t\bar{t}(2\ell)$ rejection")
            save_summary_plot("TT2L_TailEfficiency_vs_mPhi")

        # TTNuNu tail plot: shows why it remains a separate classification
        # problem even after visible ttV is included in the main DNN.
        ttnunu_tail = tail_table[tail_table["group"] == "TTNuNu"].sort_values("mPhi")
        if len(ttnunu_tail) > 0:
            plt.figure(figsize=(9, 6))
            for threshold in TAIL_THRESHOLDS:
                suffix = str(threshold).replace(".", "p")
                plt.plot(
                    ttnunu_tail["mPhi"],
                    ttnunu_tail[f"eff_gt_{suffix}"],
                    marker="o",
                    linewidth=2,
                    label=rf"DNN $>{threshold:.2f}$",
                )
            plt.xlabel(r"$m_\Phi$ [GeV]")
            plt.ylabel("TTNuNu efficiency")
            plt.ylim(0.0, 1.0)
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.title("TTNuNu survival in the main-DNN tail")
            save_summary_plot("TTNuNu_TailEfficiency_vs_mPhi")

        # High-tail ranking at representative masses.
        tightest_threshold = max(TAIL_THRESHOLDS)
        tightest_suffix = str(tightest_threshold).replace(".", "p")
        tightest_eff_column = f"eff_gt_{tightest_suffix}"
        available_tail_masses = sorted(tail_table["mPhi"].unique())
        used_masses = set()

        for requested_mass in SUMMARY_MASSES:
            mass = nearest_available_mass(requested_mass, available_tail_masses)
            if mass is None or mass in used_masses:
                continue
            used_masses.add(mass)

            ranking = (
                tail_table[np.isclose(tail_table["mPhi"], mass)]
                .set_index("group")[tightest_eff_column]
                .dropna()
                .sort_values(ascending=True)
            )

            plt.figure(figsize=(8, max(5, 0.55 * len(ranking))))
            positions = np.arange(len(ranking))
            plt.barh(positions, ranking.values)
            plt.yticks(positions, ranking.index)
            plt.xlabel(rf"Efficiency for DNN $>{tightest_threshold:.2f}$")
            plt.xscale("log")
            plt.grid(True, axis="x", which="both", alpha=0.3)
            plt.title(rf"High-score background hierarchy at $m_\Phi={mass:g}$ GeV")
            for pos, value in zip(positions, ranking.values):
                plt.text(value * 1.08, pos, f"{100.0 * value:.2f}%", va="center")
            save_summary_plot(f"HighTailHierarchy_mPhi{mass:g}")

        # Heatmap of the tightest tail efficiency.
        tail_pivot = tail_table.pivot(
            index="mPhi",
            columns="group",
            values=tightest_eff_column,
        ).sort_index()
        tail_heatmap_groups = ordered_existing_groups(tail_pivot.columns)
        tail_pivot = tail_pivot[tail_heatmap_groups]
        tail_heatmap_values = tail_pivot.T.to_numpy(dtype=float)

        plt.figure(figsize=(13, 6))
        image = plt.imshow(
            np.log10(np.clip(tail_heatmap_values, 1e-6, None)),
            aspect="auto",
            origin="lower",
        )
        colorbar = plt.colorbar(image)
        colorbar.set_label(rf"$\log_{{10}}$ efficiency for DNN $>{tightest_threshold:.2f}$")
        plt.xticks(
            np.arange(len(tail_pivot.index)),
            [f"{mass:g}" for mass in tail_pivot.index],
            rotation=45,
            ha="right",
        )
        plt.yticks(np.arange(len(tail_heatmap_groups)), tail_heatmap_groups)
        plt.xlabel(r"$m_\Phi$ [GeV]")
        plt.ylabel("Background")
        plt.title(rf"High-score survival matrix: DNN $>{tightest_threshold:.2f}$")

        for row_idx in range(tail_heatmap_values.shape[0]):
            for col_idx in range(tail_heatmap_values.shape[1]):
                value = tail_heatmap_values[row_idx, col_idx]
                if np.isfinite(value):
                    plt.text(
                        col_idx,
                        row_idx,
                        f"{100.0 * value:.2f}%",
                        ha="center",
                        va="center",
                        fontsize=7,
                    )
        save_summary_plot(f"TailEfficiency_DNNgt{tightest_suffix}_heatmap")

    print("Summary diagnostic plots saved in ./Plots")

# ============================================================
# FEATURE IMPORTANCE
# ============================================================

print("\nPlotting feature importance")

X_val = np.array(X_test_scaled)
y_val = np.array(Y_test["isSignal"])
features = var

print("Computing permutation importance")
baseline_auc = roc_auc_score(
    y_val,
    model.predict(
        X_val,
        batch_size=PREDICT_BATCH_SIZE,
        verbose=0,
    ),
)

importances = []
rng_importance = np.random.default_rng(SEED + 500)

for i in range(X_val.shape[1]):
    X_permuted = X_val.copy()
    rng_importance.shuffle(X_permuted[:, i])

    permuted_auc = roc_auc_score(
        y_val,
        model.predict(
            X_permuted,
            batch_size=PREDICT_BATCH_SIZE,
            verbose=0,
        ),
    )

    importances.append(baseline_auc - permuted_auc)

importances = np.array(importances)
idx = np.argsort(importances)[::-1]

plt.clf()
plt.figure(figsize=(10, max(6, 0.25 * len(features))))
plt.barh(np.array(features)[idx], importances[idx])
plt.xlabel("Permutation Importance (Delta AUC)")
plt.title("Feature Importance - Keras DNN")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(
    f"./Plots/FeatureImportance_Keras_Permutation_{ANALYSIS_NAME}.png",
    dpi=600,
)
plt.close()

importance_df = pd.DataFrame({
    "feature": np.array(features)[idx],
    "delta_auc": importances[idx],
})

importance_df.to_csv(
    f"./Plots/FeatureImportance_{ANALYSIS_NAME}.csv",
    index=False,
)

print("Permutation importance saved.")

# ============================================================
# OPTIONAL SHAP
# ============================================================

try:
    import shap
    HAS_SHAP = True
except ImportError:
    HAS_SHAP = False
    print("SHAP not installed — skipping SHAP plots.")

if HAS_SHAP:
    print("Computing SHAP feature importance")

    max_background = 200

    if len(X_val) > max_background:
        idx_sample = np.random.choice(len(X_val), max_background, replace=False)
        X_background = X_val[idx_sample]
    else:
        X_background = X_val

    explainer = shap.Explainer(model, X_background)
    shap_values = explainer(X_background)

    plt.clf()
    shap.summary_plot(
        shap_values,
        X_background,
        feature_names=features,
        show=False,
    )
    plt.title("Feature Importance - SHAP Values", fontsize=16)
    plt.tight_layout()
    plt.savefig(
        f"./Plots/SHAP_Feature_importance_Beeswarm_{ANALYSIS_NAME}.png",
        dpi=120,
    )
    plt.close()

    plt.clf()
    shap.summary_plot(
        shap_values,
        X_background,
        plot_type="bar",
        feature_names=features,
        show=False,
    )
    plt.title("Feature Importance - SHAP Values", fontsize=16)
    plt.tight_layout()
    plt.savefig(
        f"./Plots/SHAP_Feature_importance_Bar_{ANALYSIS_NAME}.png",
        dpi=120,
    )
    plt.close()

    print("SHAP plots saved.")

# ============================================================
# CORRELATION MATRIX HELPERS
# ============================================================


def plot_correlation_matrix(data, feature_names, outname, title):
    matrix = np.corrcoef(data, rowvar=False)

    fig, ax = plt.subplots(
        figsize=(
            max(12, 0.45 * len(feature_names)),
            max(12, 0.45 * len(feature_names)),
        )
    )
    im = ax.matshow(abs(matrix))
    plt.colorbar(im)

    for (i, j), z in np.ndenumerate(matrix):
        if len(feature_names) <= 45:
            ax.text(
                j,
                i,
                "{:0.1f}".format(abs(z)),
                ha="center",
                va="center",
                fontsize=6,
            )

    ax.set_xticks(np.arange(len(feature_names)))
    ax.set_yticks(np.arange(len(feature_names)))
    ax.set_xticklabels(feature_names)
    ax.set_yticklabels(feature_names)

    plt.setp(
        ax.get_xticklabels(),
        rotation=-45,
        ha="right",
        rotation_mode="anchor",
    )
    ax.set_title(title)

    fig.tight_layout()
    plt.savefig(outname, dpi=120)
    plt.close()


# ============================================================
# CORRELATION MATRICES
# ============================================================

print("\nPlotting correlation matrix - all training")
plot_correlation_matrix(
    X_train_scaled,
    var,
    f"./Plots/CorrelationMatrix_{ANALYSIS_NAME}.png",
    "Correlation matrix - all",
)

print("\nPlotting correlation matrix - signal")
Sig_corr = build_model_input_frame(
    Sig,
    rng=np.random.default_rng(SEED + 600),
)
plot_correlation_matrix(
    Sig_corr[var].values,
    var,
    f"./Plots/CorrelationMatrix_{ANALYSIS_NAME}_Sig.png",
    "Correlation matrix - signal",
)

print("\nPlotting correlation matrix - background")
Bkg_corr = build_model_input_frame(
    Bkg,
    rng=np.random.default_rng(SEED + 700),
)
plot_correlation_matrix(
    Bkg_corr[var].values,
    var,
    f"./Plots/CorrelationMatrix_{ANALYSIS_NAME}_Bkg.png",
    "Correlation matrix - background",
)

print("\nEND OF THE ANALYSIS")
