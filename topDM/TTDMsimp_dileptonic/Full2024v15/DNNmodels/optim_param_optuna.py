#!/usr/bin/env python3
"""
optim_pnn_architecture_optuna.py

Standalone Optuna optimization of the MAIN ttDM parametric DNN.

Designed for the current Run-3 dileptonic ttDM workflow:
  * all available analysis variables;
  * one common train / validation / untouched test split;
  * 50/50 signal-background mini-batches;
  * equal signal-mass representation in every batch;
  * analysis-aware background-group fractions;
  * background mPhi dynamically sampled from the real signal hypotheses;
  * per-mass validation AUC as the optimization target;
  * Optuna pruning and resumable SQLite storage;
  * direct search over pNN conditioning mechanisms inspired by the papers:
        concat, biasing, scaling, affine
    and conditioning positions:
        begin, all, end
  * architecture, optimization and regularization hyperparameters;
  * optional top-K multi-seed confirmation;
  * optional final training and untouched-test evaluation.

The default task includes TTNuNu directly in the MAIN classifier, because the
current analysis discussion favours one classifier over a two-DNN strategy.
Use --exclude-ttnunu to optimize the previous main task instead.

Typical usage
-------------

Broad scan, resumable:
    python3 optim_pnn_architecture_optuna.py \
        --trials 300 \
        --study-name ttdm_main_pnn_v1

Resume the same study:
    python3 optim_pnn_architecture_optuna.py \
        --trials 200 \
        --study-name ttdm_main_pnn_v1

Confirm the top 5 configurations with 3 seeds:
    python3 optim_pnn_architecture_optuna.py \
        --trials 0 \
        --study-name ttdm_main_pnn_v1 \
        --confirm-top-k 5 \
        --confirm-seeds 3

Train and save the confirmed best configuration, then evaluate the untouched
test set:
    python3 optim_pnn_architecture_optuna.py \
        --trials 0 \
        --study-name ttdm_main_pnn_v1 \
        --confirm-top-k 5 \
        --confirm-seeds 3 \
        --train-final

Important
---------
Do not optimize class balance, signal-mass balance, or arbitrary background
fractions in the same architecture scan. Those encode the physics problem.
The scan keeps them fixed and optimizes the model that solves that problem.
"""

from __future__ import annotations

import argparse
import gc
import glob
import json
import math
import os
import random
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import ROOT
import joblib
import matplotlib.pyplot as plt
import numpy as np
import optuna
import pandas as pd
import tensorflow as tf

from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from tensorflow.keras import Input
from tensorflow.keras import callbacks, initializers, regularizers
from tensorflow.keras.layers import (
    Activation,
    BatchNormalization,
    Concatenate,
    Dense,
    Dropout,
    GaussianNoise,
    Layer,
    LayerNormalization,
)


# ============================================================
# FIXED PHYSICS CONFIGURATION
# ============================================================

DEFAULT_SEED = 12345

SNAPSHOT_DIR = (
    "/afs/cern.ch/user/v/victorr/private/"
    "PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/"
    "Full2024v15/DNNmodels/files_for_training"
)

SIGNAL_PATTERN = "TTto2LDMsimpSpin0_ps_"

TRAIN_FRACTION = 0.70
VALIDATION_FRACTION = 0.15
TEST_FRACTION = 0.15

if not np.isclose(
    TRAIN_FRACTION + VALIDATION_FRACTION + TEST_FRACTION,
    1.0,
):
    raise ValueError("Train/validation/test fractions must sum to one.")

PREDICT_BATCH_SIZE = 4096
HIGH_MASS_MIN = 400.0

# The scan objective is the arithmetic mean of the per-mass validation AUCs.
# Keep these secondary metrics for robustness diagnostics.
OBJECTIVE_METRIC = "mean_mass_auc"

# Avoid models that become unnecessarily large for tabular HEP data.
MAX_MODEL_PARAMETERS = 1_500_000

# Optuna scan uses a reduced but fixed validation subset for speed.
DEFAULT_SCAN_SIGNAL_PER_MASS = 500
DEFAULT_SCAN_BACKGROUND_TOTAL = 5000

# One scan epoch can be a fraction of a full signal exposure.
DEFAULT_SCAN_STEP_FRACTION = 0.50

# Final confirmation uses the full validation set and full signal exposure.
CONFIRM_STEP_FRACTION = 1.00


# ============================================================
# FULL VARIABLE LIST
# ============================================================

OLD_VARIABLES = [
    "lep_pt1",
    "lep_pt2",
    "lep_eta1",
    "lep_eta2",

    "mll",
    "ptll",
    "drll",
    "detall",
    "dphill",
    "yll",

    "PuppiMET_pt",
    "PuppiMET_phi",
    "dphilmet",
    "dphilmet1",
    "dphilmet2",
    "dphillmet",

    "mtw1",
    "mtw2",
    "mth",
    "mTi",
    "mR",
    "mT2",
    "mTe",

    "recoil",
    "upara",
    "uperp",
    "pTWW",

    "mcoll",
    "mcollWW",
    "choiMass",

    "nbjet_jet_ratio",
    "njet",
    "ht",
    "vht_pt",
    "dphijet1met",
    "dphijet2met",
    "dphijjmet",

    "chel",
    "pdark",
    "dphi_ttbar",
    "dphi_met_llb",
]

TOPDM_EXTRA_VARIABLES = [
    "mt2blbl",

    "dphi_met_ll",
    "st",
    "met_over_sqrt_ht",
    "met_over_st",
    "dphi_min_j_met",

    "pt_llb",
    "dphi_met_llb_safe",

    "pt_llbb",
    "dphi_met_llbb",
    "m_llbb",
    "mT_llbb_met",
    "met_over_pt_llbb",

    "mt2_bell_l",
    "mbl_min",
    "mbl_max",
    "mtbl",

    "mbb",
    "drbb",
    "ptbb",

    "max_nonleading_btag",
    "pt_b2",

    "nForwardJet",
    "leadingForwardJet_pt",
    "leadingForwardJet_absEta",
    "deta_forwardJet_b",
    "dphi_forwardJet_met",

    "top1_pt_reco",
    "top2_pt_reco",
    "pdark_over_met",
]

RESTFRAME_VARIABLES = [
    "angle_ll_llbb_rf",
    "dphi_ll_llbb_rf",
    "cos_l1_llbb_rf",
    "cos_l2_llbb_rf",

    "angle_ll_llmet_rf",
    "dphi_ll_llmet_rf",
    "cos_l1_llmet_rf",
    "cos_l2_llmet_rf",
]

CANDIDATE_PHYSICS_VARIABLES = list(
    dict.fromkeys(
        OLD_VARIABLES
        + TOPDM_EXTRA_VARIABLES
        + RESTFRAME_VARIABLES
    )
)


# ============================================================
# SAMPLE DEFINITIONS
# ============================================================

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


# Previous main training.
BACKGROUND_FRACTIONS_WITHOUT_TTNUNU = {
    "tt2l": 0.60,
    "other_tt": 0.10,
    "single_top": 0.15,
    "other_ttZ": 0.05,
    "ttW": 0.05,
    "ttH": 0.05,
}

# Single-main-DNN candidate including TTNuNu.
# TTNuNu is made visible to the classifier without allowing it to dominate the
# training objective.
BACKGROUND_FRACTIONS_WITH_TTNUNU = {
    "tt2l": 0.50,
    "other_tt": 0.08,
    "single_top": 0.12,
    "other_ttZ": 0.05,
    "ttW": 0.05,
    "ttH": 0.05,
    "TTNuNu": 0.15,
}


# ============================================================
# DATA STRUCTURES
# ============================================================

@dataclass
class SplitFrames:
    train: pd.DataFrame
    validation: pd.DataFrame
    test: pd.DataFrame


@dataclass
class DatasetBundle:
    physics_features: list[str]
    signal_hypotheses: list[float]

    train_df: pd.DataFrame
    validation_df: pd.DataFrame
    test_df: pd.DataFrame

    physics_scaler: StandardScaler

    train_physics_scaled: np.ndarray
    validation_physics_scaled: np.ndarray
    test_physics_scaled: np.ndarray

    scan_validation_df: pd.DataFrame
    scan_validation_physics_scaled: np.ndarray

    background_fractions: dict[str, float]


@dataclass
class TrainingResult:
    model: tf.keras.Model
    history: pd.DataFrame
    best_mean_auc: float
    best_min_auc: float
    best_median_auc: float
    best_high_mass_auc: float
    best_epoch: int
    parameter_count: int


# ============================================================
# CLI
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Optuna optimization of vanilla and conditioned parametric DNNs "
            "for dileptonic ttDM."
        )
    )

    parser.add_argument(
        "--snapshot-dir",
        default=SNAPSHOT_DIR,
    )

    parser.add_argument(
        "--outdir",
        default="pnn_optuna_summary",
    )

    parser.add_argument(
        "--study-name",
        default="ttdm_main_pnn_architecture",
    )

    parser.add_argument(
        "--storage",
        default="pnn_optuna.db",
        help=(
            "SQLite file or full Optuna storage URL. "
            "The study is resumable."
        ),
    )

    parser.add_argument(
        "--trials",
        type=int,
        default=300,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
    )

    parser.add_argument(
        "--scan-epochs",
        type=int,
        default=80,
    )

    parser.add_argument(
        "--scan-step-fraction",
        type=float,
        default=DEFAULT_SCAN_STEP_FRACTION,
    )

    parser.add_argument(
        "--scan-signal-per-mass",
        type=int,
        default=DEFAULT_SCAN_SIGNAL_PER_MASS,
    )

    parser.add_argument(
        "--scan-background-total",
        type=int,
        default=DEFAULT_SCAN_BACKGROUND_TOTAL,
    )

    parser.add_argument(
        "--confirm-top-k",
        type=int,
        default=0,
    )

    parser.add_argument(
        "--confirm-seeds",
        type=int,
        default=3,
    )

    parser.add_argument(
        "--confirm-epochs",
        type=int,
        default=250,
    )

    parser.add_argument(
        "--train-final",
        action="store_true",
    )

    parser.add_argument(
        "--final-epochs",
        type=int,
        default=300,
    )

    parser.add_argument(
        "--include-ttnunu",
        dest="include_ttnunu",
        action="store_true",
        help="Include TTNuNu in the main classifier.",
    )

    parser.add_argument(
        "--exclude-ttnunu",
        dest="include_ttnunu",
        action="store_false",
        help="Use the previous main task without TTNuNu.",
    )

    parser.set_defaults(
        include_ttnunu=True
    )

    parser.add_argument(
        "--strict-feature-check",
        action="store_true",
    )

    parser.add_argument(
        "--enqueue-baselines",
        dest="enqueue_baselines",
        action="store_true",
    )

    parser.add_argument(
        "--no-enqueue-baselines",
        dest="enqueue_baselines",
        action="store_false",
    )

    parser.set_defaults(
        enqueue_baselines=True
    )

    return parser.parse_args()


# ============================================================
# REPRODUCIBILITY / STORAGE
# ============================================================

def set_seed(seed: int):
    np.random.seed(seed)
    random.seed(seed)
    tf.keras.utils.set_random_seed(seed)


def resolve_storage_url(storage: str) -> str:
    if "://" in storage:
        return storage

    return (
        "sqlite:///"
        + str(
            Path(storage)
            .expanduser()
            .resolve()
        )
    )


# ============================================================
# FILE DISCOVERY
# ============================================================

def matches_any(
    path: str,
    patterns: Sequence[str],
) -> bool:
    name = os.path.basename(path).lower()

    return any(
        pattern.lower() in name
        for pattern in patterns
    )


def extract_sample_tag(
    path: str,
) -> str:
    basename = os.path.basename(path)

    match = re.search(
        r"nanoLatino_(.*)__part\d+_snapshot",
        basename,
    )

    if match:
        return match.group(1)

    return (
        basename
        .replace("nanoLatino_", "")
        .replace("_snapshot.root", "")
    )


def extract_mphi(
    sample_tag: str,
) -> Optional[float]:
    match = re.search(
        r"m[Pp]hi[_-]?(\d+)",
        sample_tag,
    )

    return (
        float(match.group(1))
        if match
        else None
    )


def discover_files(
    snapshot_dir: str,
    include_ttnunu: bool,
) -> dict:
    files = sorted(
        glob.glob(
            os.path.join(
                snapshot_dir,
                "*.root",
            )
        )
    )

    if len(files) == 0:
        raise RuntimeError(
            f"No ROOT snapshots found in {snapshot_dir}"
        )

    signal_files = [
        path
        for path in files
        if SIGNAL_PATTERN in path
    ]

    ttnunu_files = [
        path
        for path in files
        if matches_any(
            path,
            TTNUNU_PATTERNS,
        )
        and path not in signal_files
    ]

    tt2l_files = [
        path
        for path in files
        if matches_any(
            path,
            TT2L_PATTERNS,
        )
        and path not in signal_files
        and path not in ttnunu_files
    ]

    other_tt_files = [
        path
        for path in files
        if matches_any(
            path,
            OTHER_TT_PATTERNS,
        )
        and path not in signal_files
        and path not in ttnunu_files
    ]

    single_top_files = [
        path
        for path in files
        if matches_any(
            path,
            SINGLE_TOP_PATTERNS,
        )
        and path not in signal_files
        and path not in ttnunu_files
    ]

    other_ttz_files = [
        path
        for path in files
        if matches_any(
            path,
            OTHER_TTZ_PATTERNS,
        )
        and path not in signal_files
        and path not in ttnunu_files
    ]

    ttw_files = [
        path
        for path in files
        if matches_any(
            path,
            TTW_PATTERNS,
        )
        and path not in signal_files
        and path not in ttnunu_files
    ]

    tth_files = [
        path
        for path in files
        if matches_any(
            path,
            TTH_PATTERNS,
        )
        and path not in signal_files
        and path not in ttnunu_files
    ]

    signal_groups: dict[str, list[str]] = {}

    for path in signal_files:
        tag = extract_sample_tag(path)

        signal_groups.setdefault(
            tag,
            [],
        ).append(path)

    background_files = {
        "tt2l": tt2l_files,
        "other_tt": other_tt_files,
        "single_top": single_top_files,
        "other_ttZ": other_ttz_files,
        "ttW": ttw_files,
        "ttH": tth_files,
    }

    if include_ttnunu:
        background_files[
            "TTNuNu"
        ] = ttnunu_files

    background_files = {
        group: group_files
        for group, group_files in (
            background_files.items()
        )
        if len(group_files) > 0
    }

    if len(signal_groups) == 0:
        raise RuntimeError(
            "No signal samples discovered."
        )

    if len(background_files) == 0:
        raise RuntimeError(
            "No background samples discovered."
        )

    print("")
    print("=" * 88)
    print("SNAPSHOT DISCOVERY")
    print("=" * 88)
    print(
        f"Signal files: {len(signal_files)}"
    )

    for group, group_files in (
        background_files.items()
    ):
        print(
            f"{group:15s}: "
            f"{len(group_files)} files"
        )

    return {
        "signal_groups": signal_groups,
        "background_files": (
            background_files
        ),
    }


# ============================================================
# FEATURE CHECKING
# ============================================================

def branches_in_file(
    path: str,
) -> set[str]:
    root_file = ROOT.TFile.Open(path)

    if (
        not root_file
        or root_file.IsZombie()
    ):
        raise RuntimeError(
            f"Could not open ROOT file: {path}"
        )

    tree = root_file.Get(
        "Events"
    )

    if tree is None:
        root_file.Close()

        raise RuntimeError(
            f"No Events tree in {path}"
        )

    branches = {
        branch.GetName()
        for branch in (
            tree.GetListOfBranches()
        )
    }

    root_file.Close()

    return branches


def resolve_physics_features(
    discovered: dict,
    strict: bool,
) -> list[str]:
    representative_files = []

    for _, files in sorted(
        discovered[
            "signal_groups"
        ].items()
    ):
        representative_files.append(
            files[0]
        )

    for _, files in sorted(
        discovered[
            "background_files"
        ].items()
    ):
        representative_files.append(
            files[0]
        )

    common_branches = None

    for path in (
        representative_files
    ):
        branches = branches_in_file(
            path
        )

        if common_branches is None:
            common_branches = (
                branches
            )
        else:
            common_branches &= (
                branches
            )

    missing = [
        feature
        for feature in (
            CANDIDATE_PHYSICS_VARIABLES
        )
        if feature not in (
            common_branches
        )
    ]

    available = [
        feature
        for feature in (
            CANDIDATE_PHYSICS_VARIABLES
        )
        if feature in (
            common_branches
        )
    ]

    print("")
    print("=" * 88)
    print("FEATURE CHECK")
    print("=" * 88)

    if missing:
        print(
            "Variables missing from at least one "
            "representative sample:"
        )

        for feature in missing:
            print(
                "  ",
                feature,
            )

        if strict:
            raise RuntimeError(
                "Missing features and "
                "--strict-feature-check was set."
            )

    print("")
    print(
        f"Using {len(available)} physics features."
    )

    for feature in available:
        print(
            "  ",
            feature,
        )

    return available


# ============================================================
# ROOT -> PANDAS
# ============================================================

def dataframe_from_files(
    files: Sequence[str],
    columns: Sequence[str],
    description: str,
) -> pd.DataFrame:
    chain = ROOT.TChain(
        "Events"
    )

    for path in files:
        chain.Add(path)

    n_entries = int(
        chain.GetEntries()
    )

    print(
        f"Reading {description}: "
        f"{len(files)} files, "
        f"{n_entries} entries"
    )

    frame = pd.DataFrame(
        ROOT.RDataFrame(
            chain
        ).AsNumpy(
            list(columns)
        )
    )

    for column in columns:
        frame[column] = (
            frame[column]
            .astype(
                np.float32,
                copy=False,
            )
        )

    frame[
        list(columns)
    ] = (
        frame[
            list(columns)
        ]
        .replace(
            [np.inf, -np.inf],
            np.nan,
        )
    )

    before = len(frame)

    frame.dropna(
        subset=list(columns),
        inplace=True,
    )

    frame.reset_index(
        drop=True,
        inplace=True,
    )

    print(
        f"  retained {len(frame)} / "
        f"{before} after numerical cleaning"
    )

    return frame


def load_data(
    discovered: dict,
    physics_features: Sequence[str],
) -> tuple[
    pd.DataFrame,
    dict[str, pd.DataFrame],
    list[float],
]:
    signal_frames = []

    for tag, files in sorted(
        discovered[
            "signal_groups"
        ].items()
    ):
        mass = extract_mphi(
            tag
        )

        if mass is None:
            raise RuntimeError(
                f"Could not infer mPhi from {tag}"
            )

        frame = (
            dataframe_from_files(
                files,
                physics_features,
                f"signal {tag}",
            )
        )

        frame["isSignal"] = 1
        frame["mPhi_true"] = (
            float(mass)
        )
        frame["bkg_group"] = (
            "signal"
        )

        signal_frames.append(
            frame
        )

    signal = pd.concat(
        signal_frames,
        ignore_index=True,
    )

    backgrounds = {}

    for group, files in sorted(
        discovered[
            "background_files"
        ].items()
    ):
        frame = (
            dataframe_from_files(
                files,
                physics_features,
                f"background {group}",
            )
        )

        frame["isSignal"] = 0
        frame["mPhi_true"] = (
            np.nan
        )
        frame["bkg_group"] = (
            group
        )

        backgrounds[
            group
        ] = frame

    signal_hypotheses = sorted(
        signal[
            "mPhi_true"
        ]
        .unique()
        .astype(float)
        .tolist()
    )

    print("")
    print("Signal entries by mass:")
    print(
        signal
        .groupby(
            "mPhi_true"
        )
        .size()
    )

    print("")
    print("Background entries:")

    for group, frame in (
        backgrounds.items()
    ):
        print(
            f"  {group:15s}: "
            f"{len(frame)}"
        )

    return (
        signal,
        backgrounds,
        signal_hypotheses,
    )


# ============================================================
# COMMON SPLIT
# ============================================================

def split_three_way(
    frame: pd.DataFrame,
    stratify_labels: pd.Series,
    seed: int,
) -> SplitFrames:
    temp_fraction = (
        VALIDATION_FRACTION
        + TEST_FRACTION
    )

    train, temp = (
        train_test_split(
            frame,
            test_size=(
                temp_fraction
            ),
            random_state=seed,
            stratify=(
                stratify_labels
            ),
        )
    )

    relative_test_fraction = (
        TEST_FRACTION
        / temp_fraction
    )

    temp_labels = (
        stratify_labels.loc[
            temp.index
        ]
    )

    validation, test = (
        train_test_split(
            temp,
            test_size=(
                relative_test_fraction
            ),
            random_state=(
                seed + 1
            ),
            stratify=(
                temp_labels
            ),
        )
    )

    return SplitFrames(
        train=(
            train
            .reset_index(
                drop=True
            )
        ),
        validation=(
            validation
            .reset_index(
                drop=True
            )
        ),
        test=(
            test
            .reset_index(
                drop=True
            )
        ),
    )


def make_common_splits(
    signal: pd.DataFrame,
    backgrounds: dict[str, pd.DataFrame],
    seed: int,
) -> tuple[
    SplitFrames,
    dict[str, SplitFrames],
]:
    signal_strata = (
        "sig_"
        + signal[
            "mPhi_true"
        ]
        .astype(int)
        .astype(str)
    )

    signal_split = (
        split_three_way(
            signal,
            signal_strata,
            seed,
        )
    )

    background_splits = {}

    for index, (
        group,
        frame,
    ) in enumerate(
        sorted(
            backgrounds.items()
        )
    ):
        labels = pd.Series(
            [group] * len(frame),
            index=frame.index,
        )

        background_splits[
            group
        ] = split_three_way(
            frame,
            labels,
            seed + 100 + index,
        )

    return (
        signal_split,
        background_splits,
    )


def combine_split(
    signal_split: SplitFrames,
    background_splits: dict[str, SplitFrames],
    split_name: str,
) -> pd.DataFrame:
    frames = [
        getattr(
            signal_split,
            split_name,
        )
    ]

    frames.extend(
        getattr(
            splits,
            split_name,
        )
        for splits in (
            background_splits.values()
        )
    )

    return pd.concat(
        frames,
        ignore_index=True,
    )


# ============================================================
# FIXED BACKGROUND FRACTIONS
# ============================================================

def normalize_background_fractions(
    requested: dict[str, float],
    available_groups: Sequence[str],
) -> dict[str, float]:
    available_groups = list(
        available_groups
    )

    fractions = {
        group: float(
            requested.get(
                group,
                0.0,
            )
        )
        for group in (
            available_groups
        )
    }

    total = sum(
        fractions.values()
    )

    if total <= 0:
        raise RuntimeError(
            "Background fractions sum to zero."
        )

    return {
        group: value / total
        for group, value in (
            fractions.items()
        )
    }


# ============================================================
# REDUCED SCAN VALIDATION SAMPLE
# ============================================================

def make_scan_validation_frame(
    validation_df: pd.DataFrame,
    signal_hypotheses: Sequence[float],
    signal_per_mass: int,
    background_total: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(
        seed
    )

    frames = []

    for mass in (
        signal_hypotheses
    ):
        mass_frame = (
            validation_df[
                (validation_df[
                    "isSignal"
                ] == 1)
                & np.isclose(
                    validation_df[
                        "mPhi_true"
                    ],
                    mass,
                )
            ]
        )

        if (
            signal_per_mass > 0
            and len(
                mass_frame
            ) > signal_per_mass
        ):
            chosen = rng.choice(
                len(mass_frame),
                size=signal_per_mass,
                replace=False,
            )

            mass_frame = (
                mass_frame.iloc[
                    chosen
                ]
            )

        frames.append(
            mass_frame
        )

    background = (
        validation_df[
            validation_df[
                "isSignal"
            ] == 0
        ]
    )

    if (
        background_total > 0
        and len(background)
        > background_total
    ):
        chosen = rng.choice(
            len(background),
            size=background_total,
            replace=False,
        )

        background = (
            background.iloc[
                chosen
            ]
        )

    frames.append(
        background
    )

    output = pd.concat(
        frames,
        ignore_index=True,
    )

    print("")
    print("=" * 88)
    print("SCAN VALIDATION SUBSET")
    print("=" * 88)
    print(
        output
        .groupby(
            [
                "isSignal",
                "bkg_group",
            ]
        )
        .size()
    )

    return output


# ============================================================
# DATASET BUNDLE / SCALING
# ============================================================

def build_dataset_bundle(
    signal: pd.DataFrame,
    backgrounds: dict[str, pd.DataFrame],
    signal_hypotheses: Sequence[float],
    physics_features: Sequence[str],
    background_fractions: dict[str, float],
    args,
) -> DatasetBundle:
    (
        signal_split,
        background_splits,
    ) = make_common_splits(
        signal,
        backgrounds,
        args.seed,
    )

    train_df = combine_split(
        signal_split,
        background_splits,
        "train",
    )

    validation_df = combine_split(
        signal_split,
        background_splits,
        "validation",
    )

    test_df = combine_split(
        signal_split,
        background_splits,
        "test",
    )

    physics_scaler = (
        StandardScaler()
    )

    physics_scaler.fit(
        train_df[
            list(
                physics_features
            )
        ]
    )

    train_physics_scaled = (
        physics_scaler.transform(
            train_df[
                list(
                    physics_features
                )
            ]
        )
        .astype(np.float32)
    )

    validation_physics_scaled = (
        physics_scaler.transform(
            validation_df[
                list(
                    physics_features
                )
            ]
        )
        .astype(np.float32)
    )

    test_physics_scaled = (
        physics_scaler.transform(
            test_df[
                list(
                    physics_features
                )
            ]
        )
        .astype(np.float32)
    )

    scan_validation_df = (
        make_scan_validation_frame(
            validation_df,
            signal_hypotheses,
            args.scan_signal_per_mass,
            args.scan_background_total,
            args.seed + 900,
        )
    )

    scan_validation_physics_scaled = (
        physics_scaler.transform(
            scan_validation_df[
                list(
                    physics_features
                )
            ]
        )
        .astype(np.float32)
    )

    print("")
    print("=" * 88)
    print("COMMON SPLIT")
    print("=" * 88)
    print(
        f"Train:      {len(train_df)}"
    )
    print(
        f"Validation: {len(validation_df)}"
    )
    print(
        f"Test:       {len(test_df)}"
    )

    return DatasetBundle(
        physics_features=list(
            physics_features
        ),
        signal_hypotheses=list(
            signal_hypotheses
        ),

        train_df=train_df,
        validation_df=(
            validation_df
        ),
        test_df=test_df,

        physics_scaler=(
            physics_scaler
        ),

        train_physics_scaled=(
            train_physics_scaled
        ),
        validation_physics_scaled=(
            validation_physics_scaled
        ),
        test_physics_scaled=(
            test_physics_scaled
        ),

        scan_validation_df=(
            scan_validation_df
        ),
        scan_validation_physics_scaled=(
            scan_validation_physics_scaled
        ),

        background_fractions=(
            background_fractions
        ),
    )


# ============================================================
# MASS ENCODING
# ============================================================

def transform_mass(
    values: np.ndarray,
    mode: str,
) -> np.ndarray:
    values = np.asarray(
        values,
        dtype=np.float32,
    )

    if mode == "raw":
        return values

    if mode == "log1p":
        return np.log1p(
            values
        ).astype(
            np.float32
        )

    raise ValueError(
        f"Unknown mass transform: {mode}"
    )


def mass_normalization_parameters(
    signal_hypotheses: Sequence[float],
    mode: str,
) -> tuple[float, float]:
    values = transform_mass(
        np.asarray(
            signal_hypotheses,
            dtype=np.float32,
        ),
        mode,
    )

    mean = float(
        np.mean(
            values
        )
    )

    std = float(
        np.std(
            values
        )
    )

    if std <= 0:
        raise RuntimeError(
            "Mass encoding has zero standard deviation."
        )

    return (
        mean,
        std,
    )


def encode_mass(
    values: np.ndarray,
    mode: str,
    signal_hypotheses: Sequence[float],
) -> np.ndarray:
    mean, std = (
        mass_normalization_parameters(
            signal_hypotheses,
            mode,
        )
    )

    transformed = (
        transform_mass(
            values,
            mode,
        )
    )

    return (
        (
            transformed
            - mean
        )
        / std
    ).astype(
        np.float32
    )


# ============================================================
# BALANCED GENERATOR
# ============================================================

def integer_partition(
    total: int,
    n_parts: int,
) -> list[int]:
    base = (
        total // n_parts
    )

    remainder = (
        total % n_parts
    )

    return [
        base
        + (
            index
            < remainder
        )
        for index in range(
            n_parts
        )
    ]


def counts_from_fractions(
    total: int,
    fractions: dict[str, float],
) -> dict[str, int]:
    groups = list(
        fractions.keys()
    )

    raw = np.asarray(
        [
            fractions[group]
            * total
            for group in groups
        ],
        dtype=float,
    )

    counts = np.floor(
        raw
    ).astype(int)

    remaining = (
        total
        - int(
            counts.sum()
        )
    )

    if remaining > 0:
        order = np.argsort(
            raw - counts
        )[::-1]

        for index in order[
            :remaining
        ]:
            counts[
                index
            ] += 1

    return {
        group: int(count)
        for group, count in zip(
            groups,
            counts,
        )
    }


def balanced_batch_generator(
    bundle: DatasetBundle,
    batch_size: int,
    mass_transform: str,
    seed: int,
):
    rng = np.random.default_rng(
        seed
    )

    signal_half = (
        batch_size // 2
    )

    background_half = (
        batch_size
        - signal_half
    )

    train_df = (
        bundle.train_df
    )

    signal_indices_by_mass = {}

    for mass in (
        bundle.signal_hypotheses
    ):
        indices = np.flatnonzero(
            (
                train_df[
                    "isSignal"
                ]
                .to_numpy(
                    dtype=int
                )
                == 1
            )
            & np.isclose(
                train_df[
                    "mPhi_true"
                ]
                .to_numpy(
                    dtype=float
                ),
                mass,
            )
        )

        if len(indices) == 0:
            raise RuntimeError(
                f"No signal training events for mPhi={mass:g}."
            )

        signal_indices_by_mass[
            mass
        ] = indices

    background_groups = sorted(
        train_df.loc[
            train_df[
                "isSignal"
            ] == 0,
            "bkg_group",
        ]
        .unique()
        .tolist()
    )

    background_indices_by_group = {}

    for group in (
        background_groups
    ):
        indices = np.flatnonzero(
            (
                train_df[
                    "isSignal"
                ]
                .to_numpy(
                    dtype=int
                )
                == 0
            )
            & (
                train_df[
                    "bkg_group"
                ]
                .to_numpy()
                == group
            )
        )

        background_indices_by_group[
            group
        ] = indices

    signal_counts = dict(
        zip(
            bundle.signal_hypotheses,
            integer_partition(
                signal_half,
                len(
                    bundle.signal_hypotheses
                ),
            ),
        )
    )

    active_fractions = (
        normalize_background_fractions(
            bundle.background_fractions,
            background_groups,
        )
    )

    background_counts = (
        counts_from_fractions(
            background_half,
            active_fractions,
        )
    )

    while True:
        selected_indices = []
        selected_masses = []

        for mass, count in (
            signal_counts.items()
        ):
            pool = (
                signal_indices_by_mass[
                    mass
                ]
            )

            chosen = rng.choice(
                pool,
                size=count,
                replace=True,
            )

            selected_indices.extend(
                chosen.tolist()
            )

            selected_masses.extend(
                [mass] * count
            )

        for group, count in (
            background_counts.items()
        ):
            pool = (
                background_indices_by_group[
                    group
                ]
            )

            chosen = rng.choice(
                pool,
                size=count,
                replace=True,
            )

            selected_indices.extend(
                chosen.tolist()
            )

            # Paper-adapted identical-sampled strategy:
            # every background draw receives one of the real signal
            # hypotheses, sampled dynamically.
            selected_masses.extend(
                rng.choice(
                    bundle.signal_hypotheses,
                    size=count,
                    replace=True,
                ).tolist()
            )

        selected_indices = np.asarray(
            selected_indices,
            dtype=int,
        )

        selected_masses = np.asarray(
            selected_masses,
            dtype=np.float32,
        )

        physics = (
            bundle.train_physics_scaled[
                selected_indices
            ]
        )

        mass_column = encode_mass(
            selected_masses,
            mass_transform,
            bundle.signal_hypotheses,
        ).reshape(
            -1,
            1,
        )

        X_batch = np.concatenate(
            [
                physics,
                mass_column,
            ],
            axis=1,
        ).astype(
            np.float32
        )

        y_batch = (
            bundle.train_df[
                "isSignal"
            ]
            .to_numpy(
                dtype=np.float32
            )[
                selected_indices
            ]
        )

        permutation = rng.permutation(
            len(
                y_batch
            )
        )

        yield (
            X_batch[
                permutation
            ],
            y_batch[
                permutation
            ],
        )


def steps_per_epoch(
    bundle: DatasetBundle,
    batch_size: int,
    fraction: float,
) -> int:
    n_signal = int(
        np.sum(
            bundle.train_df[
                "isSignal"
            ]
            .to_numpy(
                dtype=int
            )
            == 1
        )
    )

    signal_per_batch = (
        batch_size // 2
    )

    full_steps = int(
        math.ceil(
            n_signal
            / signal_per_batch
        )
    )

    return max(
        1,
        int(
            math.ceil(
                full_steps
                * fraction
            )
        ),
    )


# ============================================================
# CONDITIONING LAYERS
# ============================================================

class ConditionalBias(
    Layer
):
    def __init__(
        self,
        units: int,
        l2_value: float = 0.0,
        identity_init: bool = True,
        **kwargs,
    ):
        super().__init__(
            **kwargs
        )

        self.units = int(
            units
        )

        self.l2_value = float(
            l2_value
        )

        self.identity_init = bool(
            identity_init
        )

        zero_or_glorot = (
            "zeros"
            if identity_init
            else "glorot_uniform"
        )

        self.bias_projection = Dense(
            self.units,
            activation=None,
            kernel_initializer=(
                zero_or_glorot
            ),
            bias_initializer="zeros",
            kernel_regularizer=(
                regularizers.l2(
                    self.l2_value
                )
                if self.l2_value > 0
                else None
            ),
        )

    def call(
        self,
        inputs,
    ):
        hidden, mass = inputs

        return (
            hidden
            + self.bias_projection(
                mass
            )
        )

    def get_config(
        self,
    ):
        config = super().get_config()

        config.update(
            {
                "units": self.units,
                "l2_value": (
                    self.l2_value
                ),
                "identity_init": (
                    self.identity_init
                ),
            }
        )

        return config


class ConditionalScale(
    Layer
):
    def __init__(
        self,
        units: int,
        l2_value: float = 0.0,
        identity_init: bool = True,
        **kwargs,
    ):
        super().__init__(
            **kwargs
        )

        self.units = int(
            units
        )

        self.l2_value = float(
            l2_value
        )

        self.identity_init = bool(
            identity_init
        )

        if identity_init:
            kernel_initializer = (
                "zeros"
            )
        else:
            kernel_initializer = (
                "glorot_uniform"
            )

        self.scale_projection = Dense(
            self.units,
            activation=None,
            kernel_initializer=(
                kernel_initializer
            ),
            bias_initializer="zeros",
            kernel_regularizer=(
                regularizers.l2(
                    self.l2_value
                )
                if self.l2_value > 0
                else None
            ),
        )

    def call(
        self,
        inputs,
    ):
        hidden, mass = inputs

        raw_scale = (
            self.scale_projection(
                mass
            )
        )

        if self.identity_init:
            scale = (
                1.0
                + raw_scale
            )
        else:
            scale = (
                raw_scale
            )

        return (
            hidden
            * scale
        )

    def get_config(
        self,
    ):
        config = super().get_config()

        config.update(
            {
                "units": self.units,
                "l2_value": (
                    self.l2_value
                ),
                "identity_init": (
                    self.identity_init
                ),
            }
        )

        return config


class AffineConditioning(
    Layer
):
    def __init__(
        self,
        units: int,
        l2_value: float = 0.0,
        identity_init: bool = True,
        **kwargs,
    ):
        super().__init__(
            **kwargs
        )

        self.units = int(
            units
        )

        self.l2_value = float(
            l2_value
        )

        self.identity_init = bool(
            identity_init
        )

        if identity_init:
            scale_initializer = (
                "zeros"
            )
            bias_initializer = (
                "zeros"
            )
        else:
            scale_initializer = (
                "glorot_uniform"
            )
            bias_initializer = (
                "glorot_uniform"
            )

        kernel_regularizer = (
            regularizers.l2(
                self.l2_value
            )
            if self.l2_value > 0
            else None
        )

        self.scale_projection = Dense(
            self.units,
            activation=None,
            kernel_initializer=(
                scale_initializer
            ),
            bias_initializer="zeros",
            kernel_regularizer=(
                kernel_regularizer
            ),
        )

        self.bias_projection = Dense(
            self.units,
            activation=None,
            kernel_initializer=(
                bias_initializer
            ),
            bias_initializer="zeros",
            kernel_regularizer=(
                kernel_regularizer
            ),
        )

    def call(
        self,
        inputs,
    ):
        hidden, mass = inputs

        raw_scale = (
            self.scale_projection(
                mass
            )
        )

        bias = (
            self.bias_projection(
                mass
            )
        )

        if self.identity_init:
            scale = (
                1.0
                + raw_scale
            )
        else:
            scale = (
                raw_scale
            )

        return (
            hidden
            * scale
            + bias
        )

    def get_config(
        self,
    ):
        config = super().get_config()

        config.update(
            {
                "units": self.units,
                "l2_value": (
                    self.l2_value
                ),
                "identity_init": (
                    self.identity_init
                ),
            }
        )

        return config


# ============================================================
# SEARCH SPACE
# ============================================================

def available_optimizer_choices():
    choices = [
        "adam",
        "nadam",
        "rmsprop",
    ]

    if hasattr(
        tf.keras.optimizers,
        "AdamW",
    ):
        choices.append(
            "adamw"
        )

    return choices


def sample_trial_config(
    trial: optuna.Trial,
) -> dict:
    config = {}

    # ------------------------
    # Parametric conditioning
    # ------------------------

    config[
        "conditioning_type"
    ] = trial.suggest_categorical(
        "conditioning_type",
        [
            "concat",
            "bias",
            "scale",
            "affine",
        ],
    )

    config[
        "conditioning_position"
    ] = trial.suggest_categorical(
        "conditioning_position",
        [
            "begin",
            "all",
            "end",
        ],
    )

    config[
        "mass_transform"
    ] = trial.suggest_categorical(
        "mass_transform",
        [
            "raw",
            "log1p",
        ],
    )

    config[
        "mass_embedding_dim"
    ] = trial.suggest_categorical(
        "mass_embedding_dim",
        [
            1,
            8,
            16,
            32,
        ],
    )

    if (
        config[
            "mass_embedding_dim"
        ]
        > 1
    ):
        config[
            "mass_embedding_layers"
        ] = trial.suggest_int(
            "mass_embedding_layers",
            1,
            2,
        )
    else:
        config[
            "mass_embedding_layers"
        ] = 0

    if (
        config[
            "conditioning_type"
        ]
        in {
            "bias",
            "scale",
            "affine",
        }
    ):
        config[
            "conditioning_identity_init"
        ] = trial.suggest_categorical(
            "conditioning_identity_init",
            [
                True,
                False,
            ],
        )
    else:
        config[
            "conditioning_identity_init"
        ] = True

    # ------------------------
    # Hidden architecture
    # ------------------------

    config[
        "n_layers"
    ] = trial.suggest_int(
        "n_layers",
        2,
        6,
    )

    config[
        "base_units"
    ] = trial.suggest_categorical(
        "base_units",
        [
            64,
            96,
            128,
            192,
            256,
            384,
            512,
        ],
    )

    config[
        "width_decay"
    ] = trial.suggest_categorical(
        "width_decay",
        [
            1.0,
            1.25,
            1.5,
            2.0,
        ],
    )

    config[
        "min_units"
    ] = trial.suggest_categorical(
        "min_units",
        [
            16,
            32,
            64,
        ],
    )

    config[
        "activation"
    ] = trial.suggest_categorical(
        "activation",
        [
            "relu",
            "elu",
            "gelu",
            "swish",
        ],
    )

    config[
        "normalization"
    ] = trial.suggest_categorical(
        "normalization",
        [
            "none",
            "batch",
            "layer",
        ],
    )

    config[
        "initializer"
    ] = trial.suggest_categorical(
        "initializer",
        [
            "he_normal",
            "he_uniform",
            "glorot_uniform",
        ],
    )

    config[
        "dropout_first"
    ] = trial.suggest_float(
        "dropout_first",
        0.0,
        0.45,
        step=0.05,
    )

    config[
        "dropout_last"
    ] = trial.suggest_float(
        "dropout_last",
        0.0,
        0.30,
        step=0.05,
    )

    config[
        "use_gaussian_noise"
    ] = trial.suggest_categorical(
        "use_gaussian_noise",
        [
            False,
            True,
        ],
    )

    if config[
        "use_gaussian_noise"
    ]:
        config[
            "gaussian_noise_std"
        ] = trial.suggest_categorical(
            "gaussian_noise_std",
            [
                0.005,
                0.01,
                0.02,
                0.03,
                0.05,
            ],
        )
    else:
        config[
            "gaussian_noise_std"
        ] = 0.0

    # ------------------------
    # Optimizer / regularization
    # ------------------------

    config[
        "optimizer"
    ] = trial.suggest_categorical(
        "optimizer",
        available_optimizer_choices(),
    )

    config[
        "learning_rate"
    ] = trial.suggest_float(
        "learning_rate",
        1e-5,
        5e-3,
        log=True,
    )

    config[
        "clipnorm"
    ] = trial.suggest_categorical(
        "clipnorm",
        [
            0.0,
            1.0,
            5.0,
            10.0,
        ],
    )

    if (
        config[
            "optimizer"
        ]
        == "adamw"
    ):
        config[
            "weight_decay"
        ] = trial.suggest_float(
            "weight_decay",
            1e-7,
            1e-3,
            log=True,
        )

        config[
            "l2"
        ] = 0.0

    else:
        config[
            "weight_decay"
        ] = 0.0

        config[
            "use_l2"
        ] = trial.suggest_categorical(
            "use_l2",
            [
                False,
                True,
            ],
        )

        if config[
            "use_l2"
        ]:
            config[
                "l2"
            ] = trial.suggest_float(
                "l2",
                1e-7,
                1e-3,
                log=True,
            )
        else:
            config[
                "l2"
            ] = 0.0

    config[
        "label_smoothing"
    ] = trial.suggest_float(
        "label_smoothing",
        0.0,
        0.05,
        step=0.01,
    )

    config[
        "batch_size"
    ] = trial.suggest_categorical(
        "batch_size",
        [
            256,
            512,
            1024,
            2048,
        ],
    )

    # ------------------------
    # Learning-rate schedule / stopping
    # ------------------------

    config[
        "reduce_lr_factor"
    ] = trial.suggest_categorical(
        "reduce_lr_factor",
        [
            0.3,
            0.5,
            0.7,
        ],
    )

    config[
        "reduce_lr_patience"
    ] = trial.suggest_categorical(
        "reduce_lr_patience",
        [
            4,
            6,
            8,
            10,
        ],
    )

    config[
        "early_extra_patience"
    ] = trial.suggest_categorical(
        "early_extra_patience",
        [
            8,
            12,
            16,
            20,
        ],
    )

    config[
        "early_stopping_patience"
    ] = (
        config[
            "reduce_lr_patience"
        ]
        + config[
            "early_extra_patience"
        ]
    )

    return config


def layer_sizes_from_config(
    config: dict,
) -> list[int]:
    sizes = []

    for index in range(
        int(
            config[
                "n_layers"
            ]
        )
    ):
        value = int(
            round(
                config[
                    "base_units"
                ]
                / (
                    config[
                        "width_decay"
                    ]
                    ** index
                )
            )
        )

        value = max(
            value,
            int(
                config[
                    "min_units"
                ]
            ),
        )

        sizes.append(
            value
        )

    return sizes


def conditioning_layer_indices(
    position: str,
    n_layers: int,
) -> set[int]:
    if position == "begin":
        return {
            0
        }

    if position == "all":
        return set(
            range(
                n_layers
            )
        )

    if position == "end":
        return {
            n_layers - 1
        }

    raise ValueError(
        f"Unknown conditioning position: {position}"
    )


# ============================================================
# MODEL BUILDING
# ============================================================

def build_optimizer(
    config: dict,
):
    clipnorm = (
        None
        if config[
            "clipnorm"
        ] <= 0
        else float(
            config[
                "clipnorm"
            ]
        )
    )

    common = {
        "learning_rate": float(
            config[
                "learning_rate"
            ]
        ),
    }

    if clipnorm is not None:
        common[
            "clipnorm"
        ] = clipnorm

    name = config[
        "optimizer"
    ]

    if name == "adam":
        return (
            tf.keras.optimizers.Adam(
                **common
            )
        )

    if name == "nadam":
        return (
            tf.keras.optimizers.Nadam(
                **common
            )
        )

    if name == "rmsprop":
        return (
            tf.keras.optimizers.RMSprop(
                **common
            )
        )

    if name == "adamw":
        return (
            tf.keras.optimizers.AdamW(
                weight_decay=float(
                    config[
                        "weight_decay"
                    ]
                ),
                **common,
            )
        )

    raise ValueError(
        f"Unknown optimizer: {name}"
    )


def add_normalization(
    hidden,
    normalization: str,
    name: str,
):
    if normalization == "none":
        return hidden

    if normalization == "batch":
        return BatchNormalization(
            name=name,
        )(
            hidden
        )

    if normalization == "layer":
        return LayerNormalization(
            name=name,
        )(
            hidden
        )

    raise ValueError(
        f"Unknown normalization: {normalization}"
    )


def apply_conditioning(
    hidden,
    mass_representation,
    conditioning_type: str,
    units: int,
    l2_value: float,
    identity_init: bool,
    name: str,
):
    if conditioning_type == "bias":
        return ConditionalBias(
            units=units,
            l2_value=l2_value,
            identity_init=(
                identity_init
            ),
            name=name,
        )(
            [
                hidden,
                mass_representation,
            ]
        )

    if conditioning_type == "scale":
        return ConditionalScale(
            units=units,
            l2_value=l2_value,
            identity_init=(
                identity_init
            ),
            name=name,
        )(
            [
                hidden,
                mass_representation,
            ]
        )

    if conditioning_type == "affine":
        return AffineConditioning(
            units=units,
            l2_value=l2_value,
            identity_init=(
                identity_init
            ),
            name=name,
        )(
            [
                hidden,
                mass_representation,
            ]
        )

    raise ValueError(
        f"apply_conditioning called with "
        f"unsupported type={conditioning_type}"
    )


def build_model(
    n_physics_features: int,
    config: dict,
) -> tf.keras.Model:
    inputs = Input(
        shape=(
            n_physics_features
            + 1,
        ),
        name="inputs",
    )

    physics_inputs = (
        inputs[
            :,
            :-1
        ]
    )

    mass_input = (
        inputs[
            :,
            -1:
        ]
    )

    hidden = physics_inputs

    if (
        config[
            "gaussian_noise_std"
        ]
        > 0
    ):
        hidden = GaussianNoise(
            float(
                config[
                    "gaussian_noise_std"
                ]
            ),
            name="physics_noise",
        )(
            hidden
        )

    mass_representation = (
        mass_input
    )

    if (
        config[
            "mass_embedding_dim"
        ]
        > 1
    ):
        for index in range(
            config[
                "mass_embedding_layers"
            ]
        ):
            mass_representation = Dense(
                config[
                    "mass_embedding_dim"
                ],
                activation=(
                    config[
                        "activation"
                    ]
                ),
                name=(
                    f"mass_embedding_"
                    f"{index}"
                ),
            )(
                mass_representation
            )

    layer_sizes = (
        layer_sizes_from_config(
            config
        )
    )

    conditioned_layers = (
        conditioning_layer_indices(
            config[
                "conditioning_position"
            ],
            len(
                layer_sizes
            ),
        )
    )

    l2_value = float(
        config[
            "l2"
        ]
    )

    kernel_regularizer = (
        regularizers.l2(
            l2_value
        )
        if l2_value > 0
        else None
    )

    dropout_values = (
        np.linspace(
            float(
                config[
                    "dropout_first"
                ]
            ),
            float(
                config[
                    "dropout_last"
                ]
            ),
            len(
                layer_sizes
            ),
        )
    )

    for index, units in enumerate(
        layer_sizes
    ):
        should_condition = (
            index
            in conditioned_layers
        )

        if (
            config[
                "conditioning_type"
            ]
            == "concat"
            and should_condition
        ):
            hidden = Concatenate(
                name=(
                    f"concat_condition_"
                    f"{index}"
                )
            )(
                [
                    hidden,
                    mass_representation,
                ]
            )

        hidden = Dense(
            units,
            activation=None,
            kernel_initializer=(
                config[
                    "initializer"
                ]
            ),
            kernel_regularizer=(
                kernel_regularizer
            ),
            name=(
                f"dense_{index}"
            ),
        )(
            hidden
        )

        hidden = add_normalization(
            hidden,
            config[
                "normalization"
            ],
            name=(
                f"norm_{index}"
            ),
        )

        hidden = Activation(
            config[
                "activation"
            ],
            name=(
                f"activation_{index}"
            ),
        )(
            hidden
        )

        if (
            config[
                "conditioning_type"
            ]
            in {
                "bias",
                "scale",
                "affine",
            }
            and should_condition
        ):
            hidden = apply_conditioning(
                hidden=hidden,
                mass_representation=(
                    mass_representation
                ),
                conditioning_type=(
                    config[
                        "conditioning_type"
                    ]
                ),
                units=units,
                l2_value=l2_value,
                identity_init=(
                    config[
                        "conditioning_identity_init"
                    ]
                ),
                name=(
                    f"{config['conditioning_type']}"
                    f"_condition_{index}"
                ),
            )

        dropout = float(
            dropout_values[
                index
            ]
        )

        if dropout > 0:
            hidden = Dropout(
                dropout,
                name=(
                    f"dropout_{index}"
                ),
            )(
                hidden
            )

    output = Dense(
        1,
        activation="sigmoid",
        name="output",
    )(
        hidden
    )

    model = tf.keras.Model(
        inputs=inputs,
        outputs=output,
    )

    model.compile(
        loss=(
            tf.keras.losses.BinaryCrossentropy(
                label_smoothing=float(
                    config[
                        "label_smoothing"
                    ]
                )
            )
        ),
        optimizer=(
            build_optimizer(
                config
            )
        ),
        metrics=[
            tf.keras.metrics.AUC(
                name="auc"
            )
        ],
    )

    return model


# ============================================================
# PER-MASS VALIDATION / OPTUNA PRUNING
# ============================================================

class PerMassValidationCallback(
    callbacks.Callback
):
    def __init__(
        self,
        bundle: DatasetBundle,
        validation_df: pd.DataFrame,
        validation_physics_scaled: np.ndarray,
        mass_transform: str,
        trial: Optional[
            optuna.Trial
        ] = None,
        enable_pruning: bool = False,
    ):
        super().__init__()

        self.bundle = bundle

        self.validation_df = (
            validation_df
            .reset_index(
                drop=True
            )
        )

        self.validation_physics_scaled = (
            validation_physics_scaled
        )

        self.mass_transform = (
            mass_transform
        )

        self.trial = trial

        self.enable_pruning = (
            enable_pruning
        )

        labels = (
            self.validation_df[
                "isSignal"
            ]
            .to_numpy(
                dtype=int
            )
        )

        self.background_indices = (
            np.flatnonzero(
                labels == 0
            )
        )

        self.signal_indices_by_mass = {}

        true_masses = (
            self.validation_df[
                "mPhi_true"
            ]
            .to_numpy(
                dtype=float
            )
        )

        for mass in (
            self.bundle.signal_hypotheses
        ):
            self.signal_indices_by_mass[
                mass
            ] = np.flatnonzero(
                (labels == 1)
                & np.isclose(
                    true_masses,
                    mass,
                )
            )

        self.best_mean_auc = (
            -np.inf
        )

        self.best_min_auc = (
            np.nan
        )

        self.best_median_auc = (
            np.nan
        )

        self.best_high_mass_auc = (
            np.nan
        )

        self.best_epoch = 0

        self.rows = []

    def on_epoch_end(
        self,
        epoch,
        logs=None,
    ):
        logs = (
            logs
            if logs is not None
            else {}
        )

        auc_values = []
        high_mass_values = []

        row = {
            "epoch": (
                epoch + 1
            )
        }

        for mass in (
            self.bundle.signal_hypotheses
        ):
            indices = np.concatenate(
                [
                    self.signal_indices_by_mass[
                        mass
                    ],
                    self.background_indices,
                ]
            )

            labels = (
                self.validation_df[
                    "isSignal"
                ]
                .to_numpy(
                    dtype=int
                )[
                    indices
                ]
            )

            physics = (
                self.validation_physics_scaled[
                    indices
                ]
            )

            mass_column = encode_mass(
                np.full(
                    len(indices),
                    mass,
                    dtype=np.float32,
                ),
                self.mass_transform,
                self.bundle.signal_hypotheses,
            ).reshape(
                -1,
                1,
            )

            X = np.concatenate(
                [
                    physics,
                    mass_column,
                ],
                axis=1,
            ).astype(
                np.float32
            )

            scores = (
                self.model.predict(
                    X,
                    batch_size=(
                        PREDICT_BATCH_SIZE
                    ),
                    verbose=0,
                )
                .reshape(-1)
            )

            auc_value = float(
                roc_auc_score(
                    labels,
                    scores,
                )
            )

            auc_values.append(
                auc_value
            )

            if (
                mass
                >= HIGH_MASS_MIN
            ):
                high_mass_values.append(
                    auc_value
                )

            row[
                f"auc_mPhi_{mass:g}"
            ] = auc_value

        values = np.asarray(
            auc_values,
            dtype=float,
        )

        mean_auc = float(
            np.mean(
                values
            )
        )

        median_auc = float(
            np.median(
                values
            )
        )

        min_auc = float(
            np.min(
                values
            )
        )

        high_mass_auc = float(
            np.mean(
                high_mass_values
            )
        )

        logs[
            "val_mean_mass_auc"
        ] = mean_auc

        logs[
            "val_median_mass_auc"
        ] = median_auc

        logs[
            "val_min_mass_auc"
        ] = min_auc

        logs[
            "val_high_mass_auc"
        ] = high_mass_auc

        row.update(
            {
                "val_mean_mass_auc": (
                    mean_auc
                ),
                "val_median_mass_auc": (
                    median_auc
                ),
                "val_min_mass_auc": (
                    min_auc
                ),
                "val_high_mass_auc": (
                    high_mass_auc
                ),
            }
        )

        self.rows.append(
            row
        )

        if (
            mean_auc
            > self.best_mean_auc
        ):
            self.best_mean_auc = (
                mean_auc
            )

            self.best_median_auc = (
                median_auc
            )

            self.best_min_auc = (
                min_auc
            )

            self.best_high_mass_auc = (
                high_mass_auc
            )

            self.best_epoch = (
                epoch + 1
            )

        if (
            self.trial is not None
        ):
            self.trial.report(
                mean_auc,
                step=epoch,
            )

            self.trial.set_user_attr(
                "latest_min_mass_auc",
                min_auc,
            )

            self.trial.set_user_attr(
                "latest_median_mass_auc",
                median_auc,
            )

            self.trial.set_user_attr(
                "latest_high_mass_auc",
                high_mass_auc,
            )

            if (
                self.enable_pruning
                and self.trial.should_prune()
            ):
                raise (
                    optuna.TrialPruned(
                        f"Pruned at epoch "
                        f"{epoch + 1}: "
                        f"mean mass AUC="
                        f"{mean_auc:.6f}"
                    )
                )


# ============================================================
# TRAIN CONFIGURATION
# ============================================================

def train_configuration(
    bundle: DatasetBundle,
    config: dict,
    seed: int,
    max_epochs: int,
    step_fraction: float,
    validation_df: pd.DataFrame,
    validation_physics_scaled: np.ndarray,
    trial: Optional[
        optuna.Trial
    ] = None,
    enable_pruning: bool = False,
    verbose: int = 0,
) -> TrainingResult:
    tf.keras.backend.clear_session()

    set_seed(
        seed
    )

    model = build_model(
        len(
            bundle.physics_features
        ),
        config,
    )

    parameter_count = int(
        model.count_params()
    )

    if (
        parameter_count
        > MAX_MODEL_PARAMETERS
    ):
        raise (
            optuna.TrialPruned(
                f"Model has "
                f"{parameter_count} parameters, "
                f"above limit "
                f"{MAX_MODEL_PARAMETERS}."
            )
        )

    generator = (
        balanced_batch_generator(
            bundle=bundle,
            batch_size=int(
                config[
                    "batch_size"
                ]
            ),
            mass_transform=(
                config[
                    "mass_transform"
                ]
            ),
            seed=(
                seed + 100
            ),
        )
    )

    steps = steps_per_epoch(
        bundle,
        int(
            config[
                "batch_size"
            ]
        ),
        step_fraction,
    )

    per_mass_callback = (
        PerMassValidationCallback(
            bundle=bundle,
            validation_df=(
                validation_df
            ),
            validation_physics_scaled=(
                validation_physics_scaled
            ),
            mass_transform=(
                config[
                    "mass_transform"
                ]
            ),
            trial=trial,
            enable_pruning=(
                enable_pruning
            ),
        )
    )

    reduce_lr = (
        callbacks.ReduceLROnPlateau(
            monitor=(
                "val_mean_mass_auc"
            ),
            mode="max",
            factor=float(
                config[
                    "reduce_lr_factor"
                ]
            ),
            patience=int(
                config[
                    "reduce_lr_patience"
                ]
            ),
            min_delta=1e-4,
            min_lr=1e-6,
            verbose=0,
        )
    )

    early_stopping = (
        callbacks.EarlyStopping(
            monitor=(
                "val_mean_mass_auc"
            ),
            mode="max",
            patience=int(
                config[
                    "early_stopping_patience"
                ]
            ),
            min_delta=1e-4,
            restore_best_weights=True,
            verbose=0,
        )
    )

    history = model.fit(
        generator,
        steps_per_epoch=steps,
        epochs=max_epochs,
        callbacks=[
            # Must come first so the custom metric is inserted into logs.
            per_mass_callback,
            reduce_lr,
            early_stopping,
        ],
        verbose=verbose,
    )

    history_df = pd.DataFrame(
        history.history
    )

    history_df.insert(
        0,
        "epoch",
        np.arange(
            1,
            len(
                history_df
            )
            + 1,
        ),
    )

    return TrainingResult(
        model=model,
        history=history_df,
        best_mean_auc=float(
            per_mass_callback[
                "best_mean_auc"
            ]
            if isinstance(
                per_mass_callback,
                dict,
            )
            else per_mass_callback.best_mean_auc
        ),
        best_min_auc=float(
            per_mass_callback.best_min_auc
        ),
        best_median_auc=float(
            per_mass_callback.best_median_auc
        ),
        best_high_mass_auc=float(
            per_mass_callback.best_high_mass_auc
        ),
        best_epoch=int(
            per_mass_callback.best_epoch
        ),
        parameter_count=(
            parameter_count
        ),
    )


# ============================================================
# BASELINE TRIALS
# ============================================================

def current_baseline_trial() -> dict:
    return {
        "conditioning_type": (
            "concat"
        ),
        "conditioning_position": (
            "begin"
        ),
        "mass_transform": (
            "raw"
        ),
        "mass_embedding_dim": 1,

        "n_layers": 4,
        "base_units": 128,
        "width_decay": 2.0,
        "min_units": 16,

        "activation": "relu",
        "normalization": "batch",
        "initializer": "he_normal",

        "dropout_first": 0.30,
        "dropout_last": 0.00,

        "use_gaussian_noise": False,

        "optimizer": "adam",
        "learning_rate": 1e-3,
        "clipnorm": 0.0,

        "use_l2": True,
        "l2": 1e-4,

        "label_smoothing": 0.0,
        "batch_size": 1024,

        "reduce_lr_factor": 0.5,
        "reduce_lr_patience": 8,
        "early_extra_patience": 12,
    }


def paper_affine_trial() -> dict:
    output = (
        current_baseline_trial()
        .copy()
    )

    output[
        "conditioning_type"
    ] = "affine"

    output[
        "conditioning_position"
    ] = "all"

    output[
        "conditioning_identity_init"
    ] = True

    return output


# ============================================================
# OPTUNA STUDY
# ============================================================

def create_sampler(
    seed: int,
):
    # group=True is designed for conditional search spaces.
    try:
        return (
            optuna.samplers.TPESampler(
                seed=seed,
                n_startup_trials=30,
                multivariate=True,
                group=True,
            )
        )

    except TypeError:
        # Compatibility fallback for older Optuna installations.
        return (
            optuna.samplers.TPESampler(
                seed=seed,
                n_startup_trials=30,
                multivariate=True,
            )
        )


def create_study(
    args,
):
    study = optuna.create_study(
        study_name=(
            args.study_name
        ),
        direction="maximize",
        sampler=create_sampler(
            args.seed
        ),
        pruner=(
            optuna.pruners.MedianPruner(
                n_startup_trials=20,
                n_warmup_steps=10,
                interval_steps=1,
                n_min_trials=5,
            )
        ),
        storage=(
            resolve_storage_url(
                args.storage
            )
        ),
        load_if_exists=True,
    )

    if (
        args.enqueue_baselines
        and len(
            study.trials
        )
        == 0
    ):
        study.enqueue_trial(
            current_baseline_trial()
        )

        study.enqueue_trial(
            paper_affine_trial()
        )

    return study


def make_objective(
    bundle: DatasetBundle,
    args,
):
    def objective(
        trial: optuna.Trial,
    ) -> float:
        config = (
            sample_trial_config(
                trial
            )
        )

        trial.set_user_attr(
            "layer_sizes",
            "-".join(
                map(
                    str,
                    layer_sizes_from_config(
                        config
                    ),
                )
            ),
        )

        try:
            result = (
                train_configuration(
                    bundle=bundle,
                    config=config,
                    seed=(
                        args.seed
                        + 10000
                        + trial.number
                    ),
                    max_epochs=(
                        args.scan_epochs
                    ),
                    step_fraction=(
                        args.scan_step_fraction
                    ),
                    validation_df=(
                        bundle.scan_validation_df
                    ),
                    validation_physics_scaled=(
                        bundle.scan_validation_physics_scaled
                    ),
                    trial=trial,
                    enable_pruning=True,
                    verbose=0,
                )
            )

            trial.set_user_attr(
                "best_epoch",
                result.best_epoch,
            )

            trial.set_user_attr(
                "best_min_mass_auc",
                result.best_min_auc,
            )

            trial.set_user_attr(
                "best_median_mass_auc",
                result.best_median_auc,
            )

            trial.set_user_attr(
                "best_high_mass_auc",
                result.best_high_mass_auc,
            )

            trial.set_user_attr(
                "n_parameters",
                result.parameter_count,
            )

            return float(
                result.best_mean_auc
            )

        except (
            tf.errors.ResourceExhaustedError,
            MemoryError,
        ) as error:
            raise (
                optuna.TrialPruned(
                    f"Resource exhaustion: "
                    f"{error}"
                )
            )

        finally:
            tf.keras.backend.clear_session()
            gc.collect()

    return objective


# ============================================================
# STUDY OUTPUTS
# ============================================================

def completed_trials(
    study,
) -> list[
    optuna.trial.FrozenTrial
]:
    return [
        trial
        for trial in study.trials
        if (
            trial.state
            == optuna.trial.TrialState.COMPLETE
            and trial.value
            is not None
        )
    ]


def save_study_outputs(
    study,
    outdir: Path,
):
    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    trials_df = (
        study.trials_dataframe()
    )

    trials_df.to_csv(
        outdir
        / "optuna_results.csv",
        index=False,
    )

    completed = (
        completed_trials(
            study
        )
    )

    if len(completed) == 0:
        print(
            "No completed trials yet."
        )

        return

    ranked = sorted(
        completed,
        key=lambda trial: (
            trial.value
        ),
        reverse=True,
    )

    top_rows = []

    for rank, trial in enumerate(
        ranked[
            :50
        ],
        start=1,
    ):
        row = {
            "rank": rank,
            "trial": trial.number,
            "value": trial.value,
        }

        row.update(
            {
                f"param_{key}": value
                for key, value in (
                    trial.params.items()
                )
            }
        )

        row.update(
            {
                f"user_{key}": value
                for key, value in (
                    trial.user_attrs.items()
                )
            }
        )

        top_rows.append(
            row
        )

    pd.DataFrame(
        top_rows
    ).to_csv(
        outdir
        / "top_trials.csv",
        index=False,
    )

    best = ranked[0]

    with open(
        outdir
        / "best_trial.json",
        "w",
    ) as handle:
        json.dump(
            {
                "trial": (
                    best.number
                ),
                "value": (
                    best.value
                ),
                "params": (
                    best.params
                ),
                "user_attrs": (
                    best.user_attrs
                ),
            },
            handle,
            indent=2,
        )

    # Optimization history.
    plt.figure(
        figsize=(10, 6)
    )

    x = []
    y = []
    best_so_far = []

    running_best = (
        -np.inf
    )

    for trial in (
        completed
    ):
        x.append(
            trial.number
        )

        y.append(
            trial.value
        )

        running_best = max(
            running_best,
            trial.value,
        )

        best_so_far.append(
            running_best
        )

    order = np.argsort(
        x
    )

    x = np.asarray(
        x
    )[
        order
    ]

    y = np.asarray(
        y
    )[
        order
    ]

    best_so_far = np.asarray(
        best_so_far
    )[
        order
    ]

    plt.scatter(
        x,
        y,
        s=20,
        alpha=0.5,
        label="Completed trial",
    )

    plt.plot(
        x,
        best_so_far,
        linewidth=2,
        label="Best so far",
    )

    plt.xlabel(
        "Trial"
    )

    plt.ylabel(
        "Mean per-mass validation AUC"
    )

    plt.grid(
        True,
        alpha=0.3,
    )

    plt.legend()

    plt.savefig(
        outdir
        / "optimization_history.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()

    # Conditioning mechanism summary.
    summary_rows = []

    for trial in completed:
        summary_rows.append(
            {
                "value": trial.value,
                "conditioning_type": (
                    trial.params.get(
                        "conditioning_type"
                    )
                ),
                "conditioning_position": (
                    trial.params.get(
                        "conditioning_position"
                    )
                ),
            }
        )

    summary_df = pd.DataFrame(
        summary_rows
    )

    grouped = (
        summary_df
        .groupby(
            [
                "conditioning_type",
                "conditioning_position",
            ]
        )[
            "value"
        ]
        .agg(
            [
                "count",
                "mean",
                "median",
                "max",
            ]
        )
        .reset_index()
    )

    grouped.to_csv(
        outdir
        / "conditioning_summary.csv",
        index=False,
    )

    # Parameter importances.
    try:
        importances = (
            optuna.importance
            .get_param_importances(
                study
            )
        )

        pd.DataFrame(
            {
                "parameter": (
                    list(
                        importances.keys()
                    )
                ),
                "importance": (
                    list(
                        importances.values()
                    )
                ),
            }
        ).to_csv(
            outdir
            / "parameter_importances.csv",
            index=False,
        )

        plt.figure(
            figsize=(10, 7)
        )

        names = list(
            importances.keys()
        )[
            :20
        ]

        values = [
            importances[name]
            for name in names
        ]

        positions = np.arange(
            len(
                names
            )
        )

        plt.barh(
            positions,
            values,
        )

        plt.yticks(
            positions,
            names,
        )

        plt.gca().invert_yaxis()

        plt.xlabel(
            "Optuna importance"
        )

        plt.tight_layout()

        plt.savefig(
            outdir
            / "parameter_importances.png",
            dpi=300,
            bbox_inches="tight",
        )

        plt.close()

    except Exception as error:
        print(
            "Could not compute parameter "
            f"importances: {error}"
        )


# ============================================================
# CONFIG RECONSTRUCTION
# ============================================================

def complete_config_from_params(
    params: dict,
) -> dict:
    config = dict(
        params
    )

    if (
        config[
            "mass_embedding_dim"
        ]
        == 1
    ):
        config[
            "mass_embedding_layers"
        ] = 0

    if (
        config[
            "conditioning_type"
        ]
        not in {
            "bias",
            "scale",
            "affine",
        }
    ):
        config[
            "conditioning_identity_init"
        ] = True

    if (
        config[
            "optimizer"
        ]
        != "adamw"
    ):
        config[
            "weight_decay"
        ] = 0.0

        if not config.get(
            "use_l2",
            False,
        ):
            config[
                "l2"
            ] = 0.0

    else:
        config[
            "l2"
        ] = 0.0

    if not config.get(
        "use_gaussian_noise",
        False,
    ):
        config[
            "gaussian_noise_std"
        ] = 0.0

    config[
        "early_stopping_patience"
    ] = (
        int(
            config[
                "reduce_lr_patience"
            ]
        )
        + int(
            config[
                "early_extra_patience"
            ]
        )
    )

    return config


# ============================================================
# TOP-K MULTI-SEED CONFIRMATION
# ============================================================

def confirm_top_trials(
    study,
    bundle: DatasetBundle,
    args,
    outdir: Path,
) -> Optional[dict]:
    if (
        args.confirm_top_k
        <= 0
    ):
        return None

    completed = sorted(
        completed_trials(
            study
        ),
        key=lambda trial: (
            trial.value
        ),
        reverse=True,
    )

    if len(completed) == 0:
        raise RuntimeError(
            "No completed trials to confirm."
        )

    selected = completed[
        :args.confirm_top_k
    ]

    rows = []

    best_config = None
    best_average = (
        -np.inf
    )

    for rank, trial in enumerate(
        selected,
        start=1,
    ):
        config = (
            complete_config_from_params(
                trial.params
            )
        )

        seed_values = []

        for seed_index in range(
            args.confirm_seeds
        ):
            seed = (
                args.seed
                + 200000
                + 1000
                * trial.number
                + seed_index
            )

            print("")
            print("=" * 88)
            print(
                f"CONFIRMING RANK {rank}, "
                f"TRIAL {trial.number}, "
                f"SEED {seed}"
            )
            print("=" * 88)

            result = (
                train_configuration(
                    bundle=bundle,
                    config=config,
                    seed=seed,
                    max_epochs=(
                        args.confirm_epochs
                    ),
                    step_fraction=(
                        CONFIRM_STEP_FRACTION
                    ),
                    validation_df=(
                        bundle.validation_df
                    ),
                    validation_physics_scaled=(
                        bundle.validation_physics_scaled
                    ),
                    trial=None,
                    enable_pruning=False,
                    verbose=0,
                )
            )

            seed_values.append(
                result.best_mean_auc
            )

            rows.append(
                {
                    "rank": rank,
                    "trial": trial.number,
                    "seed": seed,
                    "scan_value": trial.value,
                    "confirmed_mean_auc": (
                        result.best_mean_auc
                    ),
                    "confirmed_min_auc": (
                        result.best_min_auc
                    ),
                    "confirmed_median_auc": (
                        result.best_median_auc
                    ),
                    "confirmed_high_mass_auc": (
                        result.best_high_mass_auc
                    ),
                    "best_epoch": (
                        result.best_epoch
                    ),
                    "n_parameters": (
                        result.parameter_count
                    ),
                }
            )

            tf.keras.backend.clear_session()
            gc.collect()

        average = float(
            np.mean(
                seed_values
            )
        )

        spread = float(
            np.std(
                seed_values
            )
        )

        print(
            f"Trial {trial.number}: "
            f"confirmed mean="
            f"{average:.6f} "
            f"+/- {spread:.6f}"
        )

        if average > best_average:
            best_average = (
                average
            )

            best_config = {
                "trial_number": (
                    trial.number
                ),
                "scan_value": (
                    trial.value
                ),
                "confirmed_mean_auc": (
                    average
                ),
                "confirmed_std_auc": (
                    spread
                ),
                "params": (
                    trial.params
                ),
                "config": config,
            }

    confirmation_df = pd.DataFrame(
        rows
    )

    confirmation_df.to_csv(
        outdir
        / "confirmation_results.csv",
        index=False,
    )

    with open(
        outdir
        / "best_confirmed_config.json",
        "w",
    ) as handle:
        json.dump(
            best_config,
            handle,
            indent=2,
        )

    return best_config


# ============================================================
# FINAL MODEL / UNTOUCHED TEST
# ============================================================

def evaluate_test_per_mass(
    model: tf.keras.Model,
    bundle: DatasetBundle,
    mass_transform: str,
) -> pd.DataFrame:
    labels_all = (
        bundle.test_df[
            "isSignal"
        ]
        .to_numpy(
            dtype=int
        )
    )

    background_indices = (
        np.flatnonzero(
            labels_all == 0
        )
    )

    true_masses = (
        bundle.test_df[
            "mPhi_true"
        ]
        .to_numpy(
            dtype=float
        )
    )

    rows = []

    for mass in (
        bundle.signal_hypotheses
    ):
        signal_indices = (
            np.flatnonzero(
                (labels_all == 1)
                & np.isclose(
                    true_masses,
                    mass,
                )
            )
        )

        indices = np.concatenate(
            [
                signal_indices,
                background_indices,
            ]
        )

        labels = (
            labels_all[
                indices
            ]
        )

        physics = (
            bundle.test_physics_scaled[
                indices
            ]
        )

        mass_column = encode_mass(
            np.full(
                len(indices),
                mass,
                dtype=np.float32,
            ),
            mass_transform,
            bundle.signal_hypotheses,
        ).reshape(
            -1,
            1,
        )

        X = np.concatenate(
            [
                physics,
                mass_column,
            ],
            axis=1,
        ).astype(
            np.float32
        )

        scores = (
            model.predict(
                X,
                batch_size=(
                    PREDICT_BATCH_SIZE
                ),
                verbose=0,
            )
            .reshape(-1)
        )

        auc_value = float(
            roc_auc_score(
                labels,
                scores,
            )
        )

        fpr, tpr, thresholds = (
            roc_curve(
                labels,
                scores,
            )
        )

        denominator = np.sqrt(
            tpr + fpr
        )

        sigma_ratio_curve = (
            np.divide(
                tpr,
                denominator,
                out=np.zeros_like(
                    tpr,
                    dtype=float,
                ),
                where=(
                    denominator > 0
                ),
            )
        )

        best_index = int(
            np.argmax(
                sigma_ratio_curve
            )
        )

        rows.append(
            {
                "mPhi": mass,
                "auc": auc_value,
                "sigma_ratio": float(
                    sigma_ratio_curve[
                        best_index
                    ]
                ),
                "best_cut": float(
                    thresholds[
                        best_index
                    ]
                ),
            }
        )

    return pd.DataFrame(
        rows
    )


def train_final_model(
    best_config: dict,
    bundle: DatasetBundle,
    args,
    outdir: Path,
):
    config = (
        best_config[
            "config"
        ]
    )

    result = (
        train_configuration(
            bundle=bundle,
            config=config,
            seed=(
                args.seed
                + 999999
            ),
            max_epochs=(
                args.final_epochs
            ),
            step_fraction=(
                CONFIRM_STEP_FRACTION
            ),
            validation_df=(
                bundle.validation_df
            ),
            validation_physics_scaled=(
                bundle.validation_physics_scaled
            ),
            trial=None,
            enable_pruning=False,
            verbose=2,
        )
    )

    model_path = (
        outdir
        / "best_pnn_model.keras"
    )

    result.model.save(
        model_path
    )

    joblib.dump(
        bundle.physics_scaler,
        outdir
        / "best_pnn_physics_scaler.pkl",
    )

    joblib.dump(
        bundle.physics_features,
        outdir
        / "best_pnn_features.pkl",
    )

    joblib.dump(
        bundle.signal_hypotheses,
        outdir
        / "best_pnn_mass_hypotheses.pkl",
    )

    result.history.to_csv(
        outdir
        / "best_pnn_training_history.csv",
        index=False,
    )

    test_metrics = (
        evaluate_test_per_mass(
            result.model,
            bundle,
            config[
                "mass_transform"
            ],
        )
    )

    test_metrics.to_csv(
        outdir
        / "best_pnn_test_metrics.csv",
        index=False,
    )

    print("")
    print("Final untouched-test metrics:")
    print(
        test_metrics.to_string(
            index=False
        )
    )


# ============================================================
# MAIN
# ============================================================

def main():
    args = parse_args()

    set_seed(
        args.seed
    )

    outdir = Path(
        args.outdir
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    ROOT.EnableImplicitMT(
        8
    )

    print("")
    print("=" * 88)
    print("OPTUNA PARAMETRIC-DNN ARCHITECTURE SEARCH")
    print("=" * 88)
    print(
        "TTNuNu in main training: "
        f"{args.include_ttnunu}"
    )

    discovered = discover_files(
        args.snapshot_dir,
        args.include_ttnunu,
    )

    physics_features = (
        resolve_physics_features(
            discovered,
            args.strict_feature_check,
        )
    )

    (
        signal,
        backgrounds,
        signal_hypotheses,
    ) = load_data(
        discovered,
        physics_features,
    )

    requested_fractions = (
        BACKGROUND_FRACTIONS_WITH_TTNUNU
        if args.include_ttnunu
        else BACKGROUND_FRACTIONS_WITHOUT_TTNUNU
    )

    active_fractions = (
        normalize_background_fractions(
            requested_fractions,
            sorted(
                backgrounds.keys()
            ),
        )
    )

    print("")
    print("Fixed background fractions:")
    for group, fraction in (
        active_fractions.items()
    ):
        print(
            f"  {group:15s}: "
            f"{fraction:.3f}"
        )

    bundle = build_dataset_bundle(
        signal=signal,
        backgrounds=backgrounds,
        signal_hypotheses=(
            signal_hypotheses
        ),
        physics_features=(
            physics_features
        ),
        background_fractions=(
            active_fractions
        ),
        args=args,
    )

    joblib.dump(
        bundle.physics_features,
        outdir
        / "scan_physics_features.pkl",
    )

    joblib.dump(
        bundle.signal_hypotheses,
        outdir
        / "scan_mass_hypotheses.pkl",
    )

    study = create_study(
        args
    )

    if args.trials > 0:
        study.optimize(
            make_objective(
                bundle,
                args,
            ),
            n_trials=(
                args.trials
            ),
            n_jobs=1,
            gc_after_trial=True,
        )

    save_study_outputs(
        study,
        outdir,
    )

    if len(
        completed_trials(
            study
        )
    ) > 0:
        print("")
        print(
            "Best scan value:",
            study.best_value,
        )

        print(
            "Best scan params:"
        )

        for key, value in (
            study.best_params.items()
        ):
            print(
                f"  {key}: {value}"
            )

    best_confirmed = (
        confirm_top_trials(
            study,
            bundle,
            args,
            outdir,
        )
    )

    if args.train_final:
        if best_confirmed is None:
            if len(
                completed_trials(
                    study
                )
            ) == 0:
                raise RuntimeError(
                    "No completed Optuna trial."
                )

            best_confirmed = {
                "trial_number": (
                    study.best_trial.number
                ),
                "scan_value": (
                    study.best_trial.value
                ),
                "confirmed_mean_auc": (
                    None
                ),
                "confirmed_std_auc": (
                    None
                ),
                "params": (
                    study.best_trial.params
                ),
                "config": (
                    complete_config_from_params(
                        study.best_trial.params
                    )
                ),
            }

        train_final_model(
            best_confirmed,
            bundle,
            args,
            outdir,
        )


if __name__ == "__main__":
    main()
