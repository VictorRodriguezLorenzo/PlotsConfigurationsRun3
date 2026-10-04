#!/usr/bin/env python3
"""Train a Run-3 parametric three-class top-DM classifier.

The three targets are background, ttDM, and tWDM.  By default snapshots from
every available ``Full*/DNNmodels/files_for_training`` campaign are combined;
``--campaign`` can be repeated to select a subset.  This file contains the
shared implementation for both b-jet categories.  The companion
``train_categorical_tWDM.py`` selects the 1b category.
"""

import argparse
import glob
import json
import re
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import tensorflow as tf
from sklearn.metrics import auc, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from tensorflow.keras import callbacks, regularizers
from tensorflow.keras.layers import Dense, Dropout, Input


CLASS_NAMES = ("background", "ttDM", "tWDM")
SIGNAL_TYPES = ("ps", "s")
SIGNAL_PREFIXES = {
    "ttDM": "TTto2LDMsimpSpin0",
    "tWDM": "TWto2LDMsimpSpin0",
}
DEFAULT_CAMPAIGNS = (
    "Full2022v12",
    "Full2022EEv12",
    "Full2023v12",
    "Full2023BPixv12",
    "Full2024v15",
)

FEATURES = list(dict.fromkeys([
    "lep_pt1", "lep_pt2", "lep_eta1", "lep_eta2", "mll", "ptll", "drll",
    "detall", "dphill", "yll", "PuppiMET_pt", "PuppiMET_phi", "dphilmet",
    "dphilmet1", "dphilmet2", "dphillmet", "mtw1", "mtw2", "mth", "mTi",
    "mR", "mT2", "mTe", "recoil", "upara", "uperp", "pTWW", "mcoll",
    "mcollWW", "choiMass", "nbjet_jet_ratio", "njet", "ht", "vht_pt",
    "dphijet1met", "dphijet2met", "dphijjmet", "chel", "pdark",
    "dphi_ttbar", "dphi_met_llb", "mt2blbl", "dphi_met_ll", "st",
    "met_over_sqrt_ht", "met_over_st", "dphi_min_j_met", "pt_llb",
    "dphi_met_llb_safe", "pt_llbb", "dphi_met_llbb", "m_llbb",
    "mT_llbb_met", "met_over_pt_llbb", "mt2_bell_l", "mbl_min", "mbl_max",
    "mtbl", "mbb", "drbb", "ptbb", "max_nonleading_btag", "pt_b2",
    "nForwardJet", "leadingForwardJet_pt", "leadingForwardJet_absEta",
    "deta_forwardJet_b", "dphi_forwardJet_met", "top1_pt_reco",
    "top2_pt_reco", "pdark_over_met", "angle_ll_llbb_rf",
    "dphi_ll_llbb_rf", "cos_l1_llbb_rf", "cos_l2_llbb_rf",
    "angle_ll_llmet_rf", "dphi_ll_llmet_rf", "cos_l1_llmet_rf",
    "cos_l2_llmet_rf",
]))

BACKGROUND_PATTERNS = {
    "tt2l": ("TTTo2L2Nu",),
    "other_tt": ("TTToSemiLeptonic",),
    "single_top": (
        "ST_t-channel_top", "ST_t-channel_antitop", "ST_s-channel_plus",
        "ST_s-channel_minus", "TWminusto2L2Nu", "TbarWplusto2L2Nu",
        "ST_tW_top", "ST_tW_antitop",
    ),
    "other_ttZ": ("TTLL_MLL-4to50", "TTLL_MLL-50", "TTZ-ZtoQQ"),
    "ttW": ("TTLNu", "TTW"),
    "ttH": ("ttH", "TTH"),
    "TTNuNu": ("TTNuNu",),
}
BACKGROUND_FRACTIONS = {
    "tt2l": 0.50,
    "other_tt": 0.08,
    "single_top": 0.12,
    "other_ttZ": 0.05,
    "ttW": 0.05,
    "ttH": 0.05,
    "TTNuNu": 0.15,
}


def parse_args(default_category="ttDM"):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--category", choices=("ttDM", "tWDM"), default=default_category)
    parser.add_argument("--campaign", action="append", choices=DEFAULT_CAMPAIGNS)
    parser.add_argument("--snapshot-dir", action="append", type=Path,
                        help="Explicit snapshot directory; may be repeated.")
    parser.add_argument("--signal-type", choices=SIGNAL_TYPES,
                        help="Train only one CP hypothesis (default: train both).")
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parent / "Models")
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--seed", type=int, default=12345)
    return parser.parse_args()


def matches(name, patterns):
    lower_name = name.lower()
    return any(pattern.lower() in lower_name for pattern in patterns)


def sample_tag(path):
    match = re.search(r"nanoLatino_(.*)__part\d+_snapshot", path.name)
    return match.group(1) if match else path.stem


def sample_mass(name):
    match = re.search(r"m[Pp]hi[_-]?(\d+)", name)
    return float(match.group(1)) if match else None


def classify_file(path, signal_type):
    name = path.name
    for mode, prefix in SIGNAL_PREFIXES.items():
        if f"{prefix}_{signal_type}_" in name:
            return mode
    # TTNuNu must be checked before the broad tt patterns.
    for group in ("TTNuNu", "tt2l", "other_tt", "single_top", "other_ttZ", "ttW", "ttH"):
        if matches(name, BACKGROUND_PATTERNS[group]):
            return group
    return None


def discover_inputs(args):
    analysis_dir = Path(__file__).resolve().parents[1]
    if args.snapshot_dir:
        directories = args.snapshot_dir
    else:
        campaigns = args.campaign or list(DEFAULT_CAMPAIGNS)
        directories = [analysis_dir / campaign / "DNNmodels" / "files_for_training"
                       for campaign in campaigns]

    inputs = []
    for directory in directories:
        campaign = directory.resolve().parents[1].name
        files = [Path(name) for name in sorted(glob.glob(str(directory / "*.root")))]
        if not files:
            raise RuntimeError(f"No ROOT snapshots found in {directory}")
        inputs.extend((campaign, path) for path in files)
    return inputs


def read_group(entries, selection, description):
    chain = ROOT.TChain("Events")
    for _, path in entries:
        chain.Add(str(path))
    rdf = ROOT.RDataFrame(chain).Filter(selection)
    count = rdf.Count().GetValue()
    print(f"{description}: {len(entries)} files, {count} selected events")
    if count == 0:
        return pd.DataFrame(columns=FEATURES)
    return pd.DataFrame(rdf.AsNumpy(FEATURES))


def build_dataframe(inputs, signal_type, selection):
    grouped = {}
    for campaign, path in inputs:
        group = classify_file(path, signal_type)
        if group is not None:
            tag = sample_tag(path)
            key = (campaign, group, tag if group in SIGNAL_PREFIXES else group)
            grouped.setdefault(key, []).append((campaign, path))

    frames = []
    for (campaign, group, tag), entries in sorted(grouped.items()):
        frame = read_group(entries, selection, f"{campaign}/{tag}")
        if frame.empty:
            continue
        is_signal = group in SIGNAL_PREFIXES
        mass = sample_mass(tag) if is_signal else np.nan
        if is_signal and mass is None:
            raise RuntimeError(f"Cannot extract mPhi from {tag}")
        frame["target"] = CLASS_NAMES.index(group) if is_signal else 0
        frame["class_name"] = group if is_signal else "background"
        frame["process_group"] = group
        frame["campaign"] = campaign
        frame["mPhi_true"] = mass
        frames.append(frame)

    if not frames:
        raise RuntimeError("No selected training events were found")
    data = pd.concat(frames, ignore_index=True)
    data[FEATURES] = data[FEATURES].replace([np.inf, -np.inf], np.nan)
    data.dropna(subset=FEATURES, inplace=True)
    for class_name in CLASS_NAMES:
        if not np.any(data["class_name"] == class_name):
            raise RuntimeError(f"Training class {class_name} is empty")
    return data.reset_index(drop=True)


def split_data(data, seed):
    labels = data["campaign"].astype(str) + ":" + data["class_name"].astype(str)
    signal = data["target"] != 0
    labels.loc[signal] += ":" + data.loc[signal, "mPhi_true"].astype(int).astype(str)
    try:
        return train_test_split(data, test_size=0.20, random_state=seed, stratify=labels)
    except ValueError as error:
        print(f"Fine-grained stratification unavailable ({error}); stratifying by target")
        return train_test_split(data, test_size=0.20, random_state=seed,
                                stratify=data["target"])


def add_mass_hypothesis(frame, masses, rng, fixed_mass=None):
    result = frame[FEATURES].copy()
    assigned = frame["mPhi_true"].to_numpy(dtype=float, copy=True)
    background = frame["target"].to_numpy() == 0
    assigned[background] = fixed_mass if fixed_mass is not None else rng.choice(
        masses, size=np.count_nonzero(background), replace=True)
    result["mPhi"] = assigned
    return result


def balanced_generator(frame, scaler, features, masses, batch_size, seed):
    rng = np.random.default_rng(seed)
    signal_pools = {}
    for target in (1, 2):
        target_frame = frame[frame["target"] == target]
        keys = target_frame["campaign"] + ":" + target_frame["mPhi_true"].astype(str)
        signal_pools[target] = [
            indices.to_numpy() for _, indices in target_frame.groupby(keys).groups.items()
        ]

    background_pools = {}
    background_frame = frame[frame["target"] == 0]
    for group, group_frame in background_frame.groupby("process_group"):
        background_pools[group] = [
            indices.to_numpy()
            for _, indices in group_frame.groupby("campaign").groups.items()
        ]

    class_counts = [batch_size // len(CLASS_NAMES)] * len(CLASS_NAMES)
    for index in range(batch_size % len(CLASS_NAMES)):
        class_counts[index] += 1

    while True:
        chosen = []
        for target, total in enumerate(class_counts):
            if target != 0:
                target_pools = signal_pools[target]
                base, remainder = divmod(total, len(target_pools))
                for index, pool in enumerate(target_pools):
                    count = base + (index < remainder)
                    chosen.extend(rng.choice(pool, count, replace=True))
                continue

            groups = list(background_pools)
            weights = np.array([BACKGROUND_FRACTIONS[group] for group in groups])
            raw_counts = total * weights / weights.sum()
            group_counts = np.floor(raw_counts).astype(int)
            for index in np.argsort(raw_counts - group_counts)[::-1][:total - group_counts.sum()]:
                group_counts[index] += 1
            for group, group_total in zip(groups, group_counts):
                campaign_pools = background_pools[group]
                base, remainder = divmod(int(group_total), len(campaign_pools))
                for index, pool in enumerate(campaign_pools):
                    count = base + (index < remainder)
                    chosen.extend(rng.choice(pool, count, replace=True))
        batch = frame.loc[chosen].copy()
        batch = batch.iloc[rng.permutation(len(batch))]
        x_batch = add_mass_hypothesis(batch, masses, rng)
        yield scaler.transform(x_batch[features]).astype(np.float32), batch["target"].to_numpy()


@tf.keras.utils.register_keras_serializable(package="topDM")
class AffineConditioning(tf.keras.layers.Layer):
    """Mass-dependent feature-wise affine transformation."""

    def __init__(self, units, **kwargs):
        super().__init__(**kwargs)
        self.units = units
        self.gamma = Dense(units, kernel_initializer="zeros", bias_initializer="ones")
        self.beta = Dense(units, kernel_initializer="zeros", bias_initializer="zeros")

    def call(self, inputs):
        hidden, mass = inputs
        return self.gamma(mass) * hidden + self.beta(mass)

    def get_config(self):
        return {**super().get_config(), "units": self.units}


def build_model(n_inputs):
    inputs = Input(shape=(n_inputs,), name="inputs")
    physics, mass = inputs[:, :-1], inputs[:, -1:]
    hidden = physics
    for units, dropout in ((128, 0.3), (64, 0.3), (32, 0.2), (16, 0.0)):
        hidden = Dense(units, activation="relu",
                       kernel_regularizer=regularizers.l2(1e-4))(hidden)
        hidden = AffineConditioning(units, name=f"affine_{units}")([hidden, mass])
        if dropout:
            hidden = Dropout(dropout)(hidden)
    output = Dense(len(CLASS_NAMES), activation="softmax", name="class_probabilities")(hidden)
    return tf.keras.Model(inputs, output)


def plot_roc(model, test, scaler, features, masses, output_path, seed):
    rng = np.random.default_rng(seed)
    x_test = scaler.transform(add_mass_hypothesis(test, masses, rng)[features])
    predictions = model.predict(x_test, batch_size=4096, verbose=0)
    plt.figure(figsize=(7, 6))
    for target, class_name in enumerate(CLASS_NAMES):
        truth = (test["target"].to_numpy() == target).astype(int)
        fpr, tpr, _ = roc_curve(truth, predictions[:, target])
        plt.plot(fpr, tpr, label=f"{class_name} (AUC={auc(fpr, tpr):.3f})")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlabel("False-positive rate")
    plt.ylabel("True-positive rate")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def train_one(args, signal_type):
    tf.keras.utils.set_random_seed(args.seed)
    np.random.seed(args.seed)
    selection = "nbjets > 1 && tt_reco" if args.category == "ttDM" else "nbjets == 1"
    inputs = discover_inputs(args)
    data = build_dataframe(inputs, signal_type, selection)
    train, test = split_data(data, args.seed)
    train, test = train.reset_index(drop=True), test.reset_index(drop=True)
    masses = sorted(data.loc[data["target"] != 0, "mPhi_true"].unique())
    features = FEATURES + ["mPhi"]

    scaler_rng = np.random.default_rng(args.seed + 1)
    scaler = StandardScaler().fit(add_mass_hypothesis(train, masses, scaler_rng)[features])
    generator = balanced_generator(train, scaler, features, masses,
                                   args.batch_size, args.seed + 2)
    validation_rng = np.random.default_rng(args.seed + 3)
    x_test = scaler.transform(add_mass_hypothesis(test, masses, validation_rng)[features])
    y_test = test["target"].to_numpy()

    model = build_model(len(features))
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
        loss="sparse_categorical_crossentropy",
        metrics=[tf.keras.metrics.SparseCategoricalAccuracy(name="accuracy")],
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"model_DNN_categorical_{args.category}_{signal_type}"
    history = model.fit(
        generator,
        steps_per_epoch=max(1, int(np.ceil(len(train) / args.batch_size))),
        epochs=args.epochs,
        validation_data=(x_test, y_test),
        callbacks=[
            callbacks.EarlyStopping(monitor="val_loss", patience=15,
                                    restore_best_weights=True),
            callbacks.ReduceLROnPlateau(monitor="val_loss", patience=8,
                                       factor=0.5, min_lr=1e-6),
        ],
        verbose=2,
    )
    model.save(args.output_dir / f"{stem}.keras")
    joblib.dump(scaler, args.output_dir / f"scaler_{stem}.pkl")
    joblib.dump(features, args.output_dir / f"features_{stem}.pkl")
    pd.DataFrame(history.history).to_csv(args.output_dir / f"history_{stem}.csv", index=False)
    with open(args.output_dir / f"config_{stem}.json", "w", encoding="utf-8") as handle:
        json.dump({
            "classes": CLASS_NAMES,
            "class_indices": {name: index for index, name in enumerate(CLASS_NAMES)},
            "category": args.category,
            "category_selection": selection,
            "signal_type": signal_type,
            "campaigns": sorted(data["campaign"].unique().tolist()),
            "masses": masses,
            "features": features,
            "output_contract": "[background, ttDM, tWDM] softmax probabilities",
        }, handle, indent=2)
    plot_roc(model, test, scaler, features, masses,
             args.output_dir / f"roc_{stem}.png", args.seed + 4)


def main(default_category="ttDM"):
    args = parse_args(default_category)
    for signal_type in (args.signal_type,) if args.signal_type else SIGNAL_TYPES:
        print(f"\nTraining {args.category} category, {signal_type} hypothesis")
        train_one(args, signal_type)


if __name__ == "__main__":
    main()
