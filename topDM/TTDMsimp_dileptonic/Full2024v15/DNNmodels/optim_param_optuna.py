#!/usr/bin/env python3
"""Hyperparameter scan with Optuna for TTDM dileptonic DNN."""

import argparse
import glob
import os
import random
import re
import multiprocessing

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ROOT
import tensorflow as tf
import seaborn as sns
import optuna

from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from tensorflow.keras import callbacks, regularizers
from tensorflow.keras.layers import BatchNormalization, Dense, Dropout, GaussianNoise, Input
from tensorflow.keras.models import Sequential

# prevent TF thread contention
tf.config.threading.set_intra_op_parallelism_threads(1)
tf.config.threading.set_inter_op_parallelism_threads(1)

FEATURES = [
    "PuppiMET_pt","mT2","dphill",
    "chel","pdark","dphi_ttbar","dphi_met_llb",
]

BKG_NAMES = [
    "TTTo2L2Nu","ST_t-channel_top","ST_t-channel_antitop","ST_s-channel_plus","ST_s-channel_minus",
    "TWminusto2L2Nu","TbarWplusto2L2Nu","ST_tW_top","ST_tW_antitop","DYto2L-2Jets_MLL-50",
    "DYto2L-2Jets_MLL-10to50",
]


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("--snapshot-glob",default="files_for_training/*.root")
    parser.add_argument("--outdir",default="hyperopt_summary")

    parser.add_argument("--epochs",type=int,default=80)
    parser.add_argument("--patience",type=int,default=12)
    parser.add_argument("--test-size",type=float,default=0.2)

    parser.add_argument("--workers",type=int,
        default=max(1,multiprocessing.cpu_count()//2))

    parser.add_argument("--trials",type=int,default=5000)

    parser.add_argument("--seed",type=int,default=42)

    return parser.parse_args()


def downsample(chain,target):

    current=chain.GetEntries()

    rdf=ROOT.RDataFrame(chain)

    if current<=target:
        return rdf

    frac=target/current

    return rdf.Define("rnd","gRandom->Rndm()").Filter(f"rnd<{frac}")


def build_dataset(snapshot_glob):

    files=glob.glob(snapshot_glob)

    files_bkg=[f for f in files if any(name in f for name in BKG_NAMES)]
    files_sig=[f for f in files if f not in files_bkg]

    sample_groups={}

    for filepath in files_sig:
        match=re.search(r"nanoLatino_(.*)__part\d+_snapshot",filepath)
        if match:
            tag=match.group(1)
            sample_groups.setdefault(tag,[]).append(filepath)

    chain_bkg=ROOT.TChain("Events")

    for f in files_bkg:
        chain_bkg.Add(f)

    sig_chains={}
    sig_counts={}

    for tag,flist in sample_groups.items():

        chain=ROOT.TChain("Events")

        for f in flist:
            chain.Add(f)

        sig_chains[tag]=chain
        sig_counts[tag]=chain.GetEntries()

    target=min(min(sig_counts.values()),chain_bkg.GetEntries())

    sig_rdfs={k:downsample(c,target) for k,c in sig_chains.items()}

    bkg_rdf=ROOT.RDataFrame(chain_bkg)

    sig_pd={}

    for key,rdf in sig_rdfs.items():

        df=pd.DataFrame(rdf.AsNumpy(FEATURES))

        df["isSignal"]=1

        sig_pd[key]=df

    bkg=pd.DataFrame(bkg_rdf.AsNumpy(FEATURES))
    bkg["isSignal"]=0

    sig=pd.concat(sig_pd.values(),ignore_index=True)

    df_all=pd.concat([sig,bkg],ignore_index=True)

    df_all[FEATURES]=df_all[FEATURES].replace([np.inf,-np.inf],np.nan)

    df_all.dropna(inplace=True)

    return df_all


def build_model(input_dim,hp):

    model=Sequential([Input(shape=(input_dim,))])

    if hp["gaussian_noise_std"]>0:
        model.add(GaussianNoise(hp["gaussian_noise_std"]))

    for size in hp["layer_sizes"]:

        model.add(Dense(
            size,
            activation=hp["activation"],
            kernel_regularizer=regularizers.l2(hp["l2"])
        ))

        if hp["batch_norm"]:
            model.add(BatchNormalization())

        if hp["dropout"]>0:
            model.add(Dropout(hp["dropout"]))

    model.add(Dense(1,activation="sigmoid"))

    # --- Optimizer choice from Optuna ---
    if hp["optimizer"] == "adam":
        opt = tf.keras.optimizers.Adam(
            learning_rate=hp["learning_rate"]
        )

    elif hp["optimizer"] == "rmsprop":
        opt = tf.keras.optimizers.RMSprop(
            learning_rate=hp["learning_rate"]
        )

    elif hp["optimizer"] == "nadam":
        opt = tf.keras.optimizers.Nadam(
            learning_rate=hp["learning_rate"]
        )

    else:
        raise ValueError(f"Unknown optimizer: {hp['optimizer']}")

    model.compile(
        loss=tf.keras.losses.BinaryCrossentropy(label_smoothing=hp["label_smoothing"]),
        optimizer=opt,
        metrics=[tf.keras.metrics.AUC(name="auc")]
    )

    return model


def balanced_batch_generator(X, y, batch_size=1024):
    y = np.array(y)
    sig_idx = np.where(y == 1)[0]
    bkg_idx = np.where(y == 0)[0]
    while True:
        sig_batch = np.random.choice(sig_idx, batch_size // 2)
        bkg_batch = np.random.choice(bkg_idx, batch_size // 2)
        idx = np.concatenate([sig_batch, bkg_batch])
        np.random.shuffle(idx)
        yield X[idx], y[idx]

def objective(trial):
    # --- Hyperparameters ---
    hp = {}
    hp["n_layers"] = trial.suggest_int("n_layers", 2, 5)
    hp["size_base"] = trial.suggest_categorical("size_base", [32,64,96,128,192,256])
    hp["size_decay_factor"] = trial.suggest_categorical("size_decay_factor", [1.0,1.5,2.0])
    hp["min_units"] = trial.suggest_categorical("min_units", [8,16])
    hp["activation"] = trial.suggest_categorical("activation", ["relu","elu","selu","gelu"])
    hp["batch_size"] = trial.suggest_categorical("batch_size", [32,64,128,256,512])
    hp["dropout"] = trial.suggest_categorical("dropout", [0.0,0.1,0.2,0.3,0.4])
    hp["gaussian_noise_std"] = trial.suggest_categorical("gaussian_noise_std", [0.0,0.01,0.03])
    hp["learning_rate"] = trial.suggest_categorical("learning_rate", [1e-2,5e-3,1e-3,5e-4,1e-4])
    hp["batch_norm"] = trial.suggest_categorical("batch_norm", [True, False])
    hp["l2"] = trial.suggest_categorical("l2", [0.0,1e-5,1e-4,1e-3])
    hp["label_smoothing"] = trial.suggest_categorical("label_smoothing", [0.0,0.01,0.03])
    hp["optimizer"] = trial.suggest_categorical("optimizer", ["adam","rmsprop","nadam"])

    hp["layer_sizes"] = [
        max(int(hp["size_base"] / (hp["size_decay_factor"] ** i)), hp["min_units"])
        for i in range(hp["n_layers"])
    ]

    model = build_model(X_train_scaled.shape[1], hp)
    
    early = callbacks.EarlyStopping(
        monitor='val_auc',
        mode='max',
        patience=args.patience,
        restore_best_weights=True
    )
    
    reduce_lr = callbacks.ReduceLROnPlateau(
        monitor='val_loss',
        factor=0.5,
        patience=8,
        min_lr=1e-6
    )
    
    train_gen = balanced_batch_generator(
        X_train_scaled,
        y_train,
        batch_size=hp["batch_size"]
    )
    
    steps_per_epoch = max(
        1,
        len(X_train_scaled) // (15 * hp["batch_size"])
    )
    
    # ------------------------------------------------------------
    # Reduced validation sample for faster Optuna trials
    # ------------------------------------------------------------
    y_test_np = np.asarray(y_test).reshape(-1)
    
    rng = np.random.default_rng(args.seed)
    
    sig_val_idx = np.where(y_test_np == 1)[0]
    bkg_val_idx = np.where(y_test_np == 0)[0]
    
    n_val_per_class = 50_000
    
    n_sig = min(n_val_per_class, len(sig_val_idx))
    n_bkg = min(n_val_per_class, len(bkg_val_idx))
    
    chosen_sig = rng.choice(sig_val_idx, size=n_sig, replace=False)
    chosen_bkg = rng.choice(bkg_val_idx, size=n_bkg, replace=False)
    
    val_idx = np.concatenate([chosen_sig, chosen_bkg])
    rng.shuffle(val_idx)
    
    X_val_small = X_test_scaled[val_idx]
    y_val_small = y_test_np[val_idx]
    
    history = model.fit(
        train_gen,
        steps_per_epoch=steps_per_epoch,
        validation_data=(X_val_small, y_val_small),
        epochs=args.epochs,
        callbacks=[early, reduce_lr],
        verbose=0
    )
    
    best_val_auc = float(np.max(history.history["val_auc"]))
    
    trial.set_user_attr(
        "layer_sizes",
        "-".join(map(str, hp["layer_sizes"]))
    )
    
    return best_val_auc


def plot_hyperparameter_vs_auc(results_df,outdir):

    os.makedirs(outdir,exist_ok=True)

    for col in results_df.columns:

        if col in ["value","number"]:
            continue

        plt.figure(figsize=(8,6))

        sns.boxplot(x=results_df[col],y=results_df["value"])

        if col=="layer_sizes":
            plt.xticks(rotation=90,fontsize=8)
        else:
            plt.xticks(rotation=45)

        plt.ylabel("Validation AUC")
        plt.xlabel(col)

        plt.tight_layout()

        plt.savefig(os.path.join(outdir,f"{col}_vs_auc.png"))

        plt.close()


def plot_feature_correlation(df,outdir):

    corr=df[FEATURES].corr()

    plt.figure(figsize=(8,6))

    sns.heatmap(corr,annot=True,cmap="coolwarm")

    plt.title("Feature Correlation Matrix")

    plt.tight_layout()

    plt.savefig(os.path.join(outdir,"feature_correlations.png"))

    plt.close()


def main():

    global args,X_train,X_test,y_train,y_test

    args=parse_args()

    os.makedirs(args.outdir,exist_ok=True)

    np.random.seed(args.seed)
    random.seed(args.seed)
    tf.keras.utils.set_random_seed(args.seed)

    df=build_dataset(args.snapshot_glob)

    plot_feature_correlation(df,args.outdir)

    X=df[FEATURES]
    y=df["isSignal"]

    X_train,X_test,y_train,y_test=train_test_split(
        X,y,test_size=args.test_size,random_state=args.seed,stratify=y
    )

    scaler=StandardScaler()

    global X_train_scaled,X_test_scaled  

    X_train_scaled=scaler.fit_transform(X_train)
    X_test_scaled=scaler.transform(X_test)

    study=optuna.create_study(direction="maximize")

    study.optimize(
        objective,
        n_trials=args.trials,
        n_jobs=args.workers
    )

    df_results=study.trials_dataframe()

    df_results["layer_sizes"]=[
        t.user_attrs.get("layer_sizes","")
        for t in study.trials
    ]

    df_results.to_csv(os.path.join(args.outdir,"optuna_results.csv"),index=False)

    plot_hyperparameter_vs_auc(df_results,args.outdir)

    print("Best AUC:",study.best_value)

    print("Best params:",study.best_params)


if __name__=="__main__":
    main()
