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
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization, Activation
from tensorflow.keras import callbacks, optimizers, regularizers
from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from sklearn import metrics

####### CONTROL PARAMETERS FOR THE WHOLE CODE #######

ANALYSIS_NAME = "dnn_model"
save_model = True

loaded_model = False

MODEL_NAME = "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2024v15/DNNmodels/model_DNN.h5"

PARAMETRIC = True

# ============================================================
# VARIABLES
# ============================================================

# Branches to load
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

if PARAMETRIC:
    var.append("mPhi")

#var = [
#    'dphill',
#    'PuppiMET_pt',
#    'mT2',
#    'pdark',
#    'chel',
#    'dphi_ttbar',
#    'dphi_met_llb',
#]

# Find all snapshot ROOT files
files = glob.glob("/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2024v15/DNNmodels/files_for_training/*.root")

# Define background substrings
bkg_names = ['TTTo2L2Nu', 'ST_t-channel_top', 'ST_t-channel_antitop', 'ST_s-channel_plus', 'ST_s-channel_minus', 'TWminusto2L2Nu', 'TbarWplusto2L2Nu', 'ST_tW_top', 'ST_tW_antitop', 'DYto2L-2Jets_MLL-50', 'DYto2L-2Jets_MLL-10to50']

# Sort files
files_bkg = [f for f in files if any(name in f for name in bkg_names)]
files_sig = [f for f in files if f not in files_bkg]

print(f"Signal files: {len(files_sig)}, Background files: {len(files_bkg)}")

# -------------------------------
# GROUP SIGNAL FILES BY SAMPLE TAG
# -------------------------------
sample_groups = {}
for f in files_sig:
    # extract the tag between 'nanoLatino_' and '__part'
    m = re.search(r"nanoLatino_(.*)__part\d+_snapshot", f)
    if m:
        tag = m.group(1)
        if tag not in sample_groups:
            sample_groups[tag] = []
        sample_groups[tag].append(f)

print("Signal sample groups:", {k: len(v) for k, v in sample_groups.items()})

def extract_mphi(sample_tag):
    match = re.search(r"m[Pp]hi[_-]?(\d+)", sample_tag)
    return float(match.group(1)) if match else None

# -------------------------------
# CREATE TCHAINS AND GET COUNTS
# -------------------------------

# Background chain
chain_bkg = ROOT.TChain("Events")
for f in files_bkg:
    chain_bkg.Add(f)
bkg_count = chain_bkg.GetEntries()

# Signal chains (one chain per sample)
sig_chains = {}
sig_counts = {}
for sample_tag, file_list in sample_groups.items():
    chain = ROOT.TChain("Events")
    for f in file_list:
        chain.Add(f)
    sig_chains[sample_tag] = chain
    sig_counts[sample_tag] = chain.GetEntries()

print("Signal counts per sample:", sig_counts)
print("Background count:", bkg_count)

# -------------------------------
# DETERMINE TARGET SIZE FOR DOWNSAMPLING
# -------------------------------
target_size = min(min(sig_counts.values()), bkg_count)
print("Target size for downsampling:", target_size)

# -------------------------------
# DOWNSAMPLE FUNCTION
# -------------------------------
def downsample(chain, target):
    current = chain.GetEntries()
    if current <= target:
        return ROOT.RDataFrame(chain)
    frac = target / current
    rdf = ROOT.RDataFrame(chain)
    rdf = rdf.Define("rnd", "gRandom->Rndm()").Filter(f"rnd < {frac}")
    return rdf

# -------------------------------
# DOWNSAMPLE SIGNALS
# -------------------------------
dataframes = {key: downsample(chain, target_size) for key, chain in sig_chains.items()}
rdf_bkg = ROOT.RDataFrame(chain_bkg)

# -------------------------------
# CONVERT TO PANDAS AND ADD LABELS
# -------------------------------
pd_dataframes = {}
for key, rdf in dataframes.items():
    features_to_read = [v for v in var if v != "mPhi"]
    pd_df = pd.DataFrame(rdf.AsNumpy(features_to_read))

    if PARAMETRIC:
        mphi = extract_mphi(key)
        if mphi is None:
            raise RuntimeError(f"Could not infer mPhi from signal tag: {key}")
        pd_df["mPhi"] = mphi
    pd_df['isSignal'] = 1
    pd_df['isBkg'] = 0
    pd_dataframes[key] = pd_df

features_to_read = [v for v in var if v != "mPhi"]
Bkg = pd.DataFrame(rdf_bkg.AsNumpy(features_to_read))

if PARAMETRIC:
    signal_hypotheses = sorted({extract_mphi(k) for k in sample_groups if extract_mphi(k) is not None})
    if not signal_hypotheses:
        raise RuntimeError("No signal mPhi hypotheses found for parametric training.")
    Bkg["mPhi"] = np.random.choice(signal_hypotheses, len(Bkg))
Bkg['isSignal'] = 0
Bkg['isBkg'] = 1

# -------------------------------
# CONCATENATE AND BALANCE
# -------------------------------
Sig = pd.concat([df for df in pd_dataframes.values()])

df_all = pd.concat([Sig, Bkg])
df_all[var] = df_all[var].replace([np.inf, -np.inf], np.nan)

print("NaNs on dataframe:\n", df_all.isna().sum())

df_all.dropna(inplace=True)

print("\nNaNs remaining:\n", df_all.isna().sum())

print("Final dataset length:", len(df_all))
print("Signal fraction:", df_all['isSignal'].mean())

X_train, X_test, Y_train, Y_test = train_test_split(
    df_all[var],
    df_all[['isSignal']],
    test_size=0.2,
)

scaler = StandardScaler()

X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled  = scaler.transform(X_test)

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

batch_size = 1024
train_generator = balanced_batch_generator(
    X_train_scaled,
    Y_train,
    batch_size=batch_size
)

# ============================================================
# MODEL 
# ============================================================

model = Sequential([
    Input(shape=(len(var),)),
    Dense(128, activation='relu', kernel_regularizer=regularizers.l2(1e-4)),
    BatchNormalization(),
    Dropout(0.3),
    Dense(64, activation='relu', kernel_regularizer=regularizers.l2(1e-4)),
    BatchNormalization(),
    Dropout(0.3),
    Dense(32, activation='relu', kernel_regularizer=regularizers.l2(1e-4)),
    Dense(8, activation='relu', kernel_regularizer=regularizers.l2(1e-4)),
    Dense(1, activation='sigmoid')
])

model.compile(
    loss='binary_crossentropy',
    optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
    metrics=[tf.keras.metrics.AUC(name="auc")]
)

early = callbacks.EarlyStopping(
    monitor='val_auc',
    mode='max',
    patience=15,
    restore_best_weights=True
)

reduce_lr = callbacks.ReduceLROnPlateau(
    monitor='val_loss',
    factor=0.5,
    patience=8,
    min_lr=1e-6,
    verbose=1
)

training = model.fit(
    train_generator,
    steps_per_epoch=len(X_train_scaled) // batch_size,
    epochs=300,
    validation_data=(X_test_scaled, Y_test),
    callbacks=[early, reduce_lr],
    verbose=2
)

# Get history
loss = training.history['loss']
val_loss = training.history['val_loss']

epochs = range(1, len(loss) + 1)

plt.figure()
plt.plot(epochs, loss)
plt.plot(epochs, val_loss)
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title('Training and Validation Loss')
plt.legend(['Train Loss', 'Validation Loss'])
plt.grid(True)
plt.yscale('log')

os.makedirs("./Plots", exist_ok=True)
plt.savefig(f'./Plots/loss_{ANALYSIS_NAME}.png', dpi=300, bbox_inches='tight')

plt.close()

if save_model and loaded_model == False:
    print("")
    print("The current used Machine Learning model is being saved:")
    print("")

    os.makedirs("./Models", exist_ok=True)
    joblib.dump(scaler, './Models/scaler_'+ANALYSIS_NAME+'.pkl')
    model.save('./Models/model_'+ANALYSIS_NAME+'.h5')

##### FOR THE DEEP NEURAL NETWORK ######
y_pred = model.predict(X_test_scaled)
y_pred_t = model.predict(X_train_scaled)

y_pred_L = y_pred
y_pred_t_L = y_pred_t

########################

#### SIGNAL category

# MODEL 1
discriminant = np.squeeze(np.asarray(y_pred_L))
true_labels = np.squeeze(np.asarray(Y_test['isSignal']))

discriminant0 = discriminant[np.array(true_labels == 0)]
discriminant1 = discriminant[np.array(true_labels == 1)]

discriminant_t = np.squeeze(np.asarray(y_pred_t_L))
true_labels_t = np.squeeze(np.asarray(Y_train['isSignal']))

discriminant0_t = discriminant_t[np.array(true_labels_t == 0)]
discriminant1_t = discriminant_t[np.array(true_labels_t == 1)]

#discriminant   = y_pred[:, 1]
#true_labels    = np.asarray(Y_test['isSignal'])
#
#discriminant0 = discriminant[true_labels == 0]
#discriminant1 = discriminant[true_labels == 1]
#
#discriminant_t = y_pred_t[:, 1]
#true_labels_t  = np.asarray(Y_train['isSignal'])
#
#discriminant0_t = discriminant_t[true_labels_t == 0]
#discriminant1_t = discriminant_t[true_labels_t == 1]

binning = np.linspace(0, 1, 51)

# Plot the discriminant distributions --------------------

print("")
print("Plotinng discriminant distribution---------------------------------------------------------------------------------------------------")
plt.clf()
plt.figure(num=None, figsize=(6, 6))
plt.subplot(111)
pdf0, bins0, patches0 = plt.hist(discriminant0, bins = binning, color = 'm', alpha = 0.0, histtype = 'stepfilled', linewidth = 1, edgecolor='r', density=True)
pdf1, bins1, patches1 = plt.hist(discriminant1, bins = binning, color = 'y', alpha = 0.0, histtype = 'stepfilled', linewidth = 1, edgecolor='b', density=True)

pdf0_t, bins0_t, patches0_t = plt.hist(discriminant0_t, bins = binning, color = 'r', alpha = 0.3, histtype = 'stepfilled', linewidth = 2, edgecolor='r', label = 'Backgrounds (train)', density=True)
pdf1_t, bins1_t, patches1_t = plt.hist(discriminant1_t, bins = binning, color = 'b', alpha = 0.3, histtype = 'stepfilled', linewidth = 2, edgecolor='b', label = 'ttDM signal (train)', density=True)

plt.scatter(bins0[:-1]+ 0.5*(bins0[1:] - bins0[:-1]), pdf0, marker='.', c='m', s=30, alpha=0.8, label = 'Backgrounds')
plt.scatter(bins1[:-1]+ 0.5*(bins1[1:] - bins1[:-1]), pdf1, marker='.', c='y', s=30, alpha=0.8, label = 'ttDM signal')

plt.legend(loc = 'upper center')
plt.ylabel('Density', fontsize = 12)
plt.xlabel('DNN discriminant', fontsize = 12)
plt.savefig('./Plots/Discriminant_distribution'+ANALYSIS_NAME+'.png', dpi = 600)

plt.clf()
plt.figure(num=None, figsize=(6, 6))
plt.subplot(111)

plt.clf()
plt.figure(num=None, figsize=(6, 6))
plt.subplot(111)
pdf0, bins0, patches0 = plt.hist(discriminant0, bins = binning, color = 'm', alpha = 0.0, histtype = 'stepfilled', linewidth = 1, edgecolor='r', density=True)
pdf1, bins1, patches1 = plt.hist(discriminant1, bins = binning, color = 'y', alpha = 0.0, histtype = 'stepfilled', linewidth = 1, edgecolor='b', density=True)

pdf0_t, bins0_t, patches0_t = plt.hist(discriminant0_t, bins = binning, color = 'r', alpha = 0.3, histtype = 'stepfilled', linewidth = 2, edgecolor='r', label = 'Backgrounds (train)', density=True)
pdf1_t, bins1_t, patches1_t = plt.hist(discriminant1_t, bins = binning, color = 'b', alpha = 0.3, histtype = 'stepfilled', linewidth = 2, edgecolor='b', label = 'ttDM signal (train)', density=True)

plt.scatter(bins0[:-1]+ 0.5*(bins0[1:] - bins0[:-1]), pdf0, marker='.', c='m', s=30, alpha=0.8, label = 'Backgrounds')
plt.scatter(bins1[:-1]+ 0.5*(bins1[1:] - bins1[:-1]), pdf1, marker='.', c='y', s=30, alpha=0.8, label = 'ttDM signal')
plt.legend(loc = 'upper center')
plt.yscale('log')
plt.ylabel('Density', fontsize = 12)
plt.xlabel('DNN discriminant (Signal)', fontsize = 12)
plt.savefig('./Plots/Log_Discriminant_distribution'+ANALYSIS_NAME+'.png', dpi = 600)

print("DONE!")
##### --------------------------------
print("")
print("Plotting ROC-------------------------------------------------------------------------------------------------------------------------")
fpr, tpr, thresholds = metrics.roc_curve(Y_test["isSignal"], y_pred_L)
#fpr, tpr, thresholds = metrics.roc_curve(Y_test["isSignal"], y_pred_L[:, 1])
auc = metrics.auc(fpr, tpr)

plt.clf()
plt.figure(num=None, figsize=(6, 6))
plt.subplot(111)
plt.plot(fpr, tpr, color = 'r', label = "ROC curve")
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label = "Random guess")
plt.legend(loc = "lower right")
plt.xlabel('False Positive rate', fontsize = 12)
plt.ylabel('True Positive rate', fontsize = 12)
plt.text(0.55, 0.3, f'AUC = {auc:.3f}', fontsize=12, color='r')
plt.axvline(x=0, color = 'black', linestyle = '--', linewidth = 0.5)
plt.axhline(y=1, color = 'black', linestyle = '--', linewidth = 0.5)
plt.savefig('./Plots/ROC'+ANALYSIS_NAME+'.png', dpi = 600)
print("")
print("The AUC of the model is: ", auc)

if PARAMETRIC and "mPhi" in X_test.columns:
    print("")
    print("Plotting individual ROCs for different mPhi hypotheses-----------------------------------------------------------------------")

    roc_df = pd.DataFrame({
        "isSignal": np.asarray(Y_test["isSignal"]).reshape(-1),
        "score": np.asarray(y_pred_L).reshape(-1),
        "mPhi": np.asarray(X_test["mPhi"]).reshape(-1),
    })

    roc_curves_by_mass = []
    unique_masses = sorted(roc_df["mPhi"].dropna().unique())

    for mass in unique_masses:
        mass_slice = roc_df[roc_df["mPhi"] == mass]

        # Need both classes to compute ROC/AUC.
        if mass_slice["isSignal"].nunique() < 2:
            print(f"Skipping mPhi={mass}: only one class present in test split")
            continue

        fpr_mass, tpr_mass, _ = metrics.roc_curve(mass_slice["isSignal"], mass_slice["score"])
        auc_mass = metrics.auc(fpr_mass, tpr_mass)
        roc_curves_by_mass.append((mass, fpr_mass, tpr_mass, auc_mass))

    if roc_curves_by_mass:
        plt.clf()
        plt.figure(num=None, figsize=(7, 6))
        plt.subplot(111)

        color_map = plt.cm.get_cmap('tab20', len(roc_curves_by_mass))

        print("AUC per mPhi hypothesis:")
        for i, (mass, fpr_mass, tpr_mass, auc_mass) in enumerate(roc_curves_by_mass):
            print(f"  mPhi={mass}: AUC={auc_mass:.4f}")
            plt.plot(fpr_mass, tpr_mass, color=color_map(i), label=rf"m_\Phi={mass:g} (AUC={auc_mass:.3f})")

        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random guess')
        plt.xlabel('False Positive rate', fontsize=12)
        plt.ylabel('True Positive rate', fontsize=12)
        plt.legend(loc='lower right', fontsize=9)
        plt.axvline(x=0, color='black', linestyle='--', linewidth=0.5)
        plt.axhline(y=1, color='black', linestyle='--', linewidth=0.5)
        plt.savefig('./Plots/ROC'+ANALYSIS_NAME+'_mPhi.png', dpi=600)
    else:
        print("No per-mPhi ROC could be produced (insufficient class mixture per mass in test split).")

print("DONE!")

##### Feature importance
print("")
print("Plotting feature importance-----------------------------------------------------------------------------------------------------------")

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

# -----------------------------
# Feature Importance for Keras DNN
# -----------------------------

X_val = np.array(X_test_scaled)
y_val = np.array(Y_test["isSignal"])
features = var  # list of feature names

# --- 1) Manual Permutation Importance ---
print("Computing permutation importance (manual, works with Keras)...")
baseline_auc = roc_auc_score(y_val, model.predict(X_val))#[:, 1])

importances = []
for i in range(X_val.shape[1]):
    X_permuted = X_val.copy()
    np.random.shuffle(X_permuted[:, i])
    permuted_auc = roc_auc_score(y_val, model.predict(X_permuted))#[:, 1])
    importances.append(baseline_auc - permuted_auc)

importances = np.array(importances)
idx = np.argsort(importances)[::-1]

plt.clf()
plt.barh(np.array(features)[idx], importances[idx])
plt.xlabel("Permutation Importance (Δ AUC)")
plt.title("Feature Importance (Manual, Keras DNN)")
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig("./Plots/FeatureImportance_Keras_Permutation.png", dpi=600)
print("Permutation importance plot saved.")

# -----------------------------
# Optional: SHAP feature importance
# -----------------------------
try:
    import shap
    HAS_SHAP = True
except ImportError:
    HAS_SHAP = False
    print("SHAP not installed — skipping SHAP plots.")

if HAS_SHAP:
    print("Computing SHAP feature importance (subsampled for speed)...")

    # Subsample to avoid huge memory
    max_background = 200
    if len(X_val) > max_background:
        idx_sample = np.random.choice(len(X_val), max_background, replace=False)
        X_background = X_val[idx_sample]
    else:
        X_background = X_val

    # Create SHAP explainer
    explainer = shap.Explainer(model, X_background)
    shap_values = explainer(X_background)

    # Beeswarm plot
    plt.clf()
    shap.summary_plot(shap_values, X_background, feature_names=features, show=False)
    plt.title('Feature Importance - SHAP Values (Beeswarm)', fontsize=16)
    plt.tight_layout()
    plt.savefig('./Plots/SHAP_Feature_importance_Beeswarm.png', dpi=120)

    # Bar plot
    plt.clf()
    shap.summary_plot(shap_values, X_background, plot_type="bar", feature_names=features, show=False)
    plt.title('Feature Importance - SHAP Values (Bar)', fontsize=16)
    plt.tight_layout()
    plt.savefig('./Plots/SHAP_Feature_importance_Bar.png', dpi=120)

    # Log-scale bar plot
    plt.clf()
    shap.summary_plot(shap_values, X_background, plot_type="bar", feature_names=features, show=False)
    plt.xscale('log')
    plt.title('Feature Importance - SHAP Values (Log Bar)', fontsize=16)
    plt.tight_layout()
    plt.savefig('./Plots/SHAP_Feature_importance_LogBar.png', dpi=120)

    print("SHAP plots saved.")

print("DONE!")

print("")
print("Finding and ploting correlation matrix-----------------------------------------------------------------------------------------------")
import matplotlib.cm as cm
m = np.corrcoef(X_train_scaled, rowvar=False) # Correlation matrix with numpy

# Name of variable in order
tickets = var
## PLOTING
fig, ax = plt.subplots(figsize=(20, 20))
im = ax.matshow(abs(m))
plt.colorbar(im)
for (i, j), z in np.ndenumerate(m):
    ax.text(j, i, '{:0.1f}'.format(abs(z)), ha='center', va='center')

ax.set_xticks(np.arange(len(tickets)))
ax.set_yticks(np.arange(len(tickets)))
ax.set_xticklabels(tickets)
ax.set_yticklabels(tickets)
plt.setp(ax.get_xticklabels(), rotation=-45, ha="right", rotation_mode="anchor")
ax.set_title("Correlation matrix")
fig.tight_layout()
plt.savefig('./Plots/CorrelationMatrix'+ANALYSIS_NAME+'.png', dpi = 60)

print("DONE!")

print("")
print("Finding and ploting correlation matrix for the signal events-------------------------------------------------------------------------")
import matplotlib.cm as cm
m = np.corrcoef(Sig[var], rowvar=False) # Correlation matrix with numpy

# Name of variable in order
tickets = var
## PLOTING
fig, ax = plt.subplots(figsize=(20, 20))
im = ax.matshow(abs(m))
plt.colorbar(im)
for (i, j), z in np.ndenumerate(m):
    ax.text(j, i, '{:0.1f}'.format(abs(z)), ha='center', va='center')

ax.set_xticks(np.arange(len(tickets)))
ax.set_yticks(np.arange(len(tickets)))
ax.set_xticklabels(tickets)
ax.set_yticklabels(tickets)
plt.setp(ax.get_xticklabels(), rotation=-45, ha="right", rotation_mode="anchor")
ax.set_title("Correlation matrix")
fig.tight_layout()
plt.savefig('./Plots/CorrelationMatrix'+ANALYSIS_NAME+'Sig.png', dpi = 60)

print("DONE!")
print("")
print("Finding and ploting correlation matrix for the background events---------------------------------------------------------------------")
import matplotlib.cm as cm
m = np.corrcoef(Bkg[var], rowvar=False) # Correlation matrix with numpy

# Name of variable in order
tickets = var
## PLOTING
fig, ax = plt.subplots(figsize=(20, 20))
im = ax.matshow(abs(m))
plt.colorbar(im)
for (i, j), z in np.ndenumerate(m):
    ax.text(j, i, '{:0.1f}'.format(abs(z)), ha='center', va='center')

ax.set_xticks(np.arange(len(tickets)))
ax.set_yticks(np.arange(len(tickets)))
ax.set_xticklabels(tickets)
ax.set_yticklabels(tickets)
plt.setp(ax.get_xticklabels(), rotation=-45, ha="right", rotation_mode="anchor")
ax.set_title("Correlation matrix")
fig.tight_layout()
plt.savefig('./Plots/CorrelationMatrix'+ANALYSIS_NAME+'Bkg.png', dpi = 60)

print("DONE!")
print("")

print("END OF THE ANALYSIS")
