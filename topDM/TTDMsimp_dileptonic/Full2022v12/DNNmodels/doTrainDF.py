import ROOT
import uproot
import pandas as pd
import numpy as np
import random
import subprocess
import joblib

import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc, classification_report, roc_auc_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier

import xgboost as xgb

import tensorflow.keras
#from keras.utils import np_utils
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout, Activation, LeakyReLU
from tensorflow.keras import regularizers, backend as K, optimizers, callbacks
from tensorflow.keras.optimizers import Adam

from sklearn import metrics

ROOT.EnableImplicitMT()

####### CONTROL PARAMETERS FOR THE WHOLE CODE #######

ANALYSIS_NAME = "dnn_model"
save_model = True

loaded_model = False 

MODEL_NAME = "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022v12/DNNmodels/model_DNN.h5"

#####################################################

var = [
  'dphill',
  'PuppiMET_pt',
  'mT2',
  'pdark',
  'chel',
  'dphi_ttbar',
  'dphi_met_llb',
        ]

print("")
print("Starting to load the events for the analysis------------------------------------------------------------------------------------------")
print("")

files_ttDM_temp = {}

from mkShapesRDF.lib.search_files import SearchFiles

searchFiles = SearchFiles()

limitFiles = -1
redirector = ""

def nanoGetSampleFiles(path, name):
    _files = searchFiles.searchFiles(path, name, redirector=redirector)
    if limitFiles != -1:
        _files = _files[:limitFiles]
    return _files

signalDirectory = '/eos/user/v/victorr/HWWNano/Summer22_130x_nAODv12_Full2022v12/MCl2loose2022v12__MCCorr2022v12JetScaling__l2tight'

mcDirectory = '/eos/cms/store/group/phys_higgs/cmshww/calderon/HWWNano/Summer22_130x_nAODv12_Full2022v12/MCl2loose2022v12__MCCorr2022v12JetScaling__sblancof__l2tight'

mPhi = ['10','50','100','150', '200', '250', '300', '350', '400', '500', '600' '700', '800', '1000']

# tt+DM dilepton scalar
for phi in mPhi:
    files_ttDM_temp[f'TTto2LDMsimpSpin0_s_mphi_{phi}'] = {
            'names': nanoGetSampleFiles(signalDirectory, f'TTto2LDMsimpSpin0_s_mphi_{phi}'),
            }

# tt+DM inclusive scalar
for phi in mPhi:
    files_ttDM_temp[f'TTDMsimpSpin0_s_mphi_{phi}'] = {
            'names': nanoGetSampleFiles(signalDirectory, f'TTDMsimpSpin0-mphi-{phi}'),
            }

# tt+DM dilepton pseudoscalar
for phi in mPhi:
    files_ttDM_temp[f'TTto2LDMsimpSpin0_ps_mphi_{phi}'] = {
            'names': nanoGetSampleFiles(signalDirectory, f'TTto2LDMsimpSpin0_s_mphi_{phi}'),
            }

# tt+DM inclusive pseudoscalar
for phi in mPhi:
    files_ttDM_temp[f'TTDMsimpSpin0_ps_mphi_{phi}'] = {
            'names': nanoGetSampleFiles(signalDirectory, f'TTDMsimpSpin0_s_mphi_{phi}'),
            }

files_DH = files_ttDM_temp.copy()
for key, value in files_ttDM_temp.items():
    if isinstance(value["names"], list) and len(value["names"]) == 0:
        del files_DH[key]

print("")
print("Files for ttDM loaded----------------------------------------------------------------------------------------------------------")
print("")

# DATA FROM RUN 3
#files_WW = nanoGetSampleFiles(mcDirectory, 'WWTo2L2Nu')

files_top = nanoGetSampleFiles(mcDirectory, 'TTTo2L2Nu') + \
    nanoGetSampleFiles(mcDirectory,'ST_t-channel_top') +  \
    nanoGetSampleFiles(mcDirectory,'ST_t-channel_antitop')

files_BKG = files_top # files_WW 
print("")
print("Background files processed-----------------------------------------------------------------------------------------------------------")
print("")

dataframes = {}
for key, df in files_DH.items():
    dataframes[key] = ROOT.RDataFrame("Events", files_DH[key]["names"])

df_bkg = ROOT.RDataFrame("Events", files_BKG)

print("")
print("Dataframes created-------------------------------------------------------------------------------------------------------------------")
print("")

ROOT.gROOT.ProcessLine(
    """
    template
    <typename container>
    float Alt(container c, int index, float alt){
        if (index < c.size()) {
            return c[index];
        }
        else{
            return alt;
        }
    }
    """
)

##################################################
########### ttDM relevant variables ##############
##################################################

ROOT.gInterpreter.Declare(
"""
    // FOR MT2 CALCULATION
    #include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/computeMT2.cc"
    // FOR JETS IN HORNS
    #include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/jet_horns.cc"

    // FOR NUNU RECONSTRUCTION
    #include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/doubleNu_producer.h"
    
    #include <cmath>
    #include <cstdlib>
    #include <iostream>
    #include "TLorentzVector.h"
    #include "TMath.h"
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

print("")
print("New variables already declared for the use within the dataframes-----------------------------------------------------------------------------")
print("")

bWP = '0.0583'
bAlgo = 'DeepFlavB'

dataframes[key] = dataframes[key].Define(
    "doubleNu",
    "doubleNu_producer(nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, PuppiMET_pt, PuppiMET_phi, Jet_btag{0}, {1})".format(bAlgo, bWP)
    )

for key, df in dataframes.items():
    dataframes[key] = dataframes[key].Define("mT2", 'computeMT2(Lepton_pt[0], Lepton_eta[0], Lepton_phi[0], Lepton_pt[1], Lepton_eta[1], Lepton_phi[1], PuppiMET_pt, PuppiMET_phi)')
    dataframes[key] = dataframes[key].Define("pdark", "doubleNu[8]")
    dataframes[key] = dataframes[key].Define("chel",  "doubleNu[6]")
    dataframes[key] = dataframes[key].Define("dphi_ttbar", "doubleNu[7]")
    dataframes[key] = dataframes[key].Define("tt_reco", "doubleNu[9]")
    dataframes[key] = dataframes[key].Define("bjet_idx", '(nCleanJet > 0 && Jet_btag{algo}[CleanJet_jetIdx[0]] > {wp}) ? 0 : (nCleanJet > 1 && Jet_btag{algo}[CleanJet_jetIdx[1]] > {wp}) ? 1 : (nCleanJet > 2 && Jet_btag{algo}[CleanJet_jetIdx[2]] > {wp}) ? 2 : -1'.format(algo=bAlgo, wp=bWP))
    dataframes[key] = dataframes[key].Define(
    "dphi_met_llb",
    "abs(DeltaPhi(PuppiMET_phi, atan2(Lepton_pt[0]*sin(Lepton_phi[0]) + Lepton_pt[1]*sin(Lepton_phi[1]) + CleanJet_pt[bjet_idx]*sin(CleanJet_phi[bjet_idx]), Lepton_pt[0]*cos(Lepton_phi[0]) + Lepton_pt[1]*cos(Lepton_phi[1]) + CleanJet_pt[bjet_idx]*cos(CleanJet_phi[bjet_idx]))))"
    )
    dataframes[key] = dataframes[key].Define("first_btag_ID", 'Alt(Jet_btag{}, Alt(CleanJet_jetIdx, 0, -1), -1)'.format(bAlgo))
    dataframes[key] = dataframes[key].Define("second_btag_ID", 'Alt(Jet_btag{}, Alt(CleanJet_jetIdx, 1, -1), -1)'.format(bAlgo))
    dataframes[key] = dataframes[key].Define("lep_eta1", "Lepton_eta[0]")
    dataframes[key] = dataframes[key].Define("lep_eta2", "Lepton_eta[1]")
    dataframes[key] = dataframes[key].Define("lep_pt1", "Lepton_pt[0]")
    dataframes[key] = dataframes[key].Define("lep_pt2", "Lepton_pt[1]")
    dataframes[key] = dataframes[key].Define("noJetInHorn", "Jet_inHorns(CleanJet_pt, CleanJet_eta)")
    dataframes[key] = dataframes[key].Define("bReq", 'Sum(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{}, CleanJet_jetIdx) > {}) >= 1'.format(bAlgo, bWP))    

df_bkg = df_bkg.Define("mT2", 'computeMT2(Lepton_pt[0], Lepton_eta[0], Lepton_phi[0], Lepton_pt[1], Lepton_eta[1], Lepton_phi[1], PuppiMET_pt, PuppiMET_phi)')
df_bkg = df_bkg.Define(
    "doubleNu",
    "doubleNu_producer(nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, PuppiMET_pt, PuppiMET_phi, Jet_btag{0}, {1})".format(bAlgo, bWP)
    )
df_bkg = df_bkg.Define("pdark", "doubleNu[8]")
df_bkg = df_bkg.Define("chel", "doubleNu[6]")
df_bkg = df_bkg.Define("dphi_ttbar", "doubleNu[7]")
df_bkg = df_bkg.Define("tt_reco", "doubleNu[9]")
df_bkg = df_bkg.Define("bjet_idx", '(nCleanJet > 0 && Jet_btag{algo}[CleanJet_jetIdx[0]] > {wp}) ? 0 : (nCleanJet > 1 && Jet_btag{algo}[CleanJet_jetIdx[1]] > {wp}) ? 1 : (nCleanJet > 2 && Jet_btag{algo}[CleanJet_jetIdx[2]] > {wp}) ? 2 : -1'.format(algo=bAlgo, wp=bWP))
df_bkg = df_bkg.Define(
    "dphi_met_llb",
    "abs(DeltaPhi(PuppiMET_phi, atan2(Lepton_pt[0]*sin(Lepton_phi[0]) + Lepton_pt[1]*sin(Lepton_phi[1]) + CleanJet_pt[bjet_idx]*sin(CleanJet_phi[bjet_idx]), Lepton_pt[0]*cos(Lepton_phi[0]) + Lepton_pt[1]*cos(Lepton_phi[1]) + CleanJet_pt[bjet_idx]*cos(CleanJet_phi[bjet_idx]))))"
)
df_bkg = df_bkg.Define("first_btag_ID", 'Alt(Jet_btag{}, Alt(CleanJet_jetIdx, 0, -1), -1)'.format(bAlgo))
df_bkg = df_bkg.Define("second_btag_ID", 'Alt(Jet_btag{}, Alt(CleanJet_jetIdx, 1, -1), -1)'.format(bAlgo))
df_bkg = df_bkg.Define("lep_eta1", "Lepton_eta[0]")
df_bkg = df_bkg.Define("lep_eta2", "Lepton_eta[1]")
df_bkg = df_bkg.Define("lep_pt1", "Lepton_pt[0]")
df_bkg = df_bkg.Define("lep_pt2", "Lepton_pt[1]")
df_bkg = df_bkg.Define("noJetInHorn", "Jet_inHorns(CleanJet_pt, CleanJet_eta)")
df_bkg = df_bkg.Define("bReq", 'Sum(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{}, CleanJet_jetIdx) > {}) >= 1'.format(bAlgo, bWP))

print("")
print("New variables loaded within the dataframes-------------------------------------------------------------------------------------------")
print("")

# -------------------------------
# FILTERING 
# -------------------------------
for key, df in dataframes.items():
    dataframes[key] = dataframes[key].Filter(
        "((abs(Lepton_pdgId[0]) == 11 || abs(Lepton_pdgId[0]) == 13) && "
        "(abs(Lepton_pdgId[1]) == 11 || abs(Lepton_pdgId[1]) == 13)) && "
        "Lepton_pt[0]>25 && Lepton_pt[1]>20 && Alt(Lepton_pt,2, 0)<10 && "
        "abs(Lepton_eta[0]) < 2.4 && abs(Lepton_eta[1]) < 2.4 && "
        "mll > 20 && noJetInHorn && bReq && tt_reco"
    )

df_bkg = df_bkg.Filter(
    "((abs(Lepton_pdgId[0]) == 11 || abs(Lepton_pdgId[0]) == 13) && "
    "(abs(Lepton_pdgId[1]) == 11 || abs(Lepton_pdgId[1]) == 13)) && "
    "Lepton_pt[0]>25 && Lepton_pt[1]>20 && Alt(Lepton_pt,2, 0)<10 && "
    "abs(Lepton_eta[0]) < 2.4 && abs(Lepton_eta[1]) < 2.4 && "
    "mll > 20 && noJetInHorn && bReq && tt_reco"
)

print("")
print("Events filtering---------------------------------------------------------------------------------------------------------------------")

columns = [
    'dphill',
    'PuppiMET_pt',
    'mT2',
    'pdark',
    'chel',
    'dphi_ttbar',
    'dphi_met_llb',
]

print("DONE!")
print("")

# -------------------------------
# DOWNSAMPLING
# -------------------------------
sig_counts = {key: df.Count() for key, df in dataframes.items()}
bkg_count = df_bkg.Count()

# Trigger event loop once
sig_counts = {key: c.GetValue() for key, c in sig_counts.items()}
bkg_count = bkg_count.GetValue()

print("Signal counts:", sig_counts)
print("Background count:", bkg_count)

target_size = min(min(sig_counts.values()), bkg_count)
print("Target size:", target_size)

import ROOT

def downsample(df, current, target):
    if current <= target:
        return df
    frac = target / current
    return (
        df
        .Define("rnd", "gRandom->Rndm()")
        .Filter(f"rnd < {frac}")
    )

for key, df in dataframes.items():
    dataframes[key] = downsample(df, sig_counts[key], target_size)

df_bkg = downsample(df_bkg, bkg_count, target_size)

# -------------------------------
# CONVERSION TO PANDAS
# -------------------------------
print("Creating pandas dataframes!----------------------------------------------------------------------------------------------------------")

numpy_dataframes = {}
for key, df in dataframes.items():
    numpy_dataframes[key] = df.AsNumpy(var)

dfBkg = df_bkg.AsNumpy(var)

pd_dataframes = {}
for key, df_sig in numpy_dataframes.items():
    pd_dataframes[key] = pd.DataFrame(df_sig)

Bkg = pd.DataFrame(dfBkg)

print("DONE!")
print("")

# -------------------------------
# CATEGORIES 
# -------------------------------
for key, df_sig in pd_dataframes.items():
    pd_dataframes[key]['isSignal'] = np.ones(len(df_sig))
    pd_dataframes[key]['isBkg'] = np.zeros(len(df_sig))

Bkg['isSignal'] = np.zeros(len(Bkg))
Bkg['isBkg'] = np.ones(len(Bkg))

print("")
print("Categories inserted in the dataframes for isSignal and isBkg-------------------------------------------------------------------------")
print("")

# -------------------------------
# CONCATENATION 
# -------------------------------
print([len(df_sig) for key, df_sig in pd_dataframes.items()])

Sig = pd.concat([df_sig for key, df_sig in pd_dataframes.items()])

print("")
print("Seeing how each of the dataframes for Signal and Background look:")
print(Sig)
print(Bkg)

print("Initial lenght of Sig:", len(Sig))
print("Initial lenght of Bkg:", len(Bkg))

# Sample the background to the length of the concatenated signal dataframe
if len(Sig) > len(Bkg):
    Sig = Sig.sample(len(Bkg))
else:
    Bkg = Bkg.sample(len(Sig))

print("")
print("Statistics after sampling")  
print("Statistics for Sig: " + str(len(Sig)))
print("Statistics for Bkg: " + str(len(Bkg)))
print("")

# Concatenate all DataFrames in concatenated_dfs with the 'Bkg' DataFrame
df_all = pd.concat([Sig, Bkg])

print("Length of dataset: " + str(len(df_all)))
df_all.dropna(inplace=True)
print("After removing NANs: " + str(len(df_all)))

X_train, X_test, Y_train, Y_test = train_test_split(df_all[var], df_all[['isSignal']], test_size=0.2, random_state=6)

print("")
print("Variables to study for the analysis:")
print(X_train)
print(Y_train)
print("")
print("DONE!")

ANALYSIS_NAME += "_DNN_"

if loaded_model:
    model= load_model(MODEL_NAME)
else:
        print("Start training the Deep Neural Network!----------------------------------------------------------------------------------------------------")
        print("")

        # MODEL 1
        model = Sequential()
        # Hidden layers 1 
        model.add(Dense(128, activation='relu', input_dim=len(var)))
        model.add(Dense(64, activation='relu'))
        model.add(Dense(32, activation='relu'))
        model.add(Dense(8, activation='relu'))
        # Output layer
        model.add(Dense(1, activation='sigmoid'))
        # Compile
        model.compile(loss='binary_crossentropy', optimizer=optimizers.RMSprop(0.00015), metrics=['accuracy'])
        
        n_epochs = 300
        n_batch = 128
        training = model.fit(X_train[var].values, Y_train, epochs=n_epochs, validation_split=0.15, batch_size=n_batch,
                                     callbacks = [callbacks.EarlyStopping(monitor='val_loss', patience=20, verbose=1)],
                                     verbose=2, shuffle= True)

#        # MODEL 2
#        model = Sequential()
#        # Hidden layers 2
#        model.add(Dense(40, input_dim=len(var)))
#        model.add(LeakyReLU(alpha=0.01))
#        model.add(Dense(40))
#        model.add(LeakyReLU(alpha=0.01))
#        model.add(Dense(30))
#        model.add(LeakyReLU(alpha=0.01))
#        # Output layer
#        model.add(Dense(2, activation='softmax'))
#
#        model.compile(optimizer=Adam(), loss='sparse_categorical_crossentropy', metrics=['accuracy'])
#
#        n_epochs = 50
#        n_batch = 100
#        training = model.fit(X_train[var].values, Y_train.values, epochs=n_epochs, batch_size=n_batch, validation_split=0.15, 
#                shuffle=True, verbose=2, callbacks=[callbacks.EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)])

        

print("DONE!")
print("Save and plot test distributions-----------------------------------------------------------------------------------------------------")

if save_model and loaded_model == False:
    print("")
    print("The current used Machine Learning model is being saved:")
    print("")
    
    model.save('./Models/model_'+ANALYSIS_NAME+'.h5')



##### FOR THE DEEP NEURAL NETWORK ######
y_pred = model.predict(X_test)
y_pred_t = model.predict(X_train)

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


# Plot the discriminant distributions ----------------------------

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
plt.savefig('Discriminant_distribution'+ANALYSIS_NAME+'.png', dpi = 600)

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
plt.savefig('Log_Discriminant_distribution'+ANALYSIS_NAME+'.png', dpi = 600)
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
plt.savefig('ROC'+ANALYSIS_NAME+'.png', dpi = 600)
print("")
print("The AUC of the model is: ", auc)

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

X_val = np.array(X_test)
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
plt.savefig("FeatureImportance_Keras_Permutation.png", dpi=600)
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
    plt.savefig('SHAP_Feature_importance_Beeswarm.png', dpi=600)

    # Bar plot
    plt.clf()
    shap.summary_plot(shap_values, X_background, plot_type="bar", feature_names=features, show=False)
    plt.title('Feature Importance - SHAP Values (Bar)', fontsize=16)
    plt.tight_layout()
    plt.savefig('SHAP_Feature_importance_Bar.png', dpi=600)

    # Log-scale bar plot
    plt.clf()
    shap.summary_plot(shap_values, X_background, plot_type="bar", feature_names=features, show=False)
    plt.xscale('log')
    plt.title('Feature Importance - SHAP Values (Log Bar)', fontsize=16)
    plt.tight_layout()
    plt.savefig('SHAP_Feature_importance_LogBar.png', dpi=1200)

    print("SHAP plots saved.")

print("DONE!")

print("")
print("Finding and ploting correlation matrix-----------------------------------------------------------------------------------------------")
import matplotlib.cm as cm
m = np.corrcoef(X_train, rowvar=False) # Correlation matrix with numpy

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
plt.savefig('CorrelationMatrix'+ANALYSIS_NAME+'.png', dpi = 600)

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
plt.savefig('CorrelationMatrix'+ANALYSIS_NAME+'Sig.png', dpi = 600)

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
plt.savefig('CorrelationMatrix'+ANALYSIS_NAME+'Bkg.png', dpi = 600)

print("DONE!")
print("")

print("END OF THE ANALYSIS------------------------------------------------------------------------------------------------------------------")
