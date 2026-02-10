import os
import copy
import inspect

configurations = os.path.realpath(inspect.getfile(inspect.currentframe())) # this file

aliases = {}
aliases = OrderedDict()

mc     = [skey for skey in samples if skey not in ('Fake', 'DATA', 'Dyemb', 'DATA_EG', 'DATA_Mu', 'DATA_EMu', 'Fake_EG', 'Fake_Mu', 'Fake_EMu')]
mc_emb = [skey for skey in samples if skey not in ('Fake', 'DATA', 'DATA_Mu', 'DATA_EMu', 'Fake_EG', 'Fake_Mu', 'Fake_EMu')]

# LepSF2l__ele_cutBased_LooseID_tthMVA_Run3__mu_cut_TightID_pfIsoLoose_HWW_tthmva_67
eleWP = 'cutBased_LooseID_tthMVA_Run3'
muWP  = 'cut_TightID_pfIsoLoose_HWW_tthmva_67'

aliases['LepWPCut'] = {
    'expr': 'LepCut2l__ele_'+eleWP+'__mu_'+muWP,
    'samples': mc + ['DATA'],
}

aliases['LepWPSF'] = {
    'expr': 'LepSF2l__ele_'+eleWP+'__mu_'+muWP,
    'samples': mc
}

aliases['CleanestJet_mask'] = {
    #'expr': '(CleanestJet_pt>50 || !(abs(CleanestJet_eta)>2.6 && abs(CleanestJet_eta)<3.1))'
    'expr': '!(CleanJet_pt < 50 && abs(CleanJet_eta) > 2.6 && abs(CleanJet_eta) < 3.1)'
}
aliases['CleanestJet_pt'] = {
    'expr': 'CleanJet_pt[CleanestJet_mask]'
}
aliases['CleanestJet_eta'] = {
    'expr': 'CleanJet_eta[CleanestJet_mask]'
}
aliases['CleanestJet_phi'] = {
    'expr': 'CleanJet_pt[CleanestJet_mask]'
}
aliases['CleanestJet_mass'] = {
    'expr': 'CleanJet_mass[CleanestJet_mask]'
}
aliases['CleanestJet_jetIdx'] = {
    'expr': 'CleanJet_jetIdx[CleanestJet_mask]'
}
aliases['nCleanestJet'] = {
    'expr': 'Sum(CleanestJet_mask)'
}
aliases['cleanestPuppiMET'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/cleanestMET.cc"'],
    'expr' : 'cleanestMET(CleanJet_pt,CleanJet_eta,CleanJet_phi,Jet_mass,CleanJet_jetIdx,CleanestJet_mask,PuppiMET_pt,PuppiMET_phi,Lepton_pt,Lepton_eta,Lepton_phi)'
}
aliases['cleanestPuppiMET_pt'] = {
    'expr' : 'cleanestPuppiMET[0]'
}
aliases['cleanestPuppiMET_phi'] = {
    'expr' : 'cleanestPuppiMET[1]'
}
aliases['cleanest_dphilmet'] = {
    'expr' : 'abs(DeltaPhi(Lepton_phi[0], cleanestPuppiMET_phi)) < abs(DeltaPhi(Lepton_phi[1], cleanestPuppiMET_phi)) ? abs(DeltaPhi(Lepton_phi[0], cleanestPuppiMET_phi)) : abs(DeltaPhi(Lepton_phi[1], cleanestPuppiMET_phi))'
}
aliases['cleanest_dphilmet1'] = {
    'expr' : 'abs(DeltaPhi(Lepton_phi[0], cleanestPuppiMET_phi))'
}
aliases['cleanest_dphilmet2'] = {
    'expr' : 'abs(DeltaPhi(Lepton_phi[1], cleanestPuppiMET_phi))'
}

aliases['cleanest_dphillmet'] = {
    'expr' : 'cleanestPuppiMET[2]'
}

aliases['cleanest_projpfmet'] = {
    'expr' : 'cleanest_dphilmet < 0.5*TMath::Pi() ? sin(cleanest_dphilmet) * cleanestPuppiMET_pt : cleanestPuppiMET_pt'
}
aliases['cleanest_mth'] = {
    'expr': 'sqrt(2. * ptll * cleanestPuppiMET_pt * (1. - cos(cleanest_dphillmet)))'
}
aliases['cleanest_mtw1'] = {
    'expr': 'sqrt(2. * Lepton_pt[0] * cleanestPuppiMET_pt * (1. - cos(cleanest_dphilmet1)))'
}
aliases['cleanest_mtw2'] = {
    'expr': 'sqrt(2. * Lepton_pt[0] * cleanestPuppiMET_pt * (1. - cos(cleanest_dphilmet2)))'
}

##### New jet variables
aliases['cleanestKinematics'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/cleanestKinematics.cc"'],
    'expr' : 'cleanestKinematics(CleanestJet_pt,CleanestJet_eta,CleanestJet_phi,Jet_mass,CleanestJet_jetIdx,cleanestPuppiMET_pt,cleanestPuppiMET_phi)'
}
aliases['cleanest_mjj'] = {
    'expr': 'cleanestKinematics[0]'
}
aliases['cleanest_dphijj'] = {
    'expr': 'cleanestKinematics[1]'
}
aliases['cleanest_drjj'] = {
    'expr': 'cleanestKinematics[2]'
}
aliases['cleanest_detajj'] = {
    'expr': 'cleanestKinematics[3]'
}
aliases['cleanest_dphijjmet'] = {
    'expr': 'cleanestKinematics[4]'
}

# gen-matching to prompt only (GenLepMatch2l matches to *any* gen lepton)
aliases['PromptGenLepMatch2l'] = {
    'expr': 'Alt(Lepton_promptgenmatched, 0, 0) * Alt(Lepton_promptgenmatched, 1, 0)',
    'samples': mc
}


# Fake leptons transfer factor --------------------------------------
aliases['fakeW'] = {
    'linesToAdd' : [f'#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.Declare('fake_rate_reader fr_reader = fake_rate_reader(\"2022EE\", \"Run3\", \"67\", \"nominal\", 2, \"std\");')"],
    'expr': 'fr_reader(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanestJet_pt, nCleanestJet)',
    'samples'    : ['Fake']
}

# And variations - already divided by central values in formulas !
aliases['fakeWEleUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.Declare('fake_rate_reader fr_reader_EleUp = fake_rate_reader(\"2022EE\", \"Run3\", \"67\", \"EleUp\", 2, \"std\");')"],
    'expr': 'fr_reader_EleUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanestJet_pt, nCleanestJet)',
    'samples': ['Fake']
}
aliases['fakeWEleDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.Declare('fake_rate_reader fr_reader_EleDown = fake_rate_reader(\"2022EE\", \"Run3\", \"67\", \"EleDown\", 2, \"std\");')"],
    'expr': 'fr_reader_EleDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanestJet_pt, nCleanestJet)',
    'samples': ['Fake']
}

aliases['fakeWMuUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.Declare('fake_rate_reader fr_reader_MuUp = fake_rate_reader(\"2022EE\", \"Run3\", \"67\", \"MuUp\", 2, \"std\");')"],
    'expr': 'fr_reader_MuUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanestJet_pt, nCleanestJet)',
    'samples': ['Fake']
}

aliases['fakeWMuDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.Declare('fake_rate_reader fr_reader_MuDown = fake_rate_reader(\"2022EE\", \"Run3\", \"67\", \"MuDown\", 2, \"std\");')"],
    'expr': 'fr_reader_MuDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanestJet_pt, nCleanestJet)',
    'samples': ['Fake']
}

aliases['fakeWStatEleUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.Declare('fake_rate_reader fr_reader_StatEleUp = fake_rate_reader(\"2022EE\", \"Run3\", \"67\", \"StatEleUp\", 2, \"std\");')"],
    'expr': 'fr_reader_StatEleUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanestJet_pt, nCleanestJet)',
    'samples': ['Fake']
}
aliases['fakeWStatEleDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.Declare('fake_rate_reader fr_reader_StatEleDown = fake_rate_reader(\"2022EE\", \"Run3\", \"67\", \"StatEleDown\", 2, \"std\");')"],
    'expr': 'fr_reader_StatEleDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanestJet_pt, nCleanestJet)',
    'samples': ['Fake']
}

aliases['fakeWStatMuUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.Declare('fake_rate_reader fr_reader_StatMuUp = fake_rate_reader(\"2022EE\", \"Run3\", \"67\", \"StatMuUp\", 2, \"std\");')"],
    'expr': 'fr_reader_StatMuUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanestJet_pt, nCleanestJet)',
    'samples': ['Fake']
}

aliases['fakeWStatMuDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.Declare('fake_rate_reader fr_reader_StatMuDown = fake_rate_reader(\"2022EE\", \"Run3\", \"67\", \"StatMuDown\", 2, \"std\");')"],
    'expr': 'fr_reader_StatMuDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanestJet_pt, nCleanestJet)',
    'samples': ['Fake']
}

###### --------------------------------------

aliases['gstarLow'] = {
    'expr': 'Gen_ZGstar_mass >0 && Gen_ZGstar_mass < 4',
    'samples': ['WZ', 'VgS', 'Vg']
}
aliases['gstarHigh'] = {
    'expr': 'Gen_ZGstar_mass <0 || Gen_ZGstar_mass > 4',
    'samples': ['WZ', 'VgS', 'Vg'],
}

# Jet bins
# using Alt(CleanestJet_pt, n, 0) instead of Sum(CleanestJet_pt >= 20) because jet pt ordering is not strictly followed in JES-varied samples

# One jet: leading jet with pt > 30 GeV
aliases['oneJet'] = {
    'expr': 'Alt(CleanestJet_pt, 0, 0) > 30.',
    'afterNuis': True
}

# Multiple jets: leading jet with pt > 30, others with pt > 20 GeV
aliases['multiJet'] = {
    'expr': 'Alt(CleanestJet_pt, 0, 0) > 30. && Alt(CleanestJet_pt, 1, 0) > 20.',
    'afterNuis': True
}

# Three jets: leading jet with pt > 30, two others with pt > 20 GeV
aliases['ThreeJet'] = {
    'expr': 'Alt(CleanestJet_pt, 0, 0) > 30. && Alt(CleanestJet_pt, 1, 0) > 20. && Alt(CleanestJet_pt, 2, 0) > 20.',
    'afterNuis': True
}

# Number of jets
aliases['njets'] = {
    'expr': 'Sum(CleanestJet_pt > 20.) '
}

# Fixing issues with jets in the horns
aliases['noJetInHorn'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/jet_horns.cc"'],
    'expr': 'Jet_inHorns(CleanestJet_pt, CleanestJet_eta)',
    'afterNuis': True
}

aliases['noJetInHorn_pT15'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/jet_horns.cc"'],
    'expr': 'Jet_inHorns(CleanestJet_pt, CleanestJet_eta, true)'
}

############################################################################
# B-Tagging WP: https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22EE/
############################################################################

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
bWP    = 'loose'     # ['loose','medium','tight','xtight','xxtight']
#bSF   = 'deepjet'

# b veto
aliases['bVeto'] = {
    'expr': 'Sum(CleanestJet_pt > 20. && abs(CleanestJet_eta) < 2.5 && Take(Jet_btag{}, CleanestJet_jetIdx) > {}) == 0'.format(bAlgo, btagging_WPs[bAlgo][bWP])
}

# At least one b-tagged jet  
aliases['bReq'] = { 
    'expr': 'Sum(CleanestJet_pt > 30. && abs(CleanestJet_eta) < 2.5 && Take(Jet_btag{}, CleanestJet_jetIdx) > {}) >= 1'.format(bAlgo, btagging_WPs[bAlgo][bWP])
}

# Number of b-jets
aliases['nbjets'] = {
    'expr': 'Sum(CleanestJet_pt > 20. && abs(CleanestJet_eta) < 2.5 && Take(Jet_btag{}, CleanestJet_jetIdx) > {})'.format(bAlgo, btagging_WPs[bAlgo][bWP])
}

year = '2022_Summer22EE' 
shifts = ['central', 'up_uncorrelated', 'down_uncorrelated', 'up_correlated', 'down_correlated']
shift_str = '{"' + '","'.join(shifts) + '"}'

for flavour in ['bc', 'light']:
    btagsf_tmp = 'btagSF_TMP' + flavour
    aliases[btagsf_tmp] = {
        'linesToProcess':[
            f'ROOT.gSystem.Load("/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/evaluate_btagSF{flavour}_cc.so","", ROOT.kTRUE)',
            f"ROOT.gInterpreter.Declare('btagSF{flavour} btag_SF{flavour} = btagSF{flavour}(\"/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/data/btag_eff/bTagEff_2022EE_ttbar_loose.root\",\"{year}\",\"_parT\");')"
        ],
        'expr': f'btag_SF{flavour}(CleanestJet_pt, CleanestJet_eta, CleanestJet_jetIdx, nCleanestJet, Jet_hadronFlavour, Jet_btag{bAlgo}, "L", {shift_str})',
        'samples' : mc,
    }
    for i in range(len(shifts)):
        btagsf = 'btagSF' + flavour
        if shifts[i] != 'central':
            btagsf += '_' + shifts[i]
        aliases[btagsf] = {
            'expr': f"{btagsf_tmp}[{i}]",
            'samples' : mc,
        }
##########################################################################
# End of b tagging
##########################################################################

# Data/MC scale factors and systematic uncertainties
aliases['SFweight'] = {
    'expr': ' * '.join(['SFweight2l', 'LepWPCut', 'LepWPSF', 'btagSFbc', 'btagSFlight']), # used to apply leptons SFs
    #'expr': ' * '.join(['SFweight2l', 'LepWPCut']), # used just for leptons WP cut
    'samples': mc
}

aliases['SFweightEleUp'] = {
    'expr': 'LepSF2l__ele_'+eleWP+'__Up',
    'samples': mc
}
aliases['SFweightEleDown'] = {
    'expr': 'LepSF2l__ele_'+eleWP+'__Down',
    'samples': mc
}
aliases['SFweightMuUp'] = {
    'expr': 'LepSF2l__mu_'+muWP+'__Up',
    'samples': mc
}
aliases['SFweightMuDown'] = {
    'expr': 'LepSF2l__mu_'+muWP+'__Down',
    'samples': mc
}

############################################################################
############### Definition of tt+DM relevant variables #####################
############################################################################
aliases['doubleNu_producer'] = {
    'linesToAdd': [f'#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/doubleNu_producer.cc"'],
    'class': 'doubleNu_producer',
    'args': 'nCleanestJet, CleanestJet_pt, CleanestJet_eta, CleanestJet_phi, CleanestJet_mass, CleanestJet_jetIdx, nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, cleanestPuppiMET_pt, cleanestPuppiMET_phi, Jet_btag{}, {}'.format(bAlgo, btagging_WPs[bAlgo][bWP]),
}

### Defining other relevant variables ###
#### mT2 variable definition
aliases['mT2'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/computeMT2.cc"'],
    'expr': 'computeMT2(Lepton_pt[0], Lepton_eta[0], Lepton_phi[0], Lepton_pt[1], Lepton_eta[1], Lepton_phi[1], cleanestPuppiMET_pt, cleanestPuppiMET_phi)',
    'afterNuis': True
}

# mpmet variable definition
aliases['mpmet'] = {
    'expr' : 'projtkmet<cleanest_projpfmet ? projtkmet : cleanest_projpfmet'
}

# ttZ indices
aliases['zLep1'] = {
    'linesToAdd': ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/ttZ_Leptons.cc"'],
    'expr': 'getZLep1Index(nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId)',
    'afterNuis': True
}

aliases['zLep2'] = {
    'expr': 'getZLep2Index(nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId)',
    'afterNuis': True
}

aliases['otherLepIndex'] = {
    'expr': 'getOtherLepIndex(nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId)',
    'afterNuis': True
}

aliases['zLep_mll'] = {
    'expr': 'getZLepMll(nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId)',
    'afterNuis': True
}

#############################################################################
####################### tt+DM regions definition ############################
#############################################################################
#Signal region
aliases['sr'] = {
    'expr': 'mT2 > 80 && (!(abs(Lepton_pdgId[0]) == abs(Lepton_pdgId[1])) || abs(91.1876 - mll) > 15) && Alt(Lepton_pt, 2, 0) < 10',
    'afterNuis': True
} 

# ttZ control region                                                                                                                                                                                       
aliases['ttZcr'] = {
    'expr': 'nLepton == 3 && Lepton_pt[2] > 20 && ThreeJet && zLep1 >=0 && zLep2 >=0 && zLep_mll > 0 && (Lepton_pdgId[zLep1] == -Lepton_pdgId[zLep2]) && abs(91.1876 - zLep_mll) < 10 && Lepton_pt[otherLepIndex] > 35',
    'afterNuis': True
}

# DY control region
aliases['dycr'] = {
    'expr': 'mT2 > 80 && (abs(Lepton_pdgId[0]) == abs(Lepton_pdgId[1])) && abs(91.1876 - mll) < 15 && Alt(Lepton_pt, 2, 0) < 10',
    'afterNuis': True
}

# Validation region
aliases['ttvr'] = {
    'expr': 'mT2 < 80 && (!(abs(Lepton_pdgId[0]) == abs(Lepton_pdgId[1])) || abs(91.1876 - mll) > 15) && Alt(Lepton_pt, 2, 0) < 10',
    'afterNuis': True
}
