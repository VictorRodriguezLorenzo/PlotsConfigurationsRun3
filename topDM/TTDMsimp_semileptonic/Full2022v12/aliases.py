import os
import copy
import inspect

configurations = os.path.realpath(inspect.getfile(inspect.currentframe())) # this file

aliases = {}
aliases = OrderedDict()

mc     = [skey for skey in samples if skey not in ('Fake', 'DATA', 'Dyemb', 'DATA_EG', 'DATA_Mu', 'DATA_EMu', 'Fake_EG', 'Fake_Mu', 'Fake_EMu')]
mc_emb = [skey for skey in samples if skey not in ('Fake', 'DATA', 'DATA_Mu', 'DATA_EMu', 'Fake_EG', 'Fake_Mu', 'Fake_EMu')]

# One lep WPs
eleWP = 'cutBased_LooseID_tthMVA_Run3'
muWP  = 'cut_TightID_pfIsoLoose_HWW_tthmva_67'

aliases['LepWPCut'] = {
    'expr': "(Lepton_isTightElectron_" + eleWP + "[0]>0.5 || Lepton_isTightMuon_" + muWP + "[0]>0.5)",
    'samples': mc + ['DATA'],
}

aliases['LepWPSF'] = {
    'expr': "(abs(Lepton_pdgId[0])==11 ? Lepton_tightElectron_" + eleWP + "_IdIsoSF[0] : 1.0) * "
            "(abs(Lepton_pdgId[0])==13 ? Lepton_tightMuon_" + muWP + "_IdIsoSF[0] : 1.0)",
    'samples': mc
}

# gen-matching to prompt only (GenLepMatch1l matches to *any* gen lepton)
aliases['PromptGenLepMatch1l'] = {
    'expr': 'Alt(Lepton_promptgenmatched, 0, 0)',
    'samples': mc
}

# Fake leptons transfer factor --------------------------------------
aliases['fakeW'] = {
    'linesToAdd' : [f'#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"nominal\", 2, \"std\");')"],
    'expr': 'fr_reader(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples'    : ['Fake']
}

# And variations - already divided by central values in formulas !
aliases['fakeWEleUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_EleUp = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"EleUp\", 2, \"std\");')"],
    'expr': 'fr_reader_EleUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}
aliases['fakeWEleDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_EleDown = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"EleDown\", 2, \"std\");')"],
    'expr': 'fr_reader_EleDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}

aliases['fakeWMuUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_MuUp = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"MuUp\", 2, \"std\");')"],
    'expr': 'fr_reader_MuUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}

aliases['fakeWMuDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_MuDown = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"MuDown\", 2, \"std\");')"],
    'expr': 'fr_reader_MuDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}

aliases['fakeWStatEleUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_StatEleUp = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"StatEleUp\", 2, \"std\");')"],
    'expr': 'fr_reader_StatEleUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}
aliases['fakeWStatEleDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_StatEleDown = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"StatEleDown\", 2, \"std\");')"],
    'expr': 'fr_reader_StatEleDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}

aliases['fakeWStatMuUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_StatMuUp = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"StatMuUp\", 2, \"std\");')"],
    'expr': 'fr_reader_StatMuUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}

aliases['fakeWStatMuDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_StatMuDown = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"StatMuDown\", 2, \"std\");')"],
    'expr': 'fr_reader_StatMuDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_LooseID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
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
# using Alt(CleanJet_pt, n, 0) instead of Sum(CleanJet_pt >= 20) because jet pt ordering is not strictly followed in JES-varied samples

# One jet: leading jet with pt > 30 GeV
aliases['oneJet'] = {
    'expr': 'Alt(CleanJet_pt, 0, 0) > 30.',
    'afterNuis': True
}

# Multiple jets: leading jet with pt > 30, others also with pt > 30 GeV
aliases['multiJet'] = {
    'expr': 'Alt(CleanJet_pt, 0, 0) > 30. && Alt(CleanJet_pt, 1, 0) > 30.',
    'afterNuis': True
}

# Number of jets
aliases['njets'] = {
    'expr': 'Sum(CleanJet_pt > 30.) '
}

aliases['nForwardjets'] = {
    'expr': 'Sum(CleanJet_pt > 30. && CleanJet_eta > 2.5) '
}

# Fixing issues with jets in the horns
aliases['noJetInHorn'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/jet_horns.cc"'],
    'expr': 'Jet_inHorns(CleanJet_pt, CleanJet_eta)',
    'afterNuis': True
}

aliases['noJetInHorn_pT15'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/jet_horns.cc"'],
    'expr': 'Jet_inHorns(CleanJet_pt, CleanJet_eta, true)'
}

############################################################################
# B-Tagging WP: https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Summer22/
############################################################################

# Algo / WP / WP cut
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

# Algorithm and WP selection
bAlgo = 'RobustParTAK4B' # ['DeepFlavB','RobustParTAK4B','PNetB'] 
bWP    = 'loose'     # ['loose','medium','tight','xtight','xxtight']
#bSF   = 'deepjet'

# b veto
aliases['bVeto'] = {
    'expr': 'Sum(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{}, CleanJet_jetIdx) > {}) == 0'.format(bAlgo, btagging_WPs[bAlgo][bWP])
}

# At least one b-tagged jet  
aliases['bReq'] = { 
    'expr': 'Sum(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{}, CleanJet_jetIdx) > {}) >= 1'.format(bAlgo, btagging_WPs[bAlgo][bWP])
}

# Number of b-jets
aliases['nbjets'] = {
    'expr': 'Sum(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{}, CleanJet_jetIdx) > {})'.format(bAlgo, btagging_WPs[bAlgo][bWP])
}

year = '2022_Summer22' 
shifts = ['central', 'up_uncorrelated', 'down_uncorrelated', 'up_correlated', 'down_correlated']
shift_str = '{"' + '","'.join(shifts) + '"}'

for flavour in ['bc', 'light']:
    btagsf_tmp = 'btagSF_TMP' + flavour
    aliases[btagsf_tmp] = {
        'linesToProcess':[
            f'ROOT.gSystem.Load("/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/extended/evaluate_btagSF{flavour}_cc.so","", ROOT.kTRUE)',
            f"ROOT.gInterpreter.ProcessLine('btagSF{flavour} btag_SF{flavour} = btagSF{flavour}(\"/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/data/btag_eff/bTagEff_2022_ttbar_loose.root\",\"{year}\",\"_parT\");')"
        ],
        'expr': f'btag_SF{flavour}(CleanJet_pt, CleanJet_eta, CleanJet_jetIdx, nCleanJet, Jet_hadronFlavour, Jet_btag{bAlgo}, "L", {shift_str})',
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
aliases['SFweight1l'] = {
    'expr': ' * '.join(['puWeight', 'TriggerEffWeight_1l', 'Lepton_RecoSF[0]']),
    'samples': mc
}

aliases['SFweight'] = {
    'expr': ' * '.join(['SFweight1l', 'LepWPCut', 'LepWPSF', 'btagSFbc', 'btagSFlight']), 
    'samples': mc
}

aliases['SFweightEleUp'] = {
    'expr': '((TMath::Abs(Lepton_pdgId[0]) == 11)*(Lepton_tightElectron_'+eleWP+'_TotSF_Up[0]/Lepton_tightElectron_'+eleWP+'_TotSF[0]) + (TMath::Abs(Lepton_pdgId[0]) == 13))',
    'samples': mc
}
aliases['SFweightEleDown'] = {
    'expr': '((abs(Lepton_pdgId[0]) == 11)*(Lepton_tightElectron_'+eleWP+'_TotSF_Down[0]/Lepton_tightElectron_'+eleWP+'_TotSF[0]) + (abs(Lepton_pdgId[0]) == 13))',
    'samples': mc
}
aliases['SFweightMuUp'] = {
    'expr': '((abs(Lepton_pdgId[0]) == 13)*(Lepton_tightMuon_'+muWP+'_TotSF_Up[0]/Lepton_tightMuon_'+muWP+'_TotSF[0]) + (abs(Lepton_pdgId[0]) == 11))',
    'samples': mc
}
aliases['SFweightMuDown'] = {
    'expr': '((abs(Lepton_pdgId[0]) == 13)*(Lepton_tightMuon_'+muWP+'_TotSF_Down[0]/Lepton_tightMuon_'+muWP+'_TotSF[0]) + (abs(Lepton_pdgId[0]) == 11))',
    'samples': mc
}

############################################################################
############### Definition of tt+DM relevant variables #####################
############################################################################
# mT2 variable definition
aliases['mT2'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/macros/computeMT2.cc"'],
    'expr': 'computeMT2(Lepton_pt[0], Lepton_eta[0], Lepton_phi[0], Lepton_pt[1], Lepton_eta[1], Lepton_phi[1], PuppiMET_pt, PuppiMET_phi)',
    'afterNuis': True
}

# mT2w variable definition
aliases['mt2w'] = {
    'linesToAdd': [f'#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/macros/mt2w_producer.cc"'],
    'class': 'mt2w_producer',
    'args': 'Lepton_pt[0], Lepton_eta[0], Lepton_phi[0], CleanJet_pt, CleanJet_eta, CleanJet_phi, Take(Jet_btag{algo}, CleanJet_jetIdx), PuppiMET_pt, PuppiMET_phi, {wp}'.format(algo=bAlgo, wp=btagging_WPs[bAlgo][bWP]),
    'afterNuis': True
}

aliases['mindphi_jet_met'] = {
    'expr': 'std::min(std::abs(DeltaPhi(CleanJet_phi[0], PuppiMET_phi)), std::abs(DeltaPhi(CleanJet_phi[1], PuppiMET_phi)))',
    'afterNuis': True
}

aliases['MTb'] = {
    'linesToAdd': ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/macros/computeMTb.cc"'],
    'expr': 'computeMTb(CleanJet_pt, CleanJet_eta, CleanJet_phi, Take(Jet_btag{algo}, CleanJet_jetIdx), {wp}, PuppiMET_pt, PuppiMET_phi)'.format(algo=bAlgo, wp=btagging_WPs[bAlgo][bWP]),
    'afterNuis': True
}

aliases['H1T'] = {
    'expr': 'CleanJet_pt[0] / Sum(CleanJet_pt[CleanJet_pt > 30])',
    'afterNuis': True
}

### Defining other relevant variables ###
# topness variable definition
aliases['topness'] = {
    'linesToAdd': ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_semileptonic/macros/topness_producer.cc"'],
    'class': 'topness_producer_minuit',
    'args': 'nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, PuppiMET_pt, PuppiMET_phi, Take(Jet_btag{algo}, CleanJet_jetIdx), {wp}'.format(algo=bAlgo, wp=btagging_WPs[bAlgo][bWP]),
}

# mpmet variable definition
aliases['mpmet'] = {
    'expr' : 'min(projtkmet, projpfmet)',
    'afterNuis': True
}

#############################################################################
####################### tt+DM regions definition ############################
#############################################################################
#Signal region
aliases['sr'] = {
    'expr': 'PuppiMET_pt >= 250 && mtw1 >= 140 && mindphi_jet_met >= 0.8 && bReq && Alt(Lepton_pt, 1, 0) < 10',
    'afterNuis': True
} 

# W+jets control region
aliases['wjetscr'] = {
    'expr': 'mT2 <= 80 && PuppiMET_pt >= 250 && mtw1 >= 140 && bVeto && Alt(Lepton_pt, 1, 0) < 10',
    'afterNuis': True
}

# Validation region
aliases['ttcr'] = {
    'expr': 'mT2 <= 80 && nLepton == 2 && PuppiMET_pt >= 250 && bReq',
    'afterNuis': True
}
