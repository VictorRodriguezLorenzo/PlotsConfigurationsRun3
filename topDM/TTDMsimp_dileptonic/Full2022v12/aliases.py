import os
import copy
import inspect
import joblib

configurations = os.path.realpath(inspect.getfile(inspect.currentframe())) # this file

aliases = {}
aliases = OrderedDict()

mc     = [skey for skey in samples if skey not in ('Fake', 'DATA', 'Dyemb', 'DATA_EG', 'DATA_Mu', 'DATA_EMu', 'Fake_EG', 'Fake_Mu', 'Fake_EMu')]
mc_emb = [skey for skey in samples if skey not in ('Fake', 'DATA', 'DATA_Mu', 'DATA_EMu', 'Fake_EG', 'Fake_Mu', 'Fake_EMu')]

# LepSF2l__ele_cutBased_MediumID_tthMVA_Run3__mu_cut_TightID_pfIsoLoose_HWW_tthmva_67
eleWP = 'cutBased_MediumID_tthMVA_Run3'
muWP  = 'cut_TightID_pfIsoLoose_HWW_tthmva_67'

aliases['LepWPCut'] = {
    'expr': 'LepCut2l__ele_'+eleWP+'__mu_'+muWP,
    'samples': mc + ['DATA'],
}

aliases['LepWPSF'] = {
    'expr': 'LepSF2l__ele_'+eleWP+'__mu_'+muWP,
    'samples': mc
}

# gen-matching to prompt only (GenLepMatch2l matches to *any* gen lepton)
aliases['PromptGenLepMatch2l'] = {
    'expr': 'Alt(Lepton_promptgenmatched, 0, 0) * Alt(Lepton_promptgenmatched, 1, 0)',
    'samples': mc
}

# Fake leptons transfer factor --------------------------------------
aliases['fakeW'] = {
    'linesToAdd' : [f'#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"nominal\", 2, \"std\");')"],
    'expr': 'fr_reader(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_MediumID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples'    : ['Fake']
}

# And variations - already divided by central values in formulas !
aliases['fakeWEleUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_EleUp = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"EleUp\", 2, \"std\");')"],
    'expr': 'fr_reader_EleUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_MediumID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}
aliases['fakeWEleDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_EleDown = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"EleDown\", 2, \"std\");')"],
    'expr': 'fr_reader_EleDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_MediumID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}

aliases['fakeWMuUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_MuUp = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"MuUp\", 2, \"std\");')"],
    'expr': 'fr_reader_MuUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_MediumID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}

aliases['fakeWMuDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_MuDown = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"MuDown\", 2, \"std\");')"],
    'expr': 'fr_reader_MuDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_MediumID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}

aliases['fakeWStatEleUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_StatEleUp = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"StatEleUp\", 2, \"std\");')"],
    'expr': 'fr_reader_StatEleUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_MediumID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}
aliases['fakeWStatEleDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_StatEleDown = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"StatEleDown\", 2, \"std\");')"],
    'expr': 'fr_reader_StatEleDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_MediumID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}

aliases['fakeWStatMuUp'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_StatMuUp = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"StatMuUp\", 2, \"std\");')"],
    'expr': 'fr_reader_StatMuUp(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_MediumID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
    'samples': ['Fake']
}

aliases['fakeWStatMuDown'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/fake_rate_reader_class.cc"'],
    'linesToProcess':["ROOT.gInterpreter.ProcessLine('fake_rate_reader fr_reader_StatMuDown = fake_rate_reader(\"2022\", \"Run3\", \"67\", \"StatMuDown\", 2, \"std\");')"],
    'expr': 'fr_reader_StatMuDown(Lepton_pdgId, Lepton_pt, Lepton_eta, Lepton_isTightMuon_cut_TightID_pfIsoLoose_HWW_tthmva_67, Lepton_isTightElectron_cutBased_MediumID_tthMVA_Run3, CleanJet_pt, nCleanJet)',
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

# Top pT reweighting ##

# NNLO + NLO EW correction derived at 13 TeV
aliases['Top_pTrw'] = {
    'expr': '(topGenPt * antitopGenPt > 0.) * (TMath::Sqrt((0.103*TMath::Exp(-0.0118*topGenPt) - 0.000134*topGenPt + 0.973) * (0.103*TMath::Exp(-0.0118*antitopGenPt) - 0.000134*antitopGenPt + 0.973))) + (topGenPt * antitopGenPt <= 0.)',
    'samples': ['TTTo2L2Nu']
}

# Extrapolation of the correction from 13 TeV to 13.6 TeV, from TOP-25-018
aliases['Top_pTrw_13To13p6'] = {
    'expr': '(topGenPt * antitopGenPt > 0.) * TMath::Sqrt((0.991 + 0.000075*topGenPt) * (0.991 + 0.000075*antitopGenPt)) + (topGenPt * antitopGenPt <= 0.)',
    'samples': ['TTTo2L2Nu']
}

# Jet bins
# using Alt(CleanJet_pt, n, 0) instead of Sum(CleanJet_pt >= 20) because jet pt ordering is not strictly followed in JES-varied samples

# One jet: leading jet with pt > 20 GeV
aliases['oneJet'] = {
    'expr': 'Alt(CleanJet_pt, 0, 0) > 20.',
    'afterNuis': True
}

# Multiple jets: leading jet with pt > 20, others with pt > 20 GeV
aliases['multiJet'] = {
    'expr': 'Alt(CleanJet_pt, 0, 0) > 20. && Alt(CleanJet_pt, 1, 0) > 20.',
    'afterNuis': True
}

# Three jets: leading jet with pt > 20, two others with pt > 20 GeV
aliases['ThreeJet'] = {
    'expr': 'Alt(CleanJet_pt, 0, 0) > 20. && Alt(CleanJet_pt, 1, 0) > 20. && Alt(CleanJet_pt, 2, 0) > 20.',
    'afterNuis': True
}

# Number of jets
aliases['njets'] = {
    'expr': 'Sum(CleanJet_pt > 20.)',
    'afterNuis': True
}

# Fixing issues with jets in the horns
aliases['noJetInHorn'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/jet_horns.cc"'],
    'expr': 'Jet_inHorns(CleanJet_pt, CleanJet_eta)',
    'afterNuis': True
}

aliases['noJetInHorn_pT20'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/jet_horns.cc"'],
    'expr': 'Jet_inHorns(CleanJet_pt, CleanJet_eta, true)',
    'afterNuis': True
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
bAlgo = 'PNetB' # ['DeepFlavB','RobustParTAK4B','PNetB'] 
bWP    = 'medium'     # ['loose','medium','tight','xtight','xxtight']
#bSF   = 'deepjet'

# No b-tagged jets
aliases['bVeto'] = {
    'expr': 'Sum(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{}, CleanJet_jetIdx) > {}) == 0'.format(
        bAlgo, btagging_WPs[bAlgo][bWP]
    )
}

# At least one b-tagged jet
aliases['bReq'] = {
    'expr': 'Sum(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{}, CleanJet_jetIdx) > {}) >= 1'.format(
        bAlgo, btagging_WPs[bAlgo][bWP]
    )
}

# Number of b-jets
aliases['nbjets'] = {
    'expr': 'Sum(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && Take(Jet_btag{}, CleanJet_jetIdx) > {})'.format(
        bAlgo, btagging_WPs[bAlgo][bWP]
    )
}

## Ratio of b-jets to selected jets
aliases['nbjet_jet_ratio'] = {
    'expr': '(njets > 0 ? nbjets/njets : 0.)'.format(bAlgo, btagging_WPs[bAlgo][bWP]),
    'afterNuis': True
}

year = '2022_Summer22' 
btag_eff_file = 'bTagEff_Full2022v12_ttbar_PNetB_medium.root'
#shifts = ['central', 'down_fsrdef', 'down_hdamp', 'down_isrdef', 'down_jer', 'down_jes', 'down_mass', 'down_statistic', 'down_tune', 'up_fsrdef', 'up_hdamp','up_isrdef', 'up_jer', 'up_jes', 'up_mass', 'up_statistic', 'up_tune']
shifts = ['central', 'up_uncorrelated', 'down_uncorrelated', 'up_correlated', 'down_correlated']
shift_str = '{"' + '","'.join(shifts) + '"}'
wp_map = {
    "loose": "L",
    "medium": "M",
    "tight": "T",
    "xtight": "XT",
    "xxtight": "XXT"
}

for flavour in ['bc', 'light']:
    btagsf_tmp = 'btagSF_TMP' + flavour
    aliases[btagsf_tmp] = {
        'linesToProcess':[
            f'ROOT.gSystem.Load("/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/extended/evaluate_btagSF{flavour}_cc.so","", ROOT.kTRUE)',
            f"ROOT.gInterpreter.ProcessLine('btagSF{flavour} btag_SF{flavour} = btagSF{flavour}(\"/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/data/btag_eff/{btag_eff_file}\",\"{year}\",\"\");')"
        ],
        'expr': f'btag_SF{flavour}(CleanJet_pt, CleanJet_eta, CleanJet_jetIdx, nCleanJet, Jet_hadronFlavour, Jet_btag{bAlgo}, "{wp_map[bWP]}", {shift_str})',
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
    'args': 'nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, PuppiMET_pt, PuppiMET_phi, Jet_btag{}, {}'.format(bAlgo, btagging_WPs[bAlgo][bWP]),
    'afterNuis': True
}

### Defining other relevant variables ###
# mT2 variable definition
aliases['mT2'] = {
    'linesToAdd' : ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/computeMT2.cc"'],
    'expr': 'computeMT2(Lepton_pt[0], Lepton_eta[0], Lepton_phi[0], Lepton_pt[1], Lepton_eta[1], Lepton_phi[1], PuppiMET_pt, PuppiMET_phi)',
    'afterNuis': True
}

aliases['mt2blbl'] = {
    'linesToAdd': ['#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/computeMT2blbl.cc"'],
    'expr': 'computeMT2blbl(nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, PuppiMET_pt, PuppiMET_phi, Jet_btag{}, {})'.format(bAlgo,btagging_WPs[bAlgo][bWP]),
    'afterNuis': True
}

# mpmet variable definition
aliases['mpmet'] = {
    'expr' : 'min(projtkmet, projpfmet)',
    'afterNuis': True
}

############################################################################
############ Extra topology variables for tt+DM / top+DM ####################
############################################################################

aliases['topDMVars'] = {
    'linesToAdd': [
        '#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/topDM_vars.cc"'
    ],
    'expr': 'topDM_DNN_vars(nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, PuppiMET_pt, PuppiMET_phi, Jet_btag{}, {})'.format(
        bAlgo,
        btagging_WPs[bAlgo][bWP]
    ),
    'afterNuis': True
}

_topdm_alias_map = {
    'dphi_met_ll'              : 0,
    'st'                       : 2,
    'met_over_sqrt_ht'         : 3,
    'met_over_st'              : 4,
    'dphi_min_j_met'           : 5,

    'pt_llb'                   : 6,
    'dphi_met_llb'             : 7,

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

for _name, _idx in _topdm_alias_map.items():
    aliases[_name] = {
        'expr': 'topDMVars[{}]'.format(_idx),
        'afterNuis': True
    }

aliases['top1_pt_reco'] = {
    'expr': 'doubleNu_producer[4]',
    'afterNuis': True
}

aliases['top2_pt_reco'] = {
    'expr': 'doubleNu_producer[5]',
    'afterNuis': True
}

aliases['pdark_over_met'] = {
    'expr': '(PuppiMET_pt > 0 && doubleNu_producer[8] > -998) ? doubleNu_producer[8]/PuppiMET_pt : -999.',
    'afterNuis': True
}

aliases['topDMRestFrameVars'] = {
    'linesToAdd': [
        '#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/macros/topDM_restFrame_vars.cc"'
    ],
    'expr': 'topDM_restFrame_vars(nCleanJet, CleanJet_pt, CleanJet_eta, CleanJet_phi, CleanJet_mass, CleanJet_jetIdx, nLepton, Lepton_pt, Lepton_eta, Lepton_phi, Lepton_pdgId, PuppiMET_pt, PuppiMET_phi, Jet_btag{}, {})'.format(
        bAlgo,
        btagging_WPs[bAlgo][bWP]
    ),
    'afterNuis': True
}

_restframe_alias_map = {
    'angle_ll_llbb_rf' : 0,
    'dphi_ll_llbb_rf'  : 1,
    'cos_l1_llbb_rf'   : 2,
    'cos_l2_llbb_rf'   : 3,

    'angle_ll_llmet_rf': 4,
    'dphi_ll_llmet_rf' : 5,
    'cos_l1_llmet_rf'  : 6,
    'cos_l2_llmet_rf'  : 7,
}

for _name, _idx in _restframe_alias_map.items():
    aliases[_name] = {
        'expr': 'topDMRestFrameVars[{}]'.format(_idx),
        'afterNuis': True
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
########################## DNN discriminator ################################
#############################################################################

# Mass hypothesis points
mPhi = ['10', '50', '100', '150', '200', '250', '300', '350', '400', '500', '600', '700', '800', '1000']
mPhiRVec = 'ROOT::VecOps::RVec<float>{' + ','.join([f'{phi}.f' for phi in mPhi]) + '}'

rdf_name = {
    'lep_pt1': 'pt1',
    'lep_pt2': 'pt2',
    'lep_eta1': 'Lepton_eta[0]',
    'lep_eta2': 'Lepton_eta[1]',
    'chel': 'doubleNu_producer[6]',
    'pdark': 'doubleNu_producer[8]',
    'dphi_ttbar': 'doubleNu_producer[7]',
    'dphi_met_llb_safe': 'dphi_met_llb',
}

# ttDM ps DNN discriminant
features_ttDM = [v for v in joblib.load(
    "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022v12/DNNmodels/Models/features_model_DNN_ttDM_ps.pkl"
) if v != 'mPhi']

featuresRVec_ttDM = 'ROOT::VecOps::RVec<float>{' + ','.join(f'static_cast<float>({rdf_name.get(v, v)})' for v in features_ttDM if v != 'mPhi') + '}'

aliases['evaluate_dnn_ttDM_ps'] = {
    'linesToAdd': [
        '#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022v12/evaluate_DNN_ttDM.cc"'
    ],
    'class': 'evaluate_dnn_ttDM',
    'args': f'{featuresRVec_ttDM}, {mPhiRVec}, "ps"',
    'afterNuis': True
}

# ttDM s DNN discriminant
features_ttDM = [v for v in joblib.load(
    "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022v12/DNNmodels/Models/features_model_DNN_ttDM_s.pkl"
) if v != 'mPhi']

featuresRVec_ttDM = 'ROOT::VecOps::RVec<float>{' + ','.join(f'static_cast<float>({rdf_name.get(v, v)})' for v in features_ttDM if v != 'mPhi') + '}'

aliases['evaluate_dnn_ttDM_s'] = {
    'linesToAdd': [
        '#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022v12/evaluate_DNN_ttDM.cc"'
    ],
    'class': 'evaluate_dnn_ttDM',
    'args': f'{featuresRVec_ttDM}, {mPhiRVec}, "s"',
    'afterNuis': True
}

# tWDM ps DNN discriminant
features_tWDM = [v for v in joblib.load(
    "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022v12/DNNmodels/Models/features_model_DNN_tWDM_ps.pkl"
) if v != 'mPhi']

featuresRVec_tWDM = 'ROOT::VecOps::RVec<float>{' + ','.join(f'static_cast<float>({rdf_name.get(v, v)})' for v in features_tWDM if v != 'mPhi') + '}'

aliases['evaluate_dnn_tWDM_ps'] = {
    'linesToAdd': [
        '#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022v12/evaluate_DNN_tWDM.cc"'
    ],
    'class': 'evaluate_dnn_tWDM',
    'args': f'{featuresRVec_tWDM}, {mPhiRVec}, "ps"',
    'afterNuis': True
}

# tWDM s DNN discriminant
features_tWDM = [v for v in joblib.load(
    "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022v12/DNNmodels/Models/features_model_DNN_tWDM_s.pkl"
) if v != 'mPhi']

featuresRVec_tWDM = 'ROOT::VecOps::RVec<float>{' + ','.join(f'static_cast<float>({rdf_name.get(v, v)})' for v in features_tWDM if v != 'mPhi') + '}'

aliases['evaluate_dnn_tWDM_s'] = {
    'linesToAdd': [
        '#include "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022v12/evaluate_DNN_tWDM.cc"'
    ],
    'class': 'evaluate_dnn_tWDM',
    'args': f'{featuresRVec_tWDM}, {mPhiRVec}, "s"',
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

## ttZ control region                                                                                                                                                                                       
aliases['ttZcr'] = {
    'expr': 'nLepton == 3 && Lepton_pt[2] > 20 && ((abs(Lepton_pdgId[2]) == 11 && abs(Lepton_eta[2]) < 2.5)|| (abs(Lepton_pdgId[2]) == 13 && abs(Lepton_eta[2]) < 2.4)) && ThreeJet && zLep1 >=0 && zLep2 >=0 && zLep_mll > 0 && (Lepton_pdgId[zLep1] == -Lepton_pdgId[zLep2]) && abs(91.1876 - zLep_mll) < 10 && Lepton_pt[otherLepIndex] > 35',
    'afterNuis': True
}

# DY control region
aliases['dycr'] = {
    'expr': 'mT2 > 80 && (abs(Lepton_pdgId[0]) == abs(Lepton_pdgId[1])) && abs(91.1876 - mll) < 15 && Alt(Lepton_pt, 2, 0) < 10',
    'afterNuis': True
}

# tt(2l) control region
aliases['ttcr'] = {
    'expr': 'mT2 < 80 && (!(abs(Lepton_pdgId[0]) == abs(Lepton_pdgId[1])) || abs(91.1876 - mll) > 15) && Alt(Lepton_pt, 2, 0) < 10',
    'afterNuis': True
}

