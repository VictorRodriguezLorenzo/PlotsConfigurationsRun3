
mcProduction = 'Summer23BPix_130x_nAODv12_Full2023BPixv12'
mcSteps      = 'MCl2loose2023BPixv12__MCCorr2023BPixv12JetScaling__sblancof__l2tight'
dataReco     = 'Run2023BPix_Prompt_nAODv12_Full2023BPixv12'
dataSteps    = 'DATAl2loose2023BPixv12__sblancof__l2loose'

treeBaseDir = '/eos/cms/store/group/phys_higgs/cmshww/calderon/HWWNano/'
limitFiles = -1

mc = [skey for skey in samples if skey not in ('Fake', 'DATA')]

redirector = ""

useXROOTD = False

def makeMCDirectory(var=''):
    _treeBaseDir = treeBaseDir + ''
    if useXROOTD:
        _treeBaseDir = redirector + treeBaseDir
    if var== '':
        return '/'.join([_treeBaseDir, mcProduction, mcSteps])
    else:
        return '/'.join([_treeBaseDir, mcProduction, mcSteps + '__' + var])



mcDirectory = makeMCDirectory()
#fakeDirectory = os.path.join(treeBaseDir, dataReco, fakeSteps)
dataDirectory = os.path.join(treeBaseDir, dataReco, dataSteps)
fakeDirectory = dataDirectory
print(treeBaseDir)

cuts0j = []
cuts1j = []
cuts2j = []
cuts_vbf = []
cuts_2j = []
total_cuts = []
for k in cuts:
    for cat in cuts[k]['categories']:
        total_cuts.append(k+'_'+cat)
        if '0j' in cat:
            cuts0j.append(k+'_'+cat)
        elif '1j' in cat: 
            cuts1j.append(k+'_'+cat)
        elif '2j' in cat and '2j_vbf' not in cat: 
            cuts2j.append(k+'_'+cat)
            cuts_2j.append(k+'_'+cat)
        elif '2j_vbf' in cat:
            cuts_vbf.append(k+'_'+cat)
            cuts_2j.append(k+'_'+cat)
        else: 
            print('WARNING: name of category does not contain either 0j,1j,2j')

nuisances = {}


################################ EXPERIMENTAL UNCERTAINTIES  #################################

#### Luminosity

# https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3
nuisances['lumi_2023'] = {
    'name'    : 'lumi_2023',
    'type'    : 'lnN',
    'samples' : dict((skey, '1.013') for skey in mc)
}

#### FAKES

nuisances['fake_syst'] = {
    'name'    : 'CMS_fake_syst',
    'type'    : 'lnN',
    'samples' : {
        'Fake' : '1.3'
    },
}
nuisances['fake_ele'] = {
    'name'    : 'CMS_fake_e_2023BPix',
    'kind'    : 'weight',
    'type'    : 'shape',
    'samples' : {
        'Fake' : ['fakeWEleUp', 'fakeWEleDown'],
    }
}
nuisances['fake_ele_stat'] = {
    'name'    : 'CMS_fake_stat_e_2023BPix',
    'kind'    : 'weight',
    'type'    : 'shape',
    'samples' : {
        'Fake' : ['fakeWStatEleUp', 'fakeWStatEleDown']
    }
}
nuisances['fake_mu'] = {
    'name'    : 'CMS_fake_m_2023BPix',
    'kind'    : 'weight',
    'type'    : 'shape',
    'samples' : {
        'Fake' : ['fakeWMuUp', 'fakeWMuDown'],
    }   
}       
nuisances['fake_mu_stat'] = {
    'name'    : 'CMS_fake_stat_m_2023BPix',
    'kind'    : 'weight',
    'type'    : 'shape',
    'samples' : {
        'Fake' : ['fakeWStatMuUp', 'fakeWStatMuDown'],
    }
}

##### B-tagger

for flavour in ['bc', 'light']:
    for corr in ['uncorrelated', 'correlated']:
        btag_syst = [f'btagSF{flavour}_up_{corr}/btagSF{flavour}', f'btagSF{flavour}_down_{corr}/btagSF{flavour}']
        if corr == 'correlated':
            name = f'CMS_btagSF{flavour}_{corr}'
        else:
            name = f'CMS_btagSF{flavour}_2023BPix'
        nuisances[f'btagSF{flavour}{corr}'] = {
            'name': name,
            'skipCMS' : 1,
            'kind': 'weight',
            'type': 'shape',
            'samples': dict((skey, btag_syst) for skey in mc),
        }

##### Trigger Scale Factors                                                                                                                                                                                

trig_syst = ['TriggerSFWeight_2l_u/TriggerSFWeight_2l', 'TriggerSFWeight_2l_d/TriggerSFWeight_2l']

nuisances['trigg'] = {
    'name': 'CMS_eff_hwwtrigger_2023BPix',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, trig_syst) for skey in mc)
}

##### Electron Efficiency and energy scale

nuisances['eff_e'] = {
    'name': 'CMS_eff_e_2023BPix',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, ['SFweightEleUp', 'SFweightEleDown']) for skey in mc),
}

##### Muon Efficiency and energy scale

nuisances['eff_m'] = {
    'name': 'CMS_eff_m_2023BPix',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, ['SFweightMuUp', 'SFweightMuDown']) for skey in mc),
}

#### Lepton scale

nuisances['leppt_scale'] = {
    'name'       : 'CMS_scale_l_2023BPix',
    'kind'       : 'suffix',
    'type'       : 'shape',
    'mapUp'      : 'leptonScaleup',
    'mapDown'    : 'leptonScaledo',
    'samples'    : dict((skey, ['1', '1']) for skey in mc),
    'folderUp'   : makeMCDirectory('leptonScaleup_suffix'),
    'folderDown' : makeMCDirectory('leptonScaledo_suffix'),
    'AsLnN'      : '0'
}

nuisances['leppt_res'] = {
    'name'       : 'CMS_resolution_l_2023BPix',
    'kind'       : 'suffix',
    'type'       : 'shape',
    'mapUp'      : 'leptonResolutionup',
    'mapDown'    : 'leptonResolutiondo',
    'samples'    : dict((skey, ['1', '1']) for skey in mc),
    'folderUp'   : makeMCDirectory('leptonResolutionup_suffix'),
    'folderDown' : makeMCDirectory('leptonResolutiondo_suffix'),
    'AsLnN'      : '0'
}


##### JES

jes_systs    = ["Absolute", "Absolute_2023BPix", "FlavorQCD", "BBEC1", "EC2", "HF", "BBEC1_2023BPix", "EC2_2023BPix", "RelativeBal", "RelativeSample_2023BPix", "HF_2023BPix"] # Reduced set of 11 uncertainties
#jes_systs = ['jesTotal']

for js in jes_systs:
    
    nuisances[js] = {
        'name'      : 'CMS_scale_j_' + js,
        'kind'      : 'suffix',
        'type'      : 'shape',
        'mapUp'     : 'jesRegroed_' + js + 'up',
        'mapDown'   : 'jesRegroed_' + js + 'do',
        'samples'   : dict((skey, ['1', '1']) for skey in mc),
        'folderUp'  : makeMCDirectory('jesRegroed_' + js + 'up_suffix'),
        'folderDown': makeMCDirectory('jesRegroed_' + js + 'do_suffix'),
        'AsLnN'     : '0'
    }

##### Jet energy resolution
nuisances['JER'] = {
    'name'      : 'CMS_res_j_2023',
    'kind'      : 'suffix',
    'type'      : 'shape',
    'mapUp'     : 'jerup',
    'mapDown'   : 'jerdo',
    'samples'   : dict((skey, ['1', '1']) for skey in mc),
    'cuts'      : [cut for cut in cuts],
    'folderUp'  : makeMCDirectory('jerup_suffix'),
    'folderDown': makeMCDirectory('jerdo_suffix'),
    'AsLnN'     : '0'
}

##### MET energy scale
nuisances['met'] = {
    'name'      : 'CMS_scale_met_2023',
    'kind'      : 'suffix',
    'type'      : 'shape',
    'mapUp'     : 'unclustEnup',
    'mapDown'   : 'unclustEndo',
    'samples'   : dict((skey, ['1', '1']) for skey in mc),
    'cuts'      : [cut for cut in cuts],
    'folderUp'  : makeMCDirectory('unclustEnup_suffix'),
    'folderDown': makeMCDirectory('unclustEndo_suffix'),
    'AsLnN'     : '0'
}

##### Pileup
nuisances['PU'] = {
    'name': 'CMS_pileup_2023',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, ['puWeightUp/puWeight', 'puWeightDown/puWeight']) for skey in mc),
    'AsLnN'   : '0',
}        

##### PS

#
# As suggested by Emmanuelle, split the nuisance parameters as a function of the number of jets; in a similar behavior as it's done with the top QCD scales.
#

for ibin in ['0j','1j','2j']:
    nuisances['PS_ISR_'+ibin]  = {
        'name'    : 'PS_hww_ISR_'+ibin,
        'kind'    : 'weight',
        'type'    : 'shape',
        'samples' : dict((skey, ['PSWeight[2]', 'PSWeight[0]']) for skey in mc),
        'cutspost' : [cut for cut in total_cuts if ibin in cut],
        'AsLnN'   : '0',
    }
    nuisances['PS_FSR_'+ibin]  = {
        'name'    : 'PS_hww_FSR_'+ibin,
        'kind'    : 'weight',
        'type'    : 'shape',
        'samples' : dict((skey, ['PSWeight[3]', 'PSWeight[1]']) for skey in mc),
        'cutspost' : [cut for cut in total_cuts if ibin in cut],
        'AsLnN'   : '0',
    }


nuisances['UE_CP5']  = {
    'name'    : 'CMS_hww_UE',
    'skipCMS' : 1,
    'type'    : 'lnN',
    'samples' : dict((skey, '1.015') for skey in mc),
}

####### Generic "cross section uncertainties"

apply_on = {
    'top': [
        '(topGenPt * antitopGenPt <= 0.) * 1.0816 + (topGenPt * antitopGenPt > 0.)',
        '(topGenPt * antitopGenPt <= 0.) * 0.9184 + (topGenPt * antitopGenPt > 0.)'
    ]
}

nuisances['singleTopToTTbar'] = {
    'name': 'singleTopToTTbar',
    'skipCMS': 1,
    'kind': 'weight',
    'type': 'shape',
    'samples': apply_on
}

nuisances['WgStar'] = {
    'name'    : 'CMS_hww_WgStarScale',
    'type'    : 'lnN',
    'samples' : {
        'WgS' : '1.25'
    }
}

###### pdf uncertainties
pdf_variations = ["LHEPdfWeight[%d]" %i for i in range(1,101)] # Float_t LHE pdf variation weights (w_var / w_nominal) for LHA IDs  320901 - 321000
nuisances['pdf_WW']  = {
    'name'  : 'CMS_hww_pdf_WW',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'WW'   : pdf_variations,
    },
}

nuisances['pdf_top']  = {
    'name'  : 'CMS_hww_pdf_top',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'top'   : pdf_variations,
    },
}

nuisances['pdf_ggH']  = {
    'name'  : 'CMS_hww_pdf_ggH',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'ggH_hww'   : pdf_variations,
    },
}

nuisances['pdf_qqH']  = {
    'name'  : 'CMS_hww_pdf_qqH',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'qqH_hww'   : pdf_variations,
    },
}

nuisances['pdf_qqbar'] = {
    'name': 'pdf_qqbar',
    'type': 'lnN',
    'samples': {
        'Vg': '1.04',
        'VZ': '1.04',  # PDF: 0.0064 / 0.1427 = 0.0448493
        'VgS': '1.04', # PDF: 0.0064 / 0.1427 = 0.0448493
    },
}

##### Renormalization & factorization scales

variations = ['Alt(LHEScaleWeight,0,1)',
              'Alt(LHEScaleWeight,1,1)',
              'Alt(LHEScaleWeight,3,1)',
              'Alt(LHEScaleWeight,nLHEScaleWeight-4,1)',
              'Alt(LHEScaleWeight,nLHEScaleWeight-2,1)',
              'Alt(LHEScaleWeight,nLHEScaleWeight-1,1)']

for ibin in ['0j','1j','2j']:
    nuisances['QCDscale_top_'+ibin]  = {
        'name'  : 'QCDscale_top_'+ibin,
        'kind'  : 'weight',
        'type'  : 'shape',
        'AsLnN': '0',
        'cutspost' : [cut for cut in total_cuts if ibin in cut],
        'samples'  : {
            'top' : variations,
        }
    }

nuisances['QCDscale_V'] = {
    'name': 'QCDscale_V',
    'skipCMS': 1,
    'kind'  : 'weight_envelope',
    'type': 'shape',
    'samples': {'DY': variations},
    'AsLnN': '1'
}
nuisances['QCDscale_VV'] = {
    'name' : 'QCDscale_VV',
    'kind' : 'weight_envelope',
    'type' : 'shape',
    'samples' : {
        'WW'  : variations,
        'WW_minnlo'  : variations,
        'WWewk' : variations,
        'Vg'  : variations,
        #'ZZ'  : variations, # Not in branches
        #'WZ'  : variations,
        #'VgS' : variations,
    }
}
nuisances['QCDscale_ggH'] = {
    'name': 'QCDscale_ggVV',
    'kind'  : 'weight_envelope',
    'type'  : 'shape',
    'samples': {
        'ggH_hww': variations,
    },
}
nuisances['QCDscale_qqH']  = {
    'name'  : 'QCDscale_qqH',
    'kind'  : 'weight_envelope',
    'type'  : 'shape',
    'samples'  : {
        'qqH_hww' : variations,
    }
}
nuisances['QCDscale_ggWW'] = {
    'name': 'QCDscale_ggWW',
    'type': 'lnN',
    'samples': {
        'ggWW': '1.10',
    },
}

nuisances['WWresum0j']  = {
    'name'  : 'CMS_hww_WWresum_0j',
    'skipCMS' : 1,
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples'  : {
        'WW'   : ['nllW_Rup/wwNLL', 'nllW_Rdown/wwNLL'],
    },
    'cutspost' : [cut for cut in total_cuts if '0j' in cut],
}
nuisances['WWqscale0j']  = {
    'name'  : 'CMS_hww_WWqscale_0j',
    'skipCMS' : 1,
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples'  : {
        'WW'   : ['nllW_Qup/wwNLL', 'nllW_Qdown/wwNLL'],
    },
    'cutspost' : [cut for cut in total_cuts if '0j' in cut],
}
nuisances['WWresum1j']  = {
    'name'  : 'CMS_hww_WWresum_1j',
    'skipCMS' : 1,
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples'  : {
        'WW'   : ['nllW_Rup/wwNLL', 'nllW_Rdown/wwNLL'],
    },
    'cutspost' : [cut for cut in total_cuts if '1j' in cut],
}
nuisances['WWqscale1j']  = {
    'name'  : 'CMS_hww_WWqscale_1j',
    'skipCMS' : 1,
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples'  : {
        'WW'   : ['nllW_Qup/wwNLL', 'nllW_Qdown/wwNLL'],
    },
    'cutspost' : [cut for cut in total_cuts if '1j' in cut],
}
nuisances['WWresum2j']  = {
    'name'  : 'CMS_hww_WWresum_2j',
    'skipCMS' : 1,
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples'  : {
        'WW'   : ['nllW_Rup/wwNLL', 'nllW_Rdown/wwNLL'],
    },
    'cutspost' : [cut for cut in total_cuts if '2j' in cut],
}
nuisances['WWqscale2j']  = {
    'name'  : 'CMS_hww_WWqscale_2j',
    'skipCMS' : 1,
    'kind'  : 'weight',
    'type'  : 'shape',
    'samples'  : {
        'WW'   : ['nllW_Qup/wwNLL', 'nllW_Qdown/wwNLL'],
    },
    'cutspost' : [cut for cut in total_cuts if '2j' in cut],
}

### Not applied cause already accounted in QCDScale
#nuisances['GGWWRew'] = {
#    'name': 'CMS_ggWW_NLO_reweighting',
#    'kind'       : 'weight',
#    'type'       : 'shape',
#    'samples'    : {
#        'ggWW': ["KFactor_ggWW_Up/KFactor_ggWW", "KFactor_ggWW_Down/KFactor_ggWW"],
#        'ggWW_si': ["KFactor_ggWW_Up/KFactor_ggWW", "KFactor_ggWW_Down/KFactor_ggWW"],
#        'ggToWW': ["KFactor_ggWW_Up/KFactor_ggWW", "KFactor_ggWW_Down/KFactor_ggWW"],
#    },
#}

########## Theory uncertainties for Higgs

# From https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWG136TeVxsec_extrap
nuisances['pdf_Higgs_ggH_ACCEPT'] = {
    'name': 'pdf_Higgs_ggH_accept',
    'samples': {
        'ggH_hww' : '1.032',
        'ggH_htt' : '1.032',
    },
    'type': 'lnN',
}

nuisances['pdf_Higgs_ttH_ACCEPT'] = {
    'name': 'pdf_Higgs_ttH_accept',
    'samples': {
        'ttH_hww': '1.035',
    },
    'type': 'lnN',
}

nuisances['pdf_Higgs_qqbar_ACCEPT'] = {
    'name': 'pdf_Higgs_qqbar_accept',
    'type': 'lnN',
    'samples': {
        'qqH_hww': '1.031/0.966',
        'qqH_htt': '1.031/0.966',
        'WH_hww': '1.018',
        'WH_htt': '1.018',
        'ZH_hww': '1.016',
        'ZH_htt': '1.016',
    },
}

nuisances['QCDScale_ggH_ACCEPT'] = {
    'name': 'QCDScale_ggH_accept',
    'samples': {
        'ggH_hww' : '1.039',
        'ggH_htt' : '1.039',
    },
    'type': 'lnN',
}

nuisances['QCDscale_qqH_ACCEPT'] = {
    'name'    : 'QCDscale_qqH_accept',
    'samples' : {
        'qqH_hww' : '1.032/0.958',
        'qqH_htt' : '1.032/0.958',
    },
    'type' : 'lnN'
}

# Uncertainty on SR/CR ratio
nuisances['CRSR_accept_DY'] = {
    'name': 'CMS_hww_CRSR_accept_DY',
    'type': 'lnN',
    'samples': {
        'DY': '1.02',
        'Dyemb': '1.02'
    },
    'cuts': [cut for cut in total_cuts if '_dytt_' in cut],
    'cutspost' : [cut for cut in total_cuts if '_dytt_' in cut],
}
# Uncertainty on SR/CR ratio
nuisances['CRSR_accept_top'] = {
    'name': 'CMS_hww_CRSR_accept_top',
    'type': 'lnN',
    'samples': {'top': '1.01'},
    'cuts': [cut for cut in total_cuts if '_top_' in cut],
    'cutspost' : [cut for cut in total_cuts if '_top_' in cut],
}
# Uncertainty on SR/CR ratio
nuisances['CRSR_accept_WW'] = {
    'name': 'CMS_hww_CRSR_accept_WW',
    'type': 'lnN',
    'samples': {'WW': '1.01'},
    'cuts': [cut for cut in total_cuts if 'TeV_ww' in cut],
    'cutspost' : [cut for cut in total_cuts if '_TeV_ww_' in cut],
}

##rate parameters
nuisances['DYembnorm0j']  = {
               'name'  : 'CMS_hww_DYttnorm0j',
               'samples'  : {
                   'DY' : '1.00',
                   },
               'type'  : 'rateParam',
               'cuts'  : cuts0j
              }
nuisances['DYembnorm1j']  = {
               'name'  : 'CMS_hww_DYttnorm1j',
               'samples'  : {
                   'DY' : '1.00',
                   },
               'type'  : 'rateParam',
               'cuts'  : cuts1j
              }
nuisances['DYembnorm2j']  = {
                 'name'  : 'CMS_hww_DYttnorm2j',
                 'samples'  : {
                   'DY' : '1.00',
                     },
                 'type'  : 'rateParam',
                 'cuts'  : cuts_2j,
                }
nuisances['WWnorm0j']  = {
               'name'  : 'CMS_hww_WWnorm0j',
               'samples'  : {
                   'WW' : '1.00',
                   #'WW_minnlo' : '1.00',
                   },
               'type'  : 'rateParam',
               'cuts'  : cuts0j
              }
nuisances['WWnorm1j']  = {
               'name'  : 'CMS_hww_WWnorm1j',
               'samples'  : {
                   'WW' : '1.00',
                   #'WW_minnlo' : '1.00',
                   },
               'type'  : 'rateParam',
               'cuts'  : cuts1j
              }
nuisances['WWnorm2j']  = {
               'name'  : 'CMS_hww_WWnorm2j',
               'samples'  : {
                   'WW' : '1.00',
                   #'WW_minnlo' : '1.00',
                   },
               'type'  : 'rateParam',
               'cuts': cuts_2j,
              }
nuisances['Topnorm0j']  = {
               'name'  : 'CMS_hww_Topnorm0j',
               'samples'  : {
                   'top' : '1.00',
                   },
               'type'  : 'rateParam',
               'cuts'  : cuts0j
              }
nuisances['Topnorm1j']  = {
               'name'  : 'CMS_hww_Topnorm1j',
               'samples'  : {
                   'top' : '1.00',
                   },
               'type'  : 'rateParam',
               'cuts'  : cuts1j
              }
nuisances['Topnorm2j']  = {
               'name'  : 'CMS_hww_Topnorm2j',
               'samples'  : {
                   'top' : '1.00',
                   },
               'type'  : 'rateParam',
               'cuts'  : cuts_2j
              }

### MC statistical uncertainty
autoStats = True
if autoStats:
    ## Use the following if you want to apply the automatic combine MC stat nuisances.
    nuisances['stat'] = {
        'type': 'auto',
        'maxPoiss': '10',
        'includeSignal': '0',
        #  nuisance ['maxPoiss'] =  Number of threshold events for Poisson modelling
        #  nuisance ['includeSignal'] =  Include MC stat nuisances on signal processes (1=True, 0=False)
        'samples': {}
    }
