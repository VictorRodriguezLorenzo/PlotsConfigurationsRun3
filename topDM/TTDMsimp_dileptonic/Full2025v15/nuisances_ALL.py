treeBaseDir = '/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/'
#signalBaseDir   = '/eos/user/v/victorr/HWWNano/'
signalBaseDir   = '/eos/user/e/emunozri/ttDM/HWWNano/'

# MC backgrounds
mcProduction = 'Summer24_150x_nAODv15_Full2024v15'
mcSteps      = 'MCl2loose2024v15__MCCorr2024v15__JERFrom23BPix__l2tight'

# Signal
signalProduction = 'Summer24_150x_nAODv15_Full2024v15'
signalSteps      = 'MCl2loose2024v15__MCCorr2024v15__JERFrom23BPix__l2tight'

# Data
dataRecoMuon     = 'Run2025_PromptCDEFG_nAODv15_Full2025v15_Muon'
dataRecoMuonEG     = 'Run2025_PromptCDEFG_nAODv15_Full2025v15_EGamma'
dataRecoEGamma     = 'Run2025_PromptCDEFG_nAODv15_Full2025v15_MuonEG'
dataSteps    = 'DATAl2loose2025v15__l2loose'

limitFiles = -1

mc = [skey for skey in samples if skey not in ('Fake', 'DATA')]

redirector = ""

useXROOTD = False

def makeDirectory(baseDir, production, steps, var='', useXROOTD=False, redirector=''):
    _treeBaseDir = baseDir + ""
    if useXROOTD:
        _treeBaseDir = redirector + baseDir
    if var == '':
        return '/'.join([_treeBaseDir, production, steps])
    else:
        return '/'.join([_treeBaseDir, production, steps + '__' + var])

def makeDirectoryForSkey(skey, var=''):
    if 'DM' in skey:
        return makeDirectory(signalBaseDir, signalProduction, signalSteps, var)
    else:
        return makeDirectory(treeBaseDir, mcProduction, mcSteps, var)

mcDirectory = makeDirectory(treeBaseDir, mcProduction, mcSteps)
signalDirectory  = makeDirectory(signalBaseDir, signalProduction, signalSteps)

dataDirectoryMuon = makeDirectory(treeBaseDir, dataRecoMuon, dataSteps)
dataDirectoryMuonEG = makeDirectory(treeBaseDir, dataRecoMuonEG, dataSteps)
dataDirectoryEGamma = makeDirectory(treeBaseDir, dataRecoEGamma, dataSteps)

fakeDirectoryMuon = dataDirectoryMuon
fakeDirectoryMuonEG = dataDirectoryMuonEG
fakeDirectoryEGamma = dataDirectoryEGamma

print("\n")
print("MC Directory:", mcDirectory)
print("Signal Directory:", signalDirectory)
print("\n")

nuisances = {}

################################ EXPERIMENTAL UNCERTAINTIES  #################################

#### Luminosity

# https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3

nuisances['lumi_2025'] = {
    'name'    : 'lumi_2025',
    'type'    : 'lnN',
    'samples' : dict((skey, '1.016') for skey in mc)
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
    'name'    : 'CMS_fake_e_2025',
    'kind'    : 'weight',
    'type'    : 'shape',
    'samples' : {
        'Fake' : ['fakeWEleUp', 'fakeWEleDown'],
    }
}

nuisances['fake_ele_stat'] = {
    'name'    : 'CMS_fake_stat_e_2025',
    'kind'    : 'weight',
    'type'    : 'shape',
    'samples' : {
        'Fake' : ['fakeWStatEleUp', 'fakeWStatEleDown']
    }
}

nuisances['fake_mu'] = {
    'name'    : 'CMS_fake_m_2025',
    'kind'    : 'weight',
    'type'    : 'shape',
    'samples' : {
        'Fake' : ['fakeWMuUp', 'fakeWMuDown'],
    }   
}       

nuisances['fake_mu_stat'] = {
    'name'    : 'CMS_fake_stat_m_2025',
    'kind'    : 'weight',
    'type'    : 'shape',
    'samples' : {
        'Fake' : ['fakeWStatMuUp', 'fakeWStatMuDown'],
    }
}

##### B-tagger

#shifts = ['fsrdef', 'hdamp', 'isrdef', 'jer', 'jes', 'mass', 'statistic', 'tune']
shifts = ['uncorrelated', 'correlated']

for flavour in ['bc', 'light']:
    for corr in shifts:
        btag_syst = [f'btagSF{flavour}_up_{corr}/btagSF{flavour}', f'btagSF{flavour}_down_{corr}/btagSF{flavour}']
        if corr == 'correlated':
            name = f'CMS_btagSF{flavour}_{corr}'
        else:
            name = f'CMS_btagSF{flavour}_2025'
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
    'name': 'CMS_eff_trigger_2025',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, trig_syst) for skey in mc)
}

##### Electron Efficiency and energy scale

nuisances['eff_e'] = {
    'name': 'CMS_eff_e_2025',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, ['SFweightEleUp', 'SFweightEleDown']) for skey in mc),
}

##### Muon Efficiency and energy scale

nuisances['eff_m'] = {
    'name': 'CMS_eff_m_2025',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, ['SFweightMuUp', 'SFweightMuDown']) for skey in mc),
}

#### Lepton scale

nuisances['leppt_scale'] = {
    'name'       : 'CMS_scale_l_2025',
    'kind'       : 'suffix',
    'type'       : 'shape',
    'mapUp'      : 'leptonScaleup',
    'mapDown'    : 'leptonScaledo',
    'samples'    : dict((skey, ['1', '1']) for skey in mc),
    'folderUp'   : {skey: makeDirectoryForSkey(skey, 'leptonScaleup_suffix') for skey in mc},
    'folderDown' : {skey: makeDirectoryForSkey(skey, 'leptonScaledo_suffix')  for skey in mc},
    'AsLnN'      : '0'
}

nuisances['leppt_res'] = {
    'name'       : 'CMS_resolution_l_2025',
    'kind'       : 'suffix',
    'type'       : 'shape',
    'mapUp'      : 'leptonResolutionup',
    'mapDown'    : 'leptonResolutiondo',
    'samples'    : dict((skey, ['1', '1']) for skey in mc),
    'folderUp'   : {skey: makeDirectoryForSkey(skey, 'leptonResolutionup_suffix') for skey in mc},
    'folderDown' : {skey: makeDirectoryForSkey(skey, 'leptonResolutiondo_suffix')  for skey in mc},
    'AsLnN'      : '0'
}

##### JES

jes_systs    = ["Absolute", "Absolute_2025", "FlavorQCD", "BBEC1", "EC2", "HF", "BBEC1_2025", "EC2_2025", "RelativeBal", "RelativeSample_2025", "HF_2025"] # Reduced set of 11 uncertainties
#jes_systs = ['jesTotal']

for js in jes_systs:
    
    nuisances[js] = {
        'name'      : 'CMS_scale_j_' + js,
        'kind'      : 'suffix',
        'type'      : 'shape',
        'mapUp'     : 'jesRegroed_' + js + 'up',
        'mapDown'   : 'jesRegroed_' + js + 'do',
        'samples'   : dict((skey, ['1', '1']) for skey in mc),
        'folderUp'  : {skey: makeDirectoryForSkey(skey, 'jesRegroed_' + js + 'up_suffix') for skey in mc},
        'folderDown': {skey: makeDirectoryForSkey(skey, 'jesRegroed_' + js + 'do_suffix') for skey in mc},
        'AsLnN'     : '0'
    }

##### Jet energy resolution

nuisances['JER'] = {
    'name'      : 'CMS_res_j_2025',
    'kind'      : 'suffix',
    'type'      : 'shape',
    'mapUp'     : 'jerup',
    'mapDown'   : 'jerdo',
    'samples'   : dict((skey, ['1', '1']) for skey in mc),
    'cuts'      : [cut for cut in cuts],
    'folderUp'  : {skey: makeDirectoryForSkey(skey, 'jerup_suffix') for skey in mc},
    'folderDown': {skey: makeDirectoryForSkey(skey, 'jerdo_suffix') for skey in mc},
    'AsLnN'     : '0'
}

##### MET energy scale

nuisances['met'] = {
    'name'      : 'CMS_scale_met_2025',
    'kind'      : 'suffix',
    'type'      : 'shape',
    'mapUp'     : 'unclustEnup',
    'mapDown'   : 'unclustEndo',
    'samples'   : dict((skey, ['1', '1']) for skey in mc),
    'cuts'      : [cut for cut in cuts],
    'folderUp'  : {skey: makeDirectoryForSkey(skey, 'unclustEnup_suffix') for skey in mc},
    'folderDown': {skey: makeDirectoryForSkey(skey, 'unclustEndo_suffix') for skey in mc},
    'AsLnN'     : '0'
}

##### Pileup

nuisances['PU'] = {
    'name': 'CMS_pileup_2025',
    'kind': 'weight',
    'type': 'shape',
    'samples': dict((skey, ['puWeightUp/puWeight', 'puWeightDown/puWeight']) for skey in mc),
    'AsLnN'   : '0',
}        

##### PS

nuisances['PS_ISR']  = {
    'name'    : 'PS_ISR',
    'kind'    : 'weight',
    'type'    : 'shape',
    'samples' : dict((skey, ['PSWeight[2]', 'PSWeight[0]']) for skey in mc),
    'AsLnN'   : '0',
}

nuisances['PS_FSR']  = {
    'name'    : 'PS_FSR',
    'kind'    : 'weight',
    'type'    : 'shape',
    'samples' : dict((skey, ['PSWeight[3]', 'PSWeight[1]']) for skey in mc),
    'AsLnN'   : '0',
}

nuisances['UE_CP5']  = {
    'name'    : 'CMS_UE',
    'skipCMS' : 1,
    'type'    : 'lnN',
    'samples' : dict((skey, '1.015') for skey in mc),
}

##### Top pT reweighting uncertainty

nuisances['TopPtRew'] = {
    'name': 'CMS_topPtRew',   # Theory uncertainty
    'kind': 'weight',
    'type': 'shape',
    'samples': {'TTTo2L2Nu': ["1.", "1./Top_pTrw"]},
    'symmetrize': True
}

####### Generic "cross section uncertainties"

apply_on = {
    'ST': [
        '(topGenPt * antitopGenPt <= 0.) * 1.0816 + (topGenPt * antitopGenPt > 0.)',
        '(topGenPt * antitopGenPt <= 0.) * 0.9184 + (topGenPt * antitopGenPt > 0.)'
    ],
    'TTTo2L2Nu': [
        '(topGenPt * antitopGenPt <= 0.) * 1.0816 + (topGenPt * antitopGenPt > 0.)',
        '(topGenPt * antitopGenPt <= 0.) * 0.9184 + (topGenPt * antitopGenPt > 0.)'
    ],
    'TTToSemiLeptonic': [
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

###### pdf uncertainties

pdf_variations = ["LHEPdfWeight[%d]" %i for i in range(1,101)] # Float_t LHE pdf variation weights (w_var / w_nominal) for LHA IDs  320901 - 321000

nuisances['pdf_V'] = {
    'name'  : 'pdf_V',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples': {
        'DY': pdf_variations
    }
}

nuisances['pdf_VV']  = {
    'name'  : 'pdf_VV',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'WW'   : pdf_variations,
#        'WZ'   : pdf_variations,
#        'ZZ'   : pdf_variations,
    },
}

nuisances['pdf_TTbar']  = {
    'name'  : 'pdf_TTbar',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'TTTo2L2Nu'       : pdf_variations,
        'TTToSemiLeptonic': pdf_variations,
    },
}

nuisances['pdf_ST']  = {
    'name'  : 'pdf_ST',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'ST'              : pdf_variations,
    },
}

nuisances['pdf_ttV']  = {
    'name'  : 'pdf_ttV',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'ttW'             : pdf_variations,
        'ttZ'             : pdf_variations,
    },
}

nuisances['pdf_ttH']  = {
    'name'  : 'pdf_ttH',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'ttH'             : pdf_variations,
    },
}

nuisances['pdf_tXX']  = {
    'name'  : 'pdf_tXX',
    'kind'  : 'weight_rms',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'THQ'             : pdf_variations,
        'THW'             : pdf_variations,
    },
}

##### Renormalization & factorization scales

variations = ['Alt(LHEScaleWeight,0,1)',
              'Alt(LHEScaleWeight,1,1)',
              'Alt(LHEScaleWeight,3,1)',
              'Alt(LHEScaleWeight,nLHEScaleWeight-4,1)',
              'Alt(LHEScaleWeight,nLHEScaleWeight-2,1)',
              'Alt(LHEScaleWeight,nLHEScaleWeight-1,1)']

nuisances['QCDscale_TTbar']  = {
    'name'  : 'QCDscale_TTbar',
    'skipCMS': 1,
    'kind'  : 'weight_envelope',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'TTTo2L2Nu': variations,
        'TTToSemiLeptonic': variations,
    }
}

nuisances['QCDscale_ST']  = {
    'name'  : 'QCDscale_ST',
    'skipCMS': 1,
    'kind'  : 'weight_envelope',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'ST': variations,
    }
}

nuisances['QCDscale_ttV']  = {
    'name'  : 'QCDscale_ttV',
    'skipCMS': 1,
    'kind'  : 'weight_envelope',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'ttW': variations,
        'ttZ': variations,
    }
}

nuisances['QCDscale_ttH']  = {
    'name'  : 'QCDscale_ttH',
    'skipCMS': 1,
    'kind'  : 'weight_envelope',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'ttH': variations,
    }
}

nuisances['QCDscale_tXX']  = {
    'name'  : 'QCDscale_tXX',
    'skipCMS': 1,
    'kind'  : 'weight_envelope',
    'type'  : 'shape',
    'AsLnN': '0',
    'samples'  : {
        'THQ': variations,
        'THW': variations,
    }
}

nuisances['QCDscale_V'] = {
    'name': 'QCDscale_V',
    'skipCMS': 1,
    'kind'  : 'weight_envelope',
    'type': 'shape',
    'AsLnN': '0',
    'samples': {'DY': variations},
}

nuisances['QCDscale_VV'] = {
    'name' : 'QCDscale_VV',
    'skipCMS': 1,
    'kind' : 'weight_envelope',
    'type' : 'shape',
    'AsLnN': '0',
    'samples' : {
        'WW'  : variations,
        #'Vg'  : variations,
        #'ZZ'  : variations, 
        #'WZ'  : variations,
        #'VgS' : variations,
    }
}

## CR rate parameters
nuisances['ttZnorm']  = {
               'name'  : 'ttZnorm_2023BPix',
               'samples'  : {
                   'ttZ' : '1.00',
                   },
               'type'  : 'rateParam',
               'cuts'  : [
                   'ttZcr_inclusive',
                   ]
              }

nuisances['DYnorm']  = {
               'name'  : 'DYnorm_2023BPix',
               'samples'  : {
                   'DY' : '1.00',
                   },
               'type'  : 'rateParam',
               'cuts'  : [
                   'dycr_1b',
                   ]
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
