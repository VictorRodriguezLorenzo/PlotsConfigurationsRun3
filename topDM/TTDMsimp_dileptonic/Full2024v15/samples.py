from mkShapesRDF.lib.search_files import SearchFiles

searchFiles = SearchFiles()

redirector = ""
useXROOTD = False

mcProduction = 'Summer24_150x_nAODv15_Full2024v15'
#mcSteps      = 'MCl2loose2024v15__MCCorr2024v15__JERFrom23BPix__l2tight'
mcSteps      = 'MCl2loose2024v15__MCCorr2024v15__JERFrom23BPix__l2tight'
dataRecoMuon     = 'Run2024_ReRecoCDE_PromptFGHI_nAODv15_Full2024v15_Muon'
dataRecoMuonEG     = 'Run2024_ReRecoCDE_PromptFGHI_nAODv15_Full2024v15_MuonEG'
dataRecoEGamma     = 'Run2024_ReRecoCDE_PromptFGHI_nAODv15_Full2024v15_EGamma'
#dataSteps    = 'DATAl2loose2024v15__sblancof__l2loose'
dataSteps    = 'DATAl2loose2024v15__l2loose'

##############################################
###### Tree base directory for the site ######
##############################################
#treeBaseDir = '/eos/cms/store/group/phys_higgs/cmshww/calderon/HWWNano/'
treeBaseDir = '/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/'

limitFiles = -1

def makeMCDirectory(var=""):
    _treeBaseDir = treeBaseDir + ""
    if redirector != "":
        _treeBaseDir = redirector + treeBaseDir
    if var == "":
        return "/".join([_treeBaseDir, mcProduction, mcSteps])
    else:
        return "/".join([_treeBaseDir, mcProduction, mcSteps + "__" + var])


mcDirectory   = makeMCDirectory()
dir_map = {
    'MuonEG': os.path.join(treeBaseDir, dataRecoMuonEG, dataSteps),
    'Muon0':  os.path.join(treeBaseDir, dataRecoMuon,   dataSteps),
    'Muon1':  os.path.join(treeBaseDir, dataRecoMuon,   dataSteps),
    'EGamma0': os.path.join(treeBaseDir, dataRecoEGamma, dataSteps),
    'EGamma1': os.path.join(treeBaseDir, dataRecoEGamma, dataSteps),
}

#signalDirectory = "/eos/user/v/victorr/HWWNano/Summer24_150x_nAODv15_Full2024v15/MCl2loose2024v15__MCCorr2024v15__JERFrom23BPix__l2tight"
signalDirectory = "/eos/user/e/emunozri/ttDM/HWWNano/Summer24_150x_nAODv15_Full2024v15/MCl2loose2024v15__MCCorr2024v15__JERFrom23BPix__l2tight"
fake_dir_map = dir_map.copy()
samples = {}


def nanoGetSampleFiles(path, name):
    _files = searchFiles.searchFiles(path, name, redirector=redirector)
    if limitFiles != -1 and len(_files) > limitFiles:
        return [(name, _files[:limitFiles])]
    else:
        return [(name, _files)]


def CombineBaseW(samples, proc, samplelist):
    _filtFiles = list(filter(lambda k: k[0] in samplelist, samples[proc]["name"]))
    _files = list(map(lambda k: k[1], _filtFiles))
    _l = list(map(lambda k: len(k), _files))
    leastFiles = _files[_l.index(min(_l))]
    dfSmall = ROOT.RDataFrame("Runs", leastFiles)
    s = dfSmall.Sum("genEventSumw").GetValue()
    f = ROOT.TFile(leastFiles[0])
    t = f.Get("Events")
    t.GetEntry(1)
    xs = t.baseW * s

    __files = []
    for f in _files:
        __files += f
    df = ROOT.RDataFrame("Runs", __files)
    s = df.Sum("genEventSumw").GetValue()
    newbaseW = str(xs / s)
    weight = newbaseW + "/baseW"

    for iSample in samplelist:
        addSampleWeight(samples, proc, iSample, weight)


def addSampleWeight(samples, sampleName, sampleNameType, weight):
    obj = list(filter(lambda k: k[0] == sampleNameType, samples[sampleName]["name"]))[0]
    samples[sampleName]["name"] = list(
        filter(lambda k: k[0] != sampleNameType, samples[sampleName]["name"])
    )
    if len(obj) > 2:
        samples[sampleName]["name"].append(
            (obj[0], obj[1], obj[2] + "*(" + weight + ")")
        )
    else:
        samples[sampleName]["name"].append((obj[0], obj[1], "(" + weight + ")"))


################################################
############ DATA DECLARATION ##################
################################################

DataRun = [
    ['C','Run2024C-ReReco-v1'],
    ['D','Run2024D-ReReco-v1'],
    ['E','Run2024E-ReReco-v1'],
    ['F','Run2024F-Prompt-v1'],
    ['G','Run2024G-Prompt-v1'],
    ['H','Run2024H-Prompt-v1'],
    ['I','Run2024I-Prompt-v1'],
]

DataSets = ['MuonEG','Muon0','Muon1','EGamma0','EGamma1']

DataTrig = {
    'MuonEG'          : 'Trigger_ElMu' ,
    'Muon0'           : '!Trigger_ElMu && (Trigger_sngMu || Trigger_dblMu)',
    'Muon1'           : '!Trigger_ElMu && (Trigger_sngMu || Trigger_dblMu)',
    'EGamma0'         : '!Trigger_ElMu && !Trigger_sngMu && !Trigger_dblMu && (Trigger_sngEl || Trigger_dblEl)',
    'EGamma1'         : '!Trigger_ElMu && !Trigger_sngMu && !Trigger_dblMu && (Trigger_sngEl || Trigger_dblEl)',
} 

#########################################
############ MC COMMON ##################
#########################################

mcCommonWeight        = '2*XSWeight*METFilter_Common*SFweight*oddEvent'
#mcCommonWeight        = '2*XSWeight*METFilter_Common*PromptGenLepMatch2l*SFweight*oddEvent'

###########################################
#############  BACKGROUNDS  ###############
###########################################

########## DY #########

ptllDYW_LO = '((0.632927+0.0456956*gen_ptll-0.00154485*gen_ptll*gen_ptll+2.64397e-05*gen_ptll*gen_ptll*gen_ptll-2.19374e-07*gen_ptll*gen_ptll*gen_ptll*gen_ptll+6.99751e-10*gen_ptll*gen_ptll*gen_ptll*gen_ptll*gen_ptll)*(gen_ptll>0)*(gen_ptll<100)+(1.41713-0.00165342*gen_ptll)*(gen_ptll>=100)*(gen_ptll<300)+1*(gen_ptll>=300))'

dy_samples = [
        'DYto2E-2Jets_MLL-50',
        'DYto2Mu-2Jets_MLL-50',
        'DYto2Tau-2Jets_MLL-50',
        'DYto2E-2Jets_MLL-10to50',
        'DYto2Mu-2Jets_MLL-10to50',
        'DYto2Tau-2Jets_MLL-10to50',
        ]

files = []

for sample in dy_samples:
    files += nanoGetSampleFiles(mcDirectory, sample)

samples['DY'] = {
        'name': files,
        'weight': mcCommonWeight,
        'FilesPerJob': 4,
        }

for sample in dy_samples:
    addSampleWeight(samples, 'DY', sample, ptllDYW_LO)

########## Single top #########
files = nanoGetSampleFiles(mcDirectory, 'ST_t-channel_top') \
        + nanoGetSampleFiles(mcDirectory, 'ST_t-channel_antitop') \
        + nanoGetSampleFiles(mcDirectory, 'ST_s-channel_plus') \
        + nanoGetSampleFiles(mcDirectory, 'ST_s-channel_minus') \
        + nanoGetSampleFiles(mcDirectory, 'TWminusto2L2Nu') \
        + nanoGetSampleFiles(mcDirectory, 'TbarWplusto2L2Nu') \
        + nanoGetSampleFiles(mcDirectory, 'ST_tW_top') \
        + nanoGetSampleFiles(mcDirectory, 'ST_tW_antitop')

samples['ST'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 10
}

########## TTTo2L2Nu #########
files = nanoGetSampleFiles(mcDirectory, 'TTTo2L2Nu')

samples['TTTo2L2Nu'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 3
}

addSampleWeight(samples,'TTTo2L2Nu','TTTo2L2Nu','Top_pTrw')
addSampleWeight(samples,'TTTo2L2Nu','TTTo2L2Nu','Top_pTrw_13To13p6')

##### TTToSemiLeptonic #####
files = nanoGetSampleFiles(mcDirectory, 'TTToSemiLeptonic')

samples['TTToSemiLeptonic'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 10
}

##### TTTo4Q #####
#files = nanoGetSampleFiles(mcDirectory, 'TTTo4Q')
#
#samples['TTTo4Q'] = {
#    'name': files,
#    'weight': mcCommonWeight,
#    'FilesPerJob': 20,
#}

########### TTX ##########
files = nanoGetSampleFiles(mcDirectory, 'TTNuNu') \
        + nanoGetSampleFiles(mcDirectory, 'TTLL_MLL-4to50') \
        + nanoGetSampleFiles(mcDirectory, 'TTLL_MLL-50') \
        + nanoGetSampleFiles(mcDirectory, 'TTZ-ZtoQQ')

samples['ttZ'] = {
        'name': files,
        'weight': mcCommonWeight,
        'FilesPerJob': 2
        }


files = nanoGetSampleFiles(mcDirectory, 'TTLNu') \
         + nanoGetSampleFiles(mcDirectory, 'TTW-WtoQQ')

samples['ttW'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 10
}

files = nanoGetSampleFiles(mcDirectory, 'TTHtoNon2B') \
        + nanoGetSampleFiles(mcDirectory, 'TTHto2B')

samples['ttH'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 10
}

######## TXX ##########
files = nanoGetSampleFiles(mcDirectory, 'THW') 

samples['tHW'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 10
}

files = nanoGetSampleFiles(mcDirectory, 'THQ') \

samples['tHQ'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 10
}

########## VV ###########
files = nanoGetSampleFiles(mcDirectory, 'WWTo2L2Nu')

samples['WW'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 5
}

files = nanoGetSampleFiles(mcDirectory, 'WZ')

samples['WZ'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 20
}

files = nanoGetSampleFiles(mcDirectory, 'ZZ')

samples['ZZ'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 10
}

######## VVV #########

files = nanoGetSampleFiles(mcDirectory, 'WWW')

samples['WWW'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 20
}

files = nanoGetSampleFiles(mcDirectory, 'WWZ')

samples['WWZ'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 20
}

files = nanoGetSampleFiles(mcDirectory, 'WZZ')

samples['WZZ'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 20
}

files = nanoGetSampleFiles(mcDirectory, 'ZZZ')

samples['ZZZ'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 20
}

######## Other ########

files = []

for label in [
        'DYGto2LG-1Jets_Bin-MLL-4to50',
        'DYGto2LG-1Jets_Bin-MLL-50',
]:
    files += nanoGetSampleFiles(mcDirectory, label)

files += nanoGetSampleFiles(mcDirectory, "TZQB-ZTo2L")

#files += nanoGetSampleFiles(mcDirectory, "TTG_PTG-10to100")
#files += nanoGetSampleFiles(mcDirectory, "TTG_PTG-100to200")
#files += nanoGetSampleFiles(mcDirectory, "TTG_PTG-200")

samples['Other'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 20,
}

###########################################
###############  SIGNALS  #################
###########################################

mPhi = ['10','50','100','150', '200', '250', '300', '350', '400', '500', '600', '700', '800', '1000']

#mPhi = ['300']

# tt+DM dilepton scalar
for phi in mPhi:
    samples[f'TTto2LDMsimpSpin0_s_mphi-{phi}'] = {
            'name': nanoGetSampleFiles(signalDirectory, f'TTto2LDMsimpSpin0_s_mphi-{phi}') + nanoGetSampleFiles(signalDirectory, f'TTto2LDMsimpSpin0_s_mphi-{phi}_ext1'),
            'weight': mcCommonWeight,
            'FilesPerJob': 5
            }

# tt+DM dilepton pseudoscalar
for phi in mPhi:
    samples[f'TTto2LDMsimpSpin0_ps_mphi-{phi}'] = {
            'name': nanoGetSampleFiles(signalDirectory, f'TTto2LDMsimpSpin0_ps_mphi-{phi}') + nanoGetSampleFiles(signalDirectory, f'TTto2LDMsimpSpin0_ps_mphi-{phi}_ext1'),
            'weight': mcCommonWeight,
            'FilesPerJob': 5
            }

# tW+DM dilepton scalar
for phi in mPhi:
    samples[f'TWto2LDMsimpSpin0_s_mphi-{phi}'] = {
            'name': nanoGetSampleFiles(signalDirectory, f'TWto2LDMsimpSpin0_s_mphi-{phi}'), 
            'weight': mcCommonWeight,
            'FilesPerJob': 5
            }

# tW+DM dilepton pseudoscalar
for phi in mPhi:
    samples[f'TWto2LDMsimpSpin0_ps_mphi-{phi}'] = {
            'name': nanoGetSampleFiles(signalDirectory, f'TWto2LDMsimpSpin0_ps_mphi-{phi}'),
            'weight': mcCommonWeight,
            'FilesPerJob': 5
            }

###########################################
################## FAKE ###################
###########################################

samples['Fake'] = {
    'name': [],
    'weight': 'METFilter_DATA*fakeW',
    'weights': [],
    'isData': ['all'],
    'FilesPerJob': 50
}

for _, sd in DataRun:
    for pd in DataSets:
        tag = f"{pd}_{sd}"
        files = nanoGetSampleFiles(fake_dir_map[pd], tag)

        samples['Fake']['name'].extend(files)
        addSampleWeight(samples, 'Fake', tag, DataTrig[pd])


###########################################
################## DATA ###################
###########################################

samples['DATA'] = {
    'name': [],
    'weight': 'LepWPCut*METFilter_DATA',
    'weights': [],
    'isData': ['all'],
    'FilesPerJob': 70
}

for _, sd in DataRun:
    for pd in DataSets:
        datatag = f"{pd}_{sd}"
        files = nanoGetSampleFiles(dir_map[pd], datatag)

        print(datatag)

        samples['DATA']['name'].extend(files)
        addSampleWeight(samples, 'DATA', datatag, DataTrig[pd])

