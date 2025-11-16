from mkShapesRDF.lib.search_files import SearchFiles

searchFiles = SearchFiles()

redirector = ""
useXROOTD = False

mcProduction = 'Summer24_150x_nAODv15_Full2024v15'
mcSteps      = 'MCl2loose2024v15__MCCorr2024v15__JERFrom23BPix__l2tight'
dataReco     = 'Run2024_ReRecoCDE_PromptFGHI_nAODv15_Full2024v15'
dataSteps    = 'DATAl2loose2024v15__sblancof__l2loose'

##############################################
###### Tree base directory for the site ######
##############################################
#treeBaseDir = '/eos/user/v/victorr'
treeBaseDir = '/eos/cms/store/group/phys_higgs/cmshww/calderon/HWWNano/'
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
# fakeDirectory = os.path.join(treeBaseDir, dataReco, fakeSteps)
dataDirectory = os.path.join(treeBaseDir, dataReco, dataSteps)
fakeDirectory = dataDirectory

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
    'MuonEG'         : ' Trigger_ElMu' ,
    #'SingleMuon'     : '!Trigger_ElMu && Trigger_sngMu' ,
    'Muon0'           : '!Trigger_ElMu && (Trigger_sngMu || Trigger_dblMu)',
    'Muon1'           : '!Trigger_ElMu && (Trigger_sngMu || Trigger_dblMu)',
    'EGamma0'         : '!Trigger_ElMu && !Trigger_sngMu && !Trigger_dblMu && (Trigger_sngEl || Trigger_dblEl)',
    'EGamma1'         : '!Trigger_ElMu && !Trigger_sngMu && !Trigger_dblMu && (Trigger_sngEl || Trigger_dblEl)',
} 

#########################################
############ MC COMMON ##################
#########################################

mcCommonWeightNoMatch = 'XSWeight*METFilter_Common*SFweight'
mcCommonWeight        = 'XSWeight*METFilter_Common*PromptGenLepMatch2l*SFweight'

###########################################
#############  BACKGROUNDS  ###############
###########################################

########## DY #########
files = []
for label in [
        'DYto2E-2Jets_MLL-50',
        'DYto2Mu-2Jets_MLL-50',
        'DYto2Tau-2Jets_MLL-50',
        'DYto2E-2Jets_MLL-10to50',
        'DYto2Mu-2Jets_MLL-10to50',
        'DYto2Tau-2Jets_MLL-10to50']:
    
    files += nanoGetSampleFiles(mcDirectory, label)

samples['DY'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 15,
}

########## Single top #########
files = nanoGetSampleFiles(mcDirectory, 'ST_t-channel_top') + nanoGetSampleFiles(mcDirectory, 'ST_t-channel_antitop')

samples['ST'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 15
}
##########MISSING SAMPLES
#nanoGetSampleFiles(mcDirectory, 'TWminusto2L2Nu') + nanoGetSampleFiles(mcDirectory, 'TbarWplusto2L2Nu')

########## TTTo2L2Nu #########
files = nanoGetSampleFiles(mcDirectory, 'TTTo2L2Nu')

samples['TTTo2L2Nu'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 15
}

###### TTToSemiLeptonic #####
#files = nanoGetSampleFiles(mcDirectory, 'TTToSemiLeptonic')
#
#samples['TTToSemiLeptonic'] = {
#    'name': files,
#    'weight': mcCommonWeight,
#    'FilesPerJob': 5
#}

########## VV ###########
files = nanoGetSampleFiles(mcDirectory, 'WWTo2L2Nu')

samples['WW'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 15
}

files = nanoGetSampleFiles(mcDirectory, 'WZ')

samples['WZ'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 15
}

files = nanoGetSampleFiles(mcDirectory, 'ZZ')

samples['ZZ'] = {
    'name': files,
    'weight': mcCommonWeight,
    'FilesPerJob': 15
}

######## Other ########
#files = nanoGetSampleFiles(mcDirectory, 'ZGToLLG')
#
#samples['Other'] = {
#    'name': files,
#    'weight': mcCommonWeight,
#    'FilesPerJob': 5,
#}

###########################################
###############  SIGNALS  #################
###########################################
'''
mPhi = ['10','50','100','150', '200', '250', '300', '350', '400', '500', '600' '700', '800', '1000']

# tt+DM dilepton scalar
for phi in mPhi:
    samples[f'TTto2LDMsimpSpin0_s_mchi_{phi}'] = {
            'name': nanoGetSampleFiles(signalDirectory, f'TTto2LDMsimpSpin0_s_mchi_{phi}'),
            'weight': mcCommonWeight,
            'FilesPerJob': 1
            }

# tt+DM inclusive scalar
for phi in mPhi:
    samples[f'TTDMsimpSpin0_s_mchi_{phi}'] = {
            'name': nanoGetSampleFiles(signalDirectory, f'TTDMsimpSpin0_s_mchi_{phi}'),
            'weight': mcCommonWeight,
            'FilesPerJob': 1
            }

# tt+DM dilepton pseudoscalar
for phi in mPhi:
    samples[f'TTto2LDMsimpSpin0_ps_mchi_{phi}'] = {
            'name': nanoGetSampleFiles(signalDirectory, f'TTto2LDMsimpSpin0_s_mchi_{phi}'),
            'weight': mcCommonWeight,
            'FilesPerJob': 1
            }

# tt+DM inclusive pseudoscalar
for phi in mPhi:
    samples[f'TTDMsimpSpin0_ps_mchi_{phi}'] = {
            'name': nanoGetSampleFiles(signalDirectory, f'TTDMsimpSpin0_s_mchi_{phi}'),
            'weight': mcCommonWeight,
            'FilesPerJob': 1
            }
'''
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
    tag = pd + '_' + sd

    files = nanoGetSampleFiles(fakeDirectory, tag)

    samples['Fake']['name'].extend(files)
    addSampleWeight(samples, 'Fake', tag, DataTrig[pd])

    
###########################################
################## DATA ###################
###########################################

samples['DATA'] = {
    'name': [],
    'weight': 'LepWPCut*METFilter_DATA',
    #'weight': 'METFilter_DATA',
    'weights': [],
    'isData': ['all'],
    'FilesPerJob': 70
}

for _, sd in DataRun:
  for pd in DataSets:
    datatag = pd + '_' + sd

    files = nanoGetSampleFiles(dataDirectory, datatag)

    print(datatag)

    samples['DATA']['name'].extend(files)
    addSampleWeight(samples, 'DATA', datatag, DataTrig[pd])
