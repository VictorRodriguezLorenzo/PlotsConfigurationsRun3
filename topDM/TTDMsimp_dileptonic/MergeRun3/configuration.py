import sys,os

#: tag used to identify the configuration folder version
tag = "ttDM_dilep_Run3"

#: file to use as runner script, default uses mkShapesRDF.shapeAnalysis.runner, otherwise specify path to script
runnerFile = "default"

#: output file name
outputFile = "mkShapes__{}.root".format(tag)

#: path to output folder
outputFolder = "/eos/user/" + os.getlogin()[0] + "/" + os.getlogin() + "/www/Run3-ttDM/rootFiles/" + tag
#outputFolder = "rootFiles/" + tag

# path to batch folder (used for condor submission)
batchFolder = "condor"

# path to configuration folder (will contain all the compiled configuration files)
configsFolder = "configs"

# luminosity to normalize to (in 1/fb), Run 3 total of 2022, 2022EE, 2023, 2023BPix, 2024 and 2025
lumi = 281.73

# file with dict of aliases to define
aliasesFile = "aliases.py"

# file with dict of variables
variablesFile = "variables.py"

# file with dict of cuts
cutsFile = "cuts.py"

# file with dict of samples
samplesFile = "samples.py"

# file with dict of samples
plotFile = "plot.py"

# file with dict of structure (used to define combine processes)
structureFile = "structure.py"

# nuisances file for mkDatacards and for mkShape
nuisancesFile = "nuisances_ALL.py"

# path to folder where to save plots
plotPath = "/eos/user/" + os.getlogin()[0] + "/" + os.getlogin() + "/www/Plots" + tag
#plotPath = "Plots" + tag

# this lines are executed right before the runner on the condor node
mountEOS = [
]

# list of imports to import when compiling the whole configuration folder, it should not contain imports used by configuration.py
imports = ["os", "glob", ("collections", "OrderedDict"), "ROOT"]

# list of files to compile
filesToExec = [
    samplesFile,
    aliasesFile,
    cutsFile,
    variablesFile,
    plotFile,
    nuisancesFile,
    structureFile,
]

# list of variables to keep in the compiled configuration folder
varsToKeep = [
    "batchVars",
    "outputFolder",
    "batchFolder",
    "configsFolder",
    "outputFile",
    "runnerFile",
    "tag",
    "samples",
    "aliases",
    "variables",
    ("cuts", {"cuts": "cuts", "preselections": "preselections"}),
    ("plot", {"plot": "plot", "groupPlot": "groupPlot", "legend": "legend"}),
    "nuisances",
    "structure",
    "lumi",
]

# list of variables to keep in the batch submission script (script.py)
batchVars = varsToKeep[varsToKeep.index("samples") :]

varsToKeep += ['plotPath']
