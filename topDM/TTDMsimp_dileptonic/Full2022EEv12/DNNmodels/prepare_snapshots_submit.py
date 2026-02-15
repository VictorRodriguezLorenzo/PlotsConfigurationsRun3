import os
import subprocess

from mkShapesRDF.lib.search_files import SearchFiles

# ============================================================
# SAME FILE DISCOVERY STRUCTURE AS ORIGINAL
# ============================================================

searchFiles = SearchFiles()

limitFiles = -1
redirector = ""

def nanoGetSampleFiles(path, name):
    _files = searchFiles.searchFiles(path, name, redirector=redirector)
    if limitFiles != -1:
        _files = _files[:limitFiles]
    return _files


signalDirectory = '/eos/user/v/victorr/HWWNano/Summer22EE_130x_nAODv12_Full2022v12/MCl2loose2022EEv12__MCCorr2022EEv12JetScaling__l2tight'

mcDirectory = '/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer22EE_130x_nAODv12_Full2022v12/MCl2loose2022EEv12__MCCorr2022EEv12JetScaling__l2tight'

mPhi = ['10','50','100','150','200','250','300','350','400','500','600','700','800','1000']

files_ttDM_temp = {}

# tt+DM dilepton scalar
for phi in mPhi:
    files_ttDM_temp[f'TTto2LDMsimpSpin0_s_mphi-{phi}'] = {
        'names': nanoGetSampleFiles(signalDirectory, f'TTto2LDMsimpSpin0_s_mphi-{phi}')
    }

# tt+DM inclusive scalar
for phi in mPhi:
    files_ttDM_temp[f'TTDMsimpSpin0_s_mphi-{phi}'] = {
        'names': nanoGetSampleFiles(signalDirectory, f'TTDMsimpSpin0_s_mphi-{phi}')
    }

# tt+DM dilepton pseudoscalar
for phi in mPhi:
    files_ttDM_temp[f'TTto2LDMsimpSpin0_ps_mphi-{phi}'] = {
        'names': nanoGetSampleFiles(signalDirectory, f'TTto2LDMsimpSpin0_s_mphi-{phi}')
    }

# tt+DM inclusive scalar
for phi in mPhi:
    files_ttDM_temp[f'TTDMsimpSpin0_ps_mphi-{phi}'] = {
        'names': nanoGetSampleFiles(signalDirectory, f'TTDMsimpSpin0_s_mphi-{phi}')
    }

files_ttDM = files_ttDM_temp.copy()
for key, value in list(files_ttDM_temp.items()):
    if isinstance(value["names"], list) and len(value["names"]) == 0:
        del files_ttDM[key]

files_top = nanoGetSampleFiles(mcDirectory, 'TTTo2L2Nu') + \
    nanoGetSampleFiles(mcDirectory,'ST_t-channel_top') + \
    nanoGetSampleFiles(mcDirectory,'ST_t-channel_antitop')

files_BKG = files_top

# ============================================================
# FLATTEN LIST
# ============================================================

all_files = []

for key, val in files_ttDM.items():
    all_files.extend(val["names"])

all_files.extend(files_BKG)

print("Submitting jobs for", len(all_files), "files")

# ============================================================
# CREATE WRAPPER SCRIPT
# ============================================================

os.makedirs("condor_logs", exist_ok=True)

with open("run_snapshot.sh", "w") as f:
    f.write("""#!/bin/bash
echo 'first source of start.sh'; source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt/setup.sh
source /afs/cern.ch/user/v/victorr/private/mkShapesRDF/myenv/bin/activate
export STARTPATH=/afs/cern.ch/user/v/victorr/private/mkShapesRDF/start.sh
export PATH=/afs/cern.ch/user/v/victorr/private/mkShapesRDF/utils/bin:$PATH
export PYTHONPATH=/afs/cern.ch/user/v/victorr/private/mkShapesRDF/myenv/lib64/python3.11/site-packages:$PYTHONPATH
export X509_USER_PROXY=$HOME/.proxy

time python prepare_training_snapshots.py $1
""")

os.chmod("run_snapshot.sh", 0o755)

# ============================================================
# CREATE SUBMIT FILE
# ============================================================

with open("submit_snapshot.sub", "w") as f:
    f.write("""
universe = vanilla
executable  = run_snapshot.sh
arguments   = $(inputfile)

should_transfer_files = YES
transfer_input_files = /afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2024v15/DNNmodels/prepare_training_snapshots.py

output      = condor_logs/job_$(Cluster)_$(Process).out
error       = condor_logs/job_$(Cluster)_$(Process).err
log         = condor_logs/job_$(Cluster).log

request_cpus = 4
requirements = (OpSysAndVer =?= "AlmaLinux9")
+JobFlavour = "nextweek"

queue inputfile from (
""")

    for file in all_files:
        f.write(file + "\n")

    f.write(")\n")

# ============================================================
# SUBMIT
# ============================================================

subprocess.run(["condor_submit", "submit_snapshot.sub"])
