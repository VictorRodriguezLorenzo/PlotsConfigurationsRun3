import sys
import argparse
import os
import subprocess

def defaultParser():
    parser = argparse.ArgumentParser(description="Submit specific jobs or all jobs if no specific jobs are provided")
    
    parser.add_argument(
        "-Sub",
        "--Submit",
        action='store_true',
        help="Submit files",
        default=False,
    )
    
    parser.add_argument(
        "-Jobs",
        "--Jobs",
        nargs='+',
        help="Specific job labels to submit",
    )

    parser.add_argument(
        "-condor_q",
        "--condor_q",
        type=str,
        help="Job ID to check running jobs on condor",
    )
    
    return parser

def run(submit=False, specific_jobs=None, condor_q=None):
    prePath = os.path.abspath(os.path.dirname(__file__))

    if "examples" in prePath:
        prePath = prePath.split("examples/")[0]   ## Assume you work in processor folder
        
    path = prePath + "/condor/EGamma1_Run2024I-Prompt-v1/"
    output_path = "/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp/Full2024v15/rootFiles/ttDM_2024/"
    jobDir = path

    cmd = "find {} -type d -name '*'".format(path)
    
    fnames = subprocess.check_output(cmd, shell=True).strip().split(b'\n')
    fnames = [fname.decode('ascii').split("EGamma1_Run2024I-Prompt-v1/")[1] for fname in fnames] 
    
    failed_jobs = []
    error_files = []
    script_files = []
    
    for fname in fnames:
        file_name = output_path + "mkShapes__ttDM_2024__ALL__" + fname + ".root"
        error_file = jobDir + fname + "/" + "err.txt"
        script_file = jobDir + fname + "/" + "script.py"

        if os.path.exists(file_name) or fname=="":
            continue
        else:
            if specific_jobs is None or fname in specific_jobs:
                print("ERROR: File does not exist in output folder")
                print("LABEL: " + fname)
                failed_jobs.append(fname)
                error_files.append(error_file)
                script_files.append(script_file)

    print("=========================")
    print("Ratio of failed jobs: " + str(len(failed_jobs)) + "/" + str(len(fnames)) + " = " + str(round(100*len(failed_jobs)/len(fnames), 2)) + "%")
    
    if submit:
        resubmit = """
universe = vanilla
executable = run.sh
arguments = $(Folder)
should_transfer_files = YES
transfer_input_files = $(Folder)/script.py, /afs/cern.ch/user/v/victorr/private/mkShapesRDF/mkShapesRDF/include/headers.hh, /afs/cern.ch/user/v/victorr/private/mkShapesRDF/mkShapesRDF/shapeAnalysis/runner.py
output = $(Folder)/out.txt
error  = $(Folder)/err.txt
log    = $(Folder)/log.txt
request_cpus   = 1
+JobFlavour = "nextweek"
queue 1 Folder in  RPLME_ALLSAMPLES"""
        
        resubmit = resubmit.replace("RPLME_ALLSAMPLES", " ".join(failed_jobs))
        
        with open(jobDir + "submit_failed.jdl", "w") as f:
            f.write(resubmit)
            
        proc = subprocess.Popen(
            f"cd {jobDir}; condor_submit submit_failed.jdl;", shell=True
        )
        
        proc.wait()

    if condor_q:
        cmd = f"condor_q -nobatch {condor_q}"
        running_jobs_output = subprocess.check_output(cmd, shell=True).decode('utf-8')
        running_jobs = [line.split()[-1] for line in running_jobs_output.split('\n') if line.startswith(condor_q)]
        # Compare failed_jobs with running_jobs to find failed jobs not running on condor
        failed_jobs_not_running = [job for job in failed_jobs if job not in running_jobs]
        print("Failed jobs not running on condor:", failed_jobs_not_running)


if __name__ == "__main__":
    parser = defaultParser()
    args = parser.parse_args()
    
    doSubmit = args.Submit
    specificJobs = args.Jobs
    condor_q = args.condor_q

    run(doSubmit, specificJobs, condor_q)

