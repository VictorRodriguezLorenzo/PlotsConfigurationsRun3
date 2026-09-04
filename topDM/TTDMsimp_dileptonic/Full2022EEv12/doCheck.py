import sys
import argparse
import os
import re
import subprocess

from configuration import tag


def defaultParser():
    parser = argparse.ArgumentParser(
        description="Submit specific jobs or all jobs if no specific jobs are provided"
    )

    parser.add_argument(
        "-Sub",
        "--Submit",
        action="store_true",
        help="Submit files",
        default=False,
    )

    parser.add_argument(
        "-Jobs",
        "--Jobs",
        nargs="+",
        help="Specific job labels to submit",
    )

    parser.add_argument(
        "-condor_q",
        "--condor_q",
        nargs="*",
        default=None,
        help="Job ID to check running jobs on condor. If used without IDs, checks all jobs.",
    )

    return parser


def run(submit=False, specific_jobs=None, condor_q=None, tag=tag):
    prePath = os.path.abspath(os.path.dirname(__file__))

    if "examples" in prePath:
        prePath = prePath.split("examples/")[0]   ## Assume you work in processor folder

    output_path = f"/eos/user/v/victorr/www/Run3-ttDM/rootFiles/{tag}/"
    path = prePath + "/condor/EGamma_Run2022G-Prompt-v1/"

    jobDir = path

    cmd = "find {} -type d -name '*'".format(path)

    fnames = subprocess.check_output(cmd, shell=True).strip().split(b"\n")
    fnames = [
        fname.decode("ascii").split("EGamma_Run2022G-Prompt-v1/")[1]
        for fname in fnames
    ]

    failed_jobs = []
    error_files = []
    script_files = []

    tag_total = {}
    tag_failed = {}

    for fname in fnames:
        if fname == "":
            continue

        tag_sample = fname.rsplit("_", 1)[0]

        tag_total[tag_sample] = tag_total.get(tag_sample, 0) + 1

        file_name = output_path + f"mkShapes__{tag}__ALL__" + fname + ".root"
        error_file = jobDir + fname + "/" + "err.txt"
        script_file = jobDir + fname + "/" + "script.py"

        if os.path.exists(file_name):
            continue

        if specific_jobs is None or fname in specific_jobs:
            print("LABEL: " + fname)
            failed_jobs.append(fname)
            error_files.append(error_file)
            script_files.append(script_file)

            tag_failed[tag_sample] = tag_failed.get(tag_sample, 0) + 1

    print("=========================")

    if len(fnames) > 0:
        failed_percentage = round(100 * len(failed_jobs) / len(fnames), 2)
    else:
        failed_percentage = 0.0

    print(
        "Ratio of failed jobs: "
        + str(len(failed_jobs))
        + "/"
        + str(len(fnames))
        + " = "
        + str(failed_percentage)
        + "%"
    )

    print("\nFailure percentage per tag:\n")

    GREEN = "\033[92m"
    RED = "\033[91m"
    RESET = "\033[0m"

    for tag_sample in sorted(tag_total):
        total = tag_total[tag_sample]
        failed = tag_failed.get(tag_sample, 0)
        percentage = 100 * failed / total

        text = f"{tag_sample}: {failed}/{total} = {percentage:.2f}%"

        if percentage == 0:
            print(f"{GREEN}{text}{RESET}")
        elif percentage == 100:
            print(f"{RED}{text}{RESET}")
        else:
            print(text + "")

    jobs_to_submit = failed_jobs

    if condor_q is not None:
        running_jobs = set()

        if len(condor_q) == 0:
            cmd = [
                "condor_q",
                "-af",
                "Args",
            ]
        else:
            cluster_ids = [qid.split(".")[0] for qid in condor_q]
            cluster_ids = list(dict.fromkeys(cluster_ids))

            constraint = " || ".join(
                ["ClusterId == {}".format(cluster_id) for cluster_id in cluster_ids]
            )

            cmd = [
                "condor_q",
                "-constraint",
                constraint,
                "-af",
                "Args",
            ]

        running_jobs_output = subprocess.check_output(cmd).decode("utf-8")

        for line in running_jobs_output.splitlines():
            job = line.strip().strip('"').strip("'")

            if job:
                running_jobs.add(job.split()[0])

        jobs_to_submit = [
            job for job in failed_jobs
            if job not in running_jobs
        ]

        failed_jobs_still_running = [
            job for job in failed_jobs
            if job in running_jobs
        ]

        print("Failed jobs not running on condor:", " ".join(jobs_to_submit))

    if submit:
        if condor_q is not None:
            print("Submitting only failed jobs that are NOT currently running on condor.")
        else:
            print("Submitting all failed jobs.")

    if submit:
        if not jobs_to_submit:
            print("No jobs to submit.")
            return

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

        resubmit = resubmit.replace("RPLME_ALLSAMPLES", " ".join(jobs_to_submit))

        with open(jobDir + "submit_failed.jdl", "w") as f:
            f.write(resubmit)

        proc = subprocess.Popen(
            f"cd {jobDir}; condor_submit submit_failed.jdl;",
            shell=True,
        )

        proc.wait()


if __name__ == "__main__":
    parser = defaultParser()
    args = parser.parse_args()

    doSubmit = args.Submit
    specificJobs = args.Jobs
    condor_q = args.condor_q

    run(doSubmit, specificJobs, condor_q)

