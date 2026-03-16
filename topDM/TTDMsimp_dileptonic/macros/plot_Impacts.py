#!/usr/bin/env python3
import os
import subprocess

os.environ["PYTHONNOUSERSITE"] = "1"

# ================================
# Configuration
# ================================
datacards_dir = "datacards_files"
root_dir = "root_files"
limit_dir = "limit_files"
parallel_jobs = 8  # number of parallel Combine fits
mphi_values = [800]  # masses to run impacts for

# ================================
# Ensure directories exist
# ================================
os.makedirs(datacards_dir, exist_ok=True)
os.makedirs(root_dir, exist_ok=True)
os.makedirs(limit_dir, exist_ok=True)

# ================================
# Function to run impacts for one mass
# ================================
def run_impacts(mphi, parallel=8):
    """
    Run full CMS nuisance parameter impacts workflow for a single mass point
    """
    # Filenames
    comb_card_filename = os.path.join(datacards_dir, f"ttDM_CombCard_mphi_{mphi}.txt")
    root_filename = os.path.join(root_dir, f"ttDM_mphi_{mphi}.root")
    impacts_json = os.path.join(limit_dir, f"impacts_mphi_{mphi}.json")
    output_plot = os.path.join(limit_dir, f"impacts_mphi_{mphi}")

    print(f"[{mphi}] === Running impacts workflow ===")

    # Step 1: Initial fit
    print(f"[{mphi}] Step 1: Initial fit")
    subprocess.run([
        "combineTool.py",
        "-M", "Impacts",
        "-d", root_filename,
        "-m", str(mphi),
        "--doInitialFit",
        "--robustFit", "1"
    ], check=True)

    # Step 2: Fit each nuisance parameter
    print(f"[{mphi}] Step 2: Per-nuisance fits (parallel={parallel})")
    subprocess.run([
        "combineTool.py",
        "-M", "Impacts",
        "-d", root_filename,
        "-m", str(mphi),
        "--doFits",
        "--robustFit", "1",
        "--parallel", str(parallel)
    ], check=True)

    # Step 3: Collect results into JSON
    print(f"[{mphi}] Step 3: Collect results into JSON -> {impacts_json}")
    subprocess.run([
        "combineTool.py",
        "-M", "Impacts",
        "-d", root_filename,
        "-m", str(mphi),
        "-o", impacts_json
    ], check=True)

    # Step 4: Plot impacts
    print(f"[{mphi}] Step 4: Plot impacts -> {output_plot}.pdf / .png")
    subprocess.run([
        "plotImpacts.py",
        "-i", impacts_json,
        "-o", output_plot
    ], check=True)

    print(f"[{mphi}] === Done ===\n")

# ================================
# Loop over all mass points
# ================================
if __name__ == "__main__":
    for mphi in mphi_values:
        run_impacts(mphi, parallel=parallel_jobs)

    print("All masses processed. Check JSONs and plots in:", limit_dir)
