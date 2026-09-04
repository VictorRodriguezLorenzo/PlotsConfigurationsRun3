#!/usr/bin/env python3
"""
Build TTDMsimp dileptonic datacards, workspaces, limits, and fit jobs.

This script is the main combine helper for one TTDMsimp dileptonic campaign at a
time.  Run the workflow once per campaign by pointing ``--work-dir`` at that campaign directory (or by running
from inside the campaign directory with ``--work-dir .``).


Prerequisites: do I need combine?
---------------------------------
Use the official CMS Combine instructions for the exact recommended release/tag.
For lxplus/EL9, the current Combine documentation recommends installing Combine
v10 inside a CMSSW 14_1_X release.  A typical setup is:

  ``source /cvmfs/cms.cern.ch/cmsset_default.sh``
  ``cmsrel CMSSW_14_1_0_pre4``
  ``cd CMSSW_14_1_0_pre4/src``
  ``cmsenv``
  ``git -c advice.detachedHead=false clone --depth 1 --branch v10.6.0 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit``
  ``cd HiggsAnalysis/CombinedLimit``
  ``scramv1 b clean; scramv1 b -j$(nproc --ignore=2)``

After that, use the CMSSW ``src`` directory as ``--cmssw-dir`` when making
Condor scripts, for example:
  ``--cmssw-dir /path/to/CMSSW_14_1_0_pre4/src``

Minimal command sequence
------------------------
From ``topDM/TTDMsimp_dileptonic`` for one campaign, for example ``Full2024v15``:
  ``CAMPAIGN=Full2024v15``

First create and submit the workspace jobs:
  ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task WS --condor --cmssw-dir /path/to/CMSSW_X_Y_Z/src``
  ``condor_submit ${CAMPAIGN}/datacards_files/HTCondor/combine_jobs.submit``

After those jobs finish and ``${CAMPAIGN}/root_files/Workspace_mphi_<M>.root``
exists, create and submit only the diagnostic you want, for example impacts:
  ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task Impacts --condor --cmssw-dir /path/to/CMSSW_X_Y_Z/src``
  ``condor_submit ${CAMPAIGN}/diagnostics_files/HTCondor/combine_jobs.submit``

Replace ``Impacts`` with ``Limits``, ``GoF``, ``AsimovFit``, ``CRonlyFit``, or
``LogLikeScan`` to run only that step.  Use ``--masses 150`` while testing one
mass, then omit ``--masses`` for the full default scan.

What is combined for each mass
------------------------------
For every mediator mass in ``DEFAULT_MASSES`` (or in ``--masses``), the script
builds one combined datacard from the same four categories inside the campaign's
input datacard directory (default: ``${CAMPAIGN}/datacards``):
  * ``SR1b``: ``datacards/ttdm_sr_2l_1b/{dnn_variable}``
  * ``SR2b``: ``datacards/ttdm_sr_2l_2b/{dnn_variable}``
  * ``DYCR``: ``datacards/dycr_1b/events``
  * ``ttZCR``: ``datacards/ttZcr_Inc/events``
  * ``ttCR``: ``datacards/ttcr_Inc/events``

This matches the campaign layout with directories such as ``datacards/dycr_1b``,
``datacards/dycr_2b``, ``datacards/ttdm_sr_2l_1b``,
``datacards/ttdm_sr_2l_2b``, ``datacards/ttcr_Inc``, and
``datacards/ttZcr_Inc``.  Only the four categories above are combined by
default; the other directories are left untouched.  If your input cards are not
under ``datacards``, override the base with ``--input-datacards-dir``.

Typical per-campaign Condor workflow
------------------------------------
Run these commands from ``topDM/TTDMsimp_dileptonic``.  Repeat the sequence for
each campaign by changing ``CAMPAIGN``.

Example campaign selection:
  ``CAMPAIGN=Full2022v12``

1. Create the workspace-production Condor scripts:
   ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task WS --condor``

   This creates ``${CAMPAIGN}/datacards_files/HTCondor/combine_job_*.sh`` and
   ``${CAMPAIGN}/datacards_files/HTCondor/combine_jobs.submit``.

2. Submit the workspace jobs and wait until they finish:
   ``condor_submit ${CAMPAIGN}/datacards_files/HTCondor/combine_jobs.submit``

   Successful jobs produce, for each mass:
   ``${CAMPAIGN}/datacards_files/ttDM_CombCard_mphi_<M>.txt`` and
   ``${CAMPAIGN}/root_files/Workspace_mphi_<M>.root``.

3. After the workspaces exist, create fit/diagnostic Condor scripts.  To prepare
   all default diagnostics in one submit file, run:
   ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task TEST --condor``

   To prepare only one diagnostic, use that diagnostic as ``--task`` instead:
   ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task Impacts --condor``
   ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task Limits --condor``
   ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task GoF --condor``

   This creates ``${CAMPAIGN}/diagnostics_files/HTCondor/combine_job_*.sh`` and
   ``${CAMPAIGN}/diagnostics_files/HTCondor/combine_jobs.submit`` for only the
   requested diagnostic(s).  ``--task TEST`` uses the comma-separated list in
   ``--tests``; by default that list is AsimovFit, CRonlyFit, GoF, Impacts,
   Limits, and LogLikeScan.

4. Submit the fit/diagnostic jobs:
   ``condor_submit ${CAMPAIGN}/diagnostics_files/HTCondor/combine_jobs.submit``

   Successful jobs write outputs below ``${CAMPAIGN}/diagnostics_files/mphi_<M>/``:
   ``AsimovFit/FitDiagnostics_Asimov_mphi_<M>.root``,
   ``CRonlyFit/FitDiagnostics_CR-only_mphi_<M>.root``, GoF ROOT files,
   ``Impacts/Impacts_Asimov_mphi_<M>_r*.json``,
   ``Limits/Limits_Expected_mphi_<M>.root`` by default, and likelihood-scan ROOT
   files.

Running only one diagnostic
---------------------------
The diagnostic tasks can be separated without editing the code:
  * Asimov fit only: ``--task AsimovFit --condor``
  * CR-only fit only: ``--task CRonlyFit --condor``
  * GoF only: ``--task GoF --condor``
  * Impacts only: ``--task Impacts --condor`` (also accepts ``IMPACTS``)
  * Condor expected-limit jobs only: ``--task Limits --condor`` (or ``LIMITS --condor``)
  * Likelihood scan only: ``--task LogLikeScan --condor``

``--task TEST --condor --tests A,B`` is still useful when you want a custom
combination of diagnostics in one submit file, for example:
  ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task TEST --condor --tests Impacts,Limits``

Quick checks and reduced jobs
-----------------------------
* Print workspace commands without running combine:
  ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task WS --dry-run``
* Create workspace scripts for one mass only while checking paths:
  ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task WS --condor --masses 150``
* Create impacts scripts for one mass only:
  ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task Impacts --condor --masses 150``

Plotting limits
---------------
The built-in ``--plot`` option is only connected to the local ``LIMITS``/``ALL``
path, not to Condor diagnostic jobs.  For example:
  ``python3 limits/mkLimits.py --work-dir ${CAMPAIGN} --task LIMITS --plot``

This runs ``combine -M AsymptoticLimits --run expected`` locally for each mass by
default, writes text logs under ``${CAMPAIGN}/limit_files/Limit_mphi_<M>.txt``,
and then calls ``brazil_band.plotUpperLimits`` on those logs.  The plotting
helper saves ``BrazilBplot.png`` in the directory where the command is executed.
If the ROOT Python bindings or ``brazil_band.py`` are not available, limits can
still be run without ``--plot`` and plotting can be done later from the produced
logs.  

IMPORTANT: Use ``--limit-run blind`` for fully blinded Asimov limits, or ``--limit-run observed`` only when you really want observed limits.
"""


from __future__ import annotations

import argparse
import ast
import importlib.util
import os
import random
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

plotUpperLimits = None
if importlib.util.find_spec("ROOT") and importlib.util.find_spec("brazil_band"):
    from brazil_band import plotUpperLimits

os.environ["PYTHONNOUSERSITE"] = "1"

DEFAULT_WORK_DIR = Path(".")
DEFAULT_INPUT_DATACARDS_DIR = "datacards"


DEFAULT_MASSES = [
    "10",
    "50",
    "100",
    "150",
    "200",
    "250",
    "300",
    "350",
    "400",
    "500",
    "600",
    "700",
    "800",
    "1000",
]

DEFAULT_MODELS = ["ps", "s"]

DEFAULT_CATEGORY_MODELS = {
    "SR1b": "TWto2LDMsimpSpin0_{model}",
    "SR2b": "TTto2LDMsimpSpin0_{model}",
}

DEFAULT_CARD_TEMPLATE = "datacard_{model}_mphi-{mass}.txt"

DEFAULT_CATEGORY_DNN_TEMPLATES = {
    "SR1b": "evaluate_dnn_tWDM_{model}_{mass}",
    "SR2b": "evaluate_dnn_ttDM_{model}_{mass}",
}

# Category paths are relative to DEFAULT_INPUT_DATACARDS_DIR/--input-datacards-dir.
# SR1b and SR2b use the mass-dependent DNN output directory, e.g.
# evaluate_dnn_mPhi_150 for mphi=150.
DEFAULT_CATEGORIES = {
    "SR1b": "ttdm_sr_2l_1b/{dnn_variable}",
    "SR2b": "ttdm_sr_2l_2b/{dnn_variable}",
    "DYCR": "dycr_1b/events",
    "ttZCR": "ttZcr_Inc/events",
    "ttCR": "ttcr_Inc/events",
}
DEFAULT_DNN_VARIABLE_TEMPLATE = "evaluate_dnn_mPhi_{mass}"
DEFAULT_SIGNAL_REGIONS = ["SR2b"] #["SR1b", "SR2b"]
DEFAULT_TESTS = ["AsimovFit", "CRonlyFit", "GoF", "Impacts", "Limits", "LogLikeScan"]
TASK_ALIASES = {
    "WS": "WS",
    "TEST": "TEST",
    "ALL": "ALL",
    "LIMITS": "LIMITS",
    "ASIMOV": "AsimovFit",
    "ASIMOVFIT": "AsimovFit",
    "CRONLY": "CRonlyFit",
    "CRONLYFIT": "CRonlyFit",
    "GOF": "GoF",
    "IMPACT": "Impacts",
    "IMPACTS": "Impacts",
    "LIMIT": "Limits",
    "LOGLIKESCAN": "LogLikeScan",
    "LLSCAN": "LogLikeScan",
    "SCAN": "LogLikeScan",
}


@dataclass(frozen=True)
class Job:
    name: str
    script: str


def quote(value: str | Path) -> str:
    return shlex.quote(str(value))


def run_command(command: str, *, dry_run: bool = False) -> None:
    print(f"\n$ {command}")
    if dry_run:
        return
    subprocess.run(command, shell=True, check=True)


def clean_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def parse_csv(values: str) -> list[str]:
    return [value.strip() for value in values.split(",") if value.strip()]


def normalize_task(task: str) -> str:
    if task in {"WS", "TEST", "LIMITS", "ALL", *DEFAULT_TESTS}:
        return task
    normalized = TASK_ALIASES.get(task.upper())
    if normalized:
        return normalized
    allowed = ["WS", "TEST", "LIMITS", "ALL", *DEFAULT_TESTS, "IMPACTS"]
    raise ValueError(f"Unsupported task {task!r}. Choose one of: {', '.join(allowed)}")


def validate_tests(tests: list[str]) -> list[str]:
    unknown = sorted(set(tests) - set(DEFAULT_TESTS))
    if unknown:
        raise ValueError(f"Unsupported --tests value(s): {', '.join(unknown)}. Choose from: {', '.join(DEFAULT_TESTS)}")
    return tests


def tests_for_task(args: argparse.Namespace) -> list[str]:
    if args.task == "TEST":
        return args.tests
    if args.task == "LIMITS":
        return ["Limits"]
    if args.task in DEFAULT_TESTS:
        return [args.task]
    return []


def dnn_variable_name(template: str, model: str, mass: str) -> str:
    return template.format(model=model, mass=mass)

def card_name(model: str, mass: str) -> str:
    return DEFAULT_CARD_TEMPLATE.format(model=model, mass=mass)

def combined_card_path(datacards_dir: Path, model: str, mass: str) -> Path:
    return datacards_dir / f"ttDM_CombCard_{model}_mphi_{mass}.txt"


def workspace_path(root_dir: Path, model: str, mass: str) -> Path:
    return root_dir / f"Workspace_{model}_mphi_{mass}.root"


def limit_log_path(limit_dir: Path, model: str, mass: str) -> Path:
    return limit_dir / f"Limit_{model}_mphi_{mass}.txt"

def format_lumi_label(lumi: float) -> str:
    return f"{lumi:g} fb^{{-1}} (13.6 TeV)"


def lumi_from_configuration(work_dir: Path) -> float | None:
    config_path = work_dir / "configuration.py"
    if not config_path.exists():
        return None
    tree = ast.parse(config_path.read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "lumi":
                    value = ast.literal_eval(node.value)
                    return float(value)
    return None


def campaign_lumi_label(work_dir: Path, override: str | None = None) -> str | None:
    if override:
        return override
    lumi = lumi_from_configuration(work_dir)
    if lumi is None:
        return None
    return format_lumi_label(lumi)

def input_datacards_path(work_dir: Path, input_datacards_dir: str) -> Path:
    input_dir = Path(input_datacards_dir)
    if input_dir.is_absolute():
        return input_dir
    return work_dir / input_dir

def combine_cards_command(
    input_dir: Path,
    datacards_dir: Path,
    model: str,
    mass: str,
    dnn_variable_template: str,
) -> str:
    parts = ["combineCards.py"]
    missing_cards = []

    for label, rel_dir_template in DEFAULT_CATEGORIES.items():

        # Category-dependent DNN variable
        category_dnn_template = DEFAULT_CATEGORY_DNN_TEMPLATES.get(
            label,
            dnn_variable_template,
        )
        dnn_variable = dnn_variable_name(
            category_dnn_template,
            model,
            mass,
        )
        rel_dir = rel_dir_template.format(
            mass=mass,
            dnn_variable=dnn_variable,
        )

        # Category-dependent signal model
        category_model = DEFAULT_CATEGORY_MODELS.get(
            label,
            "TTto2LDMsimpSpin0_{model}",
        ).format(model=model)

        card = input_dir / rel_dir / card_name(
            category_model,
            mass,
        )

        if not card.exists():
            missing_cards.append(card)

        parts.append(f"{label}={quote(card)}")

    if missing_cards:
        missing = "\n  ".join(str(card) for card in missing_cards)
        raise FileNotFoundError(
            f"Missing datacard(s) for mphi={mass}:\n  {missing}"
        )

    parts.append(
        f"> {quote(combined_card_path(datacards_dir, model, mass))}"
    )

    return " ".join(parts)

def text2workspace_command(datacards_dir: Path, root_dir: Path, model: str, mass: str) -> str:
    return " ".join(
        [
            "text2workspace.py",
            "-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel",
            "--PO verbose",
            quote("--PO"),
            quote("map=.*/(TTto2LDM|TWto2LDM).*:r[1,0,10]"),
            quote(combined_card_path(datacards_dir, model, mass)),
            "-o",
            quote(workspace_path(root_dir, model, mass)),
        ]
    )


def asymptotic_limits_command(root_dir: Path, limit_dir: Path, model: str, mass: str, limit_run: str) -> str:
    run_option = "" if limit_run == "observed" else f"--run {limit_run}"
    return " ".join(
        part
        for part in [
            "combine -M AsymptoticLimits",
            quote(workspace_path(root_dir, model, mass)),
            f"-m {mass}",
            run_option,
            "--cminDefaultMinimizerStrategy 0",
            f"&> {quote(limit_log_path(limit_dir, model, mass))}",
        ]
        if part
    )


def mask_string(signal_regions: Iterable[str]) -> str:
    regions = list(signal_regions)
    if not regions:
        return ""
    masks = ",".join(f"mask_{region}=1" for region in regions)
    return f"--setParameters {quote(masks)}"


def condor_preamble(job_id: int, cmssw_dir: str | None) -> str:
    cmssw_block = ""
    if cmssw_dir:
        cmssw_block = f"""
### CMSSW ###
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd {quote(cmssw_dir)}
cmsenv
cd -
"""
    return f"""#!/bin/bash
set -e

### Job configuration ###
echo "Processing job number {job_id} ..."
CWD=$(pwd -P)
TMPDIR=$(mktemp -d /tmp/ttdm_combine_{job_id}_XXXXXX)
cd "$TMPDIR"
{cmssw_block}
ulimit -s unlimited
"""


def condor_cleanup() -> str:
    return """
### Cleaning ###
cd "$CWD"
rm -rf "$TMPDIR"
echo "shell script has finished"
"""


def ws_job(model: str, mass: str, input_dir: Path, datacards_dir: Path, root_dir: Path, dnn_variable_template: str, job_id: int, cmssw_dir: str | None) -> Job:
    commands = [
        f'echo "Using model {model}, mphi={mass}"',
        f'echo "Using SR1b DNN variable evaluate_dnn_tWDM_{model}_{mass}"',
        f'echo "Using SR2b DNN variable evaluate_dnn_ttDM_{model}_{mass}"',
        f"mkdir -p {quote(datacards_dir)} {quote(root_dir)}",
        combine_cards_command(input_dir, datacards_dir, model, mass, dnn_variable_template),
        text2workspace_command(datacards_dir, root_dir, model, mass),
    ]
    script = condor_preamble(job_id, cmssw_dir) + "\n### Creating workspace ###\n" + "\n".join(commands) + "\n" + condor_cleanup()
    return Job(f"workspace_{model}_mphi_{mass}", script)


def asimov_fit_job(model: str, mass: str, root_dir: Path, diagnostics_dir: Path, job_id: int, cmssw_dir: str | None) -> Job:
    output_dir = diagnostics_dir / model / f"mphi_{mass}" / "AsimovFit"
    ws = workspace_path(root_dir, model, mass)
    script = condor_preamble(job_id, cmssw_dir) + f"""
### Asimov fits: r=0 and r=1 ###
cp {quote(ws)} Workspace.root

### r = 0 background-only Asimov ###
combine -M FitDiagnostics -t -1 --expectSignal=0 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL \\
  --robustFit=1 --forceRecreateNLL \\
  --rMin=-0.5 --rMax=0.5 --out=. -m {mass} -d Workspace.root -n _t0

python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py \\
  -a fitDiagnostics_t0.root -a --histogram=pulls_r0.root >> fitResults_t0

### r = 1 signal+background Asimov ###
combine -M FitDiagnostics -t -1 --expectSignal=1 --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL \\
  --robustFit=1 --forceRecreateNLL \\
  --rMin=0.5 --rMax=1.5 --out=. -m {mass} -d Workspace.root -n _t1

python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py \\
  -a fitDiagnostics_t1.root -a --histogram=pulls_r1.root >> fitResults_t1

mkdir -p {quote(output_dir)}
mv fitDiagnostics_t0.root {quote(output_dir / f"FitDiagnostics_Asimov_{model}_r0_mphi_{mass}.root")}
mv pulls_r0.root {quote(output_dir / f"Pulls_Asimov_{model}_r0_mphi_{mass}.root")}
mv fitResults_t0 {quote(output_dir / f"fitResults_Asimov_{model}_r0_mphi_{mass}.txt")}
mv fitDiagnostics_t1.root {quote(output_dir / f"FitDiagnostics_Asimov_{model}_r1_mphi_{mass}.root")}
mv pulls_r1.root {quote(output_dir / f"Pulls_Asimov_{model}_r1_mphi_{mass}.root")}
mv fitResults_t1 {quote(output_dir / f"fitResults_Asimov_{model}_r1_mphi_{mass}.txt")}
""" + condor_cleanup()
    return Job(f"asimov_{model}_mphi_{mass}", script)


def cr_only_fit_job(model: str, mass: str, root_dir: Path, diagnostics_dir: Path, signal_regions: list[str], job_id: int, cmssw_dir: str | None) -> Job:
    output_dir = diagnostics_dir / model / f"mphi_{mass}" / "CRonlyFit"
    masks = mask_string(signal_regions)
    ws = workspace_path(root_dir, model, mass)
    script = condor_preamble(job_id, cmssw_dir) + f"""
### CR-only observed fit ###
cp {quote(ws)} Workspace.root
combine -M FitDiagnostics --saveNormalizations --saveShapes --saveWithUncertainties --saveNLL \\
  --robustFit=1 --forceRecreateNLL \\
  --out=. -m {mass} -d Workspace.root -n _cronly {masks}
python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py \\
  -a fitDiagnostics_cronly.root -a --histogram=pulls_cronly.root >> fitResults_cronly
mkdir -p {quote(output_dir)}
mv fitDiagnostics_cronly.root {quote(output_dir / f"FitDiagnostics_CR-only_{model}_mphi_{mass}.root")}
mv pulls_cronly.root {quote(output_dir / f"Pulls_CR-only_{model}_mphi_{mass}.root")}
mv fitResults_cronly {quote(output_dir / f"fitResults_CR-only_{model}_mphi_{mass}.txt")}
""" + condor_cleanup()
    return Job(f"cronly_{model}_mphi_{mass}", script)


def gof_job(model: str, mass: str, root_dir: Path, diagnostics_dir: Path, signal_regions: list[str], toy_index: int | None, toys_per_job: int, job_id: int, cmssw_dir: str | None) -> Job:
    output_dir = diagnostics_dir / model / f"mphi_{mass}" / "GoF"
    masks = mask_string(signal_regions)
    if toy_index is None:
        gof_cmd = f"combine -M GoodnessOfFit --algo=saturated --fixedSignalStrength=0.0 -m {mass} -d Workspace.root {masks}"
        name = f"gof_obs_{model}_mphi_{mass}"
    else:
        seed = random.randint(1, 999999)
        gof_cmd = f"combine -M GoodnessOfFit --algo=saturated --fixedSignalStrength=0.0 -t {toys_per_job} --toysFreq -s {seed} -m {mass} -d Workspace.root {masks}"
        name = f"gof_toy{toy_index}_{model}_mphi_{mass}"
    script = condor_preamble(job_id, cmssw_dir) + f"""
### Goodness-of-fit ###
cp {quote(workspace_path(root_dir, model, mass))} Workspace.root
{gof_cmd}
mkdir -p {quote(output_dir)}
mv higgsCombineTest.GoodnessOfFit*.root {quote(output_dir)}/
""" + condor_cleanup()
    return Job(name, script)


def impacts_job(model: str, mass: str, poi: int, root_dir: Path, diagnostics_dir: Path, job_id: int, cmssw_dir: str | None) -> Job:
    output_dir = diagnostics_dir / model / f"mphi_{mass}" / "Impacts"
    script = condor_preamble(job_id, cmssw_dir) + f"""
### Impacts ###
cp {quote(workspace_path(root_dir, model, mass))} Workspace.root
combineTool.py -M Impacts --doInitialFit -t -1 --rMin -10 --rMax 10 --expectSignal={poi} \\
  --robustFit=1 -m {mass} -d Workspace.root
combineTool.py -M Impacts --doFits -t -1 --rMin -10 --rMax 10 --expectSignal={poi} \\
  --robustFit=1 -m {mass} -d Workspace.root
combineTool.py -M Impacts -o impacts_r{poi}.json -m {mass} -d Workspace.root
mkdir -p {quote(output_dir)}
mv impacts_r{poi}.json {quote(output_dir / f"Impacts_Asimov_{model}_mphi_{mass}_r{poi}.json")}
""" + condor_cleanup()
    return Job(f"impacts_{model}_r{poi}_mphi_{mass}", script)


def limits_job(model: str, mass: str, root_dir: Path, diagnostics_dir: Path, limit_run: str, job_id: int, cmssw_dir: str | None) -> Job:
    output_dir = diagnostics_dir / model / f"mphi_{mass}" / "Limits"
    run_option = "" if limit_run == "observed" else f"--run {limit_run}"
    label = limit_run.capitalize()
    script = condor_preamble(job_id, cmssw_dir) + f"""
### {label} limits ###
cp {quote(workspace_path(root_dir, model, mass))} Workspace.root
combine -M AsymptoticLimits -d Workspace.root -m {mass} {run_option}
mkdir -p {quote(output_dir)}
mv higgsCombineTest.AsymptoticLimits.mH{mass}.root {quote(output_dir / f"Limits_{label}_{model}_mphi_{mass}.root")}
""" + condor_cleanup()
    return Job(f"limits_{limit_run}_{model}_mphi_{mass}", script)


def llscan_job(model: str, mass: str, npoints: int, point: int, root_dir: Path, diagnostics_dir: Path, job_id: int, cmssw_dir: str | None) -> Job:
    output_dir = diagnostics_dir / model / f"mphi_{mass}" / "LogLikeScan"
    script = condor_preamble(job_id, cmssw_dir) + f"""
### 1D likelihood scan ###
cp {quote(workspace_path(root_dir, model, mass))} Workspace.root
combine -M MultiDimFit -m {mass} -d Workspace.root --expectSignal 0 -t -1 \\
  --algo grid --robustFit 1 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v1 --alignEdges 1 \\
  --setParameterRanges r=0.0,1.0 --points {npoints} --firstPoint {point} --lastPoint {point} -n .LogLikeScan.POINTS.{point}.{point}
mkdir -p {quote(output_dir)}
mv higgsCombine.LogLikeScan.POINTS.{point}.{point}.MultiDimFit.mH{mass}.root {quote(output_dir)}/
""" + condor_cleanup()
    return Job(f"llscan_{model}_point{point}_mphi_{mass}", script)


def condor_submit_file(run_dir: Path, njobs: int, runtime: int, memory: int, requirements: str | None) -> str:
    script_name = run_dir / "combine_job"
    lines = [
        f"+RequestRuntime       = {runtime}",
        f"RequestMemory         = {memory}",
        "universe              = vanilla",
        f"executable            = {script_name}_$(ProcId).sh",
        f"output                = {script_name}_$(ProcId).out",
        f"error                 = {script_name}_$(ProcId).err",
        f"log                   = {script_name}_$(ProcId).log",
        "transfer_executable   = True",
    ]
    if requirements:
        lines.append(f"requirements          = {requirements}")
    lines.append(f"queue {njobs}")
    return "\n".join(lines) + "\n"


def write_condor_jobs(jobs: list[Job], run_dir: Path, runtime: int, memory: int, requirements: str | None) -> None:
    clean_dir(run_dir)
    for job_id, job in enumerate(jobs):
        script_path = run_dir / f"combine_job_{job_id}.sh"
        script_path.write_text(job.script)
        script_path.chmod(0o755)
    (run_dir / "combine_jobs.submit").write_text(condor_submit_file(run_dir, len(jobs), runtime, memory, requirements))


def build_workspace_jobs(args: argparse.Namespace) -> list[Job]:
    input_dir = input_datacards_path(args.work_dir, args.input_datacards_dir)
    datacards_dir = args.work_dir / args.datacards_dir
    root_dir = args.work_dir / args.root_dir
    jobs: list[Job] = []
    job_id = 0
    for model in args.models:
        for mass in args.masses:
            jobs.append(ws_job(model, mass, input_dir, datacards_dir, root_dir, args.dnn_variable_template, job_id, args.cmssw_dir))
            job_id += 1
    return jobs


def build_fit_jobs(args: argparse.Namespace, tests_to_run: list[str] | None = None) -> list[Job]:
    root_dir = args.work_dir / args.root_dir
    diagnostics_dir = args.work_dir / args.diagnostics_dir
    jobs: list[Job] = []
    tests = set(tests_to_run or args.tests)
    job_id = 0
    for model in args.models:
        for mass in args.masses:
            if "AsimovFit" in tests:
                jobs.append(asimov_fit_job(model, mass, root_dir, diagnostics_dir, job_id, args.cmssw_dir)); job_id += 1
            if "CRonlyFit" in tests:
                jobs.append(cr_only_fit_job(model, mass, root_dir, diagnostics_dir, args.signal_regions, job_id, args.cmssw_dir)); job_id += 1
            if "GoF" in tests:
                jobs.append(gof_job(model, mass, root_dir, diagnostics_dir, args.signal_regions, None, args.toys_per_job, job_id, args.cmssw_dir)); job_id += 1
                for toy_index in range(args.ntoys):
                    jobs.append(gof_job(model, mass, root_dir, diagnostics_dir, args.signal_regions, toy_index, args.toys_per_job, job_id, args.cmssw_dir)); job_id += 1
            if "Impacts" in tests:
                for poi in (0, 1):
                    jobs.append(impacts_job(model, mass, poi, root_dir, diagnostics_dir, job_id, args.cmssw_dir)); job_id += 1
            if "Limits" in tests:
                jobs.append(limits_job(model, mass, root_dir, diagnostics_dir, args.limit_run, job_id, args.cmssw_dir)); job_id += 1
            if "LogLikeScan" in tests:
                for point in range(args.npoints):
                    jobs.append(llscan_job(model, mass, args.npoints, point, root_dir, diagnostics_dir, job_id, args.cmssw_dir)); job_id += 1
    return jobs


def run_local_workspace(args: argparse.Namespace) -> None:
    input_dir = input_datacards_path(args.work_dir, args.input_datacards_dir)
    datacards_dir = args.work_dir / args.datacards_dir
    root_dir = args.work_dir / args.root_dir
    datacards_dir.mkdir(parents=True, exist_ok=True)
    root_dir.mkdir(parents=True, exist_ok=True)
    for model in args.models:
        for mass in args.masses:
            print(f"Using model {model}, mphi={mass}")
            print(f"Using SR1b DNN variable evaluate_dnn_tWDM_{model}_{mass}")
            print(f"Using SR2b DNN variable evaluate_dnn_ttDM_{model}_{mass}")
            run_command(combine_cards_command(input_dir, datacards_dir, model, mass, args.dnn_variable_template), dry_run=args.dry_run)
            run_command(text2workspace_command(datacards_dir, root_dir, model, mass), dry_run=args.dry_run)


def run_local_limits(args: argparse.Namespace) -> None:
    limit_dir = args.work_dir / args.limit_dir
    root_dir = args.work_dir / args.root_dir
    limit_dir.mkdir(parents=True, exist_ok=True)
    for model in args.models:
        labels = []
        values = []
        for mass in args.masses:
            run_command(asymptotic_limits_command(root_dir, limit_dir, model, mass, args.limit_run), dry_run=args.dry_run)
            labels.append(str(limit_log_path(limit_dir, model, mass)))
            values.append(int(mass))
        if args.plot and labels and plotUpperLimits is not None and not args.dry_run:
            plotUpperLimits(labels, values, lumi_label=campaign_lumi_label(args.work_dir, args.plot_lumi_label), model=model, output_dir=str(limit_dir))
    if args.plot and plotUpperLimits is None:
        print("WARNING: brazil_band.plotUpperLimits could not be imported; skipping plot.", file=sys.stderr)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare and run TTDMsimp dileptonic combine workspaces, limits, and fits.")
    parser.add_argument("--task", default="ALL", help="Workflow step to run. Use WS for workspaces; TEST for all diagnostics in --tests; or a single diagnostic task such as AsimovFit, CRonlyFit, GoF, Impacts, Limits, or LogLikeScan. LIMITS keeps the local limit/plot workflow without --condor and creates only limit jobs with --condor.")
    parser.add_argument("--condor", action="store_true", help="write HTCondor scripts instead of running combine commands locally")
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR, help="campaign directory, e.g. Full2024v15")
    parser.add_argument("--input-datacards-dir", default=DEFAULT_INPUT_DATACARDS_DIR, help="input datacard directory relative to --work-dir, or an absolute path (default: datacards)")
    parser.add_argument("--datacards-dir", default="datacards_files", help="relative output directory for combined datacards")
    parser.add_argument("--root-dir", default="root_files", help="relative output directory for workspaces")
    parser.add_argument("--limit-dir", default="limit_files", help="relative output directory for local limit logs")
    parser.add_argument("--diagnostics-dir", default="diagnostics_files", help="relative output directory for fit diagnostics")
    parser.add_argument("--masses", default=",".join(DEFAULT_MASSES), help="comma-separated mediator masses")
    parser.add_argument("--models", default=",".join(DEFAULT_MODELS), help="comma-separated mediator models: ps,s")
    parser.add_argument("--dnn-variable-template", default=DEFAULT_DNN_VARIABLE_TEMPLATE, help="template for the mass-dependent DNN variable directory used by signal-region datacards; may use {mass} (default: evaluate_dnn_mPhi_{mass})")
    parser.add_argument("--signal-regions", default=",".join(DEFAULT_SIGNAL_REGIONS), help="comma-separated combine channel labels to mask for CR-only/GoF jobs")
    parser.add_argument("--tests", default=",".join(DEFAULT_TESTS), help="comma-separated TEST jobs to create")
    parser.add_argument("--ntoys", type=int, default=50, help="number of GoF toy jobs")
    parser.add_argument("--toys-per-job", type=int, default=200, help="number of GoF toys per toy job")
    parser.add_argument("--npoints", type=int, default=50, help="number of likelihood-scan grid points/jobs")
    parser.add_argument("--cmssw-dir", default=None, help="optional CMSSW src directory to cmsenv inside each Condor job")
    parser.add_argument("--condor-runtime", type=int, default=10000, help="Condor RequestRuntime")
    parser.add_argument("--condor-memory", type=int, default=2000, help="Condor RequestMemory")
    parser.add_argument("--condor-requirements", default='(OpSysAndVer =?= "AlmaLinux9")', help="Condor requirements line; pass an empty string to disable")
    parser.add_argument("--limit-run", choices=["expected", "blind", "observed"], default="expected", help="AsymptoticLimits mode for Limits jobs/local LIMITS (default: expected; use blind for fully blinded Asimov, observed only for unblinded observed limits)")
    parser.add_argument("--plot", action="store_true", help="plot expected limits with brazil_band.py after local LIMITS")
    parser.add_argument("--plot-lumi-label", default=None, help="override the luminosity text drawn on the Brazil-band plot; by default read lumi from <work-dir>/configuration.py")
    parser.add_argument("--dry-run", action="store_true", help="print local commands without executing them")
    args = parser.parse_args()
    args.task = normalize_task(args.task)
    args.work_dir = args.work_dir.resolve()
    args.masses = parse_csv(args.masses)
    args.models = parse_csv(args.models)
    unknown_models = sorted(set(args.models) - set(DEFAULT_MODELS))
    if unknown_models:
        raise ValueError(f"Unsupported model(s): {', '.join(unknown_models)}. Choose from: {', '.join(DEFAULT_MODELS)}")
    args.signal_regions = parse_csv(args.signal_regions)
    args.tests = validate_tests(parse_csv(args.tests))
    args.condor_requirements = args.condor_requirements or None
    return args


def main() -> None:
    args = parse_args()

    if args.condor:
        jobs: list[Job] = []
        if args.task == "ALL":
            raise ValueError("For Condor production, run --task WS first and submit those jobs; then run --task TEST after the workspaces exist.")
        if args.task == "WS":
            jobs.extend(build_workspace_jobs(args))
            run_dir = args.work_dir / args.datacards_dir / "HTCondor"
        elif args.task == "TEST" or args.task == "LIMITS" or args.task in DEFAULT_TESTS:
            jobs.extend(build_fit_jobs(args, tests_for_task(args)))
            run_dir = args.work_dir / args.diagnostics_dir / "HTCondor"
        else:
            raise ValueError(f"Unsupported condor task: {args.task}")
        write_condor_jobs(jobs, run_dir, args.condor_runtime, args.condor_memory, args.condor_requirements)
        print(f"\nCombine jobs are ready in {run_dir} with {len(jobs)} jobs in total.")
        print(f"Submit with: condor_submit {run_dir / 'combine_jobs.submit'}")
        return

    if args.task in {"WS", "ALL"}:
        run_local_workspace(args)
    if args.task in {"LIMITS", "ALL"}:
        run_local_limits(args)
    if args.task == "TEST" or args.task in DEFAULT_TESTS:
        print("Fit/diagnostic tasks are intentionally produced as Condor jobs. Re-run with --condor, for example --task Impacts --condor or --task Limits --condor.", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()

