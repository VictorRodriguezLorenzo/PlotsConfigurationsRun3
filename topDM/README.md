# TTDMsimp dileptonic Run 3 analysis

This directory contains the `mkShapesRDF` configurations for the dileptonic and semileptonic
t(t)+ DM simplified-model analysis.  The Run 3 data-taking periods are configured
separately and can then be merged into one Run 3 combination.

## Available campaigns

| Data-taking period | Configuration directory | Configuration tag |
| --- | --- | --- |
| 2022, pre-EE | `Full2022v12` | `ttDM_dilep_2022` |
| 2022, post-EE | `Full2022EEv12` | `ttDM_dilep_2022EE` |
| 2023, pre-BPix | `Full2023v12` | `ttDM_dilep_2023` |
| 2023, post-BPix | `Full2023BPixv12` | `ttDM_dilep_2023BPix` |
| 2024 | `Full2024v15` | `ttDM_dilep_2024` |

## Prerequisites

Install this repository next to
[`mkShapesRDF`](https://github.com/latinos/mkShapesRDF). To install it, use:

```bash
git clone https://github.com/latinos/mkShapesRDF.git
```
and then initalize the repository with:

```bash
cd mkShapesRDF
source install.sh
source start.sh
cd topDM/TTDMsimp_dileptonic
```

## Run one campaign

Use the same sequence for every complete campaign.  

```bash
cd topDM/TTDMsimp_dileptonic   # omit if already here
cd "${CAMPAIGN}"
```

### 1. Compile the configuration

Compile after every change to a configuration file:

```bash
mkShapesRDF -c 1
```

This creates the compiled configuration under `configs/`, which is also used
by the plotting and datacard scripts.

### 2. Produce and merge histograms

First, test in local:

```bash
mkShapesRDF -o 0 -f . -b 0 -l 1
```

If it works, submit the shape-production jobs to HTCondor:

```bash
mkShapesRDF -o 0 -f . -b 1
```

Check their status, resubmit failed jobs if necessary, and merge successful
job outputs:

```bash
mkShapesRDF -o 1 -f .
mkShapesRDF -o 1 -f . -r 1
mkShapesRDF -o 2 -f .
```

The merged file is written to the `outputFolder` declared in
`configuration.py`, with the name `mkShapes__<tag>.root`.

### 3. Make plots (optional)

After merging, produce the standard data/MC ratio plots with:

```bash
mkPlot --onlyPlot cratio --showIntegralLegend 1 --fileFormats png
```

The destination is the `plotPath` configured for that campaign.  The usual
`mkPlot` filters, such as `--onlyCut` and `--onlyVariable`, may be used to
restrict the output.

### 4. Make datacards (optional)

The analysis-specific datacard script reads the merged ROOT file and the most
recent compiled configuration:

```bash
python mkDatacards_topDM.py --outputDirDatacard ./datacards
```

To restrict production, pass comma-separated cut and variable names:

```bash
python mkDatacards_topDM.py \
  --outputDirDatacard ./datacards \
  --onlyCuts <cut1,cut2> \
  --onlyVariables <variable1,variable2>
```

Do this for every campaign if you want to obtain the full combination.

### NOTE: For every campaign, make sure to run previously:

- To obtain theory normalizations:

```bash
python3 mkTheoryNormalizations.py --samplesFile samples.py --year {year} --verbose 2
```

- Make sure all macros in extended/ are properly compiled. If not, run for all files:

```bash
root {file}.cc+
```

- Perform the required DNN trainings for both ttDM and tWDM before running the whole
code. For this, run first:
```bash
python3 prepare_snapshots_submit.py
```
This will create the ROOT files required for the training. Then do:
```bash
python3 train_from_snapshots_ttDM.py  /  python3 train_from_snapshots_tWDM.py
```

## Build the Run 3 combination

The combination is a histogram-level merge performed by `MergeRun3`.  It
sums nominal histograms from the campaign ROOT files and constructs the
combined nuisance variations.

### 1. Select the inputs

Edit `MergeRun3/merge_configuration.py` before running the merge.  Each entry
in `foldersToMerge` provides:

- `folder`: the campaign configuration directory;
- `tag`: the campaign's ROOT-file tag; and
- `outputFolder`: the directory containing its merged ROOT file.

Keep only campaigns whose merged input files are available (you have run before).

Also verify `lumi` in `MergeRun3/configuration.py` matches the periods selected
in `foldersToMerge`.

### 2. Compile and merge

From the TTDMsimp dileptonic directory:

```bash
cd MergeRun3
mkShapesRDF -c 1
python mkMergeYears.py --nJobs 10
```

`--nJobs` controls the number of local cut/variable merge processes.  The
combined ROOT file is written to the `outputFolder` in
`MergeRun3/configuration.py` as `mkShapes__ttDM_dilep_Run3.root`.

For a small validation job, merge only one cut and one variable:

```bash
python mkMergeYears.py --onlyCut <cut_name> --onlyVar <variable_name>
```

To prepare one HTCondor job per cut/variable instead, run:

```bash
python mkMergeYears.py --doSubmit
condor_submit condor_submit.jdl
```

After those jobs finish, merge their temporary files manually:

```bash
COMBINED_OUTPUT=/eos/user/<initial>/<user>/www/Run3-ttDM/rootFiles/ttDM_dilep_Run3
hadd -fk \
  "${COMBINED_OUTPUT}/mkShapes__ttDM_dilep_Run3.root" \
  "${COMBINED_OUTPUT}"/mkShapes__ttDM_dilep_Run3__ALL__*.root
```

Set `COMBINED_OUTPUT` to the `outputFolder` from
`MergeRun3/configuration.py`.

### 3. Plot the combination

`mkMergePlots.py` loads the compiled `MergeRun3` configuration and uses its
combined ROOT file by default:

```bash
python mkMergePlots.py \
  --outDir ./Plots \
  --onlyPlot cratio \
  --fileFormats png
```

An explicit merged file can be supplied with `--inFile`.  Use `--onlyCut` or
`--onlyVar` for quick checks.

### 4. Make combined datacards

```bash
python mkDatacards_topDM.py --outputDirDatacard ./datacards
```

As for an individual campaign, `--onlyCuts` and `--onlyVariables` accept
comma-separated lists.  The resulting cards can then be passed to the CMS
Combine workflow appropriate for the statistical result being produced.


