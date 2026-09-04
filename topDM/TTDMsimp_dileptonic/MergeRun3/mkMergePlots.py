#!/usr/bin/env python3

import argparse
import copy
import os

from mkShapesRDF.shapeAnalysis.ConfigLib import ConfigLib
import mkShapesRDF.shapeAnalysis.latinos.LatinosUtils as utils
from mkShapesRDF.shapeAnalysis.latinos.PlotFactory import PlotFactory


def parser():
    p = argparse.ArgumentParser(
        description="Plot an already-merged mkShapesRDF ROOT file using the standard PlotFactory."
    )

    p.add_argument(
        "-i", "--inFile",
        default=None,
        help="Merged ROOT file. If omitted, uses outputFolder/outputFile from configuration.py.",
    )
    p.add_argument(
        "-o", "--outDir",
        default="./plots",
        help="Output plot directory.",
    )
    p.add_argument(
        "-c", "--onlyCut",
        default=None,
        help="Only plot cuts whose name contains this string.",
    )
    p.add_argument(
        "-v", "--onlyVar",
        default=None,
        help="Only plot this variable.",
    )
    p.add_argument(
        "-p", "--onlyPlot",
        default=None,
        help="Plot type(s): c, ratio/cratio, diff/cdifference. Comma separated.",
    )

    # These are the same defaults used by mkShapesRDF's standard mkPlot.py.
    p.add_argument("--scaleToPlot", type=float, default=3.0)
    p.add_argument("--minLogC", type=float, default=0.01)
    p.add_argument("--maxLogC", type=float, default=100.0)
    p.add_argument("--minLogCratio", type=float, default=0.001)
    p.add_argument("--maxLogCratio", type=float, default=10.0)
    p.add_argument("--maxLinearScale", type=float, default=1.45)

    p.add_argument("--linearOnly", action="store_true", default=False)
    p.add_argument("--logOnly", action="store_true", default=False)

    p.add_argument(
        "--fileFormats",
        default="png,root",
        help="Comma-separated: png,pdf,root,C,eps",
    )

    p.add_argument(
        "--showIntegralLegend",
        type=float,
        default=1,
        help="1 = show yields in legend, 0 = hide yields.",
    )

    p.add_argument("--showRelativeRatio", action="store_true", default=False)
    p.add_argument("--showDataMinusBkgOnly", action="store_true", default=False)
    p.add_argument("--removeWeight", action="store_true", default=False)
    p.add_argument("--invertXY", action="store_true", default=False)
    p.add_argument("--skipMissingNuisance", action="store_true", default=False)
    p.add_argument("--removeMCStat", action="store_true", default=False)
    p.add_argument("--plotFancy", action="store_true", default=False)

    p.add_argument(
        "--postFit",
        default="n",
        choices=["n", "p", "s", "b"],
    )
    p.add_argument("--extraLegend", default=None)

    p.add_argument(
        "--plotNormalizedIncludeData",
        action="store_true",
        default=False,
    )
    p.add_argument(
        "--plotNormalizedDistributions",
        action="store_true",
        default=False,
    )
    p.add_argument(
        "--plotNormalizedDistributionsTHstack",
        action="store_true",
        default=False,
    )

    p.add_argument("--NoPreliminary", action="store_true", default=False)
    p.add_argument("--RemoveAllMC", action="store_true", default=False)

    return p


def normalize_plot_names(value):
    if value is None:
        return None

    aliases = {
        "c": "c",
        "ratio": "cratio",
        "cratio": "cratio",
        "diff": "cdifference",
        "difference": "cdifference",
        "cdifference": "cdifference",
    }

    result = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue

        if item not in aliases:
            raise ValueError(
                f"Unknown plot type '{item}'. "
                "Use c, ratio/cratio, or diff/cdifference."
            )

        mapped = aliases[item]
        if mapped not in result:
            result.append(mapped)

    return result


def main():
    args = parser().parse_args()

    # ------------------------------------------------------------------
    # Load the normal MergeRun3 compiled configuration.
    # ------------------------------------------------------------------
    cfg = {}
    ConfigLib.loadLatestPickle(os.path.abspath("configs"), cfg)

    required = [
        "samples",
        "variables",
        "cuts",
        "nuisances",
        "plot",
        "lumi",
        "tag",
        "outputFolder",
        "outputFile",
    ]

    missing = [key for key in required if key not in cfg]
    if missing:
        raise RuntimeError(
            "Missing objects in latest configuration pickle: "
            + ", ".join(missing)
        )

    samples = copy.deepcopy(cfg["samples"])
    variables = copy.deepcopy(cfg["variables"])
    cuts = copy.deepcopy(cfg["cuts"])
    nuisances = copy.deepcopy(cfg["nuisances"])
    plot_cfg = copy.deepcopy(cfg["plot"])

    lumi = cfg["lumi"]
    tag = cfg["tag"]

    if isinstance(cuts, dict) and "cuts" in cuts:
        cuts = cuts["cuts"]

    groupPlot = plot_cfg["groupPlot"]
    legend = plot_cfg["legend"]
    plot = plot_cfg["plot"]

    inputFile = (
        args.inFile
        if args.inFile is not None
        else os.path.join(cfg["outputFolder"], cfg["outputFile"])
    )

    if not os.path.isfile(inputFile):
        raise FileNotFoundError(f"Input ROOT file does not exist: {inputFile}")

    os.makedirs(args.outDir, exist_ok=True)

    # ------------------------------------------------------------------
    # This is the same preprocessing done by standard mkPlot.py.
    # ------------------------------------------------------------------
    subsamplesmap = utils.flatten_samples(samples)
    categoriesmap = utils.flatten_cuts(cuts)

    utils.update_variables_with_categories(variables, categoriesmap)
    utils.update_nuisances_with_subsamples(nuisances, subsamplesmap)
    utils.update_nuisances_with_categories(nuisances, categoriesmap)

    # ------------------------------------------------------------------
    # Optional filtering.
    # ------------------------------------------------------------------
    if args.onlyVar is not None:
        if args.onlyVar not in variables:
            raise KeyError(f"Variable '{args.onlyVar}' not found")

        variables = {
            args.onlyVar: variables[args.onlyVar]
        }

    if args.onlyCut is not None:
        cuts = {
            cutName: cutDef
            for cutName, cutDef in cuts.items()
            if args.onlyCut in cutName
        }

        if not cuts:
            raise KeyError(
                f"No cuts matched '{args.onlyCut}'"
            )

    # ------------------------------------------------------------------
    # Use PlotFactory itself.
    #
    # No uproot/matplotlib reimplementation.
    # No colorPlt.
    # The normal groupPlot['...']['color'] is used exactly as in mkPlot.
    # ------------------------------------------------------------------
    factory = PlotFactory()

    factory._tag = tag
    factory._lumi = lumi

    factory._plotNormalizedDistributions = args.plotNormalizedDistributions
    factory._plotNormalizedIncludeData = args.plotNormalizedIncludeData
    factory._plotNormalizedDistributionsTHstack = (
        args.plotNormalizedDistributionsTHstack
    )

    factory._showIntegralLegend = args.showIntegralLegend

    requested_plots = normalize_plot_names(args.onlyPlot)
    if requested_plots is not None:
        factory._plotsToWrite = requested_plots

    factory._plotLinear = args.linearOnly or not args.logOnly
    factory._plotLog = args.logOnly or not args.linearOnly

    factory._scaleToPlot = args.scaleToPlot
    factory._minLogC = args.minLogC
    factory._maxLogC = args.maxLogC

    factory._minLogCratio = args.minLogCratio
    factory._maxLogCratio = args.maxLogCratio
    factory._maxLinearScale = args.maxLinearScale

    factory._minLogCdifference = args.minLogCratio
    factory._maxLogCdifference = args.maxLogCratio

    factory._showRelativeRatio = args.showRelativeRatio
    factory._showDataMinusBkgOnly = args.showDataMinusBkgOnly

    factory._removeWeight = args.removeWeight
    factory._invertXY = args.invertXY

    factory._fileFormats = [
        item.strip()
        for item in args.fileFormats.split(",")
        if item.strip()
    ]

    factory._postFit = args.postFit
    factory._removeMCStat = args.removeMCStat
    factory._plotFancy = args.plotFancy
    factory._SkipMissingNuisance = args.skipMissingNuisance

    factory._extraLegend = args.extraLegend
    factory._preliminary = not args.NoPreliminary
    factory._removeAllMC = args.RemoveAllMC

    print("")
    print("==============================================================")
    print(" Plotting merged ROOT file with the standard PlotFactory")
    print("==============================================================")
    print(f" Input       : {inputFile}")
    print(f" Output      : {args.outDir}")
    print(f" Cuts        : {list(cuts)}")
    print(f" Variables   : {list(variables)}")
    print(f" Lumi        : {lumi}")
    print(f" groupPlot   : {len(groupPlot)} entries")
    print(" colorPlt    : NOT USED")
    print("==============================================================")
    print("")

    # PlotFactory directly understands the already-merged ROOT hierarchy:
    #
    #   cut/variable/histo_SAMPLE
    #   cut/variable/histo_SAMPLE_NUISANCEUp
    #   cut/variable/histo_SAMPLE_NUISANCEDown
    #
    # so no per-year reconstruction is needed at plotting time.
    factory.makePlot(
        inputFile,
        args.outDir,
        variables,
        cuts,
        samples,
        plot,
        nuisances,
        legend,
        groupPlot,
    )


if __name__ == "__main__":
    main()
