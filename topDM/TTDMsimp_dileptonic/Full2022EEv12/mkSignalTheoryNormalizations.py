#!/usr/bin/env python3

import argparse
import json
import ROOT
import os

EXPECTED_SCALE_WEIGHTS = 9
MIN_PDF_WEIGHTS = 101

RUN3_ERAS = [
    '2022',
    '2022EE',
    '2023',
    '2023BPix',
    '2024',
]

def getSampleFiles(sample):

    files = []

    for entry in sample['name']:
        if isinstance(entry, tuple):
            if len(entry) < 2:
                raise RuntimeError(f'Invalid mkShapesRDF sample entry: {entry}')
            entryFiles = entry[1]
            if isinstance(entryFiles, str):
                entryFiles = [entryFiles]
            files.extend(entryFiles)
        # Also allow a plain filename for robustness
        elif isinstance(entry, str):
            files.append(entry)
        else:
            raise RuntimeError(f'Unknown sample entry format: {entry}')
    return files


def cleanFileName(name):
    if 'rootd' not in name:
        name = name.replace('###', '')
    return name

def getBranchName(chain, candidates):
    branches = {b.GetName() for b in chain.GetListOfBranches()}
    for name in candidates:
        if name in branches:
            return name
    return None

def getTheoryNormalizations(sampleName, sample, verbose=0):
    chain = ROOT.TChain('Runs')
    sampleFiles = getSampleFiles(sample)
    if verbose >= 1:
        print(f'  Number of ROOT files: {len(sampleFiles)}')
    for treeName in sampleFiles:
        treeName = cleanFileName(treeName)
        if verbose >= 2:
            print(f'  Adding: {treeName}')
        chain.Add(treeName)
    nEntries = chain.GetEntries()
    if nEntries == 0:
        print(f'Error: no Runs entries found for {sampleName}')
        return {'qcdScaleStatus': 0,'pdfStatus': 0,}

    # Support both standard NanoAOD names and possible old
    # trailing-"_" branch convention.
    genBranch = getBranchName(
        chain,
        ['genEventSumw', 'genEventSumw_']
    )

    nScaleBranch = getBranchName(
        chain,
        ['nLHEScaleSumw', 'nLHEScaleSumw_']
    )

    scaleBranch = getBranchName(
        chain,
        ['LHEScaleSumw', 'LHEScaleSumw_']
    )

    nPdfBranch = getBranchName(
        chain,
        ['nLHEPdfSumw', 'nLHEPdfSumw_']
    )

    pdfBranch = getBranchName(
        chain,
        ['LHEPdfSumw', 'LHEPdfSumw_']
    )

    if genBranch is None:
        print(f'Error: {sampleName} has no genEventSumw branch')

        return {
            'qcdScaleStatus': 0,
            'pdfStatus': 0,
        }

    if verbose >= 1:
        print(f'  Runs entries    : {nEntries}')
        print(f'  genEventSumw    : {genBranch}')
        print(f'  nLHEScaleSumw   : {nScaleBranch}')
        print(f'  LHEScaleSumw    : {scaleBranch}')
        print(f'  nLHEPdfSumw     : {nPdfBranch}')
        print(f'  LHEPdfSumw      : {pdfBranch}')

    totalGenWeight = 0.0

    scaleSums = None
    pdfSums = None

    nScaleWeights = None
    nPdfWeights = None

    qcdStatus = 3
    pdfStatus = 3

    for ev in range(nEntries):

        chain.GetEntry(ev)

        genWeight = float(getattr(chain, genBranch))
        totalGenWeight += genWeight

        # ==========================================================
        # QCD muR / muF scales
        # ==========================================================

        if nScaleBranch is None or scaleBranch is None:
            qcdStatus = 0

        elif qcdStatus:

            nScale = int(getattr(chain, nScaleBranch))

            if nScaleWeights is None:
                nScaleWeights = nScale
                scaleSums = [0.0] * nScaleWeights

            elif nScale != nScaleWeights:
                print(
                    f'Error: {sampleName} has inconsistent number of '
                    f'LHEScale weights: {nScaleWeights} vs {nScale}'
                )
                qcdStatus = 0

            if qcdStatus:

                if nScale != EXPECTED_SCALE_WEIGHTS:
                    print(
                        f'Error: {sampleName} has {nScale} LHEScale weights; '
                        f'expected {EXPECTED_SCALE_WEIGHTS}'
                    )
                    qcdStatus = 0

                else:
                    values = getattr(chain, scaleBranch)

                    for i in range(nScale):
                        scaleSums[i] += (
                            genWeight * float(values[i])
                        )

        # ==========================================================
        # PDF
        # ==========================================================

        if nPdfBranch is None or pdfBranch is None:
            pdfStatus = 0

        elif pdfStatus:
            nPdf = int(getattr(chain, nPdfBranch))

            if nPdfWeights is None:
                nPdfWeights = nPdf
                pdfSums = [0.0] * nPdfWeights

            elif nPdf != nPdfWeights:
                print(
                    f'Error: {sampleName} has inconsistent number of '
                    f'LHEPdf weights: {nPdfWeights} vs {nPdf}'
                )
                pdfStatus = 0

            if pdfStatus:
                if nPdf < MIN_PDF_WEIGHTS:
                    print(
                        f'Error: {sampleName} has only {nPdf} PDF weights; '
                        f'expected at least {MIN_PDF_WEIGHTS}'
                    )
                    pdfStatus = 0
                else:
                    values = getattr(chain, pdfBranch)

                    for i in range(nPdf):
                        pdfSums[i] += (
                            genWeight * float(values[i])
                        )

    result = {
        'qcdScaleStatus': qcdStatus,
        'pdfStatus': pdfStatus,
    }

    if totalGenWeight == 0.0:
        print(
            f'Error: total genEventSumw is zero for '
            f'{sampleName}'
        )

        result['qcdScaleStatus'] = 0
        result['pdfStatus'] = 0

        return result

    # ==============================================================
    # Inclusive QCD-scale normalization
    # ==============================================================

    if qcdStatus:

        scaleNorms = [
            value / totalGenWeight
            for value in scaleSums
        ]

        result['nQCDScaleWeights'] = nScaleWeights
        result['qcdScale'] = scaleNorms

        if verbose >= 1:
            print(
                f'  QCD scale weights: '
                f'{nScaleWeights}'
            )

        if verbose >= 2:
            for i, value in enumerate(scaleNorms):
                print(
                    f'    scale[{i:3d}] = '
                    f'{value:.10f}'
                )

    else:
        print(
            f'Error: no usable QCD scale weights '
            f'for {sampleName}'
        )

    # ==============================================================
    # Inclusive PDF normalization
    # ==============================================================

    if pdfStatus:

        pdfNorms = [
            value / totalGenWeight
            for value in pdfSums
        ]

        result['nPdfWeights'] = nPdfWeights
        result['pdf'] = pdfNorms

        if nPdfWeights >= 103:
            result['alphaS'] = [
                pdfNorms[101],
                pdfNorms[102],
            ]

        if verbose >= 1:
            print(
                f'  PDF weights      : '
                f'{nPdfWeights}'
            )

        if verbose >= 2:

            print(
                f'    PDF central    = '
                f'{pdfNorms[0]:.10f}'
            )

            for i in range(1, min(101, nPdfWeights)):
                print(
                    f'    PDF[{i:3d}]      = '
                    f'{pdfNorms[i]:.10f}'
                )

            if nPdfWeights >= 103:
                print(
                    f'    alphaS down    = '
                    f'{pdfNorms[101]:.10f}'
                )
                print(
                    f'    alphaS up      = '
                    f'{pdfNorms[102]:.10f}'
                )

    else:
        print(
            f'Error: no usable PDF weights '
            f'for {sampleName}'
        )

    result['genEventSumw'] = totalGenWeight

    return result


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=(
            'Derive inclusive PDF and QCD-scale normalization '
            'factors for Run-3 signal samples.'
        )
    )

    parser.add_argument(
        '--samplesFile',
        default='./samples.py',
        help='Path to the mkShapesRDF samples.py file. Default: ./samples.py'
    )

    parser.add_argument(
        '--year',
        choices=RUN3_ERAS,
        default='2022EE',
        help='Run-3 era. Default: 2022EE'
    )

    parser.add_argument(
        '--sigset',
        default='Signals',
        help='Tag used for the output file. Default: Signals'
    )

    parser.add_argument(
        '--verbose',
        type=int,
        default=2,
        help='Verbose level. Default: 2'
    )

    args = parser.parse_args()

    # ==============================================================
    # Load mkShapesRDF samples.py
    # ==============================================================

    samples = {}

    sampleNamespace = {
        'samples': samples,
        'os': os,
    }

    with open(args.samplesFile) as f:
        exec(
            compile(
                f.read(),
                args.samplesFile,
                'exec'
            ),
            sampleNamespace
        )

    samples = sampleNamespace['samples']

    if args.verbose >= 1:
        print(
            f'\nLoaded {len(samples)} samples '
            f'from {args.samplesFile}'
        )

    theoryNormalizations = {}

    # ==============================================================
    # Process signal samples
    # ==============================================================

    for sam_k, sam_v in samples.items():
        isSignal = sam_v.get(
            'isSignal',
            'DM' in sam_k
        )

        if not isSignal:
            continue

        if args.verbose >= 1:
            print('\n' + '=' * 100)
            print(
                f'Deriving theory normalizations for '
                f'{sam_k} ({args.year})'
            )
            print('=' * 100)

        theoryNormalizations[sam_k] = (
            getTheoryNormalizations(
                sam_k,
                sam_v,
                args.verbose
            )
        )

    # ==============================================================
    # Write result
    # ==============================================================

    outputFile = (
        f'theoryNormalizations_'
        f'{args.year}_'
        f'{args.sigset}.py'
    )

    with open(outputFile, 'w') as f:

        f.write(
            'theoryNormalizations = {}\n\n'
        )

        for sample, values in theoryNormalizations.items():

            f.write(
                "theoryNormalizations['{}'] = {}\n\n".format(
                    sample,
                    json.dumps(
                        values,
                        sort_keys=True
                    )
                )
            )

    print(
        f'\nWritten: {outputFile}'
    )
