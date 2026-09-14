#!/usr/bin/python
"""
Title: Utility for MINpred
Author: Naveen Duhan
Lab: KAABiL(Kaundal Artificial Intelligence & Advanced Bioinformatics Lab)
Version:  1.0
"""
import numpy as np
import re, argparse
from Bio import SeqIO

from minpred.__version__ import __version__

import os, shutil
# Amino Acid one letter code
AA= ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
myAAorder = {
		'alphabetically': 'ACDEFGHIKLMNPQRSTVWY',
		'polarity': 'DENKRQHSGTAPYVMCWIFL',
		'sideChainVolume': 'GASDPCTNEVHQILMKRFYW',
	}
AADict={'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}
AAidx=([[ 2.000e-02, -4.200e-01, -7.700e-01, -1.040e+00,  7.700e-01,
        -1.100e+00, -1.140e+00, -8.000e-01,  2.600e-01,  1.810e+00,
         1.140e+00, -4.100e-01,  1.000e+00,  1.350e+00, -9.000e-02,
        -9.700e-01, -7.700e-01,  1.710e+00,  1.110e+00,  1.130e+00],
       [ 3.570e-01,  5.290e-01,  4.630e-01,  5.110e-01,  3.460e-01,
         4.930e-01,  4.970e-01,  5.440e-01,  3.230e-01,  4.620e-01,
         3.650e-01,  4.660e-01,  2.950e-01,  3.140e-01,  5.090e-01,
         5.070e-01,  4.440e-01,  3.050e-01,  4.200e-01,  3.860e-01],
       [ 4.600e-02,  2.910e-01,  1.340e-01,  1.050e-01,  1.280e-01,
         1.800e-01,  1.510e-01,  0.000e+00,  2.300e-01,  1.860e-01,
         1.860e-01,  2.190e-01,  2.210e-01,  2.900e-01,  1.310e-01,
         6.200e-02,  1.080e-01,  4.090e-01,  2.980e-01,  1.400e-01],
       [-3.680e-01, -1.030e+00,  0.000e+00,  2.060e+00,  4.530e+00,
         7.310e-01,  1.770e+00, -5.250e-01,  0.000e+00,  7.910e-01,
         1.070e+00,  0.000e+00,  6.560e-01,  1.060e+00, -2.240e+00,
        -5.240e-01,  0.000e+00,  1.600e+00,  4.910e+00,  4.010e-01],
       [ 1.150e+02,  2.250e+02,  1.600e+02,  1.500e+02,  1.350e+02,
         1.800e+02,  1.900e+02,  7.500e+01,  1.950e+02,  1.750e+02,
         1.700e+02,  2.000e+02,  1.850e+02,  2.100e+02,  1.450e+02,
         1.150e+02,  1.400e+02,  2.550e+02,  2.300e+02,  1.550e+02],
       [ 5.260e+01,  1.091e+02,  7.570e+01,  6.840e+01,  6.830e+01,
         8.970e+01,  8.470e+01,  3.630e+01,  9.190e+01,  1.020e+02,
         1.020e+02,  1.051e+02,  9.770e+01,  1.139e+02,  7.360e+01,
         5.490e+01,  7.120e+01,  1.354e+02,  1.162e+02,  8.510e+01],
       [ 5.200e-01,  6.800e-01,  7.600e-01,  7.600e-01,  6.200e-01,
         6.800e-01,  6.800e-01,  0.000e+00,  7.000e-01,  1.020e+00,
         9.800e-01,  6.800e-01,  7.800e-01,  7.000e-01,  3.600e-01,
         5.300e-01,  5.000e-01,  7.000e-01,  7.000e-01,  7.600e-01],
       [ 1.000e+02,  6.500e+01,  1.340e+02,  1.060e+02,  2.000e+01,
         9.300e+01,  1.020e+02,  4.900e+01,  6.600e+01,  9.600e+01,
         4.000e+01,  5.600e+01,  9.400e+01,  4.100e+01,  5.600e+01,
         1.200e+02,  9.700e+01,  1.800e+01,  4.100e+01,  7.400e+01]])
def DPC(seq):
    N = len(seq)
    dpc = []
    for i in AA:
        for j in AA:
            dp = i + j
            dp = round(float(seq.count(dp)) / (N - 1) * 100, 2)
            dpc.append(dp)
    return dpc

def CKSAAP(seq, gap=5, order='alphabetically'):
        if gap < 0:
            print('Error: the gap should be equal or greater than zero' + '\n\n')
            return 0

        if len(seq) < gap+2:
            print('Error: all the sequence length should be larger than the (gap value) + 2 = ' + str(gap+2) + '\n\n')
            return 0
        AA = myAAorder[order]
        encodings = []
        aaPairs = []
        for aa1 in AA:
            for aa2 in AA:
                aaPairs.append(aa1 + aa2)
        code =[]
        for g in range(gap+1):
            myDict = {}
            for pair in aaPairs:
                myDict[pair] = 0
            sum = 0
            for index1 in range(len(seq)):
                index2 = index1 + g + 1
                if index1 < len(seq) and index2 < len(seq) and seq[index1] in AA and seq[index2] in AA:
                    myDict[seq[index1] + seq[index2]] = myDict[seq[index1] + seq[index2]] + 1
                    sum = sum + 1
            for pair in aaPairs:
                encodings.append(myDict[pair] / sum)

        return encodings

def NMBroto(seq, nlag=30):
    pstd = np.std(AAidx, axis=1)
    pmean = np.average(AAidx, axis=1)

    for i in range(len(AAidx)):
        for j in range(len(AAidx[i])):
            AAidx[i][j] = (AAidx[i][j] - pmean[i]) / pstd[i]
    code = []
    N = len(seq)
    for prop in range(8):
        for n in range(1, nlag + 1):
            if len(seq) > nlag:
                # if key is '-', then the value is 0
                rn = sum([AAidx[prop][AADict.get(seq[j], 0)] * AAidx[prop][AADict.get(seq[j + n], 0)] for j in
                          range(len(seq) - n)]) / (N - n)
            else:
                rn = '0'
            code.append(rn)
    return code

def hybrid(seq,f1,f2):
    a = f1(seq)
    a.extend(f2(seq))
    return a

def preprocess(file,feat, classLabel,size):
    _seq=SeqIO.parse(file, 'fasta')
    tmp = []

    totalSequences = 0
    _classLabels = []
    seqID = []
    for s in _seq:
        res = ''
        totalSequences += 1
        seqS = str(s.seq)
        seqStr=re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(seqS).upper())
        res+=seqStr

        if feat=="DPC":
            tmp.append(DPC(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == "hybrid":
            tmp = tmp
            tmp.append(hybrid(res, DPC, NMBroto))
            _classLabels.append(classLabel)
            seqID.append(s.id)
        elif feat == 'CKSAAP':
            tmp = tmp
            tmp.append(CKSAAP(res))
            _classLabels.append(classLabel)
            seqID.append(s.id)

    tmp=np.array(tmp)

    aa=np.zeros((len(tmp),int(size)))
    for i,s in enumerate(tmp):
        aa[i]=list(s)
    return {'Samples':aa, 'Labels':_classLabels,'SeqID': seqID}

def fasta_process(fasta_file,result_file, wanted):
    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    with open(result_file, "w") as f:
        for seq in fasta_sequences:
            name = seq.id
            if name in wanted and len(str(name)) > 0:
                wanted.remove(name)  # Output only the first appearance for a name
                SeqIO.write([seq], f, "fasta")



def argument_parser(version=__version__):
    parser = argparse.ArgumentParser(description="""MINpred: Nitrogen mineralization enzyme class prediction tool""",
    epilog="""Written by Naveen Duhan (naveen.duhan@outlook.com),
    Kaundal Bioinformatics Lab, Utah State University,
    Released under the MIT License""",
    formatter_class=argparse.RawTextHelpFormatter )
    parser.add_argument('--version', action='version', version=f'MINpred {version}')
    parser.add_argument('-od', '--output_dir', default='MINpred_results', help="Output directory")
    parser.add_argument('-o', '--output_file', default="MINpred_results.txt", help="Output file")
    parser.add_argument('-i', '--fasta_file', required=True, help="Input FASTA file")
    parser.add_argument(
        '--sequence-type',
        choices=['protein', 'nucleotide'],
        default='protein',
        help=("Input sequence type. Nucleotide input is translated with "
              "TransDecoder before MINpred prediction (default: protein)"),
    )
    parser.add_argument(
        '-l', '--level', choices=['Phase1', 'Phase2', 'Phase3', 'Phase4'],
        default='Phase4', help="Prediction phase (default: Phase4)"
    )
    parser.add_argument(
        '-n', '--ecnumber',
        choices=['all', 'amidohydrolases', 'aminopeptidases', 'aspartic', 'cysteine',
                 'dipeptidases', 'dipeptidyl', 'metalloendopeptidases',
                 'metallopeptidases', 'omega', 'serine'],
        default='all',
        help="Phase IV enzyme class (default: all; route from Phase III)",
    )
    return parser
