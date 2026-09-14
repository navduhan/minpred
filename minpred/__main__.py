#!/usr/bin/python
"""
Title: Main Script for Nitrogen metabolism enzyme prediction
Author: Naveen Duhan
Lab: KAABiL(Kaundal Artificial Intelligence & Advanced Bioinformatics Lab)
Version: 1.0.0
"""

import os, sys, time
import pandas as pd
from minpred import utils
from minpred import prediction
from minpred import __version__
from minpred.model_store import download_models, model_status
from minpred.transdecoder import translate_nucleotide_fasta


def DNN(fasta_file, output_dir, level, pred, output):
    result = None  # Initialize result variable before the if blocks

    if level == "Phase1":
        enz1 = utils.preprocess(fasta_file, 'hybrid', [0, 0], 640)
        phase1_out = '%s/Phase_1_dnn_log.txt' % (output_dir)
        phase1, result = prediction.pred_ne_phase1(enz1, phase1_out)

    if level == "Phase2":
        enz1 = utils.preprocess(fasta_file, 'hybrid', [0, 0], 640)
        phase1_out = '%s/Phase_1_dnn_log.txt' % (output_dir)
        phase1, result = prediction.pred_ne_phase1(enz1, phase1_out)

        input2 = "%s/phase_2_input.fasta" % (output_dir)
        utils.fasta_process(fasta_file, input2, phase1)
        enz2 = utils.preprocess(input2, 'hybrid', [0, 0], 640)
        phase2_out = '%s/Phase_2_dnn_log.txt' % (output_dir)
        phase2, result = prediction.pred_ne_phase2(enz2, phase2_out)

    if level == "Phase3":
        enz1 = utils.preprocess(fasta_file, 'hybrid', [0, 0], 640)
        phase1_out = '%s/Phase_1_dnn_log.txt' % (output_dir)
        phase1, result = prediction.pred_ne_phase1(enz1, phase1_out)

        input2 = "%s/phase_2_input.fasta" % (output_dir)
        utils.fasta_process(fasta_file, input2, phase1)
        enz2 = utils.preprocess(input2, 'hybrid', [0, 0], 640)
        phase2_out = '%s/Phase_2_dnn_log.txt' % (output_dir)
        phase2, result = prediction.pred_ne_phase2(enz2, phase2_out)

        input3 = "%s/phase_3_input.fasta" % (output_dir)
        utils.fasta_process(fasta_file, input3, phase2)
        enz3 = utils.preprocess(input3, 'hybrid', [0, 0], 640)
        phase3_out = '%s/Phase_3_dnn_log.txt' % (output_dir)
        result = prediction.pred_ne_phase3(enz3, phase3_out)

    if level == "Phase4":

        enz1 = utils.preprocess(fasta_file, 'hybrid', [0, 0], 640)
        phase1_out = '%s/Phase_1_dnn_log.txt' % (output_dir)
        phase1, result = prediction.pred_ne_phase1(enz1, phase1_out)

        input2 = "%s/phase_2_input.fasta" % (output_dir)
        utils.fasta_process(fasta_file, input2, phase1)
        enz2 = utils.preprocess(input2, 'hybrid', [0, 0], 640)
        phase2_out = '%s/Phase_2_dnn_log.txt' % (output_dir)
        phase2, result = prediction.pred_ne_phase2(enz2, phase2_out)

        input3 = "%s/phase_3_input.fasta" % (output_dir)
        utils.fasta_process(fasta_file, input3, phase2)
        enz3 = utils.preprocess(input3, 'hybrid', [0, 0], 640)
        phase3_out = '%s/Phase_3_dnn_log.txt' % (output_dir)
        result = prediction.pred_ne_phase3(enz3, phase3_out)

        result = prediction.pred_ne_phase4(
            phase3_out, output_dir, fasta_file, output, predict=pred
        )

    return result


def main():
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

    if len(sys.argv) >= 2 and sys.argv[1] == "download-models":
        selection = sys.argv[2] if len(sys.argv) >= 3 else "all"
        for path in download_models(selection):
            print(f"Ready: {path}")
        return
    if len(sys.argv) >= 2 and sys.argv[1] == "models-status":
        for key, path, ready in model_status():
            print(f"{key}: {'ready' if ready else 'missing'} ({path})")
        return

    start=time.time()
    parser =utils.argument_parser(version=__version__)
    options = parser.parse_args()
    fasta_file = options.fasta_file
    output_dir = options.output_dir
    output= options.output_file
    level = options.level
    predi = options.ecnumber

    os.makedirs(output_dir, exist_ok=True)

    print(f"\n Created output directory {output_dir}\n")


    if options.sequence_type == "nucleotide":
        fasta_file = translate_nucleotide_fasta(fasta_file, output_dir)
        print(f" TransDecoder peptide FASTA: {fasta_file}\n")

    res =DNN(fasta_file,output_dir,level,predi,output)
    out ='%s/%s'% (output_dir,output)

    res.to_csv(out, sep="\t",index=False)


if __name__ == "__main__":
    main()
