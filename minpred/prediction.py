"""Prediction and report generation for the four MINpred phases.

Author: Naveen Duhan
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from minpred import utils
from minpred.inference import predict
from minpred.model_store import load_manifest

pd.options.mode.chained_assignment = None


def _labels(model_key: str) -> list[str]:
    return load_manifest()["models"][model_key]["labels"]


def _prediction_table(features: dict, model_key: str, outfile: str) -> pd.DataFrame:
    labels = _labels(model_key)
    if len(features["SeqID"]) == 0:
        table = pd.DataFrame(columns=["SampleID", "Prediction", *labels])
        table.to_csv(outfile, sep="\t", index=False)
        return table
    probabilities = predict(model_key, features["Samples"])
    if probabilities.shape[1] != len(labels):
        raise RuntimeError(
            f"Model {model_key} returned {probabilities.shape[1]} values, "
            f"but its manifest defines {len(labels)} labels."
        )
    predicted = probabilities.argmax(axis=1)
    table = pd.DataFrame(probabilities * 100.0, columns=labels).round(4)
    table.insert(0, "Prediction", [labels[index] for index in predicted])
    table.insert(0, "SampleID", features["SeqID"])
    table.to_csv(outfile, sep="\t", index=False)
    return table


def pred_ne_phase1(df, outfile):
    table = _prediction_table(df, "phase1", outfile)
    selected = table.loc[table["Prediction"] == "Enzyme", "SampleID"].tolist()
    return selected, table


def pred_ne_phase2(df, outfile):
    table = _prediction_table(df, "phase2", outfile)
    selected = table.loc[
        table["Prediction"] == "Nitrogen Mineralization", "SampleID"
    ].tolist()
    return selected, table


def pred_ne_phase3(df, outfile):
    return _prediction_table(df, "phase3", outfile)


def _phase4(df, outfile, model_key):
    return _prediction_table(df, model_key, outfile)


def amidohydrolasesPred(df, outfile):
    return _phase4(df, outfile, "phase4_amidohydrolases")


def aminopeptidasesPred(df, outfile):
    return _phase4(df, outfile, "phase4_aminopeptidases")


def asparticPred(df, outfile):
    return _phase4(df, outfile, "phase4_aspartic")


def cysteinePred(df, outfile):
    return _phase4(df, outfile, "phase4_cysteine")


def dipeptidasePred(df, outfile):
    return _phase4(df, outfile, "phase4_dipeptidases")


def dipeptidylPred(df, outfile):
    return _phase4(df, outfile, "phase4_dipeptidyl")


def metalloendoPred(df, outfile):
    return _phase4(df, outfile, "phase4_metalloendopeptidases")


def metallopeptidasesPred(df, outfile):
    return _phase4(df, outfile, "phase4_metallopeptidases")


def omegaPred(df, outfile):
    return _phase4(df, outfile, "phase4_omega")


def serinePred(df, outfile):
    return _phase4(df, outfile, "phase4_serine")


PHASE4 = {
    "amidohydrolases": ("Amidohydrolases", amidohydrolasesPred),
    "aminopeptidases": ("Aminopeptidases", aminopeptidasesPred),
    "aspartic": ("Aspartic endopeptidases", asparticPred),
    "cysteine": ("Cysteine endopeptidases", cysteinePred),
    "dipeptidases": ("Dipeptidases", dipeptidasePred),
    "dipeptidyl": ("Dipeptidyl-peptidases", dipeptidylPred),
    "metalloendopeptidases": ("Metalloendopeptidases", metalloendoPred),
    "metallopeptidases": ("Metallopeptidases", metallopeptidasesPred),
    "omega": ("Omega peptidases", omegaPred),
    "serine": ("Serine endopeptidases", serinePred),
}


def pred_ne_phase4(infile, output_dir, fasta_file, output, predict=None):
    if predict != "all" and predict not in PHASE4:
        raise ValueError(f"Unknown Phase IV class: {predict}")
    phase3 = pd.read_csv(infile, sep="\t")
    selected_classes = PHASE4.items() if predict == "all" else [(predict, PHASE4[predict])]
    outputs = []

    def normalized(label):
        return "".join(character for character in str(label).lower() if character.isalnum())

    phase3_labels = phase3["Prediction"].map(normalized)
    for class_key, (class_label, predictor) in selected_classes:
        identifiers = phase3.loc[
            phase3_labels == normalized(class_label), "SampleID"
        ].tolist()
        if not identifiers:
            continue
        input_fasta = f"{output_dir}/{output}_{class_key}_input.fasta"
        utils.fasta_process(fasta_file, input_fasta, identifiers)
        features = utils.preprocess(input_fasta, "CKSAAP", [0, 0], 2400)
        outfile = f"{output_dir}/Phase_4_{class_key}_log.txt"
        table = predictor(features, outfile)
        table.insert(1, "Enzyme_Class", class_label)
        table.to_csv(outfile, sep="\t", index=False)
        outputs.append(table)

    if outputs:
        return pd.concat(outputs, ignore_index=True)
    return pd.DataFrame(columns=["SampleID", "Enzyme_Class", "Prediction"])
