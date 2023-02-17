#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
from tensorflow import keras as k


def parse_input():
    """
    function for parsing input parameters
    :return: dictionary of parameters
    """
    parser = argparse.ArgumentParser(description='Method for prediction of potential smallRNA:target site '
                                                 'binding')
    parser.add_argument('--input', default="example.tsv", metavar='<input_tsv_filename>')
    parser.add_argument('--output', default="example_scores", metavar='<output_filename_prefix>')
    parser.add_argument('--model', required=True, metavar='<model_name>', help='choose one of Models/model_miRNA.h5, '
                                                                               'Models/model_tRNA.h5, Models/model_yRNA.h5')
    args = parser.parse_args()
    return vars(args)


def binding_encoding(df, tensor_dim=(50, 20, 1)):
    """
    fun encodes miRNAs and mRNAs in df into binding matrices
    :param df: dataframe containing 'gene' and 'noncodingRNA' columns
    :param tensor_dim: output shape of the matrix
    :return: numpy array of predictions
    """
    # alphabet for watson-crick interactions.
    alphabet = {"AT": 1., "TA": 1., "GC": 1., "CG": 1., "AU": 1., "UA": 1.}
    N = df.shape[0]  # number of samples in df
    shape_matrix_2d = (N, *tensor_dim)  # 2d matrix shape
    ohe_matrix_2d = np.zeros(shape_matrix_2d, dtype="float32")
    # make matrix with watson-crick interactions.
    for index, row in df.iterrows():
        for bind_index, bind_nt in enumerate(row.gene.upper()):
            for ncrna_index, ncrna_nt in enumerate(row.noncodingRNA.upper()):
                if ncrna_index >= tensor_dim[1]:
                    break
                base_pairs = bind_nt + ncrna_nt
                ohe_matrix_2d[index, bind_index, ncrna_index, 0] = alphabet.get(base_pairs, 0)
    return ohe_matrix_2d


def write_score(output_file, df, scores):
    """
    fun writes information about sequence and its score to the output_file
    :param output_file
    :param df: dataframe with smallRNA:target pairs
    :param scores: numpy array, predicted scores
    """
    scores = scores.flatten()
    df["score"] = pd.Series(scores, index=df.index)
    df.to_csv(output_file + '.tsv', sep='\t', index=False)


def predict_probs(df, model, output):
    """
    fun predicts the probability of miRNA:target site binding in df file
    :param df: input dataframe with sequences containing 'gene' and 'noncodingRNA' columns
    :param model: Keras model used for predicting
    :param output: output file to write probabilities to
    """
    gene_length = 50

    orig_len = len(df)
    df = df[df["gene"].str.len() == gene_length]
    processed_len = len(df)

    if orig_len != processed_len:
        print("Skipping " + str(orig_len - processed_len) + " pairs due to inappropriate target length.")

    ohe = binding_encoding(df)
    prob = model.predict(ohe)
    write_score(output, df, prob)


def main():
    arguments = parse_input()

    output = arguments["output"]

    try:
        model = k.models.load_model(arguments["model"])
    except (IOError, ImportError):
        print()
        print("Can't load the model ", arguments["model"])
        return

    print("===========================================")

    try:
        input_df = pd.read_csv(arguments["input"], names=['noncodingRNA', 'gene'], sep='\t')
    except IOError as e:
        print()
        print("Can't load file ", arguments["input"])
        print(e)
        return

    predict_probs(input_df, model, output)


main()
