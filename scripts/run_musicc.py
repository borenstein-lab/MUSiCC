#!/usr/bin/env python

import argparse

# for testing the new module, remove later!!!!!!
import sys
sys.path.append('/net/gs/vol1/home/ohadm/MUSiCC/PyCode/MUSiCC')

from musicc.core import correct_and_normalize

if __name__ == "__main__":
    # get options from user
    parser = argparse.ArgumentParser(description='MUSiCC: Metagenomic Universal Single-Copy Correction')
    parser.add_argument('input_file', help='Input abundance file to correct')
    parser.add_argument('-o', '--out', dest='output_file', help='Output destination for corrected abundance (default: MUSiCC.tab)', default='MUSiCC.tab')
    parser.add_argument('-if', '--input_format', dest='input_format', choices=['tab', 'csv', 'biom'], help='Option indicating the format of the input file (default: tab)', default='tab')
    parser.add_argument('-of', '--output_format', dest='output_format', choices=['tab', 'csv', 'biom'], help='Option indicating the format of the output file (default: tab)', default='tab')
    parser.add_argument('-n', '--normalize', dest='musicc_inter', help='Apply MUSiCC normalization (default: false)', action='store_true')
    parser.add_argument('-c', '--correct', dest='musicc_intra', choices=['use_generic', 'learn_model'], help='Correct abundance per-sample using MUSiCC (default: false)', default='None')
    parser.add_argument('-perf', '--performance', dest='compute_scores', help='Calculate model performance on various gene sets (may add to running time) (default: false)', action='store_true')
    parser.add_argument('-v', '--verbose', dest='verbose', help='Increase verbosity of module (default: false)', action='store_true')

    given_args = parser.parse_args()

    # run normalization and correction
    correct_and_normalize(vars(given_args))


