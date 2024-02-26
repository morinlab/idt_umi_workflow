"""Point at dir of fastq files and will return a list of sample names
from fastq files. Simple enough!
"""
import os
import pandas as pd

import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
    '-d',
    '--indir',
    required=True,
    type=str,
    help='Directory containing fastq files')

    parser.add_argument(
    '-d',
    '--indir',
    required=True,
    type=str,
    help='Directory containing fastq files')

    return parser.parse_args()

def main():
    args = get_args()

if __name__ == '__main__':
    main()