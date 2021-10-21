#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 11:04:17 2019

@author: jonahcullen
"""

import argparse
import os

def make_arg_parser():
    parser = argparse.ArgumentParser(
            prog="GeneratePBS.py",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
            "-d", "--data",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Path to dir containing the split dbsnp chrom pos files [required]")
    return parser

if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    with open(data + "/../ncbi_remap.sh", "w") as f:
        for filename in os.listdir(data):
            print(f"perl remap_api.pl --mode asm-asm --from GCF_000002305.2 --dest GCF_002863925.1 --allowdupes off --annotation {data}/{filename} --annot_out {data}/{filename}_remapped --report_out {data}/{filename}_report && mv {data}/{filename} {data}/../dbsnp_remapped_files/ && mv {data}/{filename}_* {data}/../dbsnp_remapped_files/", file = f)
