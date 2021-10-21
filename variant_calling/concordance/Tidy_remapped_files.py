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

    with open(data + "/../move_done_files.sh", "w") as f:
        for filename in os.listdir(data):
            a = filename.split("_")
            if "remapped" in a[-1]:
                print(f"mv {a[0]}_{a[1]}_* ../dbsnp_remapped_files", file = f)
                print(f"mv {a[0]}_{a[1]} ../dbsnp_remapped_files", file = f)
