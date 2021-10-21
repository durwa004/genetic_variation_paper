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
            help="Path to dir containing the split dbsnp  files [required]")
    return parser

if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    success = {}
    for filename in os.listdir(data):
        a = filename.split("_")
        if "remapped" in a[-1]:
            b = a[0] + "_" + a[1]
            success[b] = a[-1]
    with open(data + "/../move_failed_files.sh", "w") as f:
        for filename in os.listdir(data):
            a = filename.split("_")
            b = a[0] + "_" + a[1]
            if b not in success.keys():
                print(f"mv {b} ../dbsnp_chrom_pos_files", file = f)
