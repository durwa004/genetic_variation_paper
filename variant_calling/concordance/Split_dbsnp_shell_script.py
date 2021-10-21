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
            help="Path to dir containing the split dbsnp shell scripts [required]")
    return parser

if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    header = (
              "#!/bin/bash -l\n"
              "#PBS -l nodes=1:ppn=1,walltime=96:00:00,mem=1g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              "#PBS -o $PBS_JOBID.dbsnp_remap.out\n"
              "#PBS -e $PBS_JOBID.dbsnp_remap.err\n"
              "#PBS -N dbsnp_remap.pbs\n"
              "#PBS -q small\n"
              "source /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/activate ensembl-vep\n")
    count = 0
    for filename in os.listdir(data):
        if "ncbi_remap_x" in filename:
            count +=1
            pbs = os.path.join(os.getcwd(), f"dbsnp_remap_{count}.pbs")
            with open(pbs, "w") as f:
                print(header, file=f) 
                print(f"cd {data}/../\n", file=f) 
                print(f"sh EVD_dbsnp/{filename}", file = f)
