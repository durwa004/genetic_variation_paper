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
            help="Path to dir containing the ibio output files [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    tmp = []
    for file_name in os.listdir(data + f"/../joint_bcftools/"):
        if file_name.endswith(".genotyped.vcf.gz"):
            a = file_name.split(".genotyped.vcf.gz")
            tmp.append(a[0])

    for i in tmp:
        pbs = os.path.join(os.getcwd(), f"concordance_bcftools_gatk_" + i + ".pbs")
        header = (
              "#!/bin/bash -l\n"
              "#PBS -l nodes=1:ppn=8,walltime=12:00:00,mem=4g\n"
              "#PBS -m abe\n"
              "#PBS -M durwa004@umn.edu\n"
              "#PBS -o $PBS_JOBID.concordance_bcftools_gatk_" + i + ".out\n"
              "#PBS -e $PBS_JOBID.concordance_bcftools_gatk_" + i + ".err\n"
              "#PBS -N concordance_bcftools_gatk_" + i + ".pbs\n"
              "#PBS -q batch\n"
                 )
        with open(pbs, "w") as f:
            print(header, file=f)
            print(f"cd {data}\n", file=f)
            print("module load bcftools", file=f)
            print("mkdir concordance_" + i + "/", file=f)
            print("bcftools isec -p concordance_" + i + f" {data}/../joint_bcftools_without_Prze/" + i + f".genotyped.vcf.gz {data}/../joint_gatk_without_Prze/" + i + ".genotyped.vcf.gz", file=f)
            print("bcftools stats concordance_" + i + "/0000.vcf > concordance_" + i + "/0000.stats", file=f)
            print("bcftools stats concordance_" + i + "/0001.vcf > concordance_" + i + "/0001.stats", file=f)
            print("bcftools stats concordance_" + i + "/0002.vcf > concordance_" + i + "/0002.stats", file=f)
            print("bcftools stats concordance_" + i + "/0003.vcf > concordance_" + i + "/0003.stats", file=f)
