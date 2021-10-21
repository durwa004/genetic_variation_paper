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
            help="Path to dir containing the bcftools stats output files [required]")
    return parser

if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)

    vcf = ["0000.stats", "0001.stats", "0002.stats", "0003.stats"] 

    with open(data + "/overall_concordance/number_of_variants.txt", "w") as info_file:
        print("vcf\tN_samples\tN_records\tN_SNPs\tN_MNPs\tN_indels\tN_others\tN_multiallelic\tN_multiallelic_SNPs\tTS\tTV\tTSTV", file = info_file)
        for i in range(len(vcf)):
            sn = 0
            rn = 0 
            snp = 0
            mnp = 0
            indel = 0
            other = 0
            ma = 0
            ma_snp = 0
            ts = 0
            tv = 0
            for directory in os.listdir(data):
                for file_name in os.listdir(data + "/" + directory):
                    if file_name.endswith(vcf[i]):
                        with open(data + "/" + directory + "/" + vcf[i]) as f:
                            for line in f:
                                line = line.rstrip("\n").split("\t")
                                if "#" in line[0]:
                                    next
                                else:
                                     if line[2] == "number of samples:":
                                         sn += int(line[3])
                                     elif line[2] == "number of records:":
                                         rn += int(line[3])
                                     elif line[2] == "number of SNPs:": 
                                         snp += int(line[3])
                                     elif line[2] == "number of MNPs:":
                                         mnp += int(line[3])
                                     elif line[2] == "number of indels:":
                                         indel += int(line[3])
                                     elif line[2] == "number of others:":
                                         other += int(line[3])
                                     elif line[2] == "number of multiallelic sites:":
                                         ma += int(line[3])
                                     elif line[2] == "number of multiallelic SNP sites:":
                                         ma_snp += int(line[3])
                                     elif line[0] == "TSTV":
                                         ts += int(line[2])
                                         tv += int(line[3])
            tstv = ts/tv
            print(vcf[i], sn, rn, snp, mnp, indel, other, ma, ma_snp, ts, tv, tstv, file = info_file, sep = "\t")
