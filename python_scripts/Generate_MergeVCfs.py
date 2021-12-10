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
    parser.add_argument(
            "-c", "--chrs",
            default=argparse.SUPPRESS,
            metavar="",
            required=True,
            help="Path to directory containing group gvcfs [required]")
    return parser


if __name__ == '__main__':
     
    parser = make_arg_parser()
    args = parser.parse_args()

    data = os.path.abspath(args.data)
    chr_ids = os.path.abspath(args.chrs)
    
    java = 'java -Xmx4g -jar /home/mccuem/shared/.local/conda/envs/HorseGenomeProject/bin/picard.jar MergeVcfs'
 
    header = (
              "#!/bin/bash -l\n"  
              "#SBATCH --time=12:00:00\n"
              "#SBATCH --ntasks=1\n"
              "#SBATCH --mem=2g\n"
              "#SBATCH --mail-type=ALL\n"
              "#SBATCH --mail-user=perso208@umn.edu\n"
              "#SBATCH -o mergevcfs.out\n"
              "#SBATCH -e mergevcfs.err\n"
              "#SBATCH -p small,amdsmall\n"
             )
    
    with open(f"{chr_ids}/vcfs_for_merging.list", "w") as f:
        for directory in os.listdir(chr_ids):
            if os.path.isdir(directory) == True:
                for file_name in os.listdir(f"{chr_ids}/{directory}"):
                    file_name1 = file_name + ".vcf.gz"
                    print(f"{chr_ids}/{directory}/{file_name}/{file_name1}", file=f)
    
    pbs = os.path.join(os.getcwd(), "MergeVCFs.slurm")
    
    with open(pbs, "w") as f:
        print(header, file=f)
        print(f"cd {chr_ids}\n", file=f)
        print(f"{java} I=vcfs_for_merging.list O=merged_canine_vcf.vcf.gz", file=f)
        
