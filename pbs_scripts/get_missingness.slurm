#!/bin/bash -l
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=2g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=durwa004@umn.edu
#SBATCH -o get_missingness_ind.out
#SBATCH -e get_missingness_ind.err
#SBATCH -p small,amdsmall
module load bcftools

cd /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/

#Missingness_by_individual
/home/mccuem/shared/.local/bin/vcftools --gzvcf thesis_intersect.vcf.gz --missing-indv --out thesis_intersect_miss_by_ind
