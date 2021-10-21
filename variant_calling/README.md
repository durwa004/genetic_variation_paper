# variant_calling
Part of chapter 1, 3, and 4 of my thesis - need to split up for publication.
Code to extract variants (from ibio) and then gentoype gvcfs, and extract union/intersect

# JOINT
- Convert gvcfs to vcfs
```
$ python python_generation_scripts/Generate_SelectVariants.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/ -p gatk
$ python python_generation_scripts/Generate_SelectVariants.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/ -p bcftools
$ python ../python_generation_scripts/Generate_pbs_submission_shell.py -d /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/selectvariants_joint/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/selectvariants_joint/pbs_shell.sh
```
Download bcftools/gatk joint genotyped
```
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/s3_scripts/s3cmd_cp_joint_genotyped.pbs 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/s3_scripts/s3cmd_cp_joint_genotyped.pbs 
```
- Pull out Prze horses
```
$ python ../../python_scripts/Generate_bcftools_view_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_gatk/without_Prze/
$ python ../../python_scripts/Generate_bcftools_view_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_bcftools/without_Prze/
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d ../gatk_extract_prze/
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d ../bcftools_extract_prze/
$ python ../../python_scripts/Generate_get_horse_IDs.py -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_bcftools/
$ python ../../python_scripts/Generate_get_horse_IDs.py -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_gatk/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_extract_prze/pbs_shell.sh 
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/gatk_extract_prze/pbs_shell.sh
```
- To do for paper: 
- Get bcftools stats info for gatk and bcfools (including tstv etc) and union and intersect
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/python_scripts/Generate_bcftools_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_gatk_without_Prze/ -e .genotyped.vcf.gz
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/pbs_scripts/bcftools_stats.genotyped.vcf.gz.pbs 
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/python_scripts/Generate_bcftools_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_bcftools/ -e .genotyped.vcf.gz
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/pbs_scripts/bcftools_stats.genotyped.vcf.gz.pbs 
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/python_scripts/Extract_bcftools_stats.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_gatk_without_Prze/ -p gatk
/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_bcftools_without_Prze/ -p bcftools
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/python_scripts/Generate_bcftools_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_union/ -e _union_joint.vcf.gz
```
Move back to my laptop
```
$ scp durwa004@login02.msi.umn.edu://home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_bcftools_without_Prze/bcftools_number_of_variants.txt bcftools_stats_output
$ scp durwa004@login02.msi.umn.edu://home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_gatk_without_Prze/gatk_number_of_variants.txt bcftools_stats_output
```
- Concordance between gatk/samtools - bcftools isec
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling//python_generation_scripts/Generate_bcftools_isec.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/concordance_gatk_bcftools/
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling//python_generation_scripts/Generate_pbs_submission_shell.py -d ../concordance/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/pbs_shell.sh 
```
Tidy away files
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/delete_unnecessary_concordance_files.py
$ sh concordance_gatk_bctools_tidy.sh
```
Output:
vcf1 = bcftools
vcf2 = gatk
0000.vcf just variants in vcf1
0001.vcf just variants in vcf2
0002 variants from vcf1 shared by both
0003 variants from bcf2 shared by both
Extract the concordance information
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/Extract_concordance_bcftools_stats.py 
```
Move back to my laptop
```
$ scp durwa004@login02.msi.umn.edu://home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/concordance_gatk_bcftools/overall_concordance/number_of_variants.txt bcftools_stats_output/gatk_bcftools_concordance.txt
```

- Compare intersect to dbsnp
```
$ wget ftp://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_1/by_assembly/GCA_000002305.1/GCA_000002305.1_current_ids.vcf.gz
$ wget ftp://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_1/by_assembly/GCA_000002305.1/GCA_000002305.1_current_ids.vcf.gz.tbi       
```
# Try with Rob's tools
- Download fasta - can't seem to get the one from MSI to work
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/FASTAs/Equus_cab_nucl_wChrUn1_2.fasta /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/dbsnp/
$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/305/GCA_000002305.1_EquCab2.0/GCA_000002305.1_EquCab2.0_genomic.fna.gz
```
- Rename chromosomes from 1/2/3 to same format as fasta etc
```
$ Create_rename_chroms_list.py
$ bcftools_isec_dbsnp.pbs
$ Replace_fasta_ids.py
```
- Create the new environment and install all requirements (see minus80/locuspocus)
```
$ conda create --name lp_sda python=3
$ source activate lp_sda
```
**Install pysam first**
Need to install using: ```/home/mccuem/durwa004/.conda/envs/lp_sda/bin/python3 -m pip install minus80```
      
- Run for the dbsnp file with EquCab2 fasta
```
$ /home/mccuem/durwa004/.conda/envs/lp_sda/bin/python3 /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/MNEc2to3/scripts/01_make_probeseq.py --vcf GCA_000002305.1.chrs.vcf.gz --fasta GCA_000002305.1_EquCab2.0_genomic_renamed.fna --out dbsnp_probes.fa
```


- Remap to EquCab3 - pull out chrom/pos first - https://www.ncbi.nlm.nih.gov/genome/tools/remap (too big to upload manually)
Can only do upto 250000 lines and up to 4 jobs at a time (tried 100000, 10000 and wouldn't work!)
```
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/bcftools_query_chrom_pos.sh 
$ source activate ensembl-vep
$ chmod +x remap_api.pl
$ split -l 100 ../dbsnp_chrom_pos.txt dbsnp_
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/Split_dbsnp_file.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/dbsnp_chrom_pos_files/
$ source activate ensembl-vep
$ sh EVD_dbsnp/ncbi_remap.sh
```
So many files, need to tidy some away before carrying on
Easier to split the shell script so that I know where I am! 
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/Tidy_remapped_files.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/dbsnp_chrom_pos_files/
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/dbsnp_move_fileswq.pbs
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/Split_dbsnp_file.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/dbsnp_chrom_pos_files/
$ split -l 5000 EVD_dbsnp/ncbi_remap.sh EVD_dbsnp/ncbi_remap_sh
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/Split_dbsnp_shell_script.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/
$ for i in {1..26}; do qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/dbsnp_remap_$i.pbs; done
```
Realized that submitting too many jobs meant they all failed! So try limiting to just 50
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/Reorganize_failed_files.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp/dbsnp_remapped_files/
```
Rename chromosomes and compare to thesis intersect - bcftools isec doesn't seem to work for some reason
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/concordance/dbsnp_intersect_overlap.py
```

- Union scripts for each chromosome
```
$ python ../python_generation_scripts/Generate_union_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/
$ python ../python_generation_scripts/Generate_pbs_submission_shell.py -d /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/joint_calling_union_scripts/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/joint_calling_union_scripts/pbs_shell.sh
```


- Intersect scripts for each chromosome
```
$ python python_generation_scripts/Generate_intersect_joint.pbs -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/
$ python ../python_generation_scripts/Generate_pbs_submission_shell.py -d /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/joint_calling_intersect_scripts/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/joint_calling_intersect_scripts/pbs_shell.sh
```

- Concatentate intersect files for analysis
```
$  python ../python_generation_scripts/Generate_bcftools_concat.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/joint_calling_intersect_scripts/intersect_concat.pbs 
```




# Need to do this:
# INDIVIDUAL
- Scripts in:
```
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/variant_calling/individual_gvcfs/
```
#Generate_GATK_ind_calling_scripts.py - Individual variant calling for gatk
```
python python_generation_scripts/Generate_GATK_ind_calling_scripts.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files -i horse_ids.txt -p GenotypeGVCFs
```
#output files: gatk_individual_calling/gatk_GenotypeGVCFs_YAKU0171A.pbs etc. 

#Generate_union_PBS.py - create union of variant PBS for bcftools, gatk and platypus for the individually called vcfs.
```
python python_generation_scripts/Generate_union_PBS.py -i horse_ids.txt -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files -p bcftools
```
#output files: union_bcftools.pbs  union_gatk.pbs  union_platypus.pbs

#Get union of the 3 VCs: union_all.pbs
#Get intersect of the 3 VCs: intersect_all.pbs
