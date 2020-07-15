# genetic_variation_paper
Scripts and pipeline for the genetic variation paper in 534 horses

# Number of variants called ST/GATK and intersect/union
**bcftools_stats**
- Run bcftools for all files in a directory
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/Generate_bcftools_by_chr.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/ -e _intersect.vcf.gz
```
- Extract bcftools info
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/Extract_bcftools_stats.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_intersect/ -p bcftools
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/bcftools_stats.genotyped.vcf.gz.pbs 
```
- Now I have the concatenated intersect
```
$ bcftools stats thesis_intersect.vcf.gz > thesis_intersect.vcf.gz.stats
```
- Move stats info back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/thesis_intersect.vcf.gz.stats /Users/durwa004/Desktop
```
- Analyze stats and get tstv ratio etc
```
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/R_analysis/bcftools_stats_analysis.R
```

# Determine number of variants per individual - output file should also have breed information
- Need to do for union/intersect
```
$ python ../Generate_bcftools_by_individual.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ -ind /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt 
$ python ../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/bcftools_by_ind/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/bcftools_by_ind/pbs_shell.sh 
$ python Extract_bcftools_by_ind_stats.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ind_bcftools_stats_files/
```
- For gatk/bcftools
```
$ python ../../python_scripts/Generate_bcftools_by_chr_by_ind.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/genotyped_files/bcftools/ -e .genotyped.vcf.gz -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt 
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/pbs_scripts/ind_gatk_bcftools/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/bcftools_stats/pbs_scripts/ind_gatk_bcftools/pbs_shell.sh 
$ python Extract_bcftools_by_ind_by_chrom.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_bcftools_without_Prze/ind_bcftools_stats_files/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -p bcftools
$ python Extract_bcftools_by_ind_by_chrom.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_gatk_without_Prze/ind_bcftools_stats_files/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt -p gatk
```
Transfer to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ind_bcftools_stats_files/intersect_by_ind_number_of_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_union/ind_bcftools_stats_files/union_by_ind_number_of_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect/ind_bcftools_stats_files/*_AF_freq.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/AF_freq_files/
$ $ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_bcftools_without_Prze/ind_bcftools_stats_files/bcftools_ind_number_of_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_gatk_without_Prze/ind_bcftools_stats_files/gatk_ind_number_of_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
```
-- Breed differences in the number of variants per individual
```
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/R_analysis/bcftools_stats_analysis.R
```
-- AF of observed variants
```
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/genetic_burden/python_scripts/merge_bcfstats_AF_191015.py
```

**Input = Intersect of bcftools/gatk haplotype caller (group calling) on all chromosomes**

NB - will need to re-estimate for the individual calls as well.

**All python scripts that run from command line need source activate snakemake to get the correct python version**

- Look at variants across the region
1) Run bcftools stats on regions of 10,000 bp across the genome
```
$ qsub -I -l nodes=1:ppn=1,mem=4g,walltime=12:00:00
$ cd /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_view_regions/
$ source activate snakemake
$ python ../../python_scripts/Generate_get_bcftools_stats_region.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/split_chromosomes/ -v /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/thesis_intersect.vcf.gz -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/GCF_002863925.1_EquCab3.0_genomic/GCF_002863925.1_EquCab3.0_genomic_NC.fna.bed 
$ for i in /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_view_regions/*; do qsub $i; done
$ python ../../python_scripts/Extract_bcftools_stats_region.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/split_chromosomes/
```
2) get average varation across the genome and find regions with more than double/half of the average variation
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/split_chromosomes/regions_number_of_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ bcftools_stats_analysis.R
```
3) pull out high/low regions from snpeff vcf file
```
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/High_variation_regions.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/split_chromosomes/
$ scp /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/Low_variation_regions.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/split_chromosomes/
$  python ../../python_scripts/Generate_bcftools_view_regions.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/ -v /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/thesis_intersect_snpeff.ann.vcf.gz -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/split_chromosomes/High_variation_regions.txt
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_view_regions/Extract_variants_High_variation_regions.pbs
$  python ../../python_scripts/Generate_bcftools_view_regions.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/ -v /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/thesis_intersect_snpeff.ann.vcf.gz -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/split_chromosomes/High_variation_regions.txt
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_view_regions/Extract_variants_Low_variation_regions.pbs
$ python_scripts/Get_high_low_regions_details.py
```
4) Get genes information
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/High_low_variation_regions_all_genes.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ scp /Users/durwa004/Downloads/bioDBnet_db2db_200206165933_1087288859.txt durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/
$ ../python_scripts/Get_high_low_regions_details.py 
```
5) Move back to my laptop for analysis:
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/Low_variation_regions_regions_brief.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/High_variation_regions_regions_brief.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ sed -i -e 's/#/_/g' High_variation_regions_regions_brief.txt
$ sed -i -e 's/#/_/g' Low_variation_regions_regions_brief.txt
$ High_low_region_variation_gene_analysis.R
$ python ../../python_scripts/Get_constraint_metrics_high_low_genes.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/ -t /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/Low_variation_regions_regions_brief.txt 
$ python ../../python_scripts/Get_constraint_metrics_high_low_genes.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/ -t /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/High_variation_regions_all_brief.txt 
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/Low_variation_regions_regions_brief_constraint.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/high_low_regions/High_variation_regions_regions_brief_constraint.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ High_low_region_variation_gene_analysis.R
```
Probably want a figure of impact and or consequence of variant - plus what genes are involved

# Get variants present in all individuals
- All homozygotes
```
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_view_all_homozygous_sites.pbs 
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/view_no_homozygous_sites.pbs/Get_variant_details.py
```
- Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/all_homozygotes/all_homozygotes_details.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/all_homozygotes/lof_variants_all_homozygotes.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/all_homozygotes/thesis_intersect_all_homozygotes.vcf.gz.stats /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
```
- All heterozygous/homozygous
```
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/pbs_scripts/bcftools_view_af_over_0.5.pbs 
$ /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/genetic_burden/python_scripts/Get_variants_present_in_all.py
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip thesis_intersect_variants_present_in_all.vcf 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix thesis_intersect_variants_present_in_all.vcf.gz 
```
- Move back to my laptop
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/present_in_all/variants_present_in_all_details.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/present_in_all/variants_present_in_all_not_homozygous_in_all.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/present_in_all/lof_variants_present_in_all.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/SnpEff/present_in_all/thesis_intersect_variants_present_in_all.vcf.gz.stats /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/gb_analysis/
```
- Analyse data
```
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/genetic_burden/python_scripts/Extract_homozygous_heterozygous_details.py
$ homozygous_heterozygous_analysis.R
```

# Number of variants shared by populations
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_breed_shared_variants.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_breed_shared_variants_part2.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/overlap_of_variants_by_breed.txt
```
- Transfer to my laptop and analyze
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_breed_shared_unique_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/bcftools_stats_output/
$ bcftools_stats_analysis.R
```
- Pull out details for paper
```
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Generate_extract_shared_variants_pop.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/
$ sed -i 's$ xxxxx$$g' bcftools_query_pop_chrom_pos_shared.sh 
$ module load bcftols
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_query_breed_rare_common/bcftools_query_chrom_pos_shared.sh 
```
- Extract the breed_rare_common variants from the snpeff file
```
$ python ../../python_scripts/Generate_extract_variants.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_chrom_pos/ 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_regions_pop/Extract_variants_shared_pop.pbs
$ module load bcftools
$ bcftools stats shared_pop_snpeff.vcf.gz > shared_pop_snpeff.vcf.gz.stats
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/SnpSift_filter_shared.pbs
```
- Extract the breed details from the output vcf for breed rare/common variants
- Need to pull in info re which breed rare/common
```
$ qsub -I -l nodes=1:ppn=8,walltime=12:00:00,mem=40g
$ python  /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_shared_variant_information.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/ -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_chrom_pos/shared_pop_chrom_pos.txt
$ Get_rare_breed_common_breed_info.py
```
- Move files back to my laptop 
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/shared_variants.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
```
- Get genes details
https://biodbnet-abcc.ncifcrf.gov/db/db2db.php
```
$ split -l 2,500 breed_breed_disc_gene_list.txt breed_breed_genes_
```
https://david.ncifcrf.gov/gene2gene.jsp



# Number of variants unique to populations
Plan to split the thesis intersect by breed group
# Number of variants with big differences in frequency between populations (<3% in one and >10% in another)
- Get breed vcfs
```
$ python ../../python_scripts/Generate_bcftools_view_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_shell.sh
```
- Get common/rare variants from each breed
```
$ qsub -I -l nodes=1:ppn=8,walltime=12:00:00,mem=4g
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_rare_common_af_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/ -r 0.03 -c 0.10

$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip {breed}_rare.vcf 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix {breed}_rare.vcf.gz 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip {breed}_common.vcf 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix {breed}_common.vcf.gz 
```
- Get intersect between the files
```
$ python ../../python_scripts/Generate_bcftools_isec_rare_common_variants_by_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/joint_intersect_without_Prze/ -b WP
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d ../bcftools_isec_extract_rare_common_variants/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_isec_extract_rare_common_variants/pbs_shell.sh
```
- Creates a lot of spare files and lots of GB! So tidy up files and get bcftools stats information from them
```
$ delete_unnecessary_files.py 
$ sh breed_common_rare_variants_tidy.sh 
```
- Get chrom/pos list (no header) to extract variants that are rare in breed and common in other breed (0002.vcf.gz) these variants from the SnpEff file
```
$ python ../../python_scripts/Extract_bcftools_stats_isec.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/
$ sed -i 's$ xxxxx$$g' *
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_query_breed_rare_common/bcftools_query_chrom_pos.sh 
```
- Extract the breed_rare_common variants from the snpeff file
```
$ python ../../python_scripts/Generate_extract_variants.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_chrom_pos/ 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_regions/Extract_variants_rare_common_breed.pbs
$ module load bcftools
$ bcftools stats rare_common_breed_snpeff.vcf.gz > rare_common_breed_snpeff.vcf.gz.stats
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/SnpSift_filter_breed_rare_common.pbs
$ cat rare_common_breed_high.vcf.gz rare_common_breed_moderate.vcf.gz rare_common_breed_low.vcf.gz > rare_common_breed_coding.vcf.gz
```
- Extract the breed details from the output vcf for breed rare/common variants
- Need to pull in info re which breed rare/common
```
$ qsub -I -l nodes=1:ppn=8,walltime=12:00:00,mem=40g
$ python  /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_breed_rare_common_shared_variant_information.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_snpeff -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_chrom_pos/
```
- Move files back to my laptop 
```
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_snpeff/breed_rare_other_breed_common.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_rare_common_vcfs/breed_rare_common_snpeff/breed_differences.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ Get_breed_rare_breed_common_details.py
```
- Get genes details
https://biodbnet-abcc.ncifcrf.gov/db/db2db.php
```
$ split -l 2,500 breed_breed_disc_gene_list.txt breed_breed_genes_
```
https://david.ncifcrf.gov/gene2gene.jsp

# Number of variants that are rare in one breed and common in the general population and variants that are common in that breed and rare in the general population 

- Create vcfs without each breed
```
$ python ../../python_scripts/Generate_bcftools_view_extract_breed.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d ../bcftools_view_remove_breed/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_remove_breeds/pbs_shell.sh
```
- Create rare AF <3% and common AF >10% pop vcfs
```
$ qsub -I -l nodes=1:ppn=8,walltime=24:00:00,mem=4g -q lab-long
$ source activate snakemake
$ python /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_rare_common_breed_pop_af.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/ -r 0.03 -c 0.10

$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip without_{breed}_rare.vcf 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix without_{breed}_rare.vcf.gz 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/bgzip without_{breed}_common.vcf 
$ /home/mccuem/durwa004/.conda/envs/ensembl-vep/bin/tabix without_{breed}_common.vcf.gz 
```
- Get intersect between the files
```
$ python ../../python_scripts/Generate_bcftools_isec_rare_common_variants_breed_pop.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/ -i /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/horse_genomes_breeds_tidy.txt 
$ python ../../../../variant_calling/python_generation_scripts/Generate_pbs_submission_shell.py -d ../bcftools_isec_breed_pop_rare_common_variants/
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_isec_breed_pop_rare_common_variants/pbs_shell.sh
```
- Creates a lot of spare files and lots of GB! So tidy up files 
```
$ delete_unnecessary_files.py 
$ module load bcftools
$ sh breed_pop_common_rare_tidy.sh 
```
- Get chrom/pos list (no header) to extract variants that are rare in breed and common in other breed (0002.vcf.gz) these variants from the SnpEff file
NB - to get AFs in the general population - need to use 0003.vcf.gz
NB - to get AFs in the breed - need to use 0002.vcf.gz
```
$ python ../../python_scripts//Extract_bcftools_stats_isec_po.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/
$ sed -i 's$ xxxxx$$g' *
$ qsub -I -l nodes=1:ppn=1,walltime=12:00:00,mem=4g
$ sh /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_query_breed_rare_common/bcftools_query_pop_chrom_pos.sh 
```
- Extract the breed_rare_common variants from the snpeff file
```
$ python ../../python_scripts/Generate_extract_variants_pop.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_chrom_pos/ 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_regions_pop/Extract_variants_breed_common_rare_pop.pbs 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_regions_pop/Extract_variants_breed_rare_common_pop.pbs 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/bcftools_view_extract_regions_pop/Extract_variants_unique.pbs 
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/SnpSift_filter_breed_rare_common_pop.pbs
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/SnpSift_filter_breed_common_rare_pop.pbs
$ qsub /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/pbs_scripts/SnpSift_filter_unique.pbs
```
- Extract the breed details from the output vcfs for breed rare/common variants
```
$ qsub -I -l nodes=1:ppn=1,walltime=12:00:00,mem=4g
$ python  /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_breed_pop_shared_variant_information.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/ -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_chrom_pos/ -b unique
$ python  /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_breed_pop_shared_variant_information.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/ -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_chrom_pos/ -b rare
$ python  /home/mccuem/shared/Projects/HorseGenomeProject/scripts/EquCab3/genetic_burden_pipeline/split_by_breed/python_scripts/Extract_breed_pop_shared_variant_information.py -d /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/ -l /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_chrom_pos/ -b common
```
- Transfer to my laptop for analysis
```
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_unique_differences.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_unique_pop.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_rare_differences.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_rare_pop.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_common_differences.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ scp -r durwa004@login.msi.umn.edu:/home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/breed_intersect_files/breed_pop_rare_common_vcfs/breed_pop_rare_common_snpeff/breed_common_pop.txt /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/breed_pop_variants/
$ /Users/durwa004/Documents/PhD/Projects/1000_genomes/GB_project/genetic_burden_pipeline/genetic_burden/python_scripts/Get_breed_pop_rare_common_details.py
```
- Get genes details
https://biodbnet-abcc.ncifcrf.gov/db/db2db.php
https://david.ncifcrf.gov/gene2gene.jsp

