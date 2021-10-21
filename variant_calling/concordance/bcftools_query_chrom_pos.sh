cd /home/mccuem/shared/Projects/HorseGenomeProject/Data/ibio_EquCab3/ibio_output_files/joint_gvcf/dbsnp/EVD_dbsnp

bcftools query -f '%CHROM:%POS\n' GCA_000002305.1.refseq_chrs.vcf.gz > dbsnp_chrom_pos.txt 
