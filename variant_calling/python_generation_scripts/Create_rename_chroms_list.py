import os

#Import key/value pairs of NCBI chromosomes/ensembl chromosomes
#Can use bcftools annotate --rename-chrs to do this for me.
#Need a file of old name \t new name

chrom_dict = {}
count = 0
with open("../../GCF_002863925.1_EquCab3.0_genomic/GCF_002863925.1_EquCab3.0_genomic_NC.fna.bed", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        count +=1
        chrom_dict[line[0]] = count
chrom_dict['NC_001640.1'] = "MT"
chrom_dict['NC_009175.3'] = "X"

with open("rename_dbsnp_chrs.txt", "w") as output_file:
    for key, value in chrom_dict.items():
        print(value, key, sep = "\t", file = output_file)

#Get list of 1/2/3 etc to what is in fasta etc
with open("rename_dbsnp_chrs_with_chr.txt", "w") as output_file, open("GCA_000002305.1_EquCab2.0_assembly_report.txt", "r") as input_file:
    for line in input_file:
        line = line.rstrip("\n").split("\t")
        if "#" in line[0]:
            continue
        else:
            print(line[0], line[4], sep = "\t", file = output_file)
