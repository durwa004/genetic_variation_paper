setwd("/Users/durwa004/Documents/Research/PhD_papers_for_publication/chapt2_genetic_variation/Final/Final_final/Frontiers in genetics/Submitted/Review/new_analysis/")

#No homozygotes
data <- read.table("no_homozygotes_summary.txt", header=T)
mean(data$AC)
range(data$AC)
mean(data$AF)
range(data$AF)

t.test(data$AF[data$impact == "MODIFIER"], data$AF[data$impact != "MODIFIER"])

table(data$consequence)
length(data$consequence)

#a = data.frame(table(data$gene[data$consequence == "3_prime_UTR_variant"]))
#a = data.frame(table(data$gene[data$consequence == "5_prime_UTR_premature_start_codon_gain_variant"]))
a = data.frame(table(data$gene[data$consequence == "frameshift_variant"]))
a = data.frame(table(data$gene[data$consequence == "frameshift_variant&splice_region_variant"]))
a = data.frame(table(data$gene[data$consequence == "missense_variant"]))
a = data.frame(table(data$gene[data$consequence == "splice_acceptor_variant&intron_variant"]))
a = data.frame(table(data$gene[data$consequence == "splice_acceptor_variant&splice_donor_variant&intron_variant"]))
a = data.frame(table(data$gene[data$consequence == "stop_gained"]))
a = data.frame(table(data$gene[data$consequence == "synonymous_variant"]))
#All homozygotes
data1 <- read.table("all_homozygotes_summary.txt", header=T)
mean(data1$AC)
range(data1$AC)
mean(data1$AF)
range(data1$AF)
table(data1$impact)

t.test(data1$AF[data1$impact == "MODIFIER"], data1$AF[data1$impact != "MODIFIER"])

#Present in all - need to run this and then compare the heterozygous to the homozygous variants
data2 <- read.table("variants_present_in_all_summary.txt", header=T)
mean(data2$AC)
range(data2$AC)
mean(data2$AF)
range(data2$AF)
table(data2$impact)
table(data2$lof)

t.test(data2$AF[data2$impact == "MODIFIER"], data2$AF[data2$impact != "MODIFIER"])

#Present in all but not homozygous in all
data3 <- read.table("variants_present_in_all_not_homozygous_in_all_summary.txt", header=T)
mean(data3$AC)
range(data3$AC)
mean(data3$AF)
range(data3$AF)
table(data3$impact)

t.test(data3$AF[data3$impact == "MODIFIER"], data3$AF[data3$impact != "MODIFIER"])

#Compare all homozygotes with those that are present in all but not all homozygotes
t.test(data3$AF[data$impact == "MODIFIER"], data1$AF[data$impact == "MODIFIER"])
t.test(data3$AF[data$impact != "MODIFIER"], data1$AF[data$impact != "MODIFIER"])
