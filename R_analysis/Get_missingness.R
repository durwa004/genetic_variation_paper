module load R/3.2.1

#sites <- read.table("thesis_intersect_miss_by_ind.imiss", header=TRUE)
sites <- read.table("thesis_intersect_miss_by_site.lmiss", header=TRUE) 

print(summary(sites$F_MISS)) 

o_sites <- sites[order(-F_MISS),]

