# MAF calculation 
# Goal: concordance between common variants with MAF > 5%

# 1. SeqVarTools has a function called alleleFrequency() that takes a GDS file as an input argument. 
# Can you please compute the allele frequency for the reference allele in all 22 VCF gds files and the array fingerprint file.
# These values should be stored in the SNP annotation file for array data. You will need to create a SNP annotation file (data frame)
# for each chromosome with snp id, chromosome, position, reference allele and alt allele to which you can add allele frequency 
# and MAF as specified in the next step.

library(SeqVarTools); sessionInfo()$otherPkgs$SeqVarTools$Version 
library(Biobase); sessionInfo()$otherPkgs$Biobase$Version     

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
print(config)
#sourcing config file
source(config)

user<- try(system("echo $USER",intern = TRUE))

# compute the allele frequency for the reference allele for the fingerprint file

# open VCF sGDS file

arrayGDS    <- seqOpen(path_plink_gds)

# 2. Since the reference allele in the VCF may not be the same as the reference allele in the array data,
#  calculate the MAF from the allele frequency with the rule maf <- ifelse(allele frequency < 0.5, allele frequency, 1- allele frequency)


# Function format: alleleFrequency(gdsobj, n=0, use.names=FALSE)
afreq <- alleleFrequency(arrayGDS, n=0, use.names=FALSE)
maf <- ifelse(afreq < 0.5, afreq, 1-afreq)
head(maf)


var_annotation_fp <- data.frame(variant.id=seqGetData(arrayGDS, "variant.id"), chromosome=seqGetData(arrayGDS, "chromosome"), position=seqGetData(arrayGDS, "position"), afreq, maf)
head(var_annotation_fp)

save(var_annotation_fp, file=paste0(path_to_sample_annot_folder,"/",study,".sample.v03.MAF.",user,".RData"))

