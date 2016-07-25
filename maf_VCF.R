# Goal: concordance between common variants with MAF > 5%

# 1. SeqVarTools has a function called alleleFrequency() that takes a GDS file as an input argument. 
# Can you please compute the allele frequency for the reference allele in all 22 VCF gds files and the array fingerprint file.
# These values should be stored in the SNP annotation file for array data. You will need to create a SNP annotation file (data frame)
# for each chromosome with snp id, chromosome, position, reference allele and alt allele to which you can add allele frequency 
# and MAF as specified in the next step.

library(SeqVarTools); sessionInfo()$otherPkgs$SeqVarTools$Version  # "1.9.10"
library(Biobase); sessionInfo()$otherPkgs$Biobase$Version      # "2.30.0"

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
print(config)
#sourcing config file
source(config)

# chr <-22
args <- commandArgs(trailingOnly=TRUE)
chr <- args[2]
print(chr)

# compute the allele frequency for the reference allele in all 22 VCF gds files

# open VCF sGDS file
# fileName    <- paste0(file_name,chr,"full.gds")
# fileIn      <- paste(path_gds_vcf,fileName,sep="/")
# vcfGDS      <- seqOpen(fileIn)

fileIn <- paste0(path_to_vcf_gds,"/",suffix_gds_file_name,chr,end_gds_file_name)
print(fileIn)
vcfGDS<- seqOpen(fileIn)

# loading sample annotation file for vcf
load(paste0(path_to_vcf_annot_folder,"/vcf.sample.v01.acc.RData"))

# Subsetting VCF GDS to only sample for the current study
seqSetFilter(vcfGDS)
seqSetFilter(vcfGDS,sample.id=sample$sample.id)
# Checking if subsetting was successful
sample_id<-seqGetData(vcfGDS, "sample.id")
length(sample_id)

# 2. Since the reference allele in the VCF may not be the same as the reference allele in the array data,
#  calculate the MAF from the allele frequency with the rule maf <- ifelse(allele frequency < 0.5, allele frequency, 1- allele frequency)

# Function format: alleleFrequency(gdsobj, n=0, use.names=FALSE)
afreq <- alleleFrequency(vcfGDS, n=0, use.names=FALSE)
maf <- ifelse(afreq < 0.5, afreq, 1-afreq)
head(maf)

var_annotation <- data.frame(variant.id=seqGetData(vcfGDS, "variant.id"), chromosome=seqGetData(vcfGDS, "chromosome"), position=seqGetData(vcfGDS, "position"), afreq, maf)
save(var_annotation, file=paste0(path_to_vcf_annot_folder,"/vcf_var_annot_maf_af_chr",chr,".RData"))


# compute the allele frequency for the reference allele in all 22 VCF gds files and the array fingerprint file.





# The goal is to run the duplicateDiscordance function on only 
#  the bi-allelic SNP variants with MAF > 5%, ideally in both 
# VCF and array data. There are two parts to doing this 
# (a) identifying those variants with MAF > 5% in both and 
# (b) restricting the duplicateDiscordance function to just those SNPs.
# 
# (a)	loop through the ~10K SNPs in the array sample table
#  (maybe use a double loop, once on chromosome 
# [since you need to open the corresponding VCF SNP annotation file]
#  and once over the variants on that chromosome.
#  Load the VCF data frame for the chromosome looping over.
#  Create a column in the array SNP table snp$both.maf.lt.05.
#  Check if the SNP has MAF > 5% in both the array data and the  VCF.
#  If so, mark TRUE, else mark FALSE. If for some reason the SNP does not exist
#  in the Match by chromosome and position. Remember to only look at bi-allelic SNPs 
#  in VCF file as there may also be indels, or tri-allelic SNPs at that position. 
#  Adrienne’s duplicateDiscordance function may have some R code for identifying these cases.
#  If for some reason an array variants does not exist in VCF file set snp$both.maf.lt.05 to NA 
#  so we can find out why the variants is missing from the VCF
# 
# (b)	We need some way to specify the list of SNPs to exclude (snp$both.maf.lt.05==FALSE)
#  from the duplicateDiscordance function. I’m not sure there is a way just yet,
#  and this is something we will need to talk to Adrienne and Stephanie about.
