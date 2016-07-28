ptm <- proc.time()

library(SeqVarTools); sessionInfo()$otherPkgs$SeqVarTools$Version  # "1.9.10"
library(Biobase); sessionInfo()$otherPkgs$Biobase$Version      # "2.30.0"

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
print(config)
#sourcing config file
source(config)

chr <- args[2]
print(chr)

#load vcf dataframe with allele freq and maf
load(paste0(path_to_sample_annot_folder_vcf,"/vcf_var_annot_maf_af_chr",chr,".RData"))
dim(var_annotation)

#load array dataframe with allele freq and maf
load(paste0(path_to_sample_annot_folder_plink,"/maf_FP.aac.RData"))
dim(var_annotation_fp)

var_annotation_vcf <- var_annotation
va_fp <-subset(var_annotation_fp, var_annotation_fp$chromosome==chr)
tmp<-merge(va_fp,var_annotation_vcf,by=c("chromosome", "position"), all.y =TRUE)
dim(tmp)
tmp$maf_concordance <- FALSE

# check if maf in the array file and in vcf file are greater than 0.05
tmp$maf_concordance <- !is.na(tmp$maf.x) & tmp$maf.x >0.05 & tmp$maf.y > 0.05

# clean up the final vcf variant df
tmp2 <- tmp[,c(1,2,7:9)]
# tmp2[,c(3)]<-NULL
# tmp2[,c(3)]<-NULL
names(tmp2) <- c("chromosome","position","afreq","maf","maf_concordance")

# check number of concordant with an array data
table(tmp2$maf_concordance)
var_annotation_vcf_maf_concordance <-tmp2
save(var_annotation_vcf_maf_concordance, file=paste0(path_res_folder_dp,"/maf_concordance_vcf_array_chr",chr,".RData"))

proc.time() - ptm