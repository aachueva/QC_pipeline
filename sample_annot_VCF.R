# CONTENTS
# create sample table                        

rm(list=objects())

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
print(config)

#sourcing config file
source(config)
#####################
# Config Parameters:
# path_gds_data is path to VCF folder
#####################

###################################################
# create sample table
###################################################

library(SeqVarTools); sessionInfo()$otherPkgs$SeqVarTools$Version
library(Biobase);     sessionInfo()$otherPkgs$Biobase$Version  
date() 
user<- try(system("echo $USER",intern = TRUE))

annot<- get(load(path_sa_frz2))
vcf_subset_by_project <- pData(annot[annot$PI==PI & annot$study==study ,])

sample  <- data.frame(vcf_subset_by_project$sample.id,stringsAsFactors=FALSE);
meta    <- data.frame(labelDescription=c("NWD_ID unique scan ID within TOPMed project"), 
           stringsAsFactors=FALSE)
sample  <- AnnotatedDataFrame(sample,meta)
names(pData(sample)) <-"sample.id"
# dir.create(paste0(path_res_folder_dp,"/results"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
path_ann_df_vcf<- paste0(path_to_vcf_annot_folder,"/vcf.sample.v01.",user,".RData")

save(sample,file=path_ann_df_vcf)
head(pData(sample))
