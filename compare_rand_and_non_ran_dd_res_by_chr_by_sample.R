library(SeqVarTools); sessionInfo()$otherPkgs$SeqVarTools$Version  # "1.9.10"
library(Biobase); sessionInfo()$otherPkgs$Biobase$Version      # "2.30.0"

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
print(config)
#sourcing config file
source(config)

final <- data.frame(chr=integer(),rand_n.concordant=integer(),rand_n.variants=integer(),not_rand_n.concordant=integer(),not_rand_n.variants=integer())
not_rand_n.concordant <-0
not_rand_n.variants <-0
rand_n.concordant <-0
rand_n.variants <-0


for (chr in 1:22) {
rand_filein      <- get(load(paste0(path_res_folder_dp, "/rand_dup_disc_res_chr",chr,".RData"))); 
dim(rand_filein)
non_rand_filein <- get(load(paste0(path_res_folder_dp, "/dup_disc_res_chr",chr,".RData"))); 
dim(non_rand_filein)

not_rand_n.concordant <- not_rand_n.concordant+non_rand_filein$n.concordant
not_rand_n.variants <- not_rand_n.variants+non_rand_filein$n.variants
frac_conc <- round((not_rand_n.concordant /not_rand_n.variants)*100, digits=2)

rand_n.concordant <- rand_n.concordant+rand_filein$n.concordant
rand_n.variants <- rand_n.variants+rand_filein$n.variants
frac_conc_rand <- round((rand_n.concordant /rand_n.variants)*100, digits=2)

}

final_tmp <- data.frame( non_rand_filein$subject.id, non_rand_filein$sample.id.1,non_rand_filein$sample.id.2, not_rand_n.variants, not_rand_n.concordant,frac_conc, rand_filein$subject.id, rand_filein$sample.id.1, rand_filein$sample.id.2, rand_n.variants, rand_n.concordant,frac_conc_rand )


save(final_tmp, file= paste0(path_res_folder_dp,"/comparison_by_sample_rand_non_rand.RData"))

