library(plyr)
library(splitstackshape)
library(Biobase)
library(SeqVarTools)

# Loading working version of dataframe for samples that are in common between the array and VCF file (duplicate discordance results per sample)
args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
print(config)
#sourcing config file
source(config)
user<- try(system("echo $USER",intern = TRUE))

load(paste0(path_res_folder_dp,"/comparison_by_sample_rand_non_rand.RData"))
names(final_tmp) <- c( "array.id","NWD.id","broad.id","var.tested", "var.concordant", "frac.concordant", "array.id.rand","NWD.id.rand","broad.id.rand", "vars.tested.rand", "var.concordant.rand", "frac.concordant.rand")
head(final_tmp)
#   NWD.id array.id vars_tested var_concordant frac_concordant NWD.id_rand array.id_rand vars_tested_rand var_concordant_rand frac_concordant_rand
# 1 NWD104494  13303_1         356            356          100.00   NWD104494       12215_1              348                 191                54.89
# 2 NWD104494  13303_2        9880           9876           99.96   NWD104494       13084_3             9883                5284                53.47
# 3 NWD104582  14718_1         243            242           99.59   NWD104582         278_3              244                 130                53.28
# 4 NWD105166  15001_1         356            356          100.00   NWD105166        4379_1              356                 186                52.25

dim(final_tmp)
# 1833   10

# 2. add scan id
# 3. add missing call rate
# 4. add Het, Hom rate

###################################################
# calculate sample missing call rate
###################################################

# read sample table
sample<-get(load(paste0(path_to_sample_annot_folder_plink,"/",path_to_sample_annot_file_plink)))
all(sort(sample$scanID)==sample$scanID) #  TRUE
head(pData(sample),n=3)
#    sample.id_1 scanID sample.id    NWD_ID biosample_id
# 5886       14107 100061   14107_1 NWD815656 SAMN03985786
# 1196        3055 100245    3055_1      <NA>         <NA>
# 5323       12775 100295   12775_1 NWD827240 SAMN03987869

# open fingerprint sGDS file
gds         <- seqOpen(path_plink_gds)

# missing call rate by sample
res <- missingGenotypeRate(gds, margin="by.sample")
class(res) # numeric
sample$missing.e1 <- res
summary(sample$missing.e1)
hist(res)

###################################################
# calculate het/hom call rate
###################################################

sample$het.var <- heterozygosity(gds, margin="by.sample")
sample$homnr <- homozygosity(gds, margin="by.sample", allele="alt")
hethom <- sample$het.var/sample$homnr
# Generate histogram for Hom/Het for the fingerprint file (3305 samples)
hist(hethom, main="", xlab="Het/Hom Non-Ref")

# write annotated data frame(s)
fileOut     <- paste0(path_res_folder_dp,"/sample_annot_mcr_het_hom.RData")
save(sample,file=fileOut)


head(pData(sample))
#  sample.id.org scanID sample.id    NWD_ID biosample_id missing.e1   het.var     homnr
# 5886       14107 100061   14107_1 NWD815656 SAMN03985786  0.9638892 0.4515235 0.4570637
# 1196        3055 100245    3055_1      <NA>         <NA>  0.9639892 0.4000000 0.4888889
# 5323       12775 100295   12775_1 NWD827240 SAMN03987869  0.9638892 0.4321330 0.4570637
# 2419        5827 100328    5827_1      <NA>         <NA>  0.9638892 0.3988920 0.5096953

dim(pData(sample))
# [1] 11401    8


seqSetFilter(gds)
# of selected samples: 11401
# of selected variants: 9997

final_df <- merge(final_tmp, pData(sample), by.x="array.id",by.y="sample.id")
# final_df$NWD.id <- NULL

# Het/Hom rate for 176 samples only:
hethom <- final_df$het.var/final_df$homnr

# Generate histogram for Hom/Het for the samples that are in common (260 samples)
# hist(hethom, main="", xlab="Het/Hom Non-Ref")


# Set flag for failed samples if concordance is less than 60%
final_df$FLAG <- "PASS"
final_df$FLAG[final_df$frac.concordant<90] <- "FAIL"

table(final_df$FLAG)

# PASS 
#  260 


# Rounding numbers
final_df$missing.e1 <- round(final_df$missing.e1,digit =4)
final_df$het.var <- round(final_df$het.var,digit =4)
final_df$homnr <- round(final_df$homnr,digit =4)

head(final_df)
#      NWD.ID array.ID scanID vars_tested var_concordant frac_concordant missing.e1 het.var  homnr FLAG
# 1 NWD898460   SG0002 100612      555158         554531           99.89     0.0011  0.2595 0.6803 PASS
# 2 NWD570469   SG0004 101189      555331         554760           99.90     0.0009  0.2366 0.7057 PASS
# 3 NWD290601   SG0009 103243      555245         554563           99.88     0.0012  0.2477 0.6910 PASS
# 4 NWD668738   SG0036 111212      555452         554948           99.91     0.0006  0.2388 0.7016 PASS

save(final_df,file= paste0(path_res_folder_dp,"/final_result_df.RData"))

# Add family ID

# reading Tfam file:
fam <- read.table(fam, as.is=T)

# reading final df
df<- get(load(paste0(path_res_folder_dp,"/final_result_df.RData")))

df <- merge(df, fam, by.x="array.id", by.y="V2")
library(dplyr)
df1<-select(df, array.id:V1)
 head(df1)
#   array.id    NWD.id  broad.id var.tested var.concordant frac.concordant array.id.rand NWD.id.rand broad.id.rand vars.tested.rand var.concordant.rand frac.concordant.rand scanID sample.id missing.e1 het.var  homnr FLAG                   V1
# 1  PT-31LT NWD464643 PT-31LT_1        360            360          100.00       PT-31LT   NWD464643     PT-32FO_2             7179                3740                52.10 102137 PT-31LT_1     0.9514  0.3972 0.4750 PASS phs000284_26_PT-31LT
# 2  PT-31LV NWD464023 PT-31LV_1        359            358           99.72       PT-31LV   NWD464023     PT-4API_1              360                 195                54.17 104979 PT-31LV_1     0.9516  0.3955 0.4763 PASS phs000284_26_PT-31LV
# 3  PT-31MG NWD142436 PT-31MG_1        360            360          100.00       PT-31MG   NWD142436     PT-4AZC_1              361                 188                52.08 108568 PT-31MG_1     0.9514  0.3889 0.4833 PASS phs000284_26_PT-31MG
# 4  PT-31MI NWD991428 PT-31MI_1        359            359          1
# 
# Next split fam id on platform and phs number


# split fam.id into three fields and add to result table
library(stringr)
tmp <- as.data.frame(str_split_fixed(df1$V1, "_", 3),stringsAsFactors=FALSE); names(tmp) <- c("phs","platformID","subjectID")
dim(tmp) # 4835 3
head(tmp)
##         phs platformID subjectID
## 1 phs000282         26     10006
## 2 phs000342         71     10006
## 3 phs000342         46     10006
## 4 phs000282         26      1001
## 5 phs000342         71     10011
## 6 phs000342        116      1001
all(df1$V1==paste(df1$phs,df1$platformID,df1$subjectID,sep="_")) # TRUE
df1$phs        <- tmp$phs
df1$platformID <- tmp$platformID
df1$V1<-NULL

save(df1,file= paste0(path_res_folder_dp,"/final_result_df_v2.RData"))



