library(SeqVarTools); sessionInfo()$otherPkgs$SeqVarTools$Version  # "1.9.10"
library(Biobase); sessionInfo()$otherPkgs$Biobase$Version      # "2.30.0"


date() # "Fri Dec 18 11:03:08 2015"

user<- try(system("echo $USER",intern = TRUE))

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
print(config)
#sourcing config file
source(config)

chr <- args[2]
print(chr)
# chr <-22
#------------------------------------------------------

# read array sample table
sample      <- get(load(paste0(path_to_sample_annot_folder_plink,"/",path_to_sample_annot_file_plink))); dim(sample) # 3305  3
all(sort(sample$scanID)==sample$scanID) 
# TRUE

head(pData(sample),n=3)


set.seed(06042016)
idx <- sample(1:nrow(sample))
sample$rand.sample.id <- sample$sample.id_1[idx]



#     sample.id scanID    NWD_ID
# 152   PT-3DVG 100151 NWD977574
# 153   PT-3DVH 100189 NWD741479
# 154   PT-3DVI 100452 NWD291345

#library(splitstackshape)
#new<- as.data.frame(cSplit(pData(sample),splitCols="sample.id", sep = "_"))
#pData(sample)<- merge(pData(sample),new,by="scanID")
#sample$sample.id_2 <-NULL
#names(sample) <- c("scanID", "sample.id", "sample.id.org")
#head(pData(sample))
#   scanID 	 sample.id	   sample.id.org
# 1 102137   PT-31LT_1     PT-31LT
# 2 104979   PT-31LV_1     PT-31LV
# 3 105081   PT-31LX_1     PT-31LX
# 4 106022   PT-31M2_1     PT-31M2
# 

#save(sample, file=paste0(path_to_sample_annot_folder,"/",study,".sample.v03.",user,".RData"))
tmp <- unlist(lapply(pData(sample), function(x) sum(is.na(x)))); tmp[tmp!=0]


# open array sGDS file
arrayGDS    <- seqOpen(path_plink_gds)


ref <- refChar(arrayGDS)
table(ref)
alt <- altChar(arrayGDS)
table(alt)
 keep <- ref %in% c("A","C","G","T")
 keep_alt <- alt %in% c("A","C","G","T")
 table(keep)
 table(keep_alt)
 seqSetFilter(arrayGDS, variant.sel = keep)
# of selected variants: 1,661,619
 seqSetFilter(arrayGDS, variant.sel = keep_alt)
# of selected variants: 1,869,945
# excluding alleles that are 0s
# ref <- refChar(arrayGDS)
# table(ref)
# keep <- ref %in% c("A","C","G","T")
# seqSetFilter(arrayGDS, variant.sel = keep)

# build data object for array gds
arrayData <- SeqVarData(arrayGDS,sample)
class(arrayData) # "SeqVarData"
arrayData

#------------------------------------------------------

# read VCF sample table
vcf <- get(load( paste0(path_to_sample_annot_folder_vcf,"/vcf.sample.v01.RData"))); dim(vcf) 
head(pData(vcf),n=3)
#   sample.id
# 1 NWD101191
# 2 NWD119101
# 3 NWD120808

telemetry   <- read.table(telemetry,as.is=TRUE,header=TRUE,sep="\t"); dim(telemetry) # 4209   14
 options("width"=200)
 head(telemetry,n=3)


head(pData(vcf))
#   sample.id submitted_subject_id
# 1 NWD100395                18560
# 2 NWD100677                17461
# 3 NWD100944                 6286
# 4 NWD102011                 4215
# 5 NWD102416                20318
# 6 NWD103227                 5367

dim(pData(vcf))

# open VCF sGDS file
fileName    <- paste0(path_to_vcf_gds,"/",suffix_gds_file_name,chr,end_gds_file_name)
vcfGDS      <- seqOpen(fileName)

vcf_sample  <- data.frame(seqGetData(vcfGDS,"sample.id"),stringsAsFactors=FALSE); dim(vcf_sample) # 183  1
names(vcf_sample) <- c("sample.id")



meta    <- data.frame(labelDescription=c("NWD_ID unique scan ID within TOPMed project"), 
           stringsAsFactors=FALSE)
vcf_sample  <- AnnotatedDataFrame(vcf_sample,meta)
pData(vcf_sample) <- merge(pData(vcf_sample),telemetry, by.x="sample.id", by.y ="submitted_sample_id",all.x=TRUE)


#save(vcf_sample, file=paste0(path_res_folder_dp,'/VCF_NWD_subj_id_mapping.RData'))

# So far we have sample annot df for VCF and array data (fingerprint) and we have SNP annot (with allele frequency and maf frerquency)
# Next step is to run seqSetFilter function to filter those that are MAF < 5% and rerun dup disc function


snp<- data.frame(variant.id=seqGetData(vcfGDS, "variant.id"), chromosome=seqGetData(vcfGDS, "chromosome"), position=seqGetData(vcfGDS, "position"))

load(paste0(path_res_folder_dp,"/maf_concordance_vcf_array_chr",chr,".RData"))
 
keepem <-snp$variant.id[var_annotation_vcf_maf_concordance$maf_concordance]

seqSetFilter(vcfGDS, variant.id=keepem)


seqSetFilter(vcfGDS, sample.id=vcf$sample.id)



# build data object for VCF gds
seqData <- SeqVarData(vcfGDS,vcf_sample)
class(seqData) # "SeqVarData"
seqData

#------------------------------------------------------

# duplicate discordance function
res <- duplicateDiscordance(seqData, arrayData,
match.samples.on=c("submitted_subject_id", "rand.sample.id"),
match.variants.on="position", discordance.type="hethom",
by.variant=FALSE, verbose=TRUE)

class(res); names(res)

dim(res) # 260   5
head(res,n=3)

summary(res$n.concordant)

library(GWASTools)

save(res, file=paste0(path_res_folder_dp,'/rand_dup_disc_res_chr',chr,'.RData'))
rm(list=objects())


#99999999999999999999999999999999999999999999999999999999999999999999999999999999999999
#99999999999999999999999999999999999999999999999999999999999999999999999999999999999999
# START START START START START START START START START START START START START START 
#99999999999999999999999999999999999999999999999999999999999999999999999999999999999999
#99999999999999999999999999999999999999999999999999999999999999999999999999999999999999

