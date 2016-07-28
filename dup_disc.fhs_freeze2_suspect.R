library(SeqVarTools); sessionInfo()$otherPkgs$SeqVarTools$Version  # "1.9.10"
library(Biobase); sessionInfo()$otherPkgs$Biobase$Version      # "2.30.0"


date() # "Fri Dec 18 11:03:08 2015"

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
print(config)
#sourcing config file
source(config)

chr <- args[3]
print(chr)

suspect <- args[2]
print(suspect)

# chr <-22
#------------------------------------------------------

# read array sample table
sample      <- get(load(paste0(path_to_sample_annot_folder_dataprep,"/",study,".sample.v03.aac.RData"))); dim(sample) # 3305  3
all(sort(sample$scanID)==sample$scanID) 
# TRUE

sample$suspect <- suspect


head(pData(sample),n=3)
#     sample.id scanID    NWD_ID
# 152   PT-3DVG 100151 NWD977574
# 153   PT-3DVH 100189 NWD741479
# 154   PT-3DVI 100452 NWD291345

tmp <- unlist(lapply(pData(sample), function(x) sum(is.na(x)))); tmp[tmp!=0]


# open array sGDS file
arrayGDS    <- seqOpen(path_plink_gds)

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
vcf <- get(load( paste0(path_res_folder_dp,"/results/vcf.sample.v01.acc.RData"))); dim(vcf) 
head(pData(vcf),n=3)
#   sample.id
# 1 NWD101191
# 2 NWD119101
# 3 NWD120808

dim(pData(vcf))

# open VCF sGDS file
fileName    <- paste0(path_to_vcf_gds,"/",suffix_gds_file_name,chr,end_gds_file_name)
vcfGDS      <- seqOpen(fileName)

vcf_sample  <- data.frame(seqGetData(vcfGDS,"sample.id"),stringsAsFactors=FALSE); dim(vcf_sample) # 183  1
names(vcf_sample) <- c("sample.id")
meta    <- data.frame(labelDescription=c("NWD_ID unique scan ID within TOPMed project"), 
           stringsAsFactors=FALSE)
vcf_sample  <- AnnotatedDataFrame(vcf_sample,meta)


# So far we have sample annot df for VCF and array data (fingerprint) and we have SNP annot (with allele frequency and maf frerquency)
# Next step is to run seqSetFilter function to filter those that are MAF < 5% and rerun dup disc function


snp<- data.frame(variant.id=seqGetData(vcfGDS, "variant.id"), chromosome=seqGetData(vcfGDS, "chromosome"), position=seqGetData(vcfGDS, "position"))

load(paste0(path_res_folder_dp,"/results/maf_concordance_vcf_array_chr",chr,".RData"))
 
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
match.samples.on=c("sample.id", "suspect"),
match.variants.on="position", discordance.type="hethom",
by.variant=FALSE, verbose=TRUE)

class(res); names(res)

dim(res) # 260   5
head(res,n=3)

summary(res$n.concordant)

library(GWASTools)

save(res, file=paste0(path_res_folder_dp,"/dup_disc_res_",suspect,"_",chr,".RData"))
rm(list=objects())


#99999999999999999999999999999999999999999999999999999999999999999999999999999999999999
#99999999999999999999999999999999999999999999999999999999999999999999999999999999999999
# START START START START START START START START START START START START START START 
#99999999999999999999999999999999999999999999999999999999999999999999999999999999999999
#99999999999999999999999999999999999999999999999999999999999999999999999999999999999999

