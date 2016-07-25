# CONTENTS
# generate data frame to request scanIDs
# create sample table                                
# add scan ids to sample table                       
# use telemetry file to map PLINK ID to NWD_ID       
# count NWD_IDs in common between two sample tables
# run discordance function

rm(list=objects())
ptm <- proc.time()

library(SeqVarTools)
library(Biobase)
library(plyr)
library(Hmisc)
date()

# path to array data from the config file
args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
print(config)
#sourcing config file
source(config)

# data stamp with user name
DATE	<- try(system("echo $(date +%Y-%m-%d-%H-%M)_$USER",intern = TRUE))
user<- try(system("echo $USER",intern = TRUE))

# Retrieving path to fam/tfam
plink_path<- paste0(path_array_data,"/plink")
tfam_exist	<- list.files(plink_path, pattern= "\\.tfam$")
if (length(tfam_exist) !=0) {
	fam <- read.table(paste0(plink_path,"/",tfam_exist[1]),as.is=T)
} 

# Path to GDS array
gds <- seqOpen(path_plink_gds)

# Path to the latest scandb database
scandb <- get(load(paste0(path_scandb_folder,"/",scandb_name))); dim(scandb) 

###################################################
# generate data frame to request scanIDs
###################################################

head(fam,n=3)

# add column names
names(fam) <- c("familyID","plink.subjectID","fatherID","motherID","sex","phenotype")

head(fam,n=3)
# fam<-read.table(fam,as.is=TRUE)

# create data frame starting with subject ID used in the PLINK file
df <- data.frame(fam$plink.subjectID,stringsAsFactors=FALSE); dim(df) 
names(df) <- "arrayID"

# add study and PI from cathy's database
# tmp <- get(load("/projects/geneva/gcc-fs2/nhlbi_wgs/analysts/cclaurie/NWD_IDs/nhlbi_wgs_ids_v10.RData")); dim(tmp) 
# table(tmp$study)
# table(tmp$PI)


# if (study %nin% tmp$study) {print("Error, check if study name provided in the config file is correct")}
# if (PI %nin% tmp$PI) {print("Error, check if PI's name provided in the config file is correct")}

df$study <-study
df$PI <- PI

# add path to PLINK/GDS file for array data
df$path <- paste0(path_array_data,"plink")

# check data frame
dim(df) 

head(df,n=3)
df  <- data.frame(df,stringsAsFactors=FALSE);

# write data frame
dir.create(paste0(path_res_folder_dp,"/results"),showWarnings = FALSE)
save(df,file=paste0(path_res_folder_dp,"/results/scanID.request.RData"))

###################################################
# create sample table
###################################################

# get PLINK scanIDs 
sample.id <- seqGetData(gds, "sample.id"); length(sample.id)
head(sample.id)

# create sample annotation with sample.id
sample  <- data.frame(sample.id,stringsAsFactors=FALSE); dim(sample)
meta    <- data.frame(labelDescription=c("unique array identifier as given in PLINK fam file"), 
           stringsAsFactors=FALSE)
sample  <- AnnotatedDataFrame(sample,meta)

head(pData(sample),n=3)

# write annotated data frame (colums: array sample id )
save(sample,file=paste0(path_to_sample_annot_folder,"/",study,".sample.v01.",user,".RData"))

###################################################
# scanIDs
###################################################

table(scandb$PI)

# Making sure we don't duplicate scan'ids
if ( PI %in% unique(scandb$PI)){ print("Error: scan_ids for this PI name already exist in the database")} else {
	
# 	Latest scandb db's
	head(scandb,n=3)

# Scan ID request
	head(df,n=3)
	
	# indices in scandb for next set of scan IDs to assign
	(start <- max(1,min(which(scandb$arrayID %in% ""))))
	(end   <- start + nrow(df) - 1)
	stopifnot(length(start:end)==nrow(df))


	res <- sort(scandb$scanID[start:end],index.return=TRUE)
	names(res)  
	head(res$x)  
	head(res$ix) 
	min(res$ix)
	max(res$ix) 
	idx <- match(1:nrow(df),res$ix)

	# copy FHS data in scandb
	scandb$arrayID[start:end] <- df$arrayID[idx]
	scandb$study  [start:end] <- df$study  [idx]
	scandb$PI     [start:end] <- df$PI     [idx]
	scandb$path   [start:end] <- df$path   [idx]
	head(scandb)
	# checks
	(scanID.min <- min(scandb$scanID[start:end])) 
	(scanID.max <- max(scandb$scanID[start:end])) 
	df$arrayID[1]       
	df$arrayID[nrow(df)] 
	stopifnot(scandb$arrayID[scandb$scanID %in% scanID.min]==df$arrayID[1])
	stopifnot(scandb$arrayID[scandb$scanID %in% scanID.max]==df$arrayID[nrow(df)])

	# save scandb. Attention: check latest versipon of the existent database in this folder: dataprep/studyspecific/phase1/analysts/levine/results/ids , 
	# 	then increment version name by one. (if latest version in the folder is v4, you need to save your file as v5)
	
	newVersion <- scandb_latest_version+1
	newScandbFileName    <- paste0("topmed.identity.check.sample.id.6digits.random.v0",newVersion,".added_",study,".Rdata")
	newScandbFileOut     <- paste(path_scandb_folder,newScandbFileName,sep="/")
	save(scandb, file=newScandbFileOut)
	# 
	# rm(list=objects())
}

###################################################
# add scan ids to sample table
###################################################

# read sample table
# sampleTablePath <- paste0(path_array_data,"sample_snp_annot")
# sampleTableFileName    <- paste0(study,".sample.v01.",DATE,".RData")
# sampleTableFileOut     <- paste(sampleTablePath, sampleTableFileName,sep="/")


sample      <- get(load(paste0(path_to_sample_annot_folder,"/",study,".sample.v01.",user,".RData"))); dim(sample) # 3305  1


head(pData(sample),n=3)
#   sample.id
# 1   PT-3DVG
# 2   PT-3DVH
# 3   PT-3DVI

# read  database with scanID
# if(exists("newScandbFileOut")){
# 	newScandbFileName    <- paste0("topmed.identity.check.sample.id.6digits.random.v05.added_",study,".Rdata")
# 	newScandbFileOut     <- paste(machinePath,scandbFilePath,newScandbFileName,sep="/")
# 	scandb      <- get(load(newScandbFileOut))
# } else {
# 	scandb <- get(load(path_scandb))
# }
 
dim(scandb) # 900000  6


head(scandb,n=3)
#  idx scanID arrayID  study       PI                                                                              path
# 1   1 160278  SG0245 Samoan McGarvey phase1/mcgarvey_sas/genotype/affy6_0/source_data/2015_11_18__Dan/samoa-gwas-dbgap
# 2   2 196673  SG0400 Samoan McGarvey phase1/mcgarvey_sas/genotype/affy6_0/source_data/2015_11_18__Dan/samoa-gwas-dbgap
# 3   3 846961  SG1448 Samoan McGarvey phase1/mcgarvey_sas/genotype/affy6_0/source_data/2015_11_18__Dan/samoa-gwas-dbgap

 table(scandb$study)
# 
#           FHS    JHS Samoan 
# 882175  11401   3305   3119 

table(scandb$PI)
#                  McGarvey Ramachandran       Wilson 
#       882175         3119        11401         3305 

# keep only ID fields for for this project
colsToKeep <- c("scanID","arrayID")
idx <- which(scandb$study %in% study & scandb$PI %in% PI); length(idx) # 11401
scandb <- scandb[idx,colsToKeep]; dim(scandb) # 3305 2

# check match on arrayID id before merge
all(sort(scandb$arrayID)==sort(sample$sample.id)) # TRUE

# merge sample.num into sample table
sample$idx <- 1:nrow(sample)
pData(sample) <- merge(pData(sample),scandb,by.x="sample.id",by.y="arrayID"); dim(sample) # 11401 3
pData(sample) <- pData(sample)[order(sample$idx),]; dim(sample) # 11401 3
sample$idx <- NULL
head(pData(sample),n=3)
#  sample.id scanID
# 2005   14107_1 100061
# 8202    3055_1 100245
# 1376   12775_1 100295
tail(pData(sample),n=3)
# sample.id scanID
# 11329    9886_4 999923
# 11382    9975_4 999946
# 11401    9997_4 999996

# although GDS file is actually ordered by sample.id, should also now be ordered by increasing scanID also
all(sort(sample$scanID)==sample$scanID) # TRUE

# write annotated variant data frame
save(sample,file=paste0(path_to_sample_annot_folder,"/",study,".sample.v02.",user,".RData"))

proc.time() - ptm
