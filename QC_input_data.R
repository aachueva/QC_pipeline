#CONTENTS
######################
# FP:
# Check data type (Fingerprints or previous array)
# Check if number of columns is consistent in the TPED file
# Check # rows in the TPED file
# Exclude variants with no genotype
# Check TFAM file for duplicates
# Convert PLINK to BED format (in case of fingerprint data) and subsequently to GDS format

# Previous array data:
# Check fam file for duplicates
# Exclude variants with no genotype
# Check TFAM file for duplicates
# Convert PLINK to BED format (in case of fingerprint data) and subsequently to GDS format
#####################

rm(list=objects())

library(SeqVarTools); sessionInfo()$otherPkgs$SeqVarTools$Version
library(Biobase);     sessionInfo()$otherPkgs$Biobase$Version
date()

ptm <- proc.time()

args <- commandArgs(trailingOnly=TRUE)
config <- args[1]
print(config)

#sourcing config file
source(config)
print(path_array_data)
list.files(path_array_data)

# Create a unique signature(timastamp+username) that will be used for folder and file names in order to avoid overwritten files
DATE	<- try(system("echo $(date +%Y-%m-%d-%H-%M)_$USER",intern = TRUE))
dir.create(paste0(path_array_data,"/gds"),showWarnings = FALSE)
system(paste0("mkdir ",path_array_data,"/gds/",DATE))

###################################################
# Check data type (Fingerprints or previous array)
###################################################
data_type_check <- function(plink_path){
	# check if it is array or fingerprint data:

	bim_exist	<- list.files(plink_path, pattern= "\\.bim$")
	fam_exist	<- list.files(plink_path, pattern= "\\.fam$")
	bed_exist	<- list.files(plink_path, pattern= "\\.bed$")
	tfam_exist	<- list.files(plink_path, pattern= "\\.tfam$")
	tped_exist	<- list.files(plink_path, pattern= "\\.tped$")
 
	if (length(tfam_exist) !=0 &  length(tped_exist)!=0) {data_type <- "FP"}  else if (length(bim_exist) !=0 & length(fam_exist)!=0 & length(bed_exist)!=0) {data_type <- "prev_array"} else data_type <- "na" 
	data_type
	output<- list(data_type,bim_exist,fam_exist,bed_exist,tfam_exist,tped_exist)
	return(output)
}

###################################################
# Check if number of columns is consistent in the TPED file
###################################################
tped_NF_check 	<- function(tped) {
        #bash command to check how many fields in the tped file
    	print("#############")
		print("Function tped_NF")
		print("#############")
		
        sys_cmd <- paste0("awk '{print NF}' ",tped,"| sort -u")
        tped_NF <- try(system(sys_cmd,intern = TRUE))
        return(tped_NF)
}

###################################################
# Check # rows in the TPED file (should be =10000)
###################################################
tped_nrow_check 	<- function(tped) {

		print("#############")
		print("Function tped_nrow")
		print("#############")
		
        #bash command to check how many rows in the tped file (performs faster than reading file in R)
        sys_cmd2 	<- paste0("awk 'END {print NR}' ",tped)
        tped_nrows 	<- try(system(sys_cmd2,intern = TRUE))
        return(tped_nrows)
}

###################################################
# Convert BED to GDS
###################################################
bed_to_gds		<- function(path_array_data,data_type,data_type_output) {

		print("#############")
		print("Function bed_to_gds")
		print("#############")
		
		print("Data Type is:")
		print(data_type)
			
		if (data_type =="FP"){	
			bed.fn	<- paste0(path_array_data,"/gds/",DATE,"/",DATE,"_final.bed")
			fam.fn	<- paste0(path_array_data,"/gds/",DATE,"/",DATE,"_final.fam")
			bim.fn	<- paste0(path_array_data,"/gds/",DATE,"/",DATE,"_final.bim")
			# seqBED2GDS(bed.fn, fam.fn, bim.fn, paste0(path_array_data,"del/",DATE,"/",DATE,".nhlbi_fingerprint.gds"))
# 			seqSummary(paste0(path_array_data,"del/",DATE,"/",DATE,".nhlbi_fingerprint.gds"))
		} else {
			bim.fn	<- paste0(path_array_data,unlist(data_type_output[2]))
			fam.fn	<- paste0(path_array_data,unlist(data_type_output[3]))
			bed.fn <- paste0(path_array_data,unlist(data_type_output[4]))
			}
		seqBED2GDS(bed.fn, fam.fn, bim.fn, paste0(path_array_data,"/gds/",DATE,"/",DATE,".nhlbi_fingerprint.gds"))
		seqSummary(paste0(path_array_data,"/gds/",DATE,"/",DATE,".nhlbi_fingerprint.gds"))		
        return()
}

###################################################
# Exclude variants with no genotype (0,0)
###################################################
exclude_snps_with_no_gt <-function(tped,tfam.file,path_array_data,DATE){
		
		print("#############")
		print("Function exclude_snps_with_no_gt")
		print("#############")
		
        #  Generate a list genotyping/missingness rate statistics.This option creates two files:plink.imiss, plink.lmiss which detail missingness by individual and by SNP (locus)
        # AC. There is gotta be a faster way to filter out genotypes in plink. Need to investigate
        system(paste0("plink --tped ",tped," --tfam ",tfam.file," --missing  --recode --out ",path_array_data,"/gds/",DATE,"/",DATE,".cleaned"))

        # create bed file
        system(paste0(" plink --file ",path_array_data,"/gds/",DATE,"/",DATE,".cleaned --make-bed --out ",path_array_data,"/gds/",DATE,"/",DATE,"_plink.out"))

        # Count how many SNPs with no genotype
        system(paste0("cut -f5 ",path_array_data,"/gds/",DATE,"/",DATE,"_plink.out.bim| sort -n|uniq -c"))

        # make lmiss file tab delimited
        system(paste0("awk -v OFS='\t' '$1=$1' ",path_array_data,"/gds/",DATE,"/",DATE,".cleaned.lmiss> ",path_array_data,"/gds/",DATE,"/",DATE,".cleaned.lmiss_tab"))

        # Create list of missing snps
        system(paste0("awk '($5 >=0.999){print}' ",path_array_data,"/gds/",DATE,"/",DATE,".cleaned.lmiss_tab|cut -f2 > ",path_array_data,"/gds/",DATE,"/",DATE,"_missing_snps"))

        # make BED file
        system(paste0("plink --file ",path_array_data,"/gds/",DATE,"/",DATE,".cleaned --exclude  ",path_array_data,"/gds/",DATE,"/",DATE,"_missing_snps --make-bed --out ",path_array_data,"/gds/",DATE,"/",DATE,"_final"))

        # Check if snps's with missing genotype are not included this time
        system(paste0("cut -f5 ",path_array_data,"/gds/",DATE,"/",DATE,"_final.bim| sort -n|uniq -c"))
        return()
}

###################################################
# Check TFAM file for duplicates
###################################################
#  If TFAM file contains duplicated subjects then use first column(family.id_platform_subject.id) as subject id in the downstream analysis
check_tfam <- function(tfam){

	print("#############")
	print("Function check_tfam")
	print("#############")
	
	tfam.fn <-read.table(tfam)
	head(tfam.fn)
	tail(tfam.fn)
	# Check if subjectID is unique, if not than continue
	if (anyDuplicated(tfam.fn$V2) == 0) {final <- "0"} else {
		tfam.fn$og_idx	<- 1:nrow(tfam.fn)
		tfam.fn 		<- tfam.fn[order(tfam.fn$V2),]
		head(tfam.fn)
		tail(tfam.fn)
		index		<-sapply(1:length(tfam.fn$V2),function(x)sum(tfam.fn$V2[max(1,x-10):x]==tfam.fn$V2[x]))
		new_fam		<-cbind(tfam.fn, index)
		table(new_fam$index)
		new_fam$SUBJ_ID_idx <- paste(new_fam$V2,new_fam$index, sep="_")
		head(new_fam)
		new_fam$V2 	<-new_fam$SUBJ_ID_idx
		library(data.table)
		new			<-as.data.table(new_fam)
		new[,c("index", "SUBJ_ID_idx"):=NULL]
		final		<-new[order(new$og_idx)]
		head(final)
		tail(final)
		final$og_idx <- NULL
# 		final 		 <- final[-1,]
		write.table(final,file=paste0(path_array_data,"/gds/",DATE,"/dup_subj_idx.tfam"),sep="\t", quote=FALSE, row.names=F, col.names=F)
		final<- "1"}
	
	return(final)
}
	
# Main function
main <- function(x) {
		print("#############")
		print("Main function")
		print("#############")
		
		plink_path<- paste0(path_array_data,"/plink")
		data_type_output	<- data_type_check(plink_path)
		# data_type_check returns following parameters:	data_type,bim_exist,fam_exist,bed_exist,tfam_exist,tped_exist
		data_type<- unlist(data_type_output[1])
		
		print(data_type)
		
		if(data_type == "FP"){
			tfam	<- paste0(plink_path,"/",unlist(data_type_output[5]))
			tped	<- paste0(plink_path,"/",unlist(data_type_output[6]))
			
			print(tfam)
			print(tped)
			
			tped_NF    <- tped_NF_check(tped)
        	tped_nrows <- tped_nrow_check(tped)
        	
        	if(length(tped_NF) >1 | tped_nrows != 10000) return(system ("echo 'Error: Inconsistent number of columns or number of rows is not 10000 in the TPED file. Correct the file and resubmit the job'| mail -s 'QC_pipeline job status: Error' $USER@uw.edu"))

        	tfam_dup   <- check_tfam(tfam)
			if(tfam_dup != 0) {tfam.file <- paste0(path_array_data,"/gds/",DATE,"/dup_subj_idx.tfam")} else {tfam.file <- tfam}
			
			print(paste0("tfam.file",tfam.file))
			exclude_snps_with_no_gt(tped,tfam.file,path_array_data,DATE)
			bed_to_gds(path_array_data,data_type,data_type_output)
		} else if(data_type == "prev_array"){
			fam	<- paste0(plink_path,"/",unlist(data_type_output[3]))
			
			#check for duplicates 
			dup_exists <-try(system(paste0("cat ",fam,"|cut -f2 -d ' '|uniq -d")))
			if (dup_exists ==0){
				bed_to_gds(path_array_data,data_type,data_type_output)
			}
		} else if (data_type=="na"){
			return(system ("echo 'Error: Unable to find correct type of input data. Accepted data types: for Fingerprints: TPED, TFAM; for the previous array data: BED,FAM,BIM. Make sure that folder provided in the config file contains listed files and resubmit the job'| mail -s 'QC_pipeline job status: Error' $USER@uw.edu"))
		}
		
		print(paste0("Results are written to ",path_array_data,"/gds/",DATE))   
	    
        system (paste0("echo 'Results for QC_pipeline STEP 1: QC on the Input Data can be reviewed here: '",path_array_data,"/gds/",DATE,"   Log file for the executed script: ",getwd(),"| mail -s 'QC_pipeline job status: STEP1 is DONE' $USER@uw.edu"))
}

main(x)

proc.time() - ptm