# Enter Config Parameters:


#Path to the array data folder:
#E.g. path_array_data <- "/projects/topmed/deprecated/dataprep/studyspecific/phase1/redline_cfs/genotype/ncbi_fingerprint/Fingerprints"
path_array_data <- ""

#Folder Destination for results:
#E.g.path_res_folder_dp <-"/projects/topmed/qc/freeze.2a/analysts/aachueva/cfs"
path_res_folder_dp <-""
#Path to Fam file
#fam<-read.table("/projects/topmed/deprecated/dataprep/studyspecific/phase1/redline_cfs/genotype/ncbi_fingerprint/Fingerprints/del/2016-03-25-11-38_aachueva/2016-03-25-11-38_aachueva_final.fam", as.is =TRUE)
fam <-""

#TOPMed study name (Check in Cathy's database for correctness:/projects/geneva/gcc-fs2/nhlbi_wgs/analysts/cclaurie/NWD_IDs/nhlbi_wgs_ids_v10.RData)
#study name
#E.g. study <- "Sleep"
study <-""

#study PI
#PI <- "Redline"
PI <- ""

#Path to Sample annotation for freeze2
#E.g. path_sa_frz2 <- "/projects/topmed/deprecated/dataworking/phase1/freeze2/sample_annot/topmed_freeze2a_sample_annot_v01.RData"
path_sa_frz2 <- ""

# path to PLINK GDS file
#E.g. path_plink_gds <- "/projects/topmed/deprecated/dataprep/studyspecific/phase1/redline_cfs/genotype/ncbi_fingerprint/Fingerprints/del/2016-03-25-11-38_aachueva/2016-03-25-11-38_aachueva.nhlbi_fingerprint.gds"
path_plink_gds <- ""
 
# Path to scandb folder
#E.g. path_scandb_folder <-"/projects/topmed/deprecated/dataprep/studyspecific/phase1/analysts/levine/results/ids"
path_scandb_folder <- ""

# What is the latest version of scandb in this folder /projects/topmed/dataprep/studyspecific/phase1/analysts/levine/results/ids/ ?
#E.g. scandb_latest_version <-7
scandb_latest_version <-

#Path to Sample Annotation folder:
#E.g. path_to_sample_annot_folder_dataprep <- "/projects/topmed/deprecated/dataprep/studyspecific/phase1/redline_cfs/genotype/ncbi_fingerprint/Fingerprints/sample_snp_annot"
path_to_sample_annot_folder_dataprep <- ""

#Path to vcf gds
#E.g. path_to_vcf_gds <-"/projects/topmed/downloaded_data/IRC_freezes/freeze.2a/gds"
path_to_vcf_gds <-""

#Provide part of the GDS file name before the chromosome number
#E.g. suffix_gds_file_name <-"topmed_freeze2.chr"
suffix_gds_file_name <-""

# provide part of the GDS file name after chromosome number
#E.g. end_gds_file_name <- ".svm_pass_full.sftp-exchange-area.gds"
end_gds_file_name <- ""

# Path to most recent telemetry file for this study:
#E.g. telemetry <- "/projects/geneva/gcc-fs2/nhlbi_wgs/sample_tracking/studies/dbGaP_telemetry/2016-04-07-16-13/phs000954.v1.p1.cfs.xml.txt"
telemetry <- ""
