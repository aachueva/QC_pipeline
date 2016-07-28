# Enter Config Parameters:

#Path to the array data folder:
#E.g. path_array_data <- "/projects/topmed/qc/compare_to_arrays/array_data/fhs_ramachandran/ncbi_fingerprints/sample_snp_annot/FHS.mapped_vcf_telem_array_indexed_dups_sa.RData"
path_array_data <- ""

#Folder Destination for results:
#E.g.path_res_folder_dp <-"/projects/topmed/qc/compare_to_arrays/freeze.3a/fhs_ramachandran/analysts/aachueva/Data"
path_res_folder_dp <-""

#Path to Fam file
#E.g. fam <-"/projects/topmed/downloaded_data/prior_array_data/fhs_ramachandran/ncbi_fingerprints/plink/fp_geno_phs974.tfam"
fam <-""

#TOPMed study name (Check in Cathy's database for correctness:/projects/geneva/gcc-fs2/nhlbi_wgs/analysts/cclaurie/NWD_IDs/nhlbi_wgs_ids_v10.RData)
#study name
#E.g. study <-"FHS"
study <-""

#study PI
#PI <- "Ramachandran"
PI <- ""

# path to PLINK GDS file
#E.g. path_plink_gds <- "/projects/topmed/downloaded_data/prior_array_data/fhs_ramachandran/ncbi_fingerprints/gds/fhs_nhlbi_fingerprint.gds"
path_plink_gds <- ""
 
# Path to scandb folder
#E.g. path_scandb_folder <- "/projects/topmed/qc/compare_to_arrays/scan_db"
path_scandb_folder <- ""

# What is the latest version of scandb in this folder /projects/topmed/dataprep/studyspecific/phase1/analysts/levine/results/ids/ ?
#E.g. scandb_latest_version <-8
scandb_latest_version <-""

#Path to PLINK Sample Annotation folder:
#E.g. path_to_sample_annot_folder_plink <- "/projects/topmed/qc/compare_to_arrays/array_data/fhs_ramachandran/ncbi_fingerprints/sample_snp_annot"
path_to_sample_annot_folder_plink <- ""

#Path to PLINK Sample Annotation file:
#E.g. path_to_sample_annot_file_plink <-"FHS.mapped_vcf_telem_array_indexed_dups_sa.RData"
path_to_sample_annot_file_plink <- ""

#Path to VCF Sample Annotation folder:
#E.g. path_to_sample_annot_folder_vcf <- "/projects/topmed/qc/compare_to_arrays/freeze.3a/fhs_ramachandran/sample_snp_annot"
path_to_sample_annot_folder_vcf <- ""

# Path to the most recent telemetry file for this study:
#E.g. telemetry <- "/projects/geneva/gcc-fs2/nhlbi_wgs/sample_tracking/studies/dbGaP_telemetry/2016-07-21-09-21/phs000974.v1.p1.fhs.xml.txt"
telemetry <- ""

##########################
# Below doesnt need to be changed unless new version of freeze
#Path to Sample annotation for freeze2
#E.g. path_sa_frz2 <- "/projects/topmed/deprecated/dataworking/phase1/freeze2/sample_annot/topmed_freeze2a_sample_annot_v01.RData"
path_sa_frz2 <- "/projects/topmed/qc/freeze.3a/sample_annot/topmed_freeze3a_sample_annot_v01.RData"

#Path to vcf gds
#E.g. path_to_vcf_gds <-"/projects/topmed/downloaded_data/IRC_freezes/freeze.2a/gds"
path_to_vcf_gds <-"/projects/topmed/gds/freeze.3a/passgt.minDP0"

#Provide part of the GDS file name before the chromosome number
#E.g. suffix_gds_file_name <-"topmed_freeze2.chr"
suffix_gds_file_name <-"nhlbi.1575.sftp-exchange-area.keep.freeze3a.chr"

# provide part of the GDS file name after chromosome number
#E.g. end_gds_file_name <- ".svm_pass_full.sftp-exchange-area.gds"
end_gds_file_name <- ".pass.gtonly.gds"

