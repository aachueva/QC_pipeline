# Sudo script

#Data set : CFS Freeze 2

# Create a working directory for the project
mkdir /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/test

# copy config file from the code directory:
cp /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/config.R /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/test/


# open the config.R and filling up the blanks similar to examples in the commented lines
vi /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/test/config.R

#################################
#Run initial QC on the input data:
#################################

qsub -q bigmem.q -N QC_input_data  /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/runRscript_array.sh  /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/QC_input_data.R config.R
# Your job 1318188 ("QC_input_data") has been submitted

# Once the job is done, put a path to final gds file in the config file. 

#################################
#Run sample_annot_VCF.R 
#################################

qsub -q bigmem.q -N sample_annot_VCF /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/runRscript_array.sh /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/verified/sample_annot_VCF.R config.R
# Your job 1318295 ("sample_annot_VCF") has been submitted

# Once job is done, make sure to look in the output log for the previous step and check if sample dataframe looks like expected( expected number of rows, df is not empty, etc.). Let Anastasia know if there is a problem.

#################################
# Run  sample_annot_FP.R
#################################

qsub -q bigmem.q -N sampleAnnot_plink /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/runRscript_array.sh /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/verified/sample_annot_FP.R config.R

# results of the above script will create 2 files in the sample annotation folder:
# "Study name".sample.v01."user initials".RData - df with array.ids from the tfa/fam file
# "Study name".sample.v02."user initials".RData - df with array.ids and scan id's from the tfam/fam file
 
#################################
# Run  maf_FP.R 
#################################

qsub -q bigmem.q -N MAF_calc_FP /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/runRscript_array.sh /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/maf_FP.R config.R


#################################
# Run  maf_VCF.R.R
#################################

qsub -q bigmem.q -N MAF_calc_vcf -t 1-22 /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/runRscript_array.sh /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/maf_VCF.R.R config.R


#################################
# Run maf_VCF_FP.R
#################################

# dependency on the previous job
qsub -hold_jid 1318472 -q bigmem.q -N maf_VCF_FP.R -t 1-22 /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/runRscript_array.sh /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/maf_VCF_FP.R config.R

#################################
# Run dd_sample_id.R
#################################
# 
# Put a name for sample_annotation_file in the config file.
# sample_annotation_file <-"Sleep.sample.v02.aachueva.RData"

qsub -q bigmem.q -N dd_sample_id -t 1-22 /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/runRscript_array.sh /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/dd_sample_id.R config.R
# Your job-array 1322997.1-22:1 ("dd_sample_id") has been submitted

#################################
# Run dd_sample_id_rand.R
#################################

qsub -q bigmem.q -N rand_dd_sample_id -t 1-22 /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/runRscript_array.sh /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/dd_sample_id_rand.R config.R

#################################
# Run compare_rand_and_non_ran_dd_res_by_chr_by_sample.R
#################################

qsub -q bigmem.q -N compare_rand_and_non_ran_dd /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/runRscript_array.sh /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/compare_rand_and_non_ran_dd_res_by_chr_by_sample.R config.R


#################################
# Run HetHom.R
#################################

qsub -q bigmem.q -N het_hom /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/runRscript_array.sh /projects/topmed/analysts/aachueva/tmp/aachueva/src/tmp/Hom_Het.R config.R


#################################
# Re-running DD on failed samples
#################################



