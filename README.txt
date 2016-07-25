STEP 1

copy config file to your directory

STEP 2
Test
Test2
Test3
Test4

fill out blanks in the config file


STEP 3
a)QC on Input Data. 
#############
What it does:
#############

1. For FP and previous array data files,the script:
	checks number of rows and columns in the tped file;
	excludes variants with no genotype(they are problematic for the duplicate discordance function);
	checks for duplicated subjects in the tfam file and creates indeces if latter is true.
2. Converts PLINK to BED format (in case of fingerprint data) and subsequently to GDS format

###########
How to run:
###########

run script: qsub -q <queue_name> -N QC_input_data /projects/topmed/analysts/aachueva/src/runRscript_array.sh /projects/topmed/analysts/aachueva/tmp/aachueva/src/QC_input_data.R config.R
