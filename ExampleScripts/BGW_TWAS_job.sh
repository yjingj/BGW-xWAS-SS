#!/bin/bash

num_cores=$1
# num_cores=2

gene_name=TP63
# gene_name=$2

# N=$3
N=337 # sample size


##########
wkdir=/mnt/YangFSS/data2/jyang/testBGW/TP63_wk

BGW_dir=/home/jyang/GIT/BGW-TWAS
GeneExpFile=/mnt/YangFSS/data2/BGW-TWAS-GTEx/ExpressionFiles_TIGAR/Breast_Mammary_Tissue_GTEx_Exp.txt
geno_dir=/mnt/YangFSS/data2/GTEx/GTEx_V8/GenotypeFiles/WGS/LDdetect_EUR_Segmented_VCF
GTfield=GT

LDdir=/mnt/YangFSS/data2/jyang/testBGW/LDdir_GTExV8

### Varables for Step 2
p_thresh=0.001 # p-value threshold
max_blocks=50 # maximum blocks

### Variables for Step 3
hfile=${BGW_dir}/Example/hypval.txt
PCP_thresh=0.0001

### Variables for Step 4
#BGW_weight=/mnt/YangFSS/data2/BGW-TWAS-GTEx/weights/${gene_name}_BGW_eQTL_weights.txt
#test_geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files
#test_geno_filehead=${BGW_dir}/Example/ExampleData/test_geno_filehead.txt
#test_pheno=${BGW_dir}/Example/ExampleData/Test_pheno.txt
#GTfield_test=GT #or DS

################################################################
################################################################
# Step 1: obtain summary statistics (i.e., Score Statistics)
# Run single-variant eQTL analyses on the gene expression and genotypes data for the same training samples
################################################################
################################################################
Genome_Seg_Filehead=/mnt/YangFSS/data2/jyang/testBGW/TP63_wk/TP63_select_filehead.txt
#Genome_Seg_Filehead=/mnt/YangFSS/data2/GTEx/GTEx_V8/GenotypeFiles/WGS/LDdetect_EUR_Segmented_VCF/FileHead_EUR_1703.txt

${BGW_dir}/bin/Step1_get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneExpFile ${GeneExpFile} \
--geno_dir ${geno_dir} --LDdir ${LDdir} \
--Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--GTfield ${GTfield} --num_cores ${num_cores} --clean_output 1


