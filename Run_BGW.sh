############################################
############################################
# Step 0: Set up directories of the tool, working directory, genome block information file, and number of cores (for parallele computation).
############################################
############################################
# Tool directory
BGW_dir=/home/jyang51/jyang/GITHUB/BGW-TWAS-SS

# Working directory to write intermediate files and output files, unique per gene
wkdir=${BGW_dir}/Example/ExampleWorkDir

# File with fileheads of all genome blocks
Genome_Seg_Filehead=${BGW_dir}/Example/ExampleData/geno_block_filehead.txt

# Gene annotation file or Molecular trait file
GeneExpFile=${BGW_dir}/Example/ExampleData/Gene_Exp_example.txt

# Specify gene name or gene identifier as in the 5th column of the gene expression file.
gene_name=ABCA7

######## Directories for xQTL summary statistic data
# Parent directory of all LD files
LDdir=${BGW_dir}/Example/ExampleData/LDdir
ZScore_dir=${wkdir}/${gene_name}_Zscores # Zscore statistic files

# Number of cores/parallele_jobs to be used/implemented
num_cores=2

################################################################
################################################################
# Generate xQTL summary statistics (i.e., Score Statistics) if only individual-level training data are provided
# Run single-variant eQTL analyses with the molecular quantitative traits and genotype data for the same training samples
################################################################
################################################################
# Directory of genotype files (VCF)
geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files

# Specify the genotype field "GT" (called genotype) or "DS" (imputation dosage) to be used from the VCF files
GTfield=GT

${BGW_dir}/bin/get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneExpFile ${GeneExpFile} \
--geno_dir ${geno_dir} --LDdir ${LDdir} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--GTfield ${GTfield} --num_cores ${num_cores} --clean_output 1

################################################################
################################################################
# Prune Genome Blocks
# Select a subset of genome blocks (up to ${max_blocks}) for to run the EM-MCMC algorithm
# Cis blocks are always selected
# Trans blocks with minimum single-variant eQTL p-value < ${p_thresh} will be ranked by the smallest p-value within block (from the smallest to the largest).
# Top ranked trans blocks will be selected up to ${max_blocks}.
################################################################
################################################################
# p-value threshold, recommend 1e-5 (default) for real studies
p_thresh=0.001

# maximum blocks, recommend 50 (default) for real studies
max_blocks=50

# xQTL summary statistic file directory
Zscore_dir=${wkdir}/${gene_name}_Zscores

${BGW_dir}/bin/prune.sh --wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--ZScore_dir ${Zscore_dir} \
--p_thresh ${p_thresh} --max_blocks ${max_blocks} --clean_output 1


################################################################
################################################################
# Training xQTL weights by EM-MCMC algorithm
################################################################
################################################################
# Sample size
Nsample=499

# Initial hyper parameter values
hfile=${BGW_dir}/bin/hypval.init.txt

# PCP threshold
PCP_thresh=0.0001

# Selected filehead file
select_filehead=${wkdir}/${gene_name}_select_filehead.txt

#### Test generating make file
maf=0.01; em=2; burnin=10000; Nmcmc=10000

#### Run EM-MCMC.sh file
${BGW_dir}/bin/EM-MCMC.sh  --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --select_filehead ${select_filehead} \
--LDdir ${LDdir} --Zscore_dir ${Zscore_dir} \
--Nsample ${Nsample} --maf ${maf} --hfile ${hfile} \
--em ${em} --burnin ${burnin} --Nmcmc ${Nmcmc} \
--PCP_thresh ${PCP_thresh} --num_cores ${num_cores} \
--clean_output 0

################################################################
################################################################
# Predict GReX for test samples with individual-level GWAS data
## Need further test
################################################################
################################################################

### Variables for Step 4
BGW_weight=${wkdir}/${gene_name}_BGW_eQTL_weights.txt
test_geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files
test_geno_filehead=${BGW_dir}/Example/ExampleData/test_geno_filehead.txt
test_pheno=${BGW_dir}/Example/ExampleData/Test_pheno.txt
GTfield_test=GT #or DS

# n_sample=499
# n_sample=158
${BGW_dir}/bin/get_test_grex.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--BGW_weight ${BGW_weight} --test_geno_dir ${test_geno_dir} \
--test_geno_filehead ${test_geno_filehead} \
--GTfield ${GTfield} --test_pheno ${test_pheno} \
--num_cores ${num_cores}

# optional script for removing files
