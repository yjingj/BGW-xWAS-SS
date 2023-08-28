############################################
############################################
# Step 0: Set up directories of the tool, working directory, genome block information file, and number of cores (for parallele computation).
############################################
############################################
# Tool directory
BGW_dir=/home/jyang51/jyang/GITHUB/BGW-xWAS-SS

# Working directory to write intermediate files and output files, unique per gene
wkdir=${BGW_dir}/Example/ExampleWorkDir

# File with fileheads of all genome blocks
Genome_Seg_Filehead=${BGW_dir}/Example/ExampleData/geno_block_filehead.txt

# Gene annotation file or Molecular trait file
GeneInfo=${BGW_dir}/Example/ExampleData/Gene_Exp_example.txt

# Specify gene name or gene identifier as in the 5th column of the gene expression file.
gene_name=ABCA7

######## Directories for xQTL summary statistic data
# Parent directory of all LD files
LDdir=${BGW_dir}/Example/ExampleData/LDdir
Zscore_dir=${wkdir}/${gene_name}_Zscores # Zscore statistic files

# Number of cores/parallele_jobs to be used/implemented
num_cores=2

###################################################
###################################################
# Generate xQTL summary statistics (i.e., Score Statistics) if only individual-level training data are provided
# Run single-variant eQTL analyses with the molecular quantitative traits and genotype data for the same training samples
####################################################
####################################################
# Directory of genotype files (VCF)
geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files

# Specify the genotype field "GT" (called genotype) or "DS" (imputation dosage) to be used from the VCF files
GTfield=GT

${BGW_dir}/bin/get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneInfo ${GeneInfo} \
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

${BGW_dir}/bin/prune.sh --wkdir ${wkdir} --gene_name ${gene_name} \
--GeneInfo ${GeneInfo} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--Zscore_dir ${Zscore_dir} \
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
# CPP threshold
CPP_thresh=0.0001
# Selected filehead file
select_filehead=${wkdir}/${gene_name}_select_filehead.txt

#### Input variables for the MCMC algorithm
# MAF threshold to filter out rare variants
maf=0.01;
# number of EM iterations
em=2;
# number of burnin iterations
burnin=10000;
# number of mcmc iterations
Nmcmc=10000
# minimum prior causal probability for cis and trans xQTL
pp_cis=0.0003; pp_trans=0.0002;
# hyper shape and scale parameters in the inverse-gamma prior for xQTL effect size variance
a_gamma=1; b_gamma=2

#### Run EM-MCMC.sh file
${BGW_dir}/bin/EM-MCMC.sh  --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--GeneInfo ${GeneInfo} --select_filehead ${select_filehead} \
--LDdir ${LDdir} --Zscore_dir ${Zscore_dir} \
--Nsample ${Nsample} --maf ${maf} --hfile ${hfile} \
--em ${em} --burnin ${burnin} --Nmcmc ${Nmcmc} \
--CPP_thresh ${CPP_thresh} --num_cores ${num_cores} \
--pp_cis ${pp_cis} --pp_trans ${pp_trans} \
--a_gamma ${a_gamma} --b_gamma ${b_gamma} \
--clean_output 0

################################################################
################################################################
# Predict genetically regulated molecular traits for test samples with individual-level GWAS data
## Need further test
################################################################
################################################################

### Input variables
BGW_weight=${wkdir}/${gene_name}_BGW_xQTL_weights.txt
test_geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files
test_geno_filehead=${BGW_dir}/Example/ExampleData/test_geno_filehead.txt
test_pheno=${BGW_dir}/Example/ExampleData/test_pheno.txt
GTfield_test=GT #or DS

${BGW_dir}/bin/get_test_trait.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--BGW_weight ${BGW_weight} --test_geno_dir ${test_geno_dir} \
--test_geno_filehead ${test_geno_filehead} \
--GTfield ${GTfield_test} --test_pheno ${test_pheno} \
--quant_pheno "TRUE" \
--num_cores ${num_cores} --clean_output 0

