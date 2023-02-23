
################################################################
################################################################
# Step 0: Set up directories of the tool, working directory, genome block information file, and number of cores (parallele computation) for implementing BGW-TWAS.
################################################################
################################################################

# Tool directory
BGW_dir=/home/jyang51/jyang/GITHUB/BGW-TWAS-SS

# Working directory to write intermediate files and output files, unique per gene
wkdir=${BGW_dir}/Example/ExampleWorkDir

# Specify gene name/identifier as in the 5th column of the gene expression file.
gene_name=ABCA7

# File with fileheads of all VCF files as in ${geno_dir}/[filehead].vcf.gz of all genome blocks
Genome_Seg_Filehead=${BGW_dir}/Example/ExampleData/geno_block_filehead.txt

# Gene expression file has the following gene information in the first 5 columns:
# CHR, GeneStart (TSS), GeneEnd (TSE), TargetID/GeneID_1, GeneName/GeneID_2
# And gene expression levels from column 6 with samples in columns and genes in rows.
GeneExpFile=${BGW_dir}/Example/ExampleData/Gene_Exp_example.txt

# Number of cores/parallele_jobs to be used/implemented
num_cores=2

################################################################
################################################################
# Step 1: obtain summary statistics (i.e., Score Statistics)
# Run single-variant eQTL analyses on the gene expression and genotypes data for the same training samples
################################################################
################################################################

### Varables for Step 1
# Parent directory of genotype files (VCF)
geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files

# Specify the genotype field "GT" (called genotype) or "DS" (imputation dosage) to be used from the VCF files
GTfield=GT

# Parent directory of all LD files
LDdir=${BGW_dir}/Example/ExampleData/LDdir

${BGW_dir}/bin/Step1_get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneExpFile ${GeneExpFile} \
--geno_dir ${geno_dir} --LDdir ${LDdir} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--GTfield ${GTfield} --num_cores ${num_cores} --clean_output 1

################################################################
################################################################
# Step 2: Prune blocks
# Select a subset of genome blocks (up to ${max_blocks}) for joint model training by BGW-TWAS
# Cis blocks are always selected
# Trans blocks with minimum single-variant eQTL p-value < ${p_thresh} will be ranked by the smallest p-value within block (from the smallest to the largest).
# Top ranked trans blocks will be selected up to ${max_blocks}.
################################################################
################################################################

### Varables for Step 2
# p-value threshold, recommend 1e-5 (default) for real studies
p_thresh=0.001

# maximum blocks, recommend 50 (default) for real studies
max_blocks=50

# eQTL summary Zscore file directory
Zscore_dir=${wkdir}/${gene_name}_Zscores

${BGW_dir}/bin/Step2_prune.sh --wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--ZScore_dir ${Zscore_dir} \
--p_thresh ${p_thresh} --max_blocks ${max_blocks} --clean_output 1


################################################################
################################################################
# Step 3: Training BGW-TWAS/BVSR gene expression prediction model by EM-MCMC algorithm
################################################################
################################################################

### Variables for Step 3
# Sample size
N=499

# Initial hyper parameter values
hfile=${BGW_dir}/bin/hypval.init.txt

# PCP threshold
PCP_thresh=0.0001

# Selected filehead file
select_filehead=${wkdir}/${gene_name}_select_filehead.txt

#### Test one genome block
genome_block="Cis_geno_block"
Zscore_dir=${wkdir}/${gene_name}_Zscores

cd ${wkdir}/${gene_name}_EM_MCMC/
${BGW_dir}/bin/Estep_mcmc -inputSS \
-Zscore ${Zscore_dir}/${genome_block}.Zscore.txt.gz \
-LDcorr ${LDdir}/${genome_block}.LDcorr.txt.gz \
-target_chr 19 -start_pos 1040101 -end_pos 1065571 \
-hfile ${hfile} -maf 0.01 -n 499 -bvsrm -smin 0 -smax 10 \
-win 100 -o ${genome_block} -w 10000 -s 10000 -seed 2022

#### Test generating make file
mkfile="${wkdir}/${gene_name}_EM_MCMC/${gene_name}_BGW.mk"
target_chr=19; start_pos=1040101; end_pos=1065571
maf=0.01; em=2; burnin=10000; Nmcmc=10000

${BGW_dir}/bin/gen_mkf.pl \
--wkdir ${wkdir}/${gene_name}_EM_MCMC --BGW_dir ${BGW_dir} \
--LDdir ${LDdir} --ZScore_dir ${Zscore_dir} --filehead ${select_filehead} \
--hfile ${hfile} --Nsample ${N} --maf ${maf} \
--targ ${target_chr} --start ${start_pos} --end ${end_pos} \
--Nburnin ${burnin} --Nmcmc ${Nmcmc} \
--em ${em} --mf ${mkfile}


#### test Step3_EM-MCMC.sh file
Nsample=499
${BGW_dir}/bin/Step3_EM-MCMC.sh  --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --select_filehead ${select_filehead} \
--LDdir ${LDdir} --Zscore_dir ${Zscore_dir} \
--Nsample ${Nsample} --maf ${maf} --hfile ${hfile} \
--em 2 --burnin 100 --Nmcmc 100 \
--PCP_thresh ${PCP_thresh} --num_cores ${num_cores} \
--clean_output 0

################################################################
################################################################
# Step 4: Predict GReX for test samples with individual-level GWAS data
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
${BGW_dir}/bin/Step4_get_test_grex.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--BGW_weight ${BGW_weight} --test_geno_dir ${test_geno_dir} \
--test_geno_filehead ${test_geno_filehead} \
--GTfield ${GTfield} --test_pheno ${test_pheno} \
--num_cores ${num_cores}

# optional script for removing files
