
################################################################
################################################################
# Step 0: Set up directories of BGW-TWAS tool, input files, and working directory.
################################################################
################################################################
# Tool directory
BGW_dir=~/GIT/BGW-TWAS 

# Specify gene name/identifier as in the 5th column of the gene expression file.
gene_name=ABCA7 

# The gene expression file has the following gene information in the first 5 columns:
# CHR, GeneStart, GeneEnd, TargetID/GeneID_1, GeneName/GeneID_2
# And gene expression levels from column 6 with samples in columns and genes in rows.
GeneExpFile=${BGW_dir}/Example/ExampleData/Gene_Exp_example.txt

# Parent directory of genotype data in VCF files
geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files
# File with filehead of all VCF files as in ${geno_dir}/[filehead].vcf.gz of all genome blocks 
Genome_Seg_Filehead=${BGW_dir}/Example/ExampleData/geno_block_filehead.txt
# Specify the genotype field "GT" (genotype) or "DS" (dosage) to be used from the VCF files
GTfield=DS 

# Working directory to save all output
wkdir=${BGW_dir}/Example/ExampleWorkDir

# Parent directory of all LD files
LDdir=${BGW_dir}/Example/ExampleData/LDdir

# Number of cores/parallele_jobs to be used/implemented
num_cores=2 

################################################################
################################################################
# Step 1: obtain summary statistics (aka Score Statistics)
# Run single-variant eQTL analyses on the training gene expression and genotypes data profiled for the same samples
################################################################
################################################################

${BGW_dir}/bin/Step1_get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneExpFile ${GeneExpFile} \
--geno_dir ${geno_dir} --LDdir ${LDdir} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--GTfield ${GTfield} --num_cores ${num_cores}  

################################################################
################################################################
# Step 2: Prune blocks
# Select a subset of genome blocks (up to ${max_blocks}) for joint model training by BGW-TWAS
# Cis blocks are always selected
# Rank trans blocks with minimum single-variant eQTL p-value < ${p_thresh} by the smallest p-value within block (from the smallest to the largest), and select top ranked trans blocks up to ${max_blocks}.
################################################################
################################################################
p_thresh=0.001 # p-value threshold
max_blocks=100 # maximum blocks

${BGW_dir}/bin/Step2_prune.sh --wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --Genome_Seg_File ${Genome_Seg_File} \
--p_thresh 0.001 --max_blocks 100


select_filehead=${wkdir}/${gene_name}_scores/${gene_name}_select_segments.txt

################################################################
################################################################
# Step 3: Training BVSR gene expression prediction model by EM-MCMC algorithm
################################################################
################################################################

N=499 # sample size
${Scripts_dir}/Step3.sh ${gene} ${geneFile} ${geno_dir} ${Scripts_dir} ${Res_dir} ${LDdir} ${N} ${num_cores}

################################################################
################################################################
# Step 4: Extract genotypes for prediction sample that
# align with eQTL with non-zero effect size from the bayesian training result
################################################################
################################################################


pred_geno_dir=/mnt/YangFSS/data/AMP-AD/Mayo/Genotype/Impute2_1KG
pred_geno_filenames=/mnt/YangFSS/data/AMP-AD/Mayo/Genotype/Impute2_1KG/file_names.txt
pred_pheno_file=/mnt/YangFSS/data/AMP-AD/Mayo/Phenotype/MayoPhenoAD.txt
genotype_format=GT #or DS
tabix_mod=tabix/0.2.6

${Scripts_dir}/Step4.sh ${gene} ${Res_dir} ${pred_geno_dir} ${pred_geno_filenames} ${pred_pheno_file} ${genotype_format} ${tabix_mod}

################################################################
################################################################
# Step 5: Obtain predicted GREX
################################################################
################################################################

${Scripts_dir}/Step5.sh ${gene} ${Scripts_dir} ${Res_dir}


# optional script for removing files
