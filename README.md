# Bayesian Genome-wide (BGW) TWAS Software Usage

Tool **BGW-TWAS** tool is developed for leveraging both **cis-** and **trans-** eQTL based on a Bayesian variable selection model to predict genetically regulated gene expression (**GReX**) and then conduct **TWAS**. This tool is implemented through several command-line argument described in this manual.

Please cite our BGW-TWAS paper if you use the tool:
>[*Bayesian Genome-wide TWAS Method to Leverage both cis- and trans-eQTL Information through Summary Statistics.* 2020 AJHG.](https://www.cell.com/ajhg/pdfExtended/S0002-9297(20)30291-3)

---
- [Software Installation](#software-installation)
- [Input Files](#input-files)
	- [1. Gene Expression File](#1-gene-expression-file)
	- [2. Genotype Files for Training](#2-genotype-files-for-training)
	- [3. Training BGW-TWAS](#3-training-BGW-TWAS)
- [Example Usage](#example-usage)
	- [1. Obtain Summary Statistics](#1-obtain-summary-statistics)
	- [2. Prune Genome Segments](#2-prune-genome-segments)
---

## Software Installation

### 1. Install required C++ libraries
C++ libraries *zlib*, *gsl*, *eigen3*, *lapack*, *atlas*, *blas* are used to develop this tool. Please install these libraries to your system and include the library path `-I[path to libraries]` accordingly in the C++ compilation command line in the `Makefile`.

### 2. Compile *libStatGen* C++ library
C++ library *libStatGen* need to be compiled under your system by using the following commands:
```
cd libStatGen
make clean
make
```
After success Compiling: *libStatGen.a* will be created under `./libStatGen/`.

### 3. Compile C++ source code for *Estep_mcmc*
Compile C++ source code for the executive file that will be used to run the Estep MCMC algorithm to estimate eQTL effect sizes and the posterior probabilities to be an eQTL. Use the following commands under the `BGW-TWAS/` directory:
```
make clean
make
```
Executive file *Estep_mcmc* will be created under `./bin/` with success compilation.


## Input Files

The following information is required, stored in various text files with pre-defined formats: a gene expression file, genotype files for the training sample, a list of filenames for training genotype files, genotype files for the prediction sample, a list of filenames for prediction genotype files, and phenotype files for the prediction sample.

### 1. Gene Expression File for Training

A file containing gene expression levels for training samples, with one gene per row, and one sample per column starting from the 6th column. The first five columns are required to be gene annotation information including chromosome number, starting position, ending position, Gene ID, and Gene Name or Gene ID. An example of one gene and first 6 columns is below:

| CHROM |	GeneStart |	GeneEnd |	    TargetID     | GeneName |	 ROS20275399 |
| ----- | --------- | ------- | ---------------- | -------- | ------------ |
|  19	  |  1040101	| 1065571 |	 ENSG00000064687 |	 ABCA7  | 0.6707739044 |


### 2. VCF Genotype Files for Training

[VCF Genotype files](http://samtools.github.io/hts-specs/VCFv4.1.pdf) are required for training gene expression prediction models. The VCF genotype files should be one per genome block (variants of the same chromosome should be in the same block), sorted by position, and then zipped by `bgzip` (file names are of `[filehead].vcf.gz`), e.g., `./Example/ExampleData/genotype_data_files/Cis_geno_block.vcf.gz`.


### 3. VCF Genotype Files for GReX Prediction

[VCF Genotype files](http://samtools.github.io/hts-specs/VCFv4.1.pdf) are required for predicting GReX values of test samples, which should be one per chromosome, sorted by position, and then zipped by `bgzip`. The files should be VCF files, zipped with `bgzip`.


### 4. List of file heads of VCF genotype files
For the training procedure, users can split the genotype information into many segments based on LD information, each segment containing approximately 3,000-5,000 SNPs. To facilitate computation, users should provide a text file that is a list names for the genome segment files. These should contain a common name (e.g., the name of study or sample), plus the chromosome, start position, and end position of the genome block - all seperated by underscore `_`. E.g.:

|            #filename             |
| -------------------------------- |
| Rosmap_GWAS_19_610729_2098396    |
| Rosmap_GWAS_15_38530777_40384132 |


The prediction genotype file list should contain the names of the files for each chromosome. This is useful in the event that files include the sample or study names. E.g.:

|                #filename                |
| --------------------------------------- |
| MayoLOADGWAS_CHR1_filtered.dose.geno.gz |
| MayoLOADGWAS_CHR2_filtered.dose.geno.gz |
| MayoLOADGWAS_CHR3_filtered.dose.geno.gz |
| MayoLOADGWAS_CHR4_filtered.dose.geno.gz |


Note that the file extension suffix `.vcf.gz` should not be listed with the file names.

### 5. Phenotype File

A phenotype file is required for the prediction dataset. This file contains only two columns: the first column is a list of IDs that match the IDs in the prediction VCF files, and the second column is a  quantitative phenotype score. Column names should not be included.

For creating the $\hat{GReX}$ values, the actual phenotypes scores are not used, but this file is needed to match phenotype IDs with genotype IDs. Therefore, the second column does not have to be a meaningful variable and can be a placeholder of some kind.

### 6. Directory Paths

The entire process requires the storage of many results files for each gene (storage is optionally temporary). As the tool is only feasible in a high-performance computing environment, users should supply the following paths:

- `training_genotype` directory: a directory where all genotype segment files are stored for the training sample
- `prediction_genotype` directory: a directory where all genotype files are stored for the prediction sample (by chromosome)
- `results_directory`: a location where summary stats and training results can be stored
- `LD_directory`: a directory to store LD information corresponding to genome blocks. This step only needs to be completed once, and the LD information among the genotypes can then be used for all subsequent genes
- `Scripts_directory`: file path to the "Scripts" folder downloaded from this site, containing the software and scripts needed to implement the Bayesian $\hat{GReX}$ method


### 7. Other required input info:

- Total number of segments in the training genotype files
- Sample size of training sample
- $p$-value threshold for including genome segments in training model
- maximum number of segments to include in training model
- Genotype format (`GT` for ${0,1,2}$ Genotypes or `DS` for $[0-2]$ Dosages) for training and for prediction
- Name of `tabix` module available in your computing environment
- Number of available cores for parallel computing

## Example Usage

Set up bash variables:

```
BGW_dir=~/GIT/BGW-TWAS # tool directory
gene_name=ABCA7
GeneExpFile=${BGW_dir}/Example/ExampleData/Gene_Exp_example.txt
geno_dir=${BGW_dir}/Example/ExampleData/genotype_data_files
wkdir=${BGW_dir}/Example/ExampleWorkDir
LDdir=${BGW_dir}/Example/ExampleData/LDdir
Genome_Seg_Filehead=${BGW_dir}/Example/ExampleData/geno_block_filehead.txt
GTfield=DS # specify genotype field "GT" for genotype
num_cores=2 # number of cores to be used
```

### Step 1. Obtain Summary Statistics
Shell script `Step1_get_sum_stat.sh` will obtain single variant eQTL summary statistics (aka Score Statistics) in required formats.

#### Input arguments
- `--BGW_dir` : Specify the directory of BGW-TWAS tool
- `--wkdir` : Specify a working directory
- `-GeneExpFile` : Specify gene expression file directory
- `--gene_name` : Specify the gene name that should be the same used in `GeneExpFile`
- `--geno_dir` : Specify the directory of all genotype files
- `--LDdir` : Specify the directory of all LD files
- `--Genome_Seg_Filehead` : Specify the genome segmentation file
- `--GTfield` : Specify the genotype format in the vcf file that should be used: `GT` (default) or e.g., `DS` for dosage
- `--num_cores` : Specify the number of parallele sessions, default `1`.

#### Example command:
```
${BGW_dir}/bin/Step1_get_sumstat.sh --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} --GeneExpFile ${GeneExpFile} \
--geno_dir ${geno_dir} --LDdir ${LDdir} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--GTfield ${GTfield} --num_cores ${num_cores}
```

#### Output files
- This shell script will create a directory called `${gene}_scores/` under the specified working directory `${wkdir}/`, where results for each genome segment will be stored.

- This script will also generate LD files under `${LDdir}/` for all genome blocks. Since LD files will be the same per genome block, these files will only be generated once and used for training prediction models for all gene expression traits.

- These summary statistics files and LD files will be used for implementing the MCMC algorithm to fit the Bayesian model


### Step 2. Prune Genome Segments
Step 2 selects a subset of genome blocks (up to `${max_blocks}`) for joint model training by BGW-TWAS, where Cis blocks are always selected. Trans blocks with minimum single-variant eQTL `p-value < ${p_thresh}` will first be ranked by the smallest p-value within block (from the smallest to the largest), where top ranked trans blocks will be selected up to `${max_blocks}`.


#### Input arguments
- `--wkdir` : Specify a working directory
- `--gene_name` : Specify the gene name that should be the same used in `GeneExpFile`
- `-GeneExpFile` : Specify gene expression file directory
- `--Genome_Seg_Filehead` : Directory of the file containing a list of fileheads of segmented genotype files
- `--p_thresh` : Specify p-value threshold for pruning, default `1e-5`
- `--max_blocks` : Specify the maximum number of genome blocks for jointly training gene expression prediction models, default `100`

#### Example command:
```
${BGW_dir}/bin/Step2_prune.sh --wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --Genome_Seg_Filehead ${Genome_Seg_Filehead} \
--p_thresh 0.001 --max_blocks 100
```

#### Output files
A list of filehead of selected genome blocks is given in the file of `${wkdir}/${gene_name}_scores/${gene_name}_select_segments.txt`, which will be used in next step (model training).

### 3. Training BGW-TWAS gene expression prediction model with selected genome segments

Step 3 will use summary statistics generated from Step 1 for genome blocks pruned from Step 2 to carry out Bayesian variable selection regression with the EM-MCMC algorithm, with the goal of fitting a gene expression prediction model.

#### Example command:
```
N=499 # sample size
hfile=${BGW_dir}/Example/hypval.txt

${BGW_dir}/bin/Step3_EM-MCMC.sh  --BGW_dir ${BGW_dir} \
--wkdir ${wkdir} --gene_name ${gene_name} \
--GeneExpFile ${GeneExpFile} --LDdir ${LDdir} \
--N ${N} --hfile ${BGW_dir}/Example/hypval.txt \
--em 5 --burnin 10000 --Nmcmc 10000 \
--num_cores 2
```

#### Output files
The final iteration of the EM-MCMC algorithm will result in a file that lists all variants found to have  non-zero effect sizes (specifically, a non-zero posterior probability of inclusion, $PP_i$). The chromosome, position, rsID, reference and alternative alleles, $maf$, $cis$ or $trans$ status relative to the gene, $PP_i$, and effect size $w_i$ will all be saved in the file `${wkdir}/Eoutput/paramtemp${ems}.txt`.

[//]: <> (Numerous arguments can be used to modify the EM-MCMC algorithm in Step 3, but should be done with caution. These arguments are detailed in Yang et al. 2017 https://github.com/yjingj/bfGWAS/blob/master/bfGWAS_Manual.pdf)

### Step 4: Extract Genotypes from the Prediction Sample

In step 4, the variants found to have non-zero effect sizes are matched to the genotype information available in the prediction sample and extracted into a text file. The required arguments for this step are 1) `${gene}` 2) `${Res_dir}` 3) `${pred_geno_dir}` 4) `${pred_geno_filenames}` 5) `${pred_pheno_file}` 6) `${genotype_format}` 7) `${tabix_mod}`


```
################################################################
################################################################
### Step 4: Extract genotypes for prediction sample that
### align with eQTL with non-zero effect size from the bayesian training result
################################################################
################################################################

pred_geno_dir=/mnt/YangFSS/data/AMP-AD/Mayo/Genotype/Impute2_1KG
pred_geno_filenames=/Example/pred_geno_list.txt
pred_pheno_file=/mnt/YangFSS/data/AMP-AD/Mayo/Phenotype/MayoPhenoAD.txt
genotype_format=GT #or DS
tabix_mod=tabix/0.2.6

${Scripts_dir}/Step4.sh ${gene} ${Res_dir} ${pred_geno_dir} ${pred_geno_filenames}
${pred_pheno_file} ${genotype_format} ${tabix_mod}

```

This step will produce a file with genotype information for all of the prediction samples for those matching variants that were found to have non-zero effect size in the training step. The file is called `${gene}_grex_genotypes.geno`.

### Step 5

Step 5 uses the prediction sample genotypes and effect sizes from the training model to compute $\hat{GReX}$. The arguments required are 1) `${gene}` 2) `${Scripts_dir}` 3) `${Res_dir}`

```
################################################################
################################################################
### Step 5: Obtain predicted GREX
################################################################
################################################################

${Scripts_dir}/Step5.sh ${gene} ${Scripts_dir} ${Res_dir}
```

This script provides the two primary results of ultimate interest from the method. First, for each individual in the prediction sample, $\hat{GReX}$ for person $j$ is calculated as

$$\hat{GReX}_j = \sum_{i=1}^p x_{ij} \hat{PP}_i \hat{w}_i$$

where $x_{ij}$ is the genotype for person $j$ on variant $i$, $\hat{PP}_i$ is the estimated posterior probability for variant $i$ from the training model, and $\hat{w}_i$ is the estimated regression effect size for variant $i$ from the training model. These are stored for each individual in a file called `${Res_dir}/${gene}_GREX/${gene},_predicted_GREX.txt`.

Secondly, the sums of the posterior probabilities are also provided. These represent the expected number of eQTL for a given gene. Users are provided the total expected eQTL as well as the expected *cis*- and *trans*-eQTL. These are provided in the file `${Res_dir}/${gene}_GREX/${gene}_PPis.txt`.




