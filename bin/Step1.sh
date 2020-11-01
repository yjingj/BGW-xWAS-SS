#!/bin/sh

################################################################
################################################################
# Step 1: obtain summary statistics (aka Score Statistics)
# Runs single-variant GWAS on the training sample genotypes and available expression data.
# first extracts gene location information for target gene from geneFile.
################################################################
################################################################

gene=$1
GeneExpFile=$2
geno_dir=$3
BGW_dir=$4
wkdir=$5
LDdir=$6
Genome_Seg_File=$7
num_segments=$8
num_cores=$9
GTfield=${10}

echo ${gene}

mkdir -p ${wkdir}
mkdir -p ${LDdir}

# Set directory for single variant eQTL summary statistics (score statistics)
mkdir -p ${wkdir}/${gene}_scores
Score_dir=${wkdir}/${gene}_scores

cd ${wkdir}
echo ${wkdir}

# GeneExpFile columns are Gene Name, Chr, Pos, start, end, expr_data

# the following creates a phenotype file for target gene expression trait that includes subject IDs in the first column and gene expression levels in the second column.
head -1 ${GeneExpFile} | awk '{$1=$2=$3=$4=$5=""; print substr($0,6)}' | tr ' ' '\n' > temp_ID.txt
grep ${gene} ${GeneExpFile} | awk '{$1=$2=$3=$4=$5=""; print substr($0,6)}' | tr ' ' '\n' > exp_temp.txt
paste temp_ID.txt exp_temp.txt > ${wkdir}/${gene}_exp_trait.txt

# calculate variance of the gene expression trait
pv=$(awk '{delta=$2; sum+=$2; ssq+=(delta - avg/NR)^2} END {print ssq/(NR-1)}' ${wkdir}/${gene}_exp_trait.txt)
echo quantitative gene expression trait variance = $pv
echo -e ${gene} '\t' ${pv} > ${wkdir}/${gene}_geneExp_var.txt
rm -f temp_ID.txt exp_temp.txt

pheno=${wkdir}/${gene}_exp_trait.txt

seq 1 ${num_segments} | while read block ; do
echo $block
#  | xargs -I{block} -n 1 -P ${num_cores} sh 
${BGW_dir}/bin/get_score_stat.sh ${pheno} ${geno_dir} ${Score_dir} ${BGW_dir} ${LDdir} ${Genome_Seg_File} ${block} ${GTfield}
done

echo Step 1 complete for generating eQTL summary statistics!

exit
