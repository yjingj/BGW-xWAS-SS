#!/usr/bin/bash

#####################################################
#####################################################
# Obtain summary xQTL data (aka single variant Zscore Statistics and LD files)
# Run single-variant test with training genotypes and molecular traits.
#####################################################
#####################################################

# Variable needed for obtaining summary statistics
###
# --BGW_dir : Specify the directory of BGW-TWAS tool
# --wkdir : Specify a working directory
# --gene_name : Specify the gene name that should be the same used in `GeneExpFile`
# --GeneInfo : Specify gene annotation or molecular trait file directory
# --geno_dir : Specify the directory of all genotype files
# --LDdir : Specify the directory of all LD files
# --Genome_Seg_Filehead : Specify the file containing the fileheads of all genome segmentations
# --GTfield : Specify the genotype format in the vcf file that should be used: "GT" (default) or e.g., "DS" for dosage
# --num_cores : Specify the number of parallele sessions
# --clean_output : Remove intermediate files

#################################
VARS=`getopt -o "" -a -l \
BGW_dir:,wkdir:,gene_name:,GeneInfo:,geno_dir:,LDdir:,Genome_Seg_Filehead:,GTfield:,num_cores:,clean_output: \
-- "$@"`

if [ $? != 0 ]
then
    echo "Please provide required input arguments. Terminating....." >&2
    exit 1
fi

eval set -- "$VARS"

while true
do
    case "$1" in
        --BGW_dir|-BGW_dir) BGW_dir=$2; shift 2;;
        --wkdir|-wkdir) wkdir=$2; shift 2;;
        --gene_name|-gene_name) gene_name=$2; shift 2;;
        --GeneInfo|-GeneInfo) GeneInfo=$2; shift 2;;
        --geno_dir|-geno_dir) geno_dir=$2; shift 2;;
        --LDdir|-LDdir) LDdir=$2; shift 2;;
        --Genome_Seg_Filehead|-Genome_Seg_Filehead) Genome_Seg_Filehead=$2; shift 2;;
        --GTfield|-GTfield) GTfield=$2; shift 2;;
        --num_cores|-num_cores) num_cores=$2; shift 2;;
        --clean_output|-clean_output) clean_output=$2; shift 2;;
        --) shift;break;;
        *) echo "Wrong input arguments!";exit 1;;
        esac
done

##########################################
# Setting Default Input Argument Values
##########################################
GTfield=${GTfield:-"GT"}
num_cores=${num_cores:-1}
clean_output=${clean_output:-1}

if [ -s ${Genome_Seg_Filehead} ]; then
    num_segments=`wc -l ${Genome_Seg_Filehead} | awk '{print $1}'`
    if [ $num_segments -gt 0 ] ; then
        echo ${gene_name} with ${num_segments} genome blocks
    else
        echo "${Genome_Seg_Filehead} has no genome blocks (one per row). Expecting at least one genome block."
        exit
    fi
else
    echo "${Genome_Seg_Filehead} is empty. Please check."
    exit
fi

echo GTfield = $GTfield , Number of cores = $num_cores

#### Create work/output directory if not existed
mkdir -p ${wkdir}
mkdir -p ${LDdir}
echo LD directory ${LDdir}

# Set directory for single variant eQTL summary statistics (Zscore statistics)
Zscore_dir=${wkdir}/${gene_name}_Zscores
mkdir -p ${Zscore_dir}
echo Summary xQTL Zscore statistic directory ${Zscore_dir}
cd ${wkdir}
# echo ${wkdir}

# GeneInfo file columns are: Chr, start, end, gene ID, gene name, molecular trait one column per sample

# Creates a trait (pheno) file for target gene that includes sample IDs in the first column and molecular traits in the second column.
if [ -s ${GeneInfo} ] ; then
    head -1 ${GeneInfo} | awk '{$1=$2=$3=$4=$5=""; print substr($0,6)}' | tr ' ' '\n' > temp_ID.txt

    awk -F'[\t]' -v gene=${gene_name} '$5==gene{$1=$2=$3=$4=$5=""; print substr($0,6); exit; }' ${GeneInfo} | tr ' ' '\n' > trait_temp.txt

    paste temp_ID.txt trait_temp.txt > ${wkdir}/${gene_name}_trait.txt
    sed -i "s/\r//g"  ${wkdir}/${gene_name}_trait.txt

    gene_trait=${wkdir}/${gene_name}_trait.txt

    gene_info=$(awk -F'\t' -v gene=${gene_name} '$5==gene{print ; exit; }' ${GeneInfo})
    target_chr=$(echo ${gene_info} | awk 'FS {print $1}');
    start_pos=$(echo ${gene_info} | awk 'FS {print $2}');
    end_pos=$(echo ${gene_info} | awk 'FS {print $3}');
    echo Gene $gene_name with start position $start_pos and end position $end_pos
else
    echo ${GeneInfo} is empty. Please provide a valid gene information file.
fi

rm -f temp_ID.txt trait_temp.txt

## Run in parallele with specified number of processes by -P
seq 1 ${num_segments}  | xargs -I % -n 1 -P ${num_cores} sh ${BGW_dir}/bin/get_xqtl_sumstat.sh ${gene_trait} ${geno_dir} ${Zscore_dir} ${BGW_dir} ${LDdir} ${Genome_Seg_Filehead} % ${GTfield} ${target_chr} ${start_pos} ${end_pos}

if [ $clean_output -eq 1 ] ; then
rm -fr ${Zscore_dir}/output
fi

echo Complete generating xQTL summary statistics under ${Zscore_dir} and ${LDdir} !

exit
