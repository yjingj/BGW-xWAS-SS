#!/usr/bin/bash

# Tool Tabix is required for running this script
# R library data.table and tidyverse are required

#################################
VARS=`getopt -o "" -a -l \
BGW_dir:,wkdir:,gene_name:,BGW_weight:,test_geno_dir:,test_geno_filehead:,GTfield:,test_pheno:,quant_pheno:,num_cores:,clean_output: \
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
		--BGW_weight|-BGW_weight) BGW_weight=$2; shift 2;;
        --test_geno_dir|-test_geno_dir) test_geno_dir=$2; shift 2;;
		--test_geno_filehead|-test_geno_filehead) test_geno_filehead=$2; shift 2;;
		--GTfield|-GTfield) GTfield=$2; shift 2;;
		--test_pheno|-test_pheno) test_pheno=$2; shift 2;;
        --quant_pheno|-quant_pheno) quant_pheno=$2; shift 2;;
        --num_cores|-num_cores) num_cores=$2; shift 2;;
        --clean_output|-clean_output) clean_output=$2; shift 2;;
        --) shift;break;;
        *) echo "Wrong input arguments!";exit 1;;
        esac
done

##########################################
# Setting Default Input Argument Values
##########################################
GTfield=${Nmcmc:-"GT"}
num_cores=${num_cores:-1}
quant_pheno=${quant_pheno:-"TRUE"}
clean_output=${clean_output:-1}

mkdir -p ${wkdir}/${gene_name}_predict_trait
cd ${wkdir}/${gene_name}_predict_trait

## Check if $BGW_weight is generated and there is valid eQTL
if [ ! -s ${BGW_weight} ] ; then
    echo ${BGW_weight} is not generated.
    exit
else
    n_eqtl=`grep -v "#" ${BGW_weight} | wc -l  | awk '{print $1}'`
    if [ $n_eqtl -gt 0 ] ; then
        echo Extract genotype vcf file for ${n_eqtl} xQTL in ${BGW_weight} ...
    else
        echo There is no xQTL with non-zero xQTL weights in ${BGW_weight}.
        exit
    fi
fi

## make header row for vcf file from one of the existing genotype files
file_temp=$(head -n1 ${test_geno_filehead})
tabix ${test_geno_dir}/${file_temp}.vcf.gz -H | tail -n1 > ${wkdir}/${gene_name}_predict_trait/${gene_name}_test.vcf

# Loop through SNP ID and use tabix
tail -n+2 ${BGW_weight} | while read line ; do
chr=$(echo $line | awk -F'[\t ]' '{print $1}' )
pos=$(echo $line | awk -F'[\t ]' '{print $2}' )
ref=$(echo $line | awk -F'[\t ]' '{print $4}' )
alt=$(echo $line | awk -F'[\t ]' '{print $5}' )
#echo ${chr}:${pos}-${pos} $ref $alt

n_file=$(grep _CHR${chr}_ $test_geno_filehead | wc -l | awk '{print $1}')

if [ $n_file -eq 1 ] ; then
#echo ${file_temp}.vcf.gz
file_temp=$(grep _CHR${chr}_ $test_geno_filehead)
tabix ${test_geno_dir}/${file_temp}.vcf.gz ${chr}:${pos}-${pos} | awk -v ref=${ref} -v alt=${alt} -F"\t" '($4==ref && $5==alt) || ($4==alt && $5==ref) {print }' >> ${wkdir}/${gene_name}_predict_trait/${gene_name}_test.vcf
else
    echo There is no file head or multiple file heads with the pattern of _CHR${chr}_ in $test_geno_filehead.
    echo $file_temp
    echo Please combine all genotype in the same chromorome into one VCF file with the pattern of \"_CHR${chr}_\" in the filename and include the file head in $test_geno_filehead.
fi
done

echo ${wkdir}/${gene_name}_predict_trait/${gene_name}_test.vcf is created.

## Generate genotype file (dosages) with ${gene_name}_test.vcf
nsnp=`grep -v "#" ${wkdir}/${gene_name}_predict_trait/${gene_name}_test.vcf | wc -l  | awk '{print $1}'`
if [ $nsnp -gt 0 ] ; then
    ${BGW_dir}/bin/Estep_mcmc -vcf ${wkdir}/${gene_name}_predict_trait/${gene_name}_test.vcf -p ${test_pheno} -o ${gene_name}_pred -GTfield ${GTfield} -saveGeno -maf 0
    if [ -s ./output/${gene_name}_pred.geno.txt ] ; then
        sed -i 's/#//g' ./output/${gene_name}_pred.geno.txt
        mv -f ./output/${gene_name}_pred.geno.txt ${wkdir}/
        echo Converting test genotype VCF file to genotype dosage file ${wkdir}/${gene_name}_pred.geno.txt is success.
    else
        echo ${wkdir}/${gene_name}_predict_trait/output/${gene_name}_pred.geno.txt failed to be generated.
        exit
    fi
else
    echo There is no test SNPs with non-zero xQTL weights for gene $gene_name. Please check if you BGW weight file and test VCF files.
    exit
fi

## Remove the "#" from the header row
sed -i 's/#//' $BGW_weight
if [ -s ${wkdir}/${gene_name}_pred.geno.txt ] ; then
    echo Calculating predicted genetically regulated molecular trait ...
    Rscript ${BGW_dir}/bin/predict_molecular_trait.R ${gene_name} ${wkdir} ${test_pheno} ${quant_pheno}
    echo Predicted genetically regulated molecular trait file is generated under ${wkdir}
else
	echo Test genotype dosage file ${wkdir}/${gene_name}_pred.geno.txt file failed to be generated. Please check your input arguments.
	exit
fi

if [ $clean_output -eq 1 ] ; then
    rm -rf ${wkdir}/${gene_name}_predict_trait/
    rm -f ${wkdir}/ABCA7_pred.geno.txt
fi

exit

