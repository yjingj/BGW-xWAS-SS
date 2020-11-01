#!/bin/sh

#  get_sum_stat.sh
#  
#
#

gene_exp_trait=$1
geno_dir=$2
Score_dir=$3
BGW_dir=$4
LDdir=$5
filehead=$6
block=$7
GTfield=$8

# Score statistics directory
LDwindow=1000000
cd ${Score_dir}

echo $block

line=$(head -n $block $filehead | tail -n1)


if [ -f ${LDdir}/${line}.LDcorr.txt.gz ] ; then
#echo ${LDdir}/${line}.LDcorr.txt.gz exists!
LDwindow=1
fi

#echo LDwindow is $LDwindow

### With input genotype file in VCF format
${BGW_dir}/bin/Estep_mcmc -vcf ${geno_dir}/${line}.vcf.gz -p ${gene_exp_trait} -maf 0 -o ${line} -LDwindow ${LDwindow} -GTfield ${GTfield} -saveSS -zipSS

mv -f ${Score_dir}/output/${line}.score.txt.gz ${Score_dir}/
mv -f ${Score_dir}/output/${line}.score.txt.gz.tbi ${Score_dir}/

if [ ! -f ${LDdir}/${line}.LDcorr.txt.gz ] ; then
	mv -f ${Score_dir}/output/${line}.LDcorr.txt.gz ${LDdir}/
	mv -f ${Score_dir}/output/${line}.LDcorr.txt.gz.tbi ${LDdir}/
fi

# rm -rf ${Score_dir}/output

