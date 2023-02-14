#!/bin/sh

gene_exp_trait=$1
geno_dir=$2
ZScore_dir=$3
BGW_dir=$4
LDdir=$5
filehead=$6
block=$7
GTfield=$8

# LD LDwindow
LDwindow=1000000
cd ${ZScore_dir}

line=$(head -n $block $filehead | tail -n1)
echo Generate eQTL summary data for block ${line}

if [ -s ${LDdir}/${line}.LDcorr.txt.gz ] ; then
echo ${LDdir}/${line}.LDcorr.txt.gz exists!
LDwindow=1
fi
echo LDwindow is set as  $LDwindow

if [ ! -s ${gene_exp_trait} ] ; then
echo ${gene_exp_trait} is empty. Please check.
exit
fi

if [ -s ${geno_dir}/${line}.vcf.gz ] ; then
	### With input genotype file in VCF format
	${BGW_dir}/bin/Estep_mcmc -vcf ${geno_dir}/${line}.vcf.gz -p ${gene_exp_trait} -maf 0.01 -o ${line} -LDwindow ${LDwindow} -GTfield ${GTfield} -saveSS -zipSS
	mv -f ${ZScore_dir}/output/${line}.Zscore.txt.gz ${ZScore_dir}/
	mv -f ${ZScore_dir}/output/${line}.Zscore.txt.gz.tbi ${ZScore_dir}/
	echo ZScore statistics file ${ZScore_dir}/${line}.Zscore.txt.gz were generated.

	if [ ! -s ${LDdir}/${line}.LDcorr.txt.gz ] ; then
		mv -f ${ZScore_dir}/output/${line}.LDcorr.txt.gz ${LDdir}/
		mv -f ${ZScore_dir}/output/${line}.LDcorr.txt.gz.tbi ${LDdir}/
		echo LD file ${LDdir}/${line}.LDcorr.txt.gz were generated !
	fi
else
	echo Training genotype VCF file ${geno_dir}/${line}.vcf.gz is empty. Please check.
fi

# rm -rf ${ZScore_dir}/output

