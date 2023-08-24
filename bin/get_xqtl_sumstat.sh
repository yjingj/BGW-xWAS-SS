#!/bin/sh

gene_trait=$1
geno_dir=$2
Zscore_dir=$3
BGW_dir=$4
LDdir=$5
filehead=$6
block=$7
GTfield=$8
target_chr=$9
start_pos=${10}
end_pos=${11}

# LD LDwindow
LDwindow=1000000
cd ${Zscore_dir}

line=$(head -n $block $filehead | tail -n1)
echo Generate xQTL summary data for block ${line}

if [ -s ${LDdir}/${line}.LDcorr.txt.gz ] ; then
echo ${LDdir}/${line}.LDcorr.txt.gz exists!
LDwindow=1
fi
echo LDwindow is set as  $LDwindow

if [ ! -s ${gene_trait} ] ; then
echo ${gene_trait} is empty. Please check.
exit
fi


if [ -s ${geno_dir}/${line}.vcf.gz ] ; then
	### With input genotype file in VCF format
	${BGW_dir}/bin/Estep_mcmc -vcf ${geno_dir}/${line}.vcf.gz -p ${gene_trait} -maf 0.01 -o ${line} -LDwindow ${LDwindow} -GTfield ${GTfield} -target_chr ${target_chr} -start_pos ${start_pos} -end_pos ${end_pos} -saveSS -zipSS
	mv -f ${Zscore_dir}/output/${line}.Zscore.txt.gz ${Zscore_dir}/
	mv -f ${Zscore_dir}/output/${line}.Zscore.txt.gz.tbi ${Zscore_dir}/
	echo Zscore statistics file ${Zscore_dir}/${line}.Zscore.txt.gz were generated.

	if [ $LDwindow -ne 1 ] ; then
		mv -f ${Zscore_dir}/output/${line}.LDcorr.txt.gz ${LDdir}/
		mv -f ${Zscore_dir}/output/${line}.LDcorr.txt.gz.tbi ${LDdir}/
		echo LD file ${LDdir}/${line}.LDcorr.txt.gz were generated !
	fi
else
	echo Training genotype VCF file ${geno_dir}/${line}.vcf.gz is empty. Please check.
fi

# rm -rf ${Zscore_dir}/output

