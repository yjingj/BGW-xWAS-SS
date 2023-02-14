#!/usr/bin/bash

################################################################
################################################################
# Step 3: Training BGW-TWAS gene expression prediction model by EM-MCMC algorithm
# Employ Makefile generated by a Perl script
################################################################
################################################################

#################################
VARS=`getopt -o "" -a -l \
BGW_dir:,wkdir:,gene_name:,GeneExpFile:,select_filehead:,LDdir:,ZScore_dir:,N:,hfile:,em:,burnin:,Nmcmc:,PCP_thresh:,num_cores:,clean_output: \
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
        --GeneExpFile|-GeneExpFile) GeneExpFile=$2; shift 2;;
        --select_filehead|-select_filehead) select_filehead=$2; shift 2;;
        --LDdir|-LDdir) LDdir=$2; shift 2;;
        --ZScore_dir|-ZScore_dir) ZScore_dir=$2; shift 2;;
        --N|-N) N=$2; shift 2;;
		--hfile|-hfile) hfile=$2; shift 2;;
		--em|-em) em=$2; shift 2;;
		--burnin|-burnin) burnin=$2; shift 2;;
		--Nmcmc|-Nmcmc) Nmcmc=$2; shift 2;;
        --PCP_thresh|-PCP_thresh) PCP_thresh=$2; shift 2;;
        --num_cores|-num_cores) num_cores=$2; shift 2;;
        --clean_output|-clean_output) clean_output=$2; shift 2;;
        --) shift;break;;
        *) echo "Wrong input arguments!";exit 1;;
        esac
done

##########################################
# Setting Default Input Argument Values
##########################################
em=${em:-2}
burnin=${burnin:-10000}
Nmcmc=${Nmcmc:-10000}
PCP_thresh=${PCP_thresh:-0.0001}
num_cores=${num_cores:-1}
clean_output=${clean_output:-1}

### Extract gene info from ${GeneExpFile}
gene_info=$(awk -F'\t' -v gene=${gene_name} '$5==gene{print ; exit; }'  ${GeneExpFile})
target_chr=$(echo ${gene_info} | awk 'FS {print $1}');
start_pos=$(echo ${gene_info} | awk 'FS {print $2}');
end_pos=$(echo ${gene_info} | awk 'FS {print $3}');
echo Study gene ${gene_name} from CHR $target_chr with cis-region ranging from $start_pos to $end_pos...

### Phenotype variance
# pv=`awk 'NR==1{print $2}' ${wkdir}/${gene_name}_geneExp_var.txt`

### Fileheads of segmented genome blocks
if [ -z ${select_filehead} ] ; then
	echo "Selected filehead txt file is not provided."
    select_filehead=${wkdir}/${gene_name}_select_filehead.txt
    echo Default Selected filehead txt file : ${select_filehead}
fi

### Specify output make file directory and generate makefile
### Set up file directories
mkdir -p ${wkdir}/${gene_name}_EM_MCMC
cd ${wkdir}/${gene_name}_EM_MCMC

if [ -z ${ZScore_dir} ] ; then
    ZScore_dir=${wkdir}/${gene_name}_Zscores
    echo Default summary Zscore statistics file directory: $ZScore_dir
else
    echo Summary ZScore statistics file directory is provided as ${ZScore_dir} ;
fi

mkfile="${wkdir}/${gene_name}_EM_MCMC/${gene_name}_BGW.mk"

${BGW_dir}/bin/gen_mkf.pl \
--BGW_dir ${BGW_dir} --hyp ${hfile} \
-n ${N} -w ${wkdir}/${gene_name}_EM_MCMC --geno sumstat \
-f ${select_filehead} -l local \
--targ ${target_chr} --start ${start_pos} --end ${end_pos} \
--LDdir ${LDdir} --ZScoredir ${ZScore_dir} \
-j ${gene_name}_BGW --em ${em} -b ${burnin} -N ${Nmcmc} \
--mf ${mkfile}

### Run the makefile (./BGW-TWAS/bin/run_Estep.sh and ./BGW-TWAS/bin/run_Mstep.sh are called by make)
make -f ${mkfile} clean
echo Run make with $mkfile and j=$num_cores parallel jobs
make -k -C ${wkdir}/${gene_name}_EM_MCMC -f ${mkfile} -j ${num_cores} > ${wkdir}/${gene_name}_EM_MCMC/make.output 2> ${wkdir}/${gene_name}_EM_MCMC/make.err

if [ -s ${wkdir}/${gene_name}_EM_MCMC/Eoutput/paramtemp${em}.txt.gz ] ; then
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tMAF\tTrans\tPCP\tbeta\tSE_beta\tLRT\tpval_LRT\tRank" > ${wkdir}/${gene_name}_BGW_eQTL_weights.txt
    zcat ${wkdir}/${gene_name}_EM_MCMC/Eoutput/paramtemp${em}.txt.gz | awk -v PCP_thresh=${PCP_thresh} '$8>PCP_thresh{print }'  | sort -nk1,2 >> ${wkdir}/${gene_name}_BGW_eQTL_weights.txt
else
    echo File ${wkdir}/${gene_name}_EM_MCMC/Eoutput/paramtemp${em}.txt is not generated by EM-MCMC algorithm. Please Check Step 3.
    exit
fi

if [ $clean_output -eq 1  ]; then
    rm -rf ${wkdir}/${gene_name}_EM_MCMC
    rm -f ${wkdir}/${gene_name}_exp_trait.txt
fi

exit
