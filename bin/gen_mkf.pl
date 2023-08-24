#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

Generate bfGWAS Makefile

=head1 SYNOPSIS

./gen_mkf.pl [options] 

Options:
  -h      help messages
  -d      degub
  --mf    output make file

=head1 OPTIONS

=over 8

=item B<-help>

  --wkdir         work directory : location for all output files
  --BGW_dir         directory for the C++ executable file for Estep (running MCMC)
  --LDdir       directory of LD correlation files
  --Zscore_dir    directory of GWAS single variant Zscore statistic files
  --filehead         text file with a list of fileheads for genome blocks
  --hfile      initial hyper parameter value file
  --Nsample       sample size
  --maf      maf threshold: default 1%
  --targ    target chromosome number
  --start   target gene transcription start position
  --end     target gene transcription end position
  --Nburnin         number of burn ins: default 10,000
  --Nmcmc         number of MCMC iterations: default 10,000
  --win      window size of the neighborhood: default 100
  --smin     minimum number of variates per block in the model: default 0
  --smax     maximum number of variates per block in the model: default 10
  --em       number of EM iterations: default 2
  --mf       output make file

=back

=head1 DESCRIPTION

B<gen_mkf.pl> will generate a makefile for conducting Bayesian Functional GWASs by bfGWAS. 

=cut

# define default option variables
my $help;
my $verbose;
my $debug;
my $man;
my $launchMethod = "local";
my $makeFile = "BGW.mk";

my $BGW_dir="";
my $wkdir=getcwd();
my $filehead = "";
my $Zscore_dir = "";
my $LDdir="";
my $hfile="";
my $Nsample="1000";

my $EM=2;
my $maf="0.01";
my $target_chr;
my $start_pos;
my $end_pos;
my $smin="0";
my $smax="10";
my $win="100";
my $burnin="10000";
my $Nmcmc="10000";
my $NmcmcLast="10000";
my $pp_cis="1e-5"; # prior PCP for cis SNPs
my $pp_trans="1e-6"; # prior PCP for trans SNPs
my $a_gamma="1"; # prior PCP for cis SNPs
my $b_gamma="2"; # prior PCP for trans SNPs

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug, 'm'=>\$man,
                'wkdir:s'=>\$wkdir, 'BGW_dir:s' =>\$BGW_dir,
                'LDdir:s'=>\$LDdir, 'Zscore_dir:s'=>\$Zscore_dir,
                'filehead:s'=>\$filehead,
                'hfile:s'=>\$hfile, 'Nsample:s'=>\$Nsample, 'maf:s'=>\$maf,
                'targ:s'=>\$target_chr, 'start:s'=>\$start_pos, 'end:s'=>\$end_pos,
                'Nburnin:s'=>\$burnin, 'Nmcmc:s'=>\$Nmcmc, 'NmcmcLast:s'=>\$NmcmcLast,
                'pp_cis:s'=>\$pp_cis, 'pp_trans:s'=>\$pp_trans,
                'a_gamma:s'=>\$a_gamma, 'b_gamma:s'=>\$b_gamma,
                'win:s'=>\$win, 'smin:s'=>\$smin, 'smax:s'=>\$smax,
                'em:i'=>\$EM, 'mf:s'=>\$makeFile)
  || !defined($wkdir) || scalar(@ARGV)!=0)
{
    if ($help)
    {
        pod2usage(1);
        exit(0);
    }
    elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }
    else
    {
        pod2usage(1);
        exit(0);
    }
}

if ($help)
    {
        pod2usage(1);
        exit(0);
    }
elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }


my $toolE="${BGW_dir}/bin/Estep_mcmc";
my $rs="${BGW_dir}/bin/Mstep.r";

##############
#print options
##############
printf("Options\n");
printf("\n");
# printf("launch method : %s\n", $launchMethod);
printf("work directory : %s\n", $wkdir);
print "Estep: ", $toolE, "\n", "Rscript: ", $rs, "\n"; 
print "Zscore_dir: ", $Zscore_dir, "\n", "LDdir: ", $LDdir, "\n",
        "hfile: ", $hfile, "\nfilehead text file: ", $filehead, "\n",
        "Nsample ", $Nsample, "; maf ", $maf, "; smin ", $smin, "\n",
        "smax ", $smax, "; win ", $win, "\n", 
        #"a_gamma ", $a_gamma, "; b_gamma ", $b_gamma,"\n",
        "burnin ", $burnin, "; Nmcmc ", $Nmcmc, "\n";
printf("\n");


#arrays for storing targets, dependencies and commands
my @tgts = ();
my @deps = ();
my @cmds = ();

#temporary variables
my $tgt;
my $dep;
my @cmd;

mkpath($wkdir);

########################################
# Initial Set up before EM iterations
########################################

### prepare files before EM_MCMC
my $hypcurrent="$wkdir/hypval.current";
$tgt = "$wkdir/pre_em.OK";
$dep = "";
@cmd = "rm -f -r $wkdir/output $wkdir/Eoutput $wkdir/OUT";
push(@cmd, "mkdir -p $wkdir/output $wkdir/Eoutput $wkdir/OUT");
push(@cmd, "cp -f $hfile $hypcurrent");
push(@cmd, "> $wkdir/Eoutput/EM\_result.txt");
push(@cmd, "> $wkdir/Rout.txt");
makeJob($tgt, $dep, @cmd);  


###### EM step 0 without dependencies ###########
my $i;
my $premcmcOK="";
my @filehead;
my $line;

open(my $FILELIST, $filehead)
    or die "Can not open $filehead \!";
 while ($line = <$FILELIST>) {
    chop $line;
    push(@filehead, $line);
}
close $FILELIST;

if(@filehead == 0) {
    print STDERR "file list is empty\! \n";
    exit(1);
} 
else{ print "Total \# of genome blocks: ", scalar(@filehead), "\n \n"; }


for(my $j=0; $j< @filehead; ++$j)
    {
        $i=0;
        $line = $filehead[$j];
        $premcmcOK .= "$wkdir/OUT/$line.$i.OK ";
        $tgt = "$wkdir/OUT/$line.$i.OK";
        $dep = "$wkdir/pre_em.OK";

        @cmd = "$toolE -inputSS -Zscore ${Zscore_dir}/${line}.Zscore.txt.gz -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -hfile ${hypcurrent} -maf ${maf} -n ${Nsample} -target_chr ${target_chr} -start_pos ${start_pos} -end_pos ${end_pos} -bvsrm -smin ${smin} -smax ${smax} -win ${win} -o ${line} -w ${burnin} -s ${Nmcmc} -seed 2022 > ${wkdir}/OUT/${line}.output.txt";
        makeJob($tgt, $dep, @cmd);
    }

my $paramfile="$wkdir/Eoutput/paramtemp$i.txt";
my $hypfile="$wkdir/Eoutput/hyptemp$i.txt";
my $logfile="$wkdir/Eoutput/log$i.txt";

$tgt = "$wkdir/Eoutput/cp_param$i.OK";
$dep = "$premcmcOK";
  @cmd = "cat \`ls -d -1 $wkdir/output/**.paramtemp | grep paramtemp | head -n1 \` | head -n1 > $paramfile";
  push(@cmd, "cat \`ls -d -1 $wkdir/output/** | grep hyptemp | head -n1 \` | head -n1 > $hypfile");
  push(@cmd, "cat \`ls -d -1 $wkdir/output/**.paramtemp | grep paramtemp | sort  \` |  grep -v \"#\" | sort -nk1 -nk2 >> $paramfile");
  push(@cmd, "cat \`ls -d -1 $wkdir/output/** | grep hyptemp | sort \`| grep -v \"#\"  >> $hypfile");
  push(@cmd, "cat \`ls -d -1 $wkdir/output/** | grep log | sort \` > $logfile");
  push(@cmd, "bgzip -f $paramfile");
  push(@cmd, "tabix -p vcf -f $paramfile.gz");
makeJob($tgt, $dep, @cmd);

$tgt = "$wkdir/R$i.OK";
$dep = "$wkdir/Eoutput/cp_param$i.OK $wkdir/pre_em.OK";
@cmd = "Rscript --vanilla ${rs} $hypfile $i $pp_cis $pp_trans $a_gamma $b_gamma $Nsample $hypcurrent $wkdir/Eoutput/EM_result.txt >> $wkdir/Rout.txt 2>&1";
makeJob($tgt, $dep, @cmd);


####### With dependencies of previous output 
my $ipre="";

for $i (1..$EM){

    $ipre=$i-1; $premcmcOK="";

    for(my $j=0; $j< @filehead; ++$j){
        $line=$filehead[$j];
        $premcmcOK .= "$wkdir/OUT/$line.$i.OK ";
        $tgt = "$wkdir/OUT/$line.$i.OK";
        $dep = "$wkdir/R$ipre.OK";
        if($i < $EM){
          @cmd = "$toolE -inputSS -Zscore ${Zscore_dir}/${line}.Zscore.txt.gz -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -hfile ${hypcurrent} -n ${Nsample} -maf ${maf} -target_chr ${target_chr} -start_pos ${start_pos} -end_pos ${end_pos} -bvsrm -smin $smin -smax $smax -win $win -o ${line} -w ${burnin} -s ${Nmcmc} -seed 2022 > ${wkdir}/OUT/${line}.output.txt" ;
        } elsif ($i == $EM){
            @cmd = "$toolE -inputSS -Zscore ${Zscore_dir}/${line}.Zscore.txt.gz -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -hfile ${hypcurrent} -n ${Nsample} -maf ${maf} -target_chr ${target_chr} -start_pos ${start_pos} -end_pos ${end_pos} -bvsrm -smin $smin -smax $smax -win $win -o ${line} -w ${burnin} -s ${NmcmcLast} -seed 2022 > ${wkdir}/OUT/${line}.output.txt";
      }
        makeJob($tgt, $dep, @cmd);
    }

    $paramfile="$wkdir/Eoutput/paramtemp$i.txt";
    $hypfile="$wkdir/Eoutput/hyptemp$i.txt";
    $logfile="$wkdir/Eoutput/log$i.txt";

  $tgt = "$wkdir/Eoutput/cp_param$i.OK";
  $dep = "$premcmcOK";
  @cmd = "cat \`ls -d -1 $wkdir/output/**.paramtemp | grep paramtemp | head -n1 \` | head -n1 > $paramfile";
  push(@cmd, "cat \`ls -d -1 $wkdir/output/** | grep hyptemp | head -n1 \` | head -n1 > $hypfile");
  push(@cmd, "cat \`ls -d -1 $wkdir/output/**.paramtemp | grep paramtemp | sort \` |  grep -v \"#\" | sort -nk1 -nk2  >> $paramfile");
  push(@cmd, "cat \`ls -d -1 $wkdir/output/** | grep hyptemp | sort\`| grep -v \"#\"  >> $hypfile");
  push(@cmd, "cat \`ls -d -1 $wkdir/output/** | grep log | sort\` > $logfile");
  push(@cmd, "bgzip -f $paramfile");
  push(@cmd, "tabix -p vcf -f $paramfile.gz");
  makeJob($tgt, $dep, @cmd);

  $tgt = "$wkdir/R$i.OK";
  $dep = "$wkdir/Eoutput/cp_param$i.OK";
  @cmd = "Rscript --vanilla $rs $hypfile $i $pp_cis $pp_trans $a_gamma $b_gamma $Nsample $hypcurrent $wkdir/Eoutput/EM_result.txt >> $wkdir/Rout.txt 2>&1";
  makeJob($tgt, $dep, @cmd);

}

# $tgt = "$wkdir/clean_output_files.OK";
# $dep = "$wkdir/Eoutput/cp_param$EM.OK";
# @cmd = "rm -f $wkdir/output/** $wkdir/OUT/*.output.txt";
# makeJob($tgt, $dep, @cmd);


#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
#this tells GNU Make to delete files that were generated during a command that did not end successfully.
print MAK ".DELETE_ON_ERROR:\n\n";
#this ensures that all the targets are executed; exclude clean
print MAK "all: @tgts\n\n";

######## Create clean jobs command #######
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $wkdir/*.OK $wkdir/Eoutput/*.OK $wkdir/OUT/*.OK $wkdir/slurm_err/*.err");

for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#functions
##########

#run a job locally
sub makeJob
{    
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}


# Check empty directories
sub dir_is_empty
{
    my ($path) = @_;
    opendir DIR, $path;
    while(my $entry = readdir DIR) {
        next if($entry =~ /^\.\.?$/);
        closedir DIR;
        return 0;
    }
    closedir DIR;
    return 1;
}
