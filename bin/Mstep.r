Sys.setlocale('LC_ALL', 'C')
options(stringsAsFactors=F)

######## Need to pass args(hypfile, paramfile, k, hypcurrent_file) from bash
args <- commandArgs(TRUE)
# print(args)


hypfile=args[[1]]     # hyptemp file
k=as.numeric(args[[2]]) # EM iteration number
pp = as.numeric(args[[3]]) # for setting beta prior of CPP
a_gamma = as.numeric(args[[4]]) # Setting IG prior of sigma2
b_gamma = as.numeric(args[[5]])
n_sample = as.numeric(args[[6]]) # sample size
hypcurrent_file = args[[7]]    # hypval.current file
EM_result_file = args[[8]] # Save EM_result_file


print(c("a_gamma=", a_gamma, "b_gamma = ", b_bgamma, "sample size = ", n_sample))


###### Define functions to be used

## Calculate pi_hat and pi_se
CI_fish_pi <- function(sum_gamma, p, a, b){
  # p : number of cis or trans variants
  pi_hat = (sum_gamma + a - 1.0) / (p + a + b - 2.0)
  if(pi_hat <= 0 || pi_hat > 1){
    pi_hat = a/(a+b) # set as prior value
    se_pi = 0
  }else{
    pi_var = pi_hat * (1-pi_hat) / (p + a + b - 2.0)
    if(pi_var > 0){
      se_pi = sqrt(pi_var)
    }else{se_pi=0}
  }
  return(c(pi_hat, se_pi))
}

## Calculate sigma2_hat and se_sigma2
CI_fish_sigma2 <- function(sum_Ebeta2, sum_gamma, n_sample, a, b){
  sigma2_hat = (sum_Ebeta2 * n_sample + 2 * b) / (sum_gamma + 2 * (a + 1))
  temp_var = sum_Ebeta2 * n_sample / sigma2_hat - sum_gamma/2 - (a+1) + 2*b/sigma2_hat
  if( temp_var < 0){
    se_sigma2=0
  }else{
    se_sigma2 = sigma2_hat / sqrt(temp_var)
  }
  return(c(sigma2_hat, se_sigma2))
}

## log prior functions
logprior_sigma <- function(a, b, x){ return(-(1+a) * log(x) - b/x) }
logprior_pi <- function(a, b, x){ return((a-1) * log(x) + (b-1) * log(1 - x)) }


ptm <- proc.time()
########### Load hypfile ....
# hypcurrent_file="/home/jyang51/jyang/GITHUB/BGW-TWAS-SS/Example/ExampleWorkDir/ABCA7_EM_MCMC/hypval.current"
# hypfile="/home/jyang51/jyang/GITHUB/BGW-TWAS-SS/Example/ExampleWorkDir/ABCA7_EM_MCMC/Eoutput/hyptemp0.txt"
# k=0; pp = 1e-5; a_gamma=2; b_gamma=1; n_sample = 499

hypdata = read.table(hypfile, sep="\t", header=FALSE)
n_type = 2
print(paste(" Total Annotation categories : ", n_type))

temp_col_names <- c("genome_block_prefix", "log_post_likelihood", "r2")
for(i in 1:n_type){
  temp_col_names <- c(temp_col_names,
                      paste(c("mFunc", "sum_gamma", "sum_Ebeta2",), (i-1), sep = "_"))
}
colnames(hypdata) <-  temp_col_names

########### Update hyper parameter values
prehyp <- read.table(hypcurrent_file, header=TRUE)
print("hyper parameter values before MCMC: ")
print(prehyp)

######### Set hierarchical parameter values
n_vec = rep(0, n_type) # number of variants per category
for(i in 1:n_type){
  n_vec[i] <- hypdata[, paste("mFunc", (i-1), sep="_")]
}

#### updating hyper pi and sigma2 values for each group
# a_beta, b_beta will be set for cis- and trans- annotation
hypcurrent <- NULL
hypmat <- NULL

for(i in 1:n_type){
  # print(i)
  if(n_vec[i] > 0){
    a_beta = 2 * n_vec[i] * pp; b_beta = n_vec[i] - a_beta;
  }else{a_beta=1; b_beta = 1e5 - 1;}

  sum_gamma = hypdata[, paste("sum_gamma", (i-1), sep="_")]
  pi_temp = CI_fish_pi(sum_gamma, n_vec[i], a_beta, b_beta)

  sum_Ebeta2 = hypdata[, paste("sum_Ebeta2", (i-1), sep="_")]
  sigma2_temp = CI_fish_sigma2(sum_Ebeta2, sum_gamma, n_sample, a_gamma, b_gamma)

  hypcurrent <- c(hypcurrent, pi_temp, sigma2_temp)
  # print(cbind(pi_temp, sigma2_temp))
  hypmat <- rbind(hypmat, c(pi_temp[1], sigma2_temp[1]))
}


#### Summarize log-likelihood
loglike_total = sum(hypdata$log_post_likelihood)

for(i in 1:n_type){
  if(sum(prehyp[i, ]>0)==2){
    if(n_vec[i] > 0){
      a_beta = 2 * n_vec[i] * pp; b_beta = 2 * n_vec[i] - a_beta;
    }else{a_beta=1; b_beta = 1e5 - 1;}

    loglike_total = loglike_total +
      logprior_pi(a_beta, b_beta, prehyp[i, 1]) +
      logprior_sigma(a_gamma, b_gamma, prehyp[i, 2])
  }else{
    print("pre-hyper-parameter <= 0... ")
  }
}


########## Write out updated hyper parameter values
colnames(hypmat) <- c("pi", "sigma2")
print("hyper parameter values updates after MCMC: ")
print(hypmat)
write.table(format(hypmat, scientific=TRUE),
            file=hypcurrent_file,
            quote = FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)


########## Write out updated hyper parameter values and se to EM_result_file
# EM_result_file="/net/fantasia/home/yjingj/GIT/bfGWAS/1KG_example/Test_Wkdir/Eoutput/EM_result.txt"
hypcurrent <- format(hypcurrent, scientific = TRUE)
print("write to hypcurrent file with hyper parameter values after MCMC: ")
print(c(k, hypcurrent))

if(k==0){
  write.table(data.frame(EM_iteration = k, R2 = pve, Loglike = loglike_total,
                        pi_0__sigma2_0__pi_1__sigma2_1 = paste(hypcurrent, collapse="," ),
              file = EM_result_file,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append=FALSE)
}else{
  write.table(data.frame(EM_iteration = k, R2 = pve, Loglike = loglike_total,
                        pi_0__sigma2_0__pi_1__sigma2_1 = paste(hypcurrent, collapse="," ),
              file = EM_result_file,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
}


print("EM step time cost (in minutes) : ")
print((proc.time() - ptm)/60)








