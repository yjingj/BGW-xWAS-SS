options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)

gene_name <- args[[1]]
wkdir <- args[[2]]

print(paste( "Calculate predicted genetically regulated molecular trait for gene:", gene_name ))
print(paste("Working directory:", wkdir))

library(data.table)
library(tidyverse)
setwd(wkdir)
###########################

# gene_name="ABCA7"
### Load grex.geno and weight files
geno_pred<- fread(paste0(gene_name, "_pred.geno.txt"), header=TRUE)
geno_pred$ID <- paste(geno_pred$CHROM, geno_pred$POS, sep = ":")
setkey(geno_pred, "ID")

BGW_weights <- fread(paste0(gene_name, "_BGW_xQTL_weights.txt" ), header=TRUE)
BGW_weights <- BGW_weights[, c("CHR", "POS", "ID", "REF", "ALT",  "Trans", "CPP", "Beta")]
BGW_weights$ID <- paste(BGW_weights$CHR, BGW_weights$POS, sep = ":")
BGW_weights$w = BGW_weights$CPP * BGW_weights$Beta
setkey(BGW_weights, "ID")

### Predict genetically regulated molecular trait
# Sum cis CPP, sum trans CPP
tmp_test = NULL; w_vec = NULL; cis_CPP = 0; trans_CPP = 0;
n_snp = nrow(BGW_weights)
if( n_snp > 0 ){

	for(i in 1:n_snp){
		temp_snp_id = BGW_weights[i, "ID"]
		ref = BGW_weights[i, "REF"]
		alt = BGW_weights[i, "ALT"]

		temp_geno = geno_pred[BGW_weights[i, "ID"], ]
		if(!is.na(temp_geno$POS) & temp_geno$REF == ref & temp_geno$ALT == alt){
			tmp_test = rbind(tmp_test, c(BGW_weights[i, c("CHR", "POS", "ID", "REF", "ALT")], geno_pred[BGW_weights[i, "ID"], -c(1:5)]) )
			w_vec = c(w_vec, BGW_weights[i, ]$w)
			if(BGW_weights[i, ]$Trans){
				trans_CPP = trans_CPP + BGW_weights[i, ]$CPP
			}else{
				cis_CPP = cis_CPP + BGW_weights[i, ]$CPP
			}
		}
		else if(!is.na(temp_geno$POS) & temp_geno$REF == alt & temp_geno$ALT == ref){
			tmp_test = rbind(tmp_test, c(BGW_weights[i, c("CHR", "POS", "ID", "REF", "ALT")], 2 - geno_pred[BGW_weights[i, "ID"], -c(1:5)]) )
			w_vec = c(w_vec, BGW_weights[i, ]$w)
			if(BGW_weights[i, ]$Trans){
				trans_CPP = trans_CPP + BGW_weights[i, ]$CPP
			}else{
				cis_CPP = cis_CPP + BGW_weights[i, ]$CPP
			}
		}
		else{
			print(paste("SNP", temp_snp_id, ref, alt, "do not have test genotype!" ))
		}
	}
}else{
	print(gene_name, "has no xQTL in the weight file!")
}

X_geno <- matrix( as.numeric(t(tmp_test[, -(1:5)]) ), ncol = length(w_vec))
# is.numeric(X_geno)
pred_trait <- X_geno %*% (w_vec)
write_tsv(data.frame(sampleID=colnames(tmp_test)[-(1:5)], pred_trait=pred_trait), file=paste0(gene_name, "_pred_trait.txt"))

total_CPP <- cis_CPP + trans_CPP
write_tsv(data.frame(Gene=gene_name, total_CPP=total_CPP, cis_CPP=cis_CPP, trans_CPP=trans_CPP), file=paste0(gene_name, "_sumCPP.txt"))

