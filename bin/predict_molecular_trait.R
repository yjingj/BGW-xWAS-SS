options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)

gene_name <- args[[1]]
wkdir <- args[[2]]
pheno_file <- args[[3]]
quant_pheno <- as.logical(args[[4]])

print(paste( "Calculate predicted genetically regulated molecular trait for gene:", gene_name ))
print(paste("Working directory:", wkdir))
print(paste("Test phenotype directory:", pheno_file))
print(paste("quant_pheno:", quant_pheno))


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
pred_trait_dt <- data.frame(SampleID=colnames(tmp_test)[-(1:5)], Pred_trait=pred_trait)
write_tsv(pred_trait_dt, file=paste0(gene_name, "_pred_trait.txt"))

pheno_data <- read.table(file = pheno_file, header = FALSE, col.names = c("SampleID", "Pheno"))

xWAS_data <- merge(pred_trait_dt, pheno_data, by = "SampleID", sort = FALSE)
Pvalue = NA; Test_stat = NA;
if(quant_pheno){
	temp_fit = lm(Pheno ~ Pred_trait, data = xWAS_data)
}else{
	temp_fit = glm(factor(Pheno) ~ Pred_trait, family = "binomial", data = xWAS_data)
}
Test_stat = summary(temp_fit)$coefficients[2, 3] %>% format(format = "e", digits = 4)
Pvalue = summary(temp_fit)$coefficients[2, 4] %>% format(format = "e", digits = 4)

# formatC(x, format = "e", digits = 2)

total_CPP <- cis_CPP + trans_CPP
write_tsv(data.frame(Gene = gene_name, Test_stat = Test_stat,
	Pvalue = Pvalue,
	Total_CPP= format(total_CPP, format = "e", digits = 4),
	Cis_CPP = format(cis_CPP, format = "e", digits = 4),
	Trans_CPP = format(trans_CPP, format = "e", digits = 4)),
	file = paste0(gene_name, "_xWAS_results.txt"))

