


library(stringr)
library(MASS)
library(tictoc)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(randomForest)
library(knockoff)
library(glmnet)
library(mboost)
library(ggpubr) 
library(mvtnorm)
library(doParallel)
library(xgboost)
library(reticulate)


# List of target FDR for the knockoffs
list_target_fdr <- seq(0,1,0.01)

#list of lambda for the lasso penalized logistic regression
list_lambda <- seq(0, 0.40, 0.004)

library(readr) # open faster csv files

# open data

RNA.BC.HGNC <- read_csv(file = "/home/julie/Documents/Paper_codes/Data/BC/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv", col_names = TRUE)

# get HGNC names
HGNC_names <- RNA.BC.HGNC$...1

####################################### Conversion between ENS and HGNC names #######################################################################################

library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
# 
ensdb <- EnsDb.Hsapiens.v75
keytypes(ensdb)

# Check that the genome version corresponds to the version used in the paper 

metadata_v75 <- metadata(ensdb)

# conversion

ENS_names_list <- mapIds(ensdb, keys = HGNC_names, keytype = "SYMBOL",       # Input type
                         column = "GENEID",        # Desired output
                         multiVals = "list"       # keep all ENS id's when on hgnc genes has different ENS equivalent
)


save(ENS_names_list, file = '/home/julie/Documents/Paper_codes/Data/BC/ENS_names_filter.R')

df_BC_names <- do.call(rbind, lapply(names(ENS_names_list), function(hgnc_genes) {
  ensembl_ids <- ENS_names_list[[hgnc_genes]]
  data.frame(HGNC_symbol = rep(hgnc_genes, length(ensembl_ids)),
             Ensembl_gene_id = ensembl_ids, 
             stringsAsFactors = FALSE)
}))

# check whether there are duplicated ENS id 

sum(duplicated(df_BC_names$Ensembl_gene_id))

### Create X_BC ###

RNA.data <- read.table(file = "/home/julie/Documents/Paper_codes/Data/CRUKPAP/dat_vst_NS_no_Inel_selected.txt", head = T) 

names.CRUKPAP <- colnames(RNA.data) # Exctarc names of the 749 selected genes 

names.BC.ENS <- df_BC_names$Ensembl_gene_id[df_BC_names$Ensembl_gene_id %in% names.CRUKPAP] # Extract genes names from the BC data matrix

names.BC.hgnc <- df_BC_names$HGNC_symbol[df_BC_names$Ensembl_gene_id %in% names.CRUKPAP] # get the correspondance with HGNC names 

RNA.BC.HGNC.part <- RNA.BC.HGNC[RNA.BC.HGNC$...1 %in% names.BC.hgnc, ] # Extract the transcriptomic data
X_BC = t(as.matrix(RNA.BC.HGNC.part[, sample(2:3274, 500, replace = FALSE)])) # choose randomly among the patients (without the replicate)
colnames(X_BC) <- names.BC.ENS

n <- nrow(X_BC)
p <- ncol(X_BC)


save(X_BC, file =  "/home/julie/Documents/Paper_codes/Data/BC/X_BC_exp.R")
