rm(list = ls())

library(readODS)
library(tidyr)
library(readr)
library(dplyr)
library(GEOquery)




#### Create X_BC ####

## Create selected genes subset ##

# open files

df_annotated_genes <- read_ods(path = "/home/julie/Documents/Paper_codes/Data/BC/ds_po.17.00135-3_list_genes_selected.ods", col_names = TRUE)

# remove genes with null coefficients for all biomarkers

df_selected_genes <- df_annotated_genes[rowSums(df_annotated_genes[, (ncol(df_annotated_genes)-10):ncol(df_annotated_genes)] != 0) > 0, ]

hgnc.selected.names <- df_selected_genes$Gene.symbol

## Get the transcriptomic matrix for the set of selected genes ##

# open data

RNA.BC.HGNC <- read_csv(file = "/home/julie/Documents/Paper_codes/Data/BC/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv")
RNA.BC.HGNC <- as.data.frame(RNA.BC.HGNC)
rownames(RNA.BC.HGNC) <- RNA.BC.HGNC[[1]]  # Set first column as row names
RNA.BC.HGNC <- RNA.BC.HGNC[, -1]

# Vanished genes identification

vanished.gene <- which(hgnc.selected.names %in% rownames(RNA.BC.HGNC) == FALSE)
hgnc.selected.names[vanished.gene]

# equivalence with GeneCards (name + function)
hgnc.selected.names <- c(hgnc.selected.names, "ADIRF", "TICRR", "MISP", "DNAAF3", "AUNIP", "CIART", "EFCAB12", "CCDC170", "ERICH5", "SAPCD2", "EFCC1", "CTSV", "CT83", "EVA1B", "MTFR2", "C5AR2")

#### Create real outcomes vector ####
# select rows and remove replicates

df.BC <- RNA.BC.HGNC[rownames(RNA.BC.HGNC) %in% hgnc.selected.names,]
df.BC <- df.BC[, which(grepl("repl", colnames(df.BC)) == FALSE)]
X_BC <- t(as.matrix(df.BC))


gse=getGEO(filename="/home/julie/Documents/Paper_codes/Data/BC/GSE96058-GPL11154_series_matrix.txt.gz")
gse_add=getGEO(filename="/home/julie/Documents/Paper_codes/Data/BC/GSE96058-GPL18573_series_matrix.txt.gz")

df.clinical_data_BC <- rbind(data.frame(gse@phenoData@data), data.frame(gse_add@phenoData@data))
df.clinical_data_BC <- df.clinical_data_BC[which(grepl("repl", df.clinical_data_BC$title) == FALSE), ]

# ensure that patients are the same in both data 

sum(df.clinical_data_BC$title %in% colnames(df.BC))

# reorder the name of the data frame to match X_BC

df.clinical_data_BC <- df.clinical_data_BC[match(colnames(df.BC), df.clinical_data_BC$title), ]

# select important features

df.clinical_data_BC.tronc <- data.frame(df.clinical_data_BC$er.status.ch1, df.clinical_data_BC$her2.status.ch1, df.clinical_data_BC$pgr.status.ch1, df.clinical_data_BC$ki67.status.ch1)
rownames(df.clinical_data_BC.tronc) <- df.clinical_data_BC$title
df.clinical_data_BC.tronc <-  mutate_all(df.clinical_data_BC.tronc, as.numeric)

# remove rows with NA 

df.omit.all = na.omit(df.clinical_data_BC.tronc)

df.omit.er.her.pgr <- df.clinical_data_BC.tronc[complete.cases(df.clinical_data_BC.tronc[, -4]), ]

y_BC.er = df.omit.er.her.pgr$df.clinical_data_BC.er.status.ch1
y_BC.hr2 = df.omit.er.her.pgr$df.clinical_data_BC.her2.status.ch1
y_BC.pgr = df.omit.er.her.pgr$df.clinical_data_BC.pgr.status.ch1

# remove rows in X_BC
X_BC <- X_BC[rownames(X_BC) %in% rownames(df.omit.er.her.pgr), ]


save(X_BC, y_BC.er, y_BC.hr2, y_BC.pgr, file = "/home/julie/Documents/Paper_codes/Selection_with_true_y/R_files/BC_data.R")



