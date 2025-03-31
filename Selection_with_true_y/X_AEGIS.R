library(readODS)
library(GEOquery)



# get the subset of differentially expressed genes

df.DE.names <- read_ods(path = "/home/julie/Documents/Paper_codes/Data/AEGIS/djw327_supplemental_file_1.ods")
df.DE.names$leading_edge_genes <- grepl("\\*", df.DE.names$Probeset) # add a column to save leading edge genes information
df.DE.names$Probeset <- gsub("\\*", "", df.DE.names$Probeset)

# get the corresponding transcriptomic matrix

gse=getGEO(filename="/home/julie/Documents/Paper_codes/Data/AEGIS/GSE80796_series_matrix.txt")

AEGIS.expression.data <- gse@assayData$exprs

X_AEGIS = t(AEGIS.expression.data[rownames(AEGIS.expression.data) %in% df.DE.names$Probeset, ])

# get the cancer status

Clinic_data_perez_tot <- read.table(file = '/home/julie/Bureau/Data/perezRodgers2017/clindat_w_scores.tsv', head = TRUE, sep = "\t")

y_AEGIS = Clinic_data_perez_tot$cancer_status
y_AEGIS[which(y_AEGIS == "Cancer")] = 1
y_AEGIS[which(y_AEGIS == "NoCancer")] = 0
y_AEGIS = as.numeric(y_AEGIS)

save(X_AEGIS, y_AEGIS, file = "/home/julie/Bureau/Methods/Knockoffs/New_cohorts/TrueY/R_files/AEGIS_data.R")


