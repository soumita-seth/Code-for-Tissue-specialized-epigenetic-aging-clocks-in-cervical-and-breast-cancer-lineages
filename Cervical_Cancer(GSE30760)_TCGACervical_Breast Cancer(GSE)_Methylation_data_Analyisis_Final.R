setwd("D:/PhD Work of Soumita Seth/DNA Methylation Project with Roger/New Clock with Cervical Breast Merge Samples")
############Read GSE30760 27K Methylation Data Cervical Cancer######
FullMethdata<-read.table("GSE30760meth_filt2.csv",row.names=1, header=TRUE,sep=",",check.names = FALSE)
dim(FullMethdata)##[27578, 215]###normal 152 and treated 63
FullMethdata.mat<-as.matrix(apply(FullMethdata, 2, as.numeric))
dim(FullMethdata.mat)##[27578, 215]
rownames(FullMethdata.mat)<-rownames(FullMethdata)
############# Remove rows which have missing value #######################
FullMethdata.mat.omitNA <- FullMethdata.mat[complete.cases(FullMethdata.mat), ] 
dim(FullMethdata.mat.omitNA)         ####[25553   215]

############### Quantile Normalization ###############################
install.packages("BiocManager")
BiocManager::install("preprocessCore")  # Only once
library(preprocessCore)
# Add your user library to the front of the path list
#.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

#install.packages("BiocManager", lib = Sys.getenv("R_LIBS_USER"))
#BiocManager::install("preprocessCore", lib = Sys.getenv("R_LIBS_USER"))

FullMethdata.mat.omitNA.norm <- normalize.quantiles(FullMethdata.mat.omitNA)
dim(FullMethdata.mat.omitNA.norm)
# Preserve original row and column names
dimnames(FullMethdata.mat.omitNA.norm) <- dimnames(FullMethdata.mat.omitNA)
write.table(FullMethdata.mat.omitNA.norm, "GSE30760_qnorm.csv", sep = ",", quote = FALSE,row.names = T, col.names = T)

###########Boxplot Before and After Quantile Normalization#################
##subset_index <- 1:1000
par(mfrow = c(1, 2))  # side-by-side plots
boxplot(FullMethdata.mat.omitNA, main = "Before Normalization", outline = FALSE, col = "skyblue", xaxt = "n")
boxplot(FullMethdata.mat.omitNA.norm, main = "After Normalization", outline = FALSE, col = "lightgreen", xaxt = "n")

# Mean-based filtering (remove extreme methylation sites) (CpG sites with average beta values between 0.1 and 0.9 were retained to remove extremely methylated/unmethylated sites)
mean_beta <- rowMeans(FullMethdata.mat.omitNA.norm, na.rm = TRUE)
FullMethdata.mat.omitNA.norm.filtered1 <- FullMethdata.mat.omitNA.norm[mean_beta > 0.1 & mean_beta < 0.9, ]
dim(FullMethdata.mat.omitNA.norm.filtered1) ###[13179   215]
# Variance-based filtering (remove non-variable sites)(Sites with low variance (e.g., < 0.005) were removed to retain only those with meaningful variability across samples)
var_beta <- apply(FullMethdata.mat.omitNA.norm.filtered1, 1, var, na.rm = TRUE)
FullMethdata.mat.omitNA.norm.filtered2 <- FullMethdata.mat.omitNA.norm.filtered1[var_beta > 0.005, ] 
dim(FullMethdata.mat.omitNA.norm.filtered2)    ###[7013 215]
write.table(FullMethdata.mat.omitNA.norm.filtered2,file="GSE30760_Normalized_Filtered_Probe.csv",quote=F,sep=",",row.names=TRUE,col.name=T)

################Read TCGA 450K Methylation Cervical Cancer #####################
TCGA_Cervical_Methylation450K<- read.table("TCGA_Cervical_Cancer_HumanMethylation450.csv",row.names=1, header=TRUE,sep=",",check.names = FALSE)
dim(TCGA_Cervical_Methylation450K) ##[485577, 312]
TCGA_Cervical_Methylation450K.mat<-as.matrix(apply(TCGA_Cervical_Methylation450K, 2, as.numeric))
dim(TCGA_Cervical_Methylation450K.mat)##[485577, 312]
rownames(TCGA_Cervical_Methylation450K.mat)<-rownames(TCGA_Cervical_Methylation450K)
TCGA_Cervical_Methylation450K.mat.omitNA <- TCGA_Cervical_Methylation450K.mat[complete.cases(TCGA_Cervical_Methylation450K.mat), ] 
dim(TCGA_Cervical_Methylation450K.mat.omitNA)  ##[372137,312]
TCGA_Cervical_Methylation450K.mat.omitNA.norm <- normalize.quantiles(TCGA_Cervical_Methylation450K.mat.omitNA)
dim(TCGA_Cervical_Methylation450K.mat.omitNA.norm) ##[372137    312]
# Preserve original row and column names
dimnames(TCGA_Cervical_Methylation450K.mat.omitNA.norm) <- dimnames(TCGA_Cervical_Methylation450K.mat.omitNA)
View(TCGA_Cervical_Methylation450K.mat.omitNA.norm)
write.table(TCGA_Cervical_Methylation450K.mat.omitNA.norm, "TCGA_Cervical_Methylation450K_qnorm.csv", sep = ",", quote = FALSE,row.names = T, col.names = T)
###########Boxplot Before and After Quantile Normalization#################
##subset_index <- 1:1000
par(mfrow = c(1, 2))  # side-by-side plots
boxplot(TCGA_Cervical_Methylation450K.mat.omitNA, main = "Before Normalization", outline = FALSE, col = "skyblue", xaxt = "n")
boxplot(TCGA_Cervical_Methylation450K.mat.omitNA.norm, main = "After Normalization", outline = FALSE, col = "lightgreen", xaxt = "n")
############# Matched Probes GSE30760 Filtered and TCGA ################
matched_probes <- intersect(rownames(FullMethdata.mat.omitNA.norm.filtered2), rownames(TCGA_Cervical_Methylation450K.mat.omitNA.norm))
print(matched_probes)
num_matched <- length(matched_probes)
print(num_matched) #####5594
write.table(matched_probes, 
            file = "GSE30760Filtered_TCGA_Cervical450K_matched_probes_Omit_NA.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

GSE30760_matched_probeTCGA450K_Methylation_mat <- FullMethdata.mat.omitNA.norm.filtered2[match(matched_probes, rownames(FullMethdata.mat.omitNA.norm.filtered2)), , drop = FALSE]
dim(GSE30760_matched_probeTCGA450K_Methylation_mat)  ####[5594 215]
View(GSE30760_matched_probeTCGA450K_Methylation_mat)
################Read GSE32393 Breast Cancer 27K Methylation Cervical Cancer #####################
BreastCancer_GSE32393_Methdata<-read.table("GSE32393_Breast_Cancer_Probe_Data.csv",row.names=1, header=TRUE,sep=",",check.names = FALSE)
dim(BreastCancer_GSE32393_Methdata)##[27578, 137]###normal 152 and treated 63
View(BreastCancer_GSE32393_Methdata)
BreastCancer_GSE32393_Methdata.mat<-as.matrix(apply(BreastCancer_GSE32393_Methdata,2,as.numeric))
############# Remove rows which have missing value #######################
BreastCancer_GSE32393_Methdata.mat.omitNA <- BreastCancer_GSE32393_Methdata[complete.cases(BreastCancer_GSE32393_Methdata.mat), ] 
dim(BreastCancer_GSE32393_Methdata.mat.omitNA)         ####[26403   137]
BreastCancer_GSE32393_Methdata.mat.omitNA.norm <- normalize.quantiles(BreastCancer_GSE32393_Methdata.mat.omitNA)
dim(BreastCancer_GSE32393_Methdata.mat.omitNA.norm) ###[26403   137]
# Preserve original row and column names
dimnames(BreastCancer_GSE32393_Methdata.mat.omitNA.norm) <- dimnames(BreastCancer_GSE32393_Methdata.mat.omitNA)
write.table(BreastCancer_GSE32393_Methdata.mat.omitNA.norm, "BreastCancer_GSE32393_Methdata_omitNA_qnorm.csv", sep = ",", quote = FALSE,row.names = T, col.names = T)
par(mfrow = c(1, 2))  # side-by-side plots
boxplot(BreastCancer_GSE32393_Methdata.mat.omitNA, main = "Before Normalization", outline = FALSE, col = "skyblue", xaxt = "n")
boxplot(BreastCancer_GSE32393_Methdata.mat.omitNA.norm, main = "After Normalization", outline = FALSE, col = "lightgreen", xaxt = "n")
Breast_matched_probes <- intersect(matched_probes, rownames(BreastCancer_GSE32393_Methdata.mat.omitNA.norm))
print(Breast_matched_probes)
num_matched <- length(Breast_matched_probes)
print(num_matched) #####5360
write.table(Breast_matched_probes, 
            file = "GSE30760Filtered_TCGA_Cervical450K_BreastGSE32393_matched_probes_Omit_NA.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

GSE32393_matched_probeTCGA450K_GSE30760_Methylation_mat <- BreastCancer_GSE32393_Methdata.mat.omitNA.norm[match(Breast_matched_probes, rownames(BreastCancer_GSE32393_Methdata.mat.omitNA.norm)), , drop = FALSE]
dim(GSE32393_matched_probeTCGA450K_GSE30760_Methylation_mat)  ####[5360  137]
View(GSE32393_matched_probeTCGA450K_GSE30760_Methylation_mat)
GSE30760_matched_probeTCGA450K_GSE32393_27K_Methylation_mat<-GSE30760_matched_probeTCGA450K_Methylation_mat[match(Breast_matched_probes, rownames(GSE30760_matched_probeTCGA450K_Methylation_mat)), , drop = FALSE]
dim(GSE30760_matched_probeTCGA450K_GSE32393_27K_Methylation_mat) ###[5360  215]
GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata<-cbind(GSE30760_matched_probeTCGA450K_GSE32393_27K_Methylation_mat,GSE32393_matched_probeTCGA450K_GSE30760_Methylation_mat)
dim(GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata) ###[5360  352]
write.table(GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata, "GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata.csv", sep = ",", quote = FALSE,row.names = T, col.names = T)
GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata.norm <- normalize.quantiles(GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata)
dimnames(GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata.norm) <- dimnames(GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata)
dim(GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata.norm)
write.table(GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata.norm, "GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata.norm.csv", sep = ",", quote = FALSE,row.names = T, col.names = T)
par(mfrow = c(1, 2))  # side-by-side plots
boxplot(GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata, main = "Before Normalization", outline = FALSE, col = "skyblue", xaxt = "n")
boxplot(GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata.norm, main = "After Normalization", outline = FALSE, col = "lightgreen", xaxt = "n")
TCGA_Qnormed_matched_probe_GSE30760_GSE32392_Methylation_mat <- TCGA_Cervical_Methylation450K.mat.omitNA.norm[match(Breast_matched_probes, rownames(TCGA_Cervical_Methylation450K.mat.omitNA.norm)), , drop = FALSE]
dim(TCGA_Qnormed_matched_probe_GSE30760_GSE32392_Methylation_mat)  ####[5360  312]
class(TCGA_Qnormed_matched_probe_GSE30760_GSE32392_Methylation_mat)
write.table(TCGA_Qnormed_matched_probe_GSE30760_GSE32392_Methylation_mat, "TCGA_Qnormed_matched_probe_GSE30760_GSE32392_Methylation_mat.csv", sep = ",", quote = FALSE,row.names = T, col.names = T)
########## Read GSE31979 Breast Cancer 27K Methylation data ########################################
BreastCancer_GSE31979_Methdata<-read.table("GSE31979_Breast_Cancer_Probe_Data.csv",row.names=1, header=TRUE,sep=",",check.names = FALSE)
dim(BreastCancer_GSE31979_Methdata)##[27578, 118]###normal 152 and treated 63
View(BreastCancer_GSE31979_Methdata)
BreastCancer_GSE31979_Methdata.mat<-as.matrix(apply(BreastCancer_GSE31979_Methdata, 2, as.numeric))
dim(BreastCancer_GSE31979_Methdata.mat)##[27578, 118]
rownames(BreastCancer_GSE31979_Methdata.mat)<-rownames(BreastCancer_GSE31979_Methdata)
#BreastCancer_GSE31979_Methdata.mat<-as.matrix(BreastCancer_GSE31979_Methdata)
############# Remove rows which have missing value #######################
BreastCancer_GSE31979_Methdata.mat.omitNA <- BreastCancer_GSE31979_Methdata[complete.cases(BreastCancer_GSE31979_Methdata.mat), ] 
dim(BreastCancer_GSE31979_Methdata.mat.omitNA)         ####[27495   118]
BreastCancer_GSE31979_Methdata.mat.omitNA.mat<-as.matrix(apply(BreastCancer_GSE31979_Methdata.mat.omitNA,2,as.numeric))
BreastCancer_GSE31979_Methdata.mat.omitNA.mat.norm <- normalize.quantiles(BreastCancer_GSE31979_Methdata.mat.omitNA.mat)
dim(BreastCancer_GSE31979_Methdata.mat.omitNA.mat.norm) ###[27495   118]
# Preserve original row and column names
dimnames(BreastCancer_GSE31979_Methdata.mat.omitNA.mat.norm) <- dimnames(BreastCancer_GSE31979_Methdata.mat.omitNA)
View(BreastCancer_GSE31979_Methdata.mat.omitNA.mat.norm)
write.table(BreastCancer_GSE31979_Methdata.mat.omitNA.mat.norm, "BreastCancer_GSE31979_Methdata_omitNA_qnorm.csv", sep = ",", quote = FALSE,row.names = T, col.names = T)
par(mfrow = c(1, 2))  # side-by-side plots
boxplot(BreastCancer_GSE31979_Methdata.mat.omitNA.mat, main = "Before Normalization", outline = FALSE, col = "skyblue", xaxt = "n")
boxplot(BreastCancer_GSE31979_Methdata.mat.omitNA.mat.norm, main = "After Normalization", outline = FALSE, col = "lightgreen", xaxt = "n")
GSE31979_matched_probes <- intersect(Breast_matched_probes, rownames(BreastCancer_GSE31979_Methdata.mat.omitNA.mat.norm))
print(GSE31979_matched_probes)
GSE31979_num_matched <- length(GSE31979_matched_probes)
print(GSE31979_num_matched) #####5348
GSE31979_Qnormed_matched_probe_GSE30760_GSE32392_TCGA_Methylation_mat <- BreastCancer_GSE31979_Methdata.mat.omitNA.mat.norm[match(GSE31979_matched_probes, rownames(BreastCancer_GSE31979_Methdata.mat.omitNA.mat.norm)), , drop = FALSE]
dim(GSE31979_Qnormed_matched_probe_GSE30760_GSE32392_TCGA_Methylation_mat) ####[5348  118]
View(GSE31979_Qnormed_matched_probe_GSE30760_GSE32392_TCGA_Methylation_mat)
write.table(GSE31979_Qnormed_matched_probe_GSE30760_GSE32392_TCGA_Methylation_mat, "GSE31979_Qnormed_matched_probe_GSE30760_GSE32392_TCGA_Methylation_mat.csv", sep = ",", quote = FALSE,row.names = T, col.names = T)
##########GDC TCGA Liver Cancer ########################################
library(data.table)
TCGA_Liver_Methylation450K_dt <- fread("TCGA-LIHC.methylation450.csv")
dim(TCGA_Liver_Methylation450K_dt)        ##[486427 431]
class(TCGA_Liver_Methylation450K_dt)
TCGA_Liver_Methylation450K_dt.omitNA <- TCGA_Liver_Methylation450K_dt[complete.cases(TCGA_Liver_Methylation450K_dt), ] 
dim(TCGA_Liver_Methylation450K_dt.omitNA)  ##[300641    431]
probes <- TCGA_Liver_Methylation450K_dt.omitNA[[1]]
TCGA_Liver_Methylation450K_dt.omitNA[[1]] <- NULL
setDF(TCGA_Liver_Methylation450K_dt.omitNA)
rownames(TCGA_Liver_Methylation450K_dt.omitNA) <- probes
dim(TCGA_Liver_Methylation450K_dt.omitNA)  ###[300641    430]
class(TCGA_Liver_Methylation450K_dt.omitNA)
View(TCGA_Liver_Methylation450K_dt.omitNA[1:10,1:10])
Liver_matched_probes <- intersect(matched_probes, rownames(TCGA_Liver_Methylation450K_dt.omitNA))
print(Liver_matched_probes)
Liver_num_matched <- length(Liver_matched_probes)
print(Liver_num_matched)  ###4402
Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation <- TCGA_Liver_Methylation450K_dt.omitNA[match(Liver_matched_probes, rownames(TCGA_Liver_Methylation450K_dt.omitNA)), , drop = FALSE]
dim(Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation) ####[4402  430]
Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation.mat<-as.matrix(apply(Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation, 2, as.numeric))
Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation.mat.norm <- normalize.quantiles(Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation.mat)
dim(Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation.mat.norm) ###[4402  430]
# Preserve original row and column names
dimnames(Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation.mat.norm) <- dimnames(Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation)
View(Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation.mat.norm)
write.table(Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation.mat.norm, "TCGALiver_matched_probe_GSE30760_GSE32392_TCGA_Methylation_qnorm.csv", sep = ",", quote = FALSE,row.names = T, col.names = T)
par(mfrow = c(1, 2))  # side-by-side plots
boxplot(Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation.mat, main = "Before Normalization", outline = FALSE, col = "skyblue", xaxt = "n")
boxplot(Liver_matched_probe_GSE30760_GSE32392_TCGA_Methylation.mat.norm, main = "After Normalization", outline = FALSE, col = "lightgreen", xaxt = "n")
