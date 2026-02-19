setwd("D:/PhD Work of Soumita Seth/DNA Methylation Project with Roger/New Clock with Cervical Breast Merge Samples")
rm(list=ls())
bioPkgs<-c("pracma","biomaRt","GenomicFeatures","glmnet","doMC","REdaS","topsis")
BiocManager::install(bioPkgs)
install.packages("doMC", repos="http://R-Forge.R-project.org") #Not directly available in bioconductor
##install & Load libraries
library(pracma)##for time measurement (tic/toc)
library(biomaRt)##for gene name conversion
library(GenomicFeatures)
library(glmnet)##for cross-validation
library(ggplot2)
library(doMC)
#library(doParallel) ##Same function with doMC
library(REdaS)
library(topsis)##for multi-objective optimization
library('dplyr')##for "%>%" function

######garbage cleaning ##########
gc() 

####increasing the memory limit#############
memory.limit()###16180
memory.limit(size=41000)###41000

###core number for parallel processing###########
no_cores<-41

###Read Sample Age and Probes###
GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Type_Age_Data <- read.table("GSE30760_GSE32392_Merged_Sample Type and Age.csv", header=TRUE, sep=",")
dim(GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Type_Age_Data)###[352 3]
GSE30760_GSE32392_Merged_Cervical_Breast_data_AllSamples<-read.table("GSE30760Cervical_GSE32392Breast_27KMerged_TCGA450K_CommonProbes_Methdata.norm.csv", header=TRUE, sep=",", row.names = 1, stringsAsFactors=FALSE)
dim(GSE30760_GSE32392_Merged_Cervical_Breast_data_AllSamples) ###[5360  352]
GSE30760_GSE32392_Merged_Normal_Samples_Data<-GSE30760_GSE32392_Merged_Cervical_Breast_data_AllSamples[,GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Type_Age_Data$Sample_Type%in% c("Normal", "Healthy")]
dim(GSE30760_GSE32392_Merged_Normal_Samples_Data)  ###[5360 100]
View(GSE30760_GSE32392_Merged_Normal_Samples_Data) 
GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Age_Data<-cbind(GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Type_Age_Data$Sample_title,GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Type_Age_Data$Age)
dim(GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Age_Data) ###[352   2]
View(GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Age_Data)
colnames(GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Age_Data)<-c("sampleID","Age")
GSE30760_GSE32392_Merged_Normal_Sample_Age_Data<-GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Age_Data[GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Type_Age_Data$Sample_Type %in% c("Normal", "Healthy"),]
dim(GSE30760_GSE32392_Merged_Normal_Sample_Age_Data)  ###[100 2]
View(GSE30760_GSE32392_Merged_Normal_Sample_Age_Data)
class(GSE30760_GSE32392_Merged_Normal_Sample_Age_Data)
GSE30760_GSE32392_Merged_Normal_Sample_Age_Data.df<-as.data.frame(GSE30760_GSE32392_Merged_Normal_Sample_Age_Data)

#######Correcting Row Names for matching with Age Data Sample############
colnames(GSE30760_GSE32392_Merged_Normal_Samples_Data)<-GSE30760_GSE32392_Merged_Normal_Sample_Age_Data.df$sampleID
View(GSE30760_GSE32392_Merged_Normal_Samples_Data)
ratio<-GSE30760_GSE32392_Merged_Normal_Samples_Data
dim(ratio)  ##[5360  100]
View(ratio)
class(ratio)
info<- GSE30760_GSE32392_Merged_Normal_Sample_Age_Data.df
dim(info)  ##[100 2]
info$Age<-as.numeric(info$Age)
class(info)
View(info)
loop_no<-nrow(info)

##Function call "General_3obj_myfile_MultiSamples" with input parameters (raw data, age info, number of loops/iterations used for building the models, weight/preference score of each objective, criteria whether objective is maximization denoting "+" or minimization "-") and store output (biological age)##
weights<-c(1,1,1)##change weight/preference of each objective (Rho,MAE and theta) according to requirement
criteriaMaxMin2<-c("+","-","-") ###c(Rho,MAE,theta)
tic()
source("General_3obj_myfile_MultiSamples_63_37_Split.R")
Output_PredictedAge_plus_3objective_optimization = General_3obj_myfile_MultiSamples(ratio,info,loop_no,weights,criteriaMaxMin2)
Output_PredictedAge_plus_3objective_optimization[[1]]##Predicted age for samples for each model for Topsis multi-objective optimization
Output_PredictedAge_plus_3objective_optimization[[2]]##Topsis multi-objective optimized score for each iteration
Output_PredictedAge_plus_3objective_optimization[[3]]##model co-efficient 
Output_PredictedAge_plus_3objective_optimization[[4]] ##Predicted Age of All samples per iteration
Output_PredictedAge_plus_3objective_optimization[[5]]##Predicted Age of Training samples per iteration
Output_PredictedAge_plus_3objective_optimization[[6]]##Predicted Age of Test samples per iteration
View(Output_PredictedAge_plus_3objective_optimization[[2]])
write.table(as.data.frame(Output_PredictedAge_plus_3objective_optimization[[2]]),file="GSE30760_Cervical_GSE32393_Breast_MergedQNormalized_trainingdata2_models_alliterations_after_optimization.csv",sep=",",quote=F,row.names=F,col.names=T)
aMyef_cervical3<-Output_PredictedAge_plus_3objective_optimization[[3]]
View(aMyef_cervical3)
save(aMyef_cervical3, file = "GSE30760_cervical_GSE32393_breast_trainingdata3_model_Coefficient.rda")###store Myef for all iterations##
toc()##elapsed time

##################Normal Training Model Predicted Age Correlation Plot by Test sample ranking#####################
Normal_Training_Optimized_Models<-as.data.frame(Output_PredictedAge_plus_3objective_optimization[[2]])
###################Top1 Model#########################################
Normal_Test_Top1_index_value <-Normal_Training_Optimized_Models$Index[which(Normal_Training_Optimized_Models$`Topsis optimal rank(Test Samp)`==1)]
#####Plot for 63 Training Samples
Normal_Training_Age_Input1 <-Output_PredictedAge_plus_3objective_optimization[[5]][[Normal_Test_Top1_index_value]]
write.table(Normal_Training_Age_Input1,file="63Normal_Training_Rank1_predictedAge.csv",sep=",",quote=F,row.names=F,col.names=T)
View(Normal_Training_Age_Input1)
Normal_Training_Top1_MAE <-Normal_Training_Optimized_Models$`MAE(Train Samp)`[which(Normal_Training_Optimized_Models$`Topsis optimal rank(Test Samp)`==1)]
chronological_age <-  Normal_Training_Age_Input1$Age
predicted_age <- Normal_Training_Age_Input1$Predicted
# Create a data frame for plotting
age_data <- data.frame(
  Chronological = chronological_age,
  Predicted = predicted_age
)
View(age_data)
# Calculate Spearman Correlation 
spearman_corr <- cor(chronological_age, predicted_age, method = "spearman")
# Plot with Spearman Correlation of Normal Samples
ggplot(age_data, aes(x = Chronological, y = Predicted)) +
  geom_point(color = "blue", size = 3) +                  # Scatter plot
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "pink4") + # Dashed y=x line
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Linear regression line
  ggtitle(paste("63 Normal Training Samples by Freez Ranked Model 1:\nSpearman Correlation: ", round(spearman_corr, 3),"\nMAE:",round(Normal_Training_Top1_MAE,2))) +
  xlab("Chronological Age") + ylab("Predicted Age") +
  coord_cartesian(xlim = c(0, max(age_data$Chronological)), 
                  ylim = c(0, max(age_data$Predicted)))+
  ##theme_minimal()
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "White"),  # Or any light color
    panel.grid.major = element_line(color = "orange", linetype = "dashed", linewidth = 0.6),  # Bold major grid lines
    panel.grid.minor = element_line(color = "green",linetype = "dashed", linewidth = 0.3),    # Clear minor grid lines
    plot.title = element_text(face = "bold", size = 14),  # Bold title
    axis.title.x = element_text(face = "bold", size = 12),  # Bold X-axis label
    axis.title.y = element_text(face = "bold", size = 12)   # Bold Y-axis label
  )

#####Plot for 37 Test Samples
Normal_Test_Age_Input1 <-Output_PredictedAge_plus_3objective_optimization[[6]][[Normal_Test_Top1_index_value]]
write.table(Normal_Test_Age_Input1,file="37Normal_Test_Rank1_predictedAge.csv",sep=",",quote=F,row.names=F,col.names=T)
View(Normal_Test_Age_Input1)
Normal_Test_Top1_MAE <-Normal_Training_Optimized_Models$`MAE(Test Samp)`[which(Normal_Training_Optimized_Models$`Topsis optimal rank(Test Samp)`==1)]
chronological_age <-  Normal_Test_Age_Input1$Age
predicted_age <- Normal_Test_Age_Input1$Predicted
# Create a data frame for plotting
age_data <- data.frame(
  Chronological = chronological_age,
  Predicted = predicted_age
)
View(age_data)
# Calculate Spearman Correlation 
spearman_corr <- cor(chronological_age, predicted_age, method = "spearman")
# Plot with Spearman Correlation of Normal Samples
ggplot(age_data, aes(x = Chronological, y = Predicted)) +
  geom_point(color = "blue", size = 3) +                  # Scatter plot
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "pink4") + # Dashed y=x line
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Linear regression line
  ggtitle(paste("37 Normal Test Samples by Freez Ranked Model 1:\nSpearman Correlation: ", round(spearman_corr, 3),"\nMAE:",round(Normal_Test_Top1_MAE,2))) +
  xlab("Chronological Age") + ylab("Predicted Age") +
  coord_cartesian(xlim = c(0, max(age_data$Chronological)), 
                  ylim = c(0, max(age_data$Predicted)))+
  ##theme_minimal()
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "White"),  # Or any light color
    panel.grid.major = element_line(color = "orange", linetype = "dashed", linewidth = 0.6),  # Bold major grid lines
    panel.grid.minor = element_line(color = "green",linetype = "dashed", linewidth = 0.3),    # Clear minor grid lines
    plot.title = element_text(face = "bold", size = 14),  # Bold title
    axis.title.x = element_text(face = "bold", size = 12),  # Bold X-axis label
    axis.title.y = element_text(face = "bold", size = 12)   # Bold Y-axis label
  )

###### Predicting 48 Diseased Sample Age using evolved models from 77 normal sample data of Cervical Cancer ####
GSE30760_GSE32392_Merged_Cervical_Cancer_Samples_Data<-GSE30760_GSE32392_Merged_Cervical_Breast_data_AllSamples[,GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Type_Age_Data$Sample_Type=="Cervical Cancer"]
dim(GSE30760_GSE32392_Merged_Cervical_Cancer_Samples_Data)  ###[5360 48]
View(GSE30760_GSE32392_Merged_Cervical_Cancer_Samples_Data)

GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data<-GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Age_Data[GSE30760_GSE32392_Merged_Cervical_Breast_Cancer_Type_Age_Data$Sample_Type=="Cervical Cancer",]
dim(GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data)  ###[48 2]
View(GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data)
class(GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data)
GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data.df<-as.data.frame(GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data)

#######Correcting Row Names for matching with Age Data Sample############
colnames(GSE30760_GSE32392_Merged_Cervical_Cancer_Samples_Data)<-GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data.df$sampleID
View(GSE30760_GSE32392_Merged_Cervical_Cancer_Samples_Data)
GSE30760_GSE32392_Merged_Cervical_Cancer_Samples_Data_t<-t(GSE30760_GSE32392_Merged_Cervical_Cancer_Samples_Data)
dim(GSE30760_GSE32392_Merged_Cervical_Cancer_Samples_Data_t)  ###[48 5360]
GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data.df$Age<-as.numeric(GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data.df$Age)
Freez_Training_optimization<-Output_PredictedAge_plus_3objective_optimization[[2]][Output_PredictedAge_plus_3objective_optimization[[2]]$Index %in% c(Normal_Test_Top1_index_value), ]
View(Freez_Training_optimization)
aAgePred_CV=aAE_CV=aMAE=aRHO=aP=atheta_tan=atheta_cos=Count=c()
Cervical_Cancer_Sample_predictedAge_list <- list()
nModel<-length(aMyef_cervical3)
Myef1<- aMyef_cervical3[[Normal_Test_Top1_index_value]]
class(Myef1)
names(Myef1)
Cervical_Cancer_Sample_Data_selected = GSE30760_GSE32392_Merged_Cervical_Cancer_Samples_Data_t[,!is.na(match(colnames(GSE30760_GSE32392_Merged_Cervical_Cancer_Samples_Data_t),names(Myef1)))]###select the subdata of methylation test data considering only those filtered CpGs (columns here)
# Get the sorted column names of Cervical_Cancer_Sample_Data_selected while preserving all names in Myef
sorted_names <- sort(names(Myef1))
length(sorted_names)
View(Cervical_Cancer_Sample_Data_selected)
dim(Cervical_Cancer_Sample_Data_selected)  ####[48 30]
sorted_names_existing <- sorted_names[sorted_names %in% colnames(Cervical_Cancer_Sample_Data_selected)]
# Reorder Cervical_Cancer_Sample_Data_selected based on sorted column names
Cervical_Cancer_Sample_Data_selected <-Cervical_Cancer_Sample_Data_selected[, sorted_names_existing, drop = FALSE]
dim(Cervical_Cancer_Sample_Data_selected)
# Reorder Myef based on the same sorted names
Myef_intercept<-Myef1[1]
Myef1 <- Myef1[sorted_names_existing]
AgePred_CV = exp(as.matrix(Cervical_Cancer_Sample_Data_selected) %*% Myef1 + Myef_intercept)###exponential(meth*coeff each for CpG + intercept)
length(AgePred_CV)
#View(AgePred_CV)
aAgePred_CV<-c(aAgePred_CV,AgePred_CV)
AE_CV= abs(AgePred_CV - GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data.df$Age)
aAE_CV=c(aAE_CV,AE_CV)   
NewInfo <- GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data.df ##find the age of the RNA samples
dim(NewInfo)##[48 2]
NewInfo$Predicted =AgePred_CV
dim(NewInfo)##[48 3]
View(NewInfo)
length(AgePred_CV)
Cervical_Cancer_Sample_predictedAge_list[[paste0("Model ", Normal_Test_Top1_index_value)]] <- NewInfo
MAE = median(abs(AgePred_CV - GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data.df$Age))##find median value of the absolute difference between predicted age and scaled chronological age for each sample
aMAE = c(aMAE,MAE) ##merge the current MAE score for this iteration i to the MAE values of previous iterations
Count = c(Count,length(Myef1)-1)##merge the current coef score for this iteration i to the coef values of previous iterations
CORR = cor.test(AgePred_CV,GSE30760_GSE32392_Merged_Cervical_Cancer_Sample_Age_Data.df$Age,method='s',exact=FALSE)##Spearman's rank correlation between the vector of predicted age and the vector of scaled chronological age
aRHO = c(aRHO,CORR$estimate)##merge the current corr score for this iteration i to the corr values of previous iterations
aP = c(aP,CORR$p.value)##merge the current p-value of corr test for this iteration i to the p-values of corr tests of previous iterationsaFeaturenum=c(aFeaturenum,length(Myef)-1)##store number of features in the coefficient matrix for each iteration with thesame with previous iterations
#############Compute Theta###############################
sp<-ggplot(NewInfo,aes(x=Age,y=Predicted)) + geom_point(size = 0.5)+geom_abline(slope=1,linetype="dashed", linewidth=0.5,color="grey")+
  guides(color=guide_legend(title=NULL,override.aes = list(size=1.6)))+
  ylab("Predicted age")+ xlab("Chronological Age")+
  theme_bw()+
  theme(axis.text = element_text(size=unit(6,"pt")),
        axis.title.y=element_text(margin=margin(r=1,l=-2)), 
        axis.title.x=element_text(margin=margin(t=3,b=-3)),
        axis.title = element_text(size=unit(8,"pt")),
        legend.position = c(0,1.06),
        legend.background=element_blank(),
        legend.justification = c(0,1),
        legend.text = element_text(size=5),
        legend.key = element_blank(),
        #        legend.spacing.y = unit(5,"pt"),
        legend.key.size = unit(9,"pt"),
        plot.title=element_text(size=unit(9,"pt"),hjust=0.5,margin=margin(t=0,b=0)),
        axis.text.x = element_text(size=unit(6,"pt")))

##Note that, the function stat_smooth() can be used for fitting smooth models to data###
sp1<-sp + stat_smooth(method="lm", formula = y ~ x, se=FALSE)
##dev.off()

##print(sp1)

###obtain the data values of regression line#####
##https://stackoverflow.com/questions/9789871/method-to-extract-stat-smooth-line-fit##
dat_allsmooth_reglines<-ggplot_build(sp1)$data[[3]]
class(dat_allsmooth_reglines)##data.frame
dim(dat_allsmooth_reglines)##[80, 11]

reglnsize<-nrow(dat_allsmooth_reglines)##80
twpts_regln<-dat_allsmooth_reglines[c(1,reglnsize),]

slope_regln<-(twpts_regln$y[2]-twpts_regln$y[1])/(twpts_regln$x[2]-twpts_regln$x[1])##m=(y2-y1)/(x2-x1)=0.001670675
intercept_regln<-twpts_regln$y[1] - (slope_regln * twpts_regln$x[1])##c=y-mx=y1-mx1=1.780789

###optimal 1:1 line; 45 degree; gray dotted line vector##
slope_optiln<-1##m=(y2-y1)/(x2-x1)
intercept_optiln<-0 ##c=y-mx=y1-mx1

##find angel between two skewed (non-parallel and non-intersecting) straight lines###
##theta=tan^(-1)[(m1-m2)/(1+(m1*m2))]
##https://math.stackexchange.com/questions/1269050/finding-the-angle-between-two-line-equations##
theta1_tan<- abs(atan(((slope_optiln-slope_regln)/(1+(slope_optiln * slope_regln))))*180/3.14) 

if(theta1_tan>90)     
{
  theta1_tan<-180-theta1_tan
}

atheta_tan = c(atheta_tan,theta1_tan)##merge the current cosine btw two black and blue regression lines for this iteration i to cosine of them of previous iterations
################theta computed by cos theta####
#### from data3, pick up two points(x1,y1) and (x2,y2); representing line vector V1=(x2-x1,y2-y1); for 1:1 line V2=(10-5,10-5)
## see this http://home.cc.umanitoba.ca/~thomas/Courses/space.pdf, page 32, example angle between the skew lines

reglnsize<-nrow(dat_allsmooth_reglines)##80
# Select two random points from the regression line
twpts_regln<-dat_allsmooth_reglines[c(1,reglnsize),]
p1 <- c(twpts_regln$x[1],twpts_regln$y[1])
p2 <- c(twpts_regln$x[2],twpts_regln$y[2])
# Compute v1 (vector along the regression line)
v1 <- p2 - p1
# Select two random points from the line y = x
min_x<-min(dat_allsmooth_reglines$x)
max_x<-max(dat_allsmooth_reglines$x)
min_y<-min(dat_allsmooth_reglines$y)
max_y<-max(dat_allsmooth_reglines$y)
lower_limit<-max(min_x,min_y)
upper_limit<-min(max_x,max_y)

# Generate values for y = x within the given range
valid_vals <- seq(lower_limit, upper_limit, length.out = reglnsize)

# Select two random points from the y = x line within the given range
indices2 <- sample(1:reglnsize, 2)
q1 <- c(valid_vals[indices2[1]], valid_vals[indices2[1]])  # y = x
q2 <- c(valid_vals[indices2[2]], valid_vals[indices2[2]])  # y = x

# Compute v2 (vector along the y = x line)
v2 <- q2 - q1

# Compute dot product
dot_product <- sum(v1 * v2)

# Compute magnitudes
magnitude_v1 <- sqrt(sum(v1^2))
magnitude_v2 <- sqrt(sum(v2^2))

# Compute cosine of the angle
cos_theta <- dot_product / (magnitude_v1 * magnitude_v2)

# Compute angle in degrees
theta1_cos <- abs(acos(cos_theta) * (180 / pi))
if(theta1_cos>90)     
{
  theta1_cos<-180-theta1_cos
}
# Print results
#list(Cos_theta = theta1_cos, v1 = v1, v2 = v2)
atheta_cos = c(atheta_cos,theta1_cos)
print(aMAE)
print(atheta_tan)
print(atheta_cos)
print(aRHO)
aMAE_df<-as.data.frame(aMAE)
aRHo_df<-as.data.frame(aRHO)
atheta_tan_df<-as.data.frame(atheta_tan)
atheta_cos_df<-as.data.frame(atheta_cos)
aAgePred_CV_df<-as.data.frame(aAgePred_CV)
Count_df<-as.data.frame(Count)
##### 3 objective optimization with all samples ####
Cervical_Cancer_dataout_objs<-data.frame(aRHo_df,aMAE_df,atheta_tan_df)
colnames(aRHo_df)<-"RHO( Cervical Cancer TestData)"
colnames(aMAE_df)<-"MAE(Cervical Cancer TestData)"
colnames(atheta_tan_df)<-"Responsiveness Angle(Cervical Cancer TestData)(tan theta)"
colnames(atheta_cos_df)<-"Responsiveness Angle(Cervical Cancer TestData)(cos theta)"
temp1_Cancer<-cbind(Freez_Training_optimization,aRHo_df)
temp2_Cancer<-cbind(temp1_Cancer,aMAE_df)
temp3_Cancer<-cbind(temp2_Cancer,atheta_tan_df)
Freez_Model_3objective_optimization_TestData<-cbind(temp3_Cancer,atheta_cos_df)
View(Freez_Model_3objective_optimization_TestData)
write.table(Freez_Model_3objective_optimization_TestData,file="GSE30760_GSE32392_Merged_Cervical_Breast_100Normaltrainingdata_common_probesTCGA_48CervicalCancertest_Top1FreezModel_after_3optimization.csv",sep=",",quote=F,row.names=F,col.names=T)
##################Normal and Cervical Cancer Validation Model Age Prediction Plot using Freez Ranked model with Test sample Top Ranking#####################
###################Top1 Model#########################################
Cervical_Cancer_Test_Age_Input1 = Cervical_Cancer_Sample_predictedAge_list[[paste0("Model ", as.character(Normal_Test_Top1_index_value))]]
write.table(Cervical_Cancer_Test_Age_Input1,file="48CervicalCancer_Test_common_probesTCGA_Rank1_predictedAge.csv",sep=",",quote=F,row.names=F,col.names=T)
View(Cervical_Cancer_Test_Age_Input1)
Cervical_Cancer_Test_Top1_MAE <- aMAE
chronological_age <-  Cervical_Cancer_Test_Age_Input1$Age
predicted_age <- Cervical_Cancer_Test_Age_Input1$Predicted
# Create a data frame for plotting
age_data <- data.frame(
  Chronological = chronological_age,
  Predicted = predicted_age
)
View(age_data)
# Calculate Spearman Correlation 
spearman_corr <- cor(chronological_age, predicted_age, method = "spearman")
# Plot with Spearman Correlation of Normal Samples
min_CAge<-min(age_data$Chronological)
min_PAge<-min(age_data$Predicted)
minlim<-min(min_CAge,min_PAge)
ggplot(age_data, aes(x = Chronological, y = Predicted)) +
  geom_point(color = "blue", size = 3) +                  # Scatter plot
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "pink4") + # Dashed y=x line
  geom_smooth(method = "lm", color = "red", se = TRUE) +  # Linear regression line
  ggtitle(paste("48 Cervical Cancer Test Samples by Freez Ranked Model 1:\nSpearman Correlation: ", round(spearman_corr, 3),"\nMAE:",round(Cervical_Cancer_Test_Top1_MAE,2))) +
  xlab("Chronological Age") + ylab("Predicted Age") +
  coord_cartesian(xlim = c(minlim-5, max(age_data$Chronological)), 
                  ylim = c(minlim-5, max(age_data$Predicted)))+
  ##theme_minimal()
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "White"),  # Or any light color
    panel.grid.major = element_line(color = "orange", linetype = "dashed", linewidth = 0.6),  # Bold major grid lines
    panel.grid.minor = element_line(color = "green",linetype = "dashed", linewidth = 0.3),    # Clear minor grid lines
    plot.title = element_text(face = "bold", size = 14),  # Bold title
    axis.title.x = element_text(face = "bold", size = 12),  # Bold X-axis label
    axis.title.y = element_text(face = "bold", size = 12)   # Bold Y-axis label
  )



