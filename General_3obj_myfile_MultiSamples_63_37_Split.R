General_3obj_myfile_MultiSamples<- function(ratio,info,loop_no,weights,criteriaMaxMin2) {
#General_2obj_myfile <- function(GSE30760_ratio,GSE30760_info,loop_no,weights,criteriaMinMax1) {  
  # Reading sample age data
  info.1=info######select those sample IDS from rDNA samples, no whole genome data samples
  #info.1=GSE30760_info######select those sample IDS from rDNA samples, no whole genome data samples
  dim(info.1) #[ 100   2]
  
  ID.2=info
  dim(ID.2)
  
  # Reading rDNA methylation data
  ratioCG=ratio[,colnames(ratio) %in% ID.2$sampleID]###pick up common samples between rDNA methylation, no whole genome data 
  dim(ratioCG)##row=CpG and column=common samples##[5360  100]  
  colnames(ratioCG)=ID.2$sampleID#ID of the rDNA methylation samples, no whole genome data
  ratioCG_whole=ratioCG####rDNA methylation (no whole genome methylation data) for the 255 samples
  dim(ratioCG_whole)##rows=CpGs, column=common samples##[5360  100] 
  
  Methyall.RNA = (t(ratioCG_whole[,which(colnames(ratioCG_whole) %in% info.1$sampleID)]))##rows=common samples, columns=CpGs 
  dim(Methyall.RNA)##[100 5360]
  #View(Methyall.RNA)
  Ageall.RNA = info.1$Age##find the age of the samples
  length(Ageall.RNA)##100
  Ageall.RNA=as.numeric(Ageall.RNA)
  ratioCG_11=t(Methyall.RNA)
  dim(ratioCG_11)##[5360  100]
  
  ############New permutation method using TOPSIS########################
  # Here is the permutation process
  
  # Define empty variables
  aMAE = Count = aRHO =aP = atheta= aFeaturenum =i_noskip= aMAE_Train= aMAE_Test= aRHO_Train= aRHO_Test=aP_Train=aP_Test=atheta_train= atheta_test=c()
  atheta_tan= atheta_cos= atheta_train_tan=atheta_train_cos=atheta_test_tan=atheta_test_cos=c()
  aAgePred_loov_RNA<-c()
  best_Myef<-vector(mode="numeric", length=0) ##initialization
  aMyef<-list()
  aAE_loov<-c()
  aAgePred_loov<-c()
  respective_Count<-0 ##initialization
  respective_RHO<-0 ##initialization
  respective_P<-0 ##initialization
  respective_atheta<-0 ##initialization
  sampidvec <- numeric(0)
  sampidvec <- c(sampidvec, 1:loop_no)
  length(sampidvec)##100
  
  # Loops
  #
  #
  #i=2
  #
  ##Newly added
  #Training_Sample_List <- list()  # Initializes an empty list
  #Test_Sample_List <- list()  # Initializes an empty list
  Cervical_Breast_predictedAge_list_AllSample <- list()
  Cervical_Breast_predictedAge_list_TrainingSample <- list()
  Cervical_Breast_predictedAge_list_TestSample <- list()
  #
  d1_rows <- 1:77
  d2_rows <- 78:100
  #i=1
  #for (i in 1:loop_no){  
  #for (i in 1:10){  
  for(i in 1:100){  ##Newly Updated
  #for(i in 1:4){
  #for(i in seq(1, 4, by = 2)){      
    print(paste0("i=",i))  #####i=1#######
    
    # splitting the training and test datasets
    set.seed(i)###use initial seed value to keep the same result###
    ##Newly Added
    train_d1 <- sample(d1_rows, 50)
    train_d2 <- sample(d2_rows, 13)
    index_training <- c(train_d1, train_d2)  # Randomly 63 sample for training, 50 from cervical, 13 from breast
    index_test <- setdiff(sampidvec, index_training)  # Remaining 37 for testing, 27 from cervical, 10 from breast
    ##Model 2 (80% Training Samples and 20% Test Samples)####
    MethyTraining <- Methyall.RNA[index_training, ]  # Training methylation data
    dim(MethyTraining)  # [63, 5360]
    class(MethyTraining)
    #View(MethyTraining)
    row.names(MethyTraining)
    ##Newly Added
    #Training_Sample_List[[i]]<-row.names(MethyTraining) #Store Training Samples of each iteartion
    #class(Training_Sample_List[[i]])
    #View(as.matrix(Training_Sample_List[[i]]))
    MethyTest <- Methyall.RNA[index_test, ]  # Test methylation data
    dim(MethyTest)  # [37 5360]
    ##Newly Added
    #Test_Sample_List[[i]]<-row.names(MethyTest) #Store Training Samples of each iteartion
    #View(MethyTest) 
    row.names(MethyTest)
    any(is.na(MethyTraining ))
    AgeTraining = Ageall.RNA[index_training]
    #AgeTraining=as.numeric(AgeTraining)
    length(AgeTraining)##57 ##16
    any(is.na(AgeTraining))
    AgeTest = Ageall.RNA[index_test]
    length(AgeTest)##1 ##37
    
    # modeling
    # Obtain a lambda from 10-fold cross validation
    registerDoMC(cores = 2)
    #View(MethyTraining)
    #View(AgeTraining)
    custom_lambda <- 10^seq(-6, -1, length = 100)    ##custom lamda setting definition https://glmnet.stanford.edu/
    #custom_lambda <- 10^seq(-6, 3, length = 100)
    #custom_lambda <- 10^seq(-7, -2, length = 100)
    #custom_lambda <- 10^seq(-4, -1, length = 100)
    #glmnet.Training.CV = cv.glmnet(as.matrix(MethyTraining), log(AgeTraining), nfolds=nrow(MethyTraining),alpha=0.5,family="gaussian",nlambda=500,lambda.min.ratio=0.0001, parallel = T)
    glmnet.Training.CV = cv.glmnet(as.matrix(MethyTraining), log(AgeTraining), nfolds=nrow(MethyTraining),alpha=0.5,family="gaussian", lambda=custom_lambda, lambda.min.ratio=0.0001)
    sum(AgeTraining <= 0)
    lambda.glmnet.Training = glmnet.Training.CV$lambda.1se
    
    # Select CpGs and get their coefficients based on the selected lambda
    #glmnet.Training = glmnet(MethyTraining, log(AgeTraining), family="gaussian", alpha=0.5, nlambda=500, lambda.min.ratio=0.0001)
    glmnet.Training = glmnet(MethyTraining, log(AgeTraining), family="gaussian", alpha=0.5, lambda=custom_lambda, lambda.min.ratio=0.0001)
    Myef_sparse = coef(glmnet.Training,s=lambda.glmnet.Training)
    dim(Myef_sparse)##[33671, 1]
    Myef = Myef_sparse[Myef_sparse[,1]>0|Myef_sparse[,1]<0,]###find the CpGs whose coef is non-zero to remove sparse CpGs
    length(Myef)##207 ##103 ##25
    #Myef=as.matrix(Myef)
    #class(Myef)
    #
    #View(Myef)
    #####newly added####skipping the iteration i if #features in that model is zero###
    if((length(Myef)-1)==0) next # skip 3rd iteration and go to next iteration
    cat(i)
    ###########
    
    print(paste0("i=",i, "noskip"))
    
    # Test the model 
    
    MethyTest_selected = MethyTest[,!is.na(match(colnames(MethyTest),names(Myef)))]###select the subdata of methylation test data considering only those 83 filtered CpGs (columns here)
    #MethyTest_selected = MethyTest[,!is.na(match(colnames(MethyTest),rownames(Myef)))]###select the subdata of methylation test data considering only those 83 filtered CpGs (columns here)
    dim(MethyTest_selected)##row=sample, column=CpG [17,206] ##[17 102] ##[17 24]
    
    AgePred_loov = exp(as.matrix(MethyTest_selected) %*% Myef[-1] + Myef[1])###exponential(meth*coeff each for CpG + intercept)
    length(AgePred_loov)
    aAgePred_loov<-c(aAgePred_loov,AgePred_loov)
    
    AE_loov= abs(AgePred_loov - AgeTest)
    aAE_loov=c(aAE_loov,AE_loov)
    
    ###store Myef for each iteration with Myef of previous iterations as a list##
    aMyef<-c(aMyef,list(Myef))
    
    NewInfo = info.1##find the age of the RNA samples
    class(NewInfo)
    
    dim(NewInfo)##[33,3]
    
    ratioCG_t_11=Methyall.RNA
    dim(ratioCG_t_11)##[33, 33670]
    
    # Read all samples in the model
    AllSamples_selected = ratioCG_t_11[,!is.na(match(colnames(ratioCG_t_11),names(Myef)))]###for 255 samples
    #AllSamples_selected = ratioCG_t_11[,!is.na(match(colnames(ratioCG_t_11),rownames(Myef)))]
    dim(AllSamples_selected)##[33, 206] ##[33 24]
    
    pred_allsamples_RNA = exp(AllSamples_selected %*% Myef[-1] + Myef[1])###exponential(meth*coeff each for rDNA CpG + intercept) 
    dim(pred_allsamples_RNA)##[33, 1]
    NewInfo$Predicted =pred_allsamples_RNA
    Cervical_Breast_predictedAge_list_AllSample[[paste0("iteration", i)]] <- NewInfo
    Cervical_Breast_predictedAge_list_TrainingSample[[paste0("iteration", i)]] <- NewInfo[index_training,]
    Cervical_Breast_predictedAge_list_TestSample[[paste0("iteration", i)]] <- NewInfo[index_test,]
    dim(NewInfo)##[33, 4]
    length(pred_allsamples_RNA)
    MAE = median(abs(pred_allsamples_RNA - Ageall.RNA))##find median value of the absolute difference between predicted age and scaled chronological age for each sample
    aMAE = c(aMAE,MAE) ##merge the current MAE score for this iteration i to the MAE values of previous iterations
    Count = c(Count,length(Myef)-1)##merge the current coef score for this iteration i to the coef values of previous iterations
    CORR = cor.test(pred_allsamples_RNA,Ageall.RNA,method='s',exact=FALSE)##Spearman's rank correlation between the vector of predicted age and the vector of scaled chronological age
    aRHO = c(aRHO,CORR$estimate)##merge the current corr score for this iteration i to the corr values of previous iterations
    aP = c(aP,CORR$p.value)##merge the current p-value of corr test for this iteration i to the p-values of corr tests of previous iterations
    aFeaturenum=c(aFeaturenum,length(Myef)-1)##store number of features in the coefficient matrix for each iteration with thesame with previous iterations
    MAE_Train = median(abs(pred_allsamples_RNA[index_training] - AgeTraining))
    aMAE_Train=c(aMAE_Train,MAE_Train)
    MAE_Test = median(abs(pred_allsamples_RNA[index_test] - AgeTest))
    aMAE_Test=c(aMAE_Test,MAE_Test)
    CORR_Train= cor.test(pred_allsamples_RNA[index_training],AgeTraining,method='s',exact=FALSE)##Spearman's rank correlation between the vector of predicted age and the vector of scaled chronological age
    aRHO_Train= c(aRHO_Train,CORR_Train$estimate)##merge the current corr score for this iteration i to the corr values of previous iterations
    aP_Train = c(aP_Train,CORR_Train$p.value)##merge the current p-value of corr test for this iteration i to the p-values of corr tests of previous iterations
    CORR_Test= cor.test(pred_allsamples_RNA[index_test],AgeTest,method='s',exact=FALSE)##Spearman's rank correlation between the vector of predicted age and the vector of scaled chronological age
    aRHO_Test= c(aRHO_Test,CORR_Test$estimate)##merge the current corr score for this iteration i to the corr values of previous iterations
    aP_Test = c(aP_Test,CORR_Test$p.value)##merge the current p-value of corr test for this iteration i to the p-values of corr tests of previous iterations
    #pdf(file="demo.pdf",width=4.25/2.54,height=4.3/2.54)
    #datout = data.frame(Age=c(NewInfo[],NewInfo[]), rDNAage=c(NewInfo[],NewInfo[]), )
    
    ###notably, for adding more independent variables, https://cran.r-project.org/web/packages/ggiraphExtra/vignettes/ggPredict.html###
    ###http://www.sthda.com/english/wiki/ggplot2-add-straight-lines-to-a-plot-horizontal-vertical-and-regression-lines###
    
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
    
    ##find angel between two skewed (non-parallel and non-intersecting) straight lines### Updated formula###
    ##theta=tan^(-1)[(m1-m2)/(1+(m1*m2))]
    ##https://math.stackexchange.com/questions/1269050/finding-the-angle-between-two-line-equations####
    theta1_tan<- abs(atan(((slope_optiln-slope_regln))/(1+(slope_optiln * slope_regln)))*180/3.14) ### atan(0.998/1.000167)*180/3.14 = atan(0.99783)*180/3.14 
    ### atan((-0.998)/1.000167)*180/3.14 
    
    #if(theta1<0)         ##########https://math.stackexchange.com/questions/1269050/finding-the-angle-between-two-line-equations##
    #{
    # theta1<-180+theta1
    #}
    
    if(theta1_tan>90)     
    {
    theta1__tan<-180-theta1_tan
    }
    #theta1_tan<-theta1
    atheta_tan = c(atheta_tan, theta1_tan)
    
    
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
    
    #atheta = c(atheta,theta1)##merge the current cosine btw two black and blue regression lines for this iteration i to cosine of them of previous iterations
    #i_noskip<-c(i_noskip,i)
    
    
    ####theta for Training Samples####
    NewInfo_Train<-NewInfo[index_training,]
    dim(NewInfo_Train) ##[16,4]
    sp_train<-ggplot(NewInfo_Train,aes(x=Age,y=Predicted)) + geom_point(size = 0.5)+geom_abline(slope=1,linetype="dashed", linewidth=0.5,color="grey")+
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
    sp1_train<-sp_train + stat_smooth(method="lm", formula = y ~ x, se=FALSE)
    ##dev.off()
    
    ##print(sp1)
    
    ###obtain the data values of regression line#####
    ##https://stackoverflow.com/questions/9789871/method-to-extract-stat-smooth-line-fit##
    dat_allsmooth_reglines_train<-ggplot_build(sp1_train)$data[[3]]
    class(dat_allsmooth_reglines_train)##data.frame
    dim(dat_allsmooth_reglines_train)##[80, 11]
    
    reglnsize_train<-nrow(dat_allsmooth_reglines_train)##80
    twpts_regln_train<-dat_allsmooth_reglines_train[c(1,reglnsize_train),]
    slope_regln_train<-(twpts_regln_train$y[2]-twpts_regln_train$y[1])/(twpts_regln_train$x[2]-twpts_regln_train$x[1])##m=(y2-y1)/(x2-x1)=0.001670675
    intercept_regln_train<-twpts_regln_train$y[1] - (slope_regln_train * twpts_regln_train$x[1])##c=y-mx=y1-mx1=1.780789
    
    ###optimal 1:1 line; 45 degree; gray dotted line vector##
    slope_optiln_train<-1##m=(y2-y1)/(x2-x1)
    intercept_optiln_train<-0 ##c=y-mx=y1-mx1
    
    ##find angel between two skewed (non-parallel and non-intersecting) straight lines###
    ##theta=tan^(-1)[(m1-m2)/(1+(m1*m2))]
    ##https://math.stackexchange.com/questions/1269050/finding-the-angle-between-two-line-equations##
    theta1_train_tan<- abs(atan(((slope_optiln_train - slope_regln_train)/(1+(slope_optiln_train * slope_regln_train))))*180/3.14) 
    
   # if(theta1_train<0)      #####https://math.stackexchange.com/questions/1269050/finding-the-angle-between-two-line-equations##
    #{
     # theta1_train<-180+theta1_train
    #}
    if(theta1_train_tan>90)     
    {
      theta1_train_tan<-180-theta1_train_tan
    }
  
    atheta_train_tan = c(atheta_train_tan,theta1_train_tan)##merge the current cosine btw two black and blue regression lines for this iteration i to cosine of them of previous iterations
    
    ################theta_train computed by cos theta####
    #### from data3, pick up two points(x1,y1) and (x2,y2); representing line vector V1=(x2-x1,y2-y1); for 1:1 line V2=(10-5,10-5)
    ## see this http://home.cc.umanitoba.ca/~thomas/Courses/space.pdf, page 32, example angle between the skew lines
    
    reglnsize_train<-nrow(dat_allsmooth_reglines_train)##80
    # Select two random points from the regression line
    twpts_regln_train<-dat_allsmooth_reglines_train[c(1,reglnsize_train),]
    p1_train <- c(twpts_regln_train$x[1],twpts_regln_train$y[1])
    p2_train <- c(twpts_regln_train$x[2],twpts_regln_train$y[2])
    # Compute v1 (vector along the regression line)
    v1_train <- p2_train - p1_train
    # Select two random points from the line y = x
    min_x_train<-min(dat_allsmooth_reglines_train$x)
    max_x_train<-max(dat_allsmooth_reglines_train$x)
    min_y_train<-min(dat_allsmooth_reglines_train$y)
    max_y_train<-max(dat_allsmooth_reglines_train$y)
    lower_limit_train<-max(min_x_train,min_y_train)
    upper_limit_train<-min(max_x_train,max_y_train)
    
    # Generate values for y = x within the given range
    valid_vals_train <- seq(lower_limit_train, upper_limit_train, length.out = reglnsize_train)
    
    # Select two random points from the y = x line within the given range
    indices2_train <- sample(1:reglnsize_train, 2)
    q1_train <- c(valid_vals_train[indices2_train[1]], valid_vals_train[indices2_train[1]])  # y = x
    q2_train <- c(valid_vals_train[indices2_train[2]], valid_vals_train[indices2_train[2]])  # y = x
    
    # Compute v2 (vector along the y = x line)
    v2_train <- q2_train - q1_train
    
    # Compute dot product
    dot_product_train <- sum(v1_train * v2_train)
    
    # Compute magnitudes
    magnitude_v1_train <- sqrt(sum(v1_train^2))
    magnitude_v2_train <- sqrt(sum(v2_train^2))
    
    # Compute cosine of the angle
    cos_theta_train <- dot_product_train / (magnitude_v1_train * magnitude_v2_train)
    
    # Compute angle in degrees
    theta1_train_cos <- abs(acos(cos_theta_train) * (180 / pi))
    
    if(theta1_train_cos>90)     
    {
      theta1_train_cos<-180-theta1_train_cos
    }
    # Print results
    #list(Cos_theta = theta1_cos, v1 = v1, v2 = v2)
    atheta_train_cos = c(atheta_train_cos,theta1_train_cos)
    ##End of theta for Training Samples
    
    
    ##theta for Test Samples
    NewInfo_Test<-NewInfo[index_test,]
    dim(NewInfo_Test) ##[16,4]
    sp_test<-ggplot(NewInfo_Test,aes(x=Age,y=Predicted)) + geom_point(size = 0.5)+geom_abline(slope=1,linetype="dashed", linewidth=0.5,color="grey")+
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
    sp1_test<-sp_test + stat_smooth(method="lm", formula = y ~ x, se=FALSE)
    ##dev.off()
    
    ##print(sp1)
    
    ###obtain the data values of regression line#####
    ##https://stackoverflow.com/questions/9789871/method-to-extract-stat-smooth-line-fit##
    dat_allsmooth_reglines_test<-ggplot_build(sp1_test)$data[[3]]
    class(dat_allsmooth_reglines_test)##data.frame
    dim(dat_allsmooth_reglines_test)##[80, 11]
    
    reglnsize_test<-nrow(dat_allsmooth_reglines_test)##80
    twpts_regln_test<-dat_allsmooth_reglines_test[c(1,reglnsize_test),]
    slope_regln_test<-(twpts_regln_test$y[2]-twpts_regln_test$y[1])/(twpts_regln_test$x[2]-twpts_regln_test$x[1])##m=(y2-y1)/(x2-x1)=0.001670675
    intercept_regln_test<-twpts_regln_test$y[1] - (slope_regln_test * twpts_regln_test$x[1])##c=y-mx=y1-mx1=1.780789
    
    ###optimal 1:1 line; 45 degree; gray dotted line vector##
    slope_optiln_test<-1##m=(y2-y1)/(x2-x1)
    intercept_optiln_test<-0 ##c=y-mx=y1-mx1
    
    ##find angel between two skewed (non-parallel and non-intersecting) straight lines###
    ##theta=tan^(-1)[(m1-m2)/(1+(m1*m2))]
    ##https://math.stackexchange.com/questions/1269050/finding-the-angle-between-two-line-equations##
    theta1_test_tan<- abs(atan(((slope_optiln_test - slope_regln_test)/(1+(slope_optiln_test * slope_regln_test))))*180/3.14) 
    
    #if(theta1_test<0) #####https://math.stackexchange.com/questions/1269050/finding-the-angle-between-two-line-equations##
    #{
     # theta1_test<-180+theta1_test
    #}
    if(theta1_test_tan>90)     
    {
      theta1_test_tan<-180-theta1_test_tan
     }
    
    atheta_test_tan = c(atheta_test_tan,theta1_test_tan)##merge the current cosine btw two black and blue regression lines for this iteration i to cosine of them of previous iterations
    ################theta_test computed by cos theta####
    #### from data3, pick up two points(x1,y1) and (x2,y2); representing line vector V1=(x2-x1,y2-y1); for 1:1 line V2=(10-5,10-5)
    ## see this http://home.cc.umanitoba.ca/~thomas/Courses/space.pdf, page 32, example angle between the skew lines
    
    reglnsize_test<-nrow(dat_allsmooth_reglines_test)##80
    # Select two random points from the regression line
    twpts_regln_test<-dat_allsmooth_reglines_test[c(1,reglnsize_test),]
    p1_test <- c(twpts_regln_test$x[1],twpts_regln_test$y[1])
    p2_test <- c(twpts_regln_test$x[2],twpts_regln_test$y[2])
    # Compute v1 (vector along the regression line)
    v1_test <- p2_test - p1_test
    # Select two random points from the line y = x
    min_x_test<-min(dat_allsmooth_reglines_test$x)
    max_x_test<-max(dat_allsmooth_reglines_test$x)
    min_y_test<-min(dat_allsmooth_reglines_test$y)
    max_y_test<-max(dat_allsmooth_reglines_test$y)
    lower_limit_test<-max(min_x_test,min_y_test)
    upper_limit_test<-min(max_x_test,max_y_test)
    
    # Generate values for y = x within the given range
    valid_vals_test <- seq(lower_limit_test, upper_limit_test, length.out = reglnsize_test)
    
    # Select two random points from the y = x line within the given range
    indices2_test <- sample(1:reglnsize_test, 2)
    q1_test <- c(valid_vals_test[indices2_test[1]], valid_vals_test[indices2_test[1]])  # y = x
    q2_test <- c(valid_vals_test[indices2_test[2]], valid_vals_test[indices2_test[2]])  # y = x
    
    # Compute v2 (vector along the y = x line)
    v2_test <- q2_test - q1_test
    
    # Compute dot product
    dot_product_test <- sum(v1_test * v2_test)
    
    # Compute magnitudes
    magnitude_v1_test <- sqrt(sum(v1_test^2))
    magnitude_v2_test <- sqrt(sum(v2_test^2))
    
    # Compute cosine of the angle
    cos_theta_test <- dot_product_test / (magnitude_v1_test * magnitude_v2_test)
    
    # Compute angle in degrees
    theta1_test_cos <- abs(acos(cos_theta_test) * (180 / pi))
    
    if(theta1_test_cos>90)     
    {
      theta1_test_cos<-180-theta1_test_cos
    }
    # Print results
    #list(Cos_theta = theta1_cos, v1 = v1, v2 = v2)
    atheta_test_cos = c(atheta_test_cos,theta1_test_cos)
    ##End of theta for Test Samples
    i_noskip<-c(i_noskip,i)
    
  }
  ###predicted biological rNA age for each model###
  datout_f_loov_for_PredAge_rNA = as.matrix(pred_allsamples_RNA)
  length(datout_f_loov_for_PredAge_rNA)##58
  #id12<-seq(1,length(i_noskip),by=2)
  colnames(datout_f_loov_for_PredAge_rNA)<-"Predicted_age"
  
  # Output the performance of all the models length(i_noskip)#
  aFeaturenum=Count##newly added
  datout = data.frame(Index=i_noskip,Count,aFeaturenum,aRHO,aMAE,atheta_tan,atheta_cos,aRHO_Train,aMAE_Train,atheta_train_tan,atheta_train_cos,aRHO_Test,aMAE_Test,atheta_test_tan,atheta_test_cos)
  
  ##### 3 objective optimization with all samples ####
  datout_objs<-data.frame(datout$aRHO,datout$aMAE,datout$atheta_tan)
  
  #######topsis: multiple-criteria decision making ########
  result.topsis1<-topsis(as.matrix(datout_objs), weights, criteriaMaxMin2)
  temp23<-result.topsis1$score
  cs<-seq(1,length(i_noskip),by=1)
  names(temp23)<-paste0("#Iteration=",cs)
  
  ##
  temp23_sorted<-temp23[order(temp23,decreasing = TRUE)]
  temp23_sorted
  temp23_sorted_top10<-temp23_sorted[1:10]
  temp23_sorted_top10
  cs_top10<-seq(1,(10),by=1) 
  ##
  
  if(loop_no>=10){
    setEPS()
    postscript("Cervical_GSE30670_Breast_GSE32393_3-objective optimal scores for top 10 best models(All Samples).eps")
    #postscript("HSPC BM50 2-objective optimal scores for top 10 best models_4_1_foldCV_cts_11March(All Samples).eps")
    #postscript("HSPC BM 2-objective optimal scores for top 10 best models_80_20percent_New_Merged.eps")
    #postscript("HSPC BM 2-objective optimal scores for top 10 best models_50percent_New.eps")
    #postscript("GSE30760_2-objective optimal scores for top 10 best models_50percent.eps")
    par(mar=c(12.1,10.1,8.1,8.1))
    original.parameters<- par(no.readonly = TRUE )##to disble x-axis position labels
    par(xaxt="n")##to disble x-axis position labels
    plot(temp23_sorted_top10, type="b",col = "red", xlab="", ylab="TOPSIS optimal score", main = "TOPSIS optimal score for each case study")
    axis(1, at=cs_top10)
    text(cs_top10, par("usr")[3]-0.085, srt=90, pos=1, labels=names(temp23_sorted_top10), xpd=TRUE)
    dev.off()
  }
  
  file.remove(list.files(pattern = "demo.pdf"))
  
  
  ##### 2 objective optimization with training samples ####
  
  datout_objs_train<-data.frame(datout$aRHO_Train,datout$aMAE_Train,datout$atheta_train_tan)
  
  #######topsis: multiple-criteria decision making ########
  result.topsis1.train<-topsis(as.matrix(datout_objs_train), weights, criteriaMaxMin2)
  temp23.train<-result.topsis1.train$score
  cs<-seq(1,length(i_noskip),by=1)
  names(temp23.train)<-paste0("#Iteration=",cs)
  
  ##
  temp23.train_sorted<-temp23.train[order(temp23.train,decreasing = TRUE)]
  temp23.train_sorted
  temp23.train_sorted_top10<-temp23.train_sorted[1:10]
  temp23.train_sorted_top10
  cs_top10_train<-seq(1,(10),by=1) 
  ##
  
  if(loop_no>=10){
    setEPS()
    postscript("Cervical_GSE30670_Breast_GSE32393_3-objective optimal scores for top 10 best models(Training Sample).eps")
    #postscript("HSPC BM50 2-objective optimal scores for top 10 best models_4_1_foldCV_cts_11March(Training Sample).eps")
    #postscript("HSPC BM 2-objective optimal scores for top 10 best models_80_20percent(Training Samples)_New_Merged.eps")
    #postscript("GSE30760_2-objective optimal scores for top 10 best models_50percent.eps")
    par(mar=c(12.1,10.1,8.1,8.1))
    original.parameters<- par(no.readonly = TRUE )##to disble x-axis position labels
    par(xaxt="n")##to disble x-axis position labels
    plot(temp23.train_sorted_top10, type="b",col = "red", xlab="", ylab="TOPSIS optimal score", main = "TOPSIS optimal score for each case study")
    axis(1, at=cs_top10_train)
    text(cs_top10_train, par("usr")[3]-0.085, srt=90, pos=1, labels=names(temp23.train_sorted_top10), xpd=TRUE)
    dev.off()
  }
  
  file.remove(list.files(pattern = "demo.pdf"))
  
  
  ##### 2 objective optimization with test samples ####
  
  datout_objs_test<-data.frame(datout$aRHO_Test,datout$aMAE_Test,datout$atheta_test_tan)
  
  #######topsis: multiple-criteria decision making ########
  result.topsis1.test<-topsis(as.matrix(datout_objs_test), weights, criteriaMaxMin2)
  temp23.test<-result.topsis1.test$score
  cs<-seq(1,length(i_noskip),by=1)
  names(temp23.test)<-paste0("#Iteration=",cs)
  
  ##
  temp23.test_sorted<-temp23.test[order(temp23.test,decreasing = TRUE)]
  temp23.test_sorted
  temp23.test_sorted_top10<-temp23.test_sorted[1:10]
  temp23.test_sorted_top10
  cs_top10_test<-seq(1,(10),by=1) 
  ##
  
  if(loop_no>=10){
    setEPS()
    postscript("Cervical_GSE30670_Breast_GSE32393_3-objective optimal scores for top 10 best models(Test Samples).eps")
    #postscript("HSPC BM50 2-objective optimal scores for top 10 best models_4_1_foldCV_cts_11March(Test Samples).eps")
    #postscript("HSPC BM 2-objective optimal scores for top 10 best models_80_20percent(Test Samples)_New_merged.eps")
    #postscript("GSE30760_2-objective optimal scores for top 10 best models_50percent.eps")
    par(mar=c(12.1,10.1,8.1,8.1))
    original.parameters<- par(no.readonly = TRUE )##to disble x-axis position labels
    par(xaxt="n")##to disble x-axis position labels
    plot(temp23.test_sorted_top10, type="b",col = "red", xlab="", ylab="TOPSIS optimal score", main = "TOPSIS optimal score for each case study")
    axis(1, at=cs_top10_test)
    text(cs_top10_test, par("usr")[3]-0.085, srt=90, pos=1, labels=names(temp23.test_sorted_top10), xpd=TRUE)
    dev.off()
  }
  
  file.remove(list.files(pattern = "demo.pdf"))
  
  
  ##combining objectives and topsis result for each iteration##
  #temp24<-datout[,1:6]
  temp24<-cbind(datout[,1:7],result.topsis1[,2:3])
  temp25<-cbind(datout[,8:11],result.topsis1.train[,2:3])
  temp26<-cbind(temp24,temp25)
  temp27<-cbind(datout[,12:15],result.topsis1.test[,2:3])
  dataout_withtopsis<-cbind(temp26,temp27)
  #dataout_withtopsis<-cbind(datout,result.topsis1)
  #rownames(dataout_withtopsis)<-names(temp23)
  colnames(dataout_withtopsis)[3]<-"#Features in model"
  colnames(dataout_withtopsis)[4]<-"RHO(All Samp)"
  colnames(dataout_withtopsis)[5]<-"MAE(All Samp)"
  colnames(dataout_withtopsis)[6]<-"Responsiveness angle(tan theta)(All Samp)"
  colnames(dataout_withtopsis)[7]<-"Responsiveness angle(cos theta)(All Samp)"
  #colnames(dataout_withtopsis)[]<-"Iteration ID"
  colnames(dataout_withtopsis)[8]<-"Topsis optimal score(All Samp)"
  colnames(dataout_withtopsis)[9]<-"Topsis optimal rank(All Samp)"
  colnames(dataout_withtopsis)[10]<-"RHO(Train Samp)"
  colnames(dataout_withtopsis)[11]<-"MAE(Train Samp)"
  colnames(dataout_withtopsis)[12]<-"Responsiveness angle(tan theta)(Train Samp)"
  colnames(dataout_withtopsis)[13]<-"Responsiveness angle(cos theta)(Train Samp)"
  colnames(dataout_withtopsis)[14]<-"Topsis optimal score(Train Samp)"
  colnames(dataout_withtopsis)[15]<-"Topsis optimal rank(Train Samp)"
  colnames(dataout_withtopsis)[16]<-"RHO(Test Samp)"
  colnames(dataout_withtopsis)[17]<-"MAE(Test Samp)"
  colnames(dataout_withtopsis)[18]<-"Responsiveness angle(tan theta)(Test Samp)"
  colnames(dataout_withtopsis)[19]<-"Responsiveness angle(cos theta)(Test Samp)"
  colnames(dataout_withtopsis)[20]<-"Topsis optimal score(Test Samp)"
  colnames(dataout_withtopsis)[21]<-"Topsis optimal rank(Test Samp)"
  #Output_PredictedAge_plus_2objective_optimization<-list(datout_f_loov_for_PredAge_rNA,dataout_withtopsis,aMyef)
  #return(list(datout_f_loov_for_PredAge_rNA,dataout_withtopsis,aMyef,Training_Sample_List,Test_Sample_List,BM_predictedAge_list_AllSample,BM_predictedAge_list_TrainingSample,BM_predictedAge_list_TestSample)) ##Newly Updated
  return(list(datout_f_loov_for_PredAge_rNA,dataout_withtopsis,aMyef,Cervical_Breast_predictedAge_list_AllSample,Cervical_Breast_predictedAge_list_TrainingSample,Cervical_Breast_predictedAge_list_TestSample)) ##Newly Updated
}

