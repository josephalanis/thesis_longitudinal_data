library (BoostMLR)
library(missMethods)
library(bnlearn)
library(mice)
library(caret)
library(ggplot2)
library("Metrics")
library(survey)
library(readr)
library(simstudy)
library(dplyr)
library(nlme)
library(splines)
library(ggplot2)
library(lme4)
library(imputeTS)

#Step 1: Creating Simulated Data

tdef = defData(varname = "Age", dist="normal", formula=50, variance=25)
tdef = defData(tdef, varname = "Sex", dist = "binary", formula = 0.5)
tdef = defData(tdef, varname = "Trt", dist = "binary", formula = 0.5)
tdef = defData(tdef, varname = "Error", dist="normal", formula=0, variance=4)
tdef = defData(tdef, varname = "BP0", dist = "normal", formula = "Age*2.95 + 2.5*Sex + Error", variance = 4)
tdef = defData(tdef, varname = "BP1", dist = "normal", formula = "BP0 - 2 - 0.25 * Trt", variance = 4)
tdef = defData(tdef, varname = "BP2", dist = "normal", formula = "BP0 - 4 - 0.75 * Trt", variance = 4)

set.seed(12345)

dtTrial <- genData(70000, tdef)
dtTrial1 = dtTrial[which(BP0 > 139)]
dtTrial2=dtTrial1[1:50000,]

#Rounding down age and blood pressure
dtTrail3<-floor(dtTrial2)

#fixing id issue

dtTrail3$index<-1:nrow(dtTrail3)

#Blank Matrix to put results in 
b1<-matrix(nrow=500,ncol=1) #bias for method 1
b2<-matrix(nrow=500,ncol=1) #bias for method 1
b3<-matrix(nrow=500,ncol=1) #bias for method 2
b4<-matrix(nrow=500,ncol=1) #bias for method 2
b5<-matrix(nrow=500,ncol=1) #bias for method 3
b6<-matrix(nrow=500,ncol=1) ##bias for method 3
b7<-matrix(nrow=500,ncol=1) ##bias for method 4
b8<-matrix(nrow=500,ncol=1) ##bias for method 4
b9<-matrix(nrow=500,ncol=1) ##bias for method 5
b10<-matrix(nrow=500,ncol=1) ##bias for method 5
b11<-matrix(nrow=500,ncol=1) ##bias for method 6
b12<-matrix(nrow=500,ncol=1) ##bias for method 6
#RMSE
r1<-matrix(nrow=500,ncol=1) #RMSE for Method 1
r2<-matrix(nrow=500,ncol=1) #RMSE for Method 1
r3<-matrix(nrow=500,ncol=1) #RMSE for Method 2
r4<-matrix(nrow=500,ncol=1) #RMSE for Method 2
r5<-matrix(nrow=500,ncol=1) #RMSE for Method 3
r6<-matrix(nrow=500,ncol=1) #RMSE for Method 3
r7<-matrix(nrow=500,ncol=1) #RMSE for Method 4
r8<-matrix(nrow=500,ncol=1) #RMSE for Method 4
r9<-matrix(nrow=500,ncol=1) #RMSE for Method 5
r10<-matrix(nrow=500,ncol=1) #RMSE for Method 5
r11<-matrix(nrow=500,ncol=1) #RMSE for Method 6
r12<-matrix(nrow=500,ncol=1) #RMSE for Method 6

pb <- txtProgressBar(max=500, style=3, width=50)

#Running Simulation

for (i in 1:500)  {
  
  # Next getting a sample of 100
  SRS<-sample(nrow(dtTrial2),500,replace=F) #This is a new sample of 100 random data
  newSRS=round(dtTrial2[SRS,],2) #Makes a new data frame named newSRS
  
  
  #Step 2 creating missing values for scenario 
  
  m1<-delete_MAR_censoring(newSRS,0.3,"BP1",cols_ctrl = "BP0")
  m2<-delete_MAR_censoring(newSRS,0.5,"BP2",cols_ctrl = "BP0")
  
  BP1<-m1$BP1
  BP2<-m2$BP2
  
  newSRS1<-as.data.frame(select(newSRS,-BP1,-BP2))
  
  #Creating missing value subset with 
  s1<-as.data.frame(cbind(newSRS1,BP1,BP2))
  head(s1)
  
  s2<-round(s1,digit=2)
  
  #Step3 - Imputing with Last One Carried Forward
  LOCF<-as.data.frame(na_locf(BP1))
  colnames(LOCF)[1]<-"BP1"
  LOCF1<-as.data.frame(na_locf(BP2))
  colnames(LOCF1)[1]<-"BP2"
  
  
  
  #Step 4 - Imputing with Mice Linear Regression
  mla <- mice(s2, m=1, maxit=0, seed = 1234, print= FALSE)
  
  imp0=mice::complete(mla,1)
  
  # Step 5- Imputing data with 2lonly.pmm
  
  #Getting Predictor Matrix
  micelong0 <- mice(s2, maxit = 0)
  meth_micelong <- micelong0$method
  pred_micelong <- micelong0$predictorMatrix
  
  
  meth_micelong[c("BP1", "BP2")] <- "2lonly.pmm"
  
  pred_micelong[, "id"] <- -2
  pred_micelong[, "BP0"] <- 2
  
  # Next step
  micelong <- mice(s2, meth = meth_micelong, pred = pred_micelong,
                   maxit = 20, seed = 2019, printFlag = FALSE)
  
  imp=mice::complete(micelong,1)
  
  #Step 6 - Bayesian Network (Method 2)
  
  # Need subset with only variables listed in dag1
  dm2<-subset.data.frame(select(s2,-id))
  dag1 = model2network("[Age][Trt][Sex][Error][BP0|Age:Sex:Error][BP1|BP0:Trt][BP2|BP0:Trt]")
  dfitted=bn.fit(dag1,dm2)
  
  incomplete=dm2
  
  missing=matrix(c(sample(nrow(incomplete),500),sample(ncol(incomplete),500,replace=TRUE)),ncol=2)
  
  incomplete[missing] = NA
  #head(incomplete, n = 100 )
  
  completed=impute(dfitted,data=incomplete,method = "parents")
  all(complete.cases(completed))
  
  all(complete.cases(completed))
  
  #Step 7 - Method 3:Bayesian Network with Imputing with Monte Carlo Posterior Inference
  completed1=impute(dfitted,data = incomplete,method = "bayes-lw")
  
  #Step 8 - Method 4: Bayesian Network Imputing with exact inference
  completed2=impute(dfitted,data = incomplete,method = "exact")
  

  #Bias for all methods
  #Method 1
  b1[i]<-round(bias(newSRS$BP1,LOCF$BP1),3)
  b2[i]<-round(bias(newSRS$BP2,LOCF1$BP2),3)
  #Method 2
  b3[i]<-round(bias(newSRS$BP1,imp0$BP1),3)
  b4[i]<-round(bias(newSRS$BP2,imp0$BP2),3)
  #Method 3
  b5[i]<-round(bias(newSRS$BP1,imp$BP1),3)
  b6[i]<-round(bias(newSRS$BP2,imp$BP2),3)
  #Method 4
  b7[i]<-round(bias(newSRS$BP1,completed$BP1),3)
  b8[i]<-round(bias(newSRS$BP2,completed$BP2),3)
  #Method 5
  b9[i]<-round(bias(newSRS$BP1,completed1$BP1),3)
  b10[i]<-round(bias(newSRS$BP2,completed1$BP2),3)
  
  #Method 6
  b11[i]<-round(bias(newSRS$BP1,completed2$BP1),3)
  b12[i]<-round(bias(newSRS$BP2,completed2$BP2),3)
  
  #RMSE for all methods
  #Method 1
  r1[i]<-round(RMSE(newSRS$BP1,LOCF$BP1),3)
  r2[i]<-round(RMSE(newSRS$BP2,LOCF1$BP2),3)
  #Method 2
  r3[i]<-round(RMSE(newSRS$BP1,imp0$BP1),3)
  r4[i]<-round(RMSE(newSRS$BP2,imp0$BP2),3)
  #Method 3
  r5[i]<-round(RMSE(newSRS$BP1,imp$BP1),3)
  r6[i]<-round(RMSE(newSRS$BP2,imp$BP2),3)
  #Method 4
  r7[i]<-round(RMSE(newSRS$BP1,completed$BP1),3)
  r8[i]<-round(RMSE(newSRS$BP2,completed$BP2),3)
  #Method 5
  r9[i]<-round(RMSE(newSRS$BP1,completed1$BP1),3)
  r10[i]<-round(RMSE(newSRS$BP2,completed1$BP2),3)
  
  #Method 6
  r11[i]<-round(RMSE(newSRS$BP1,completed2$BP1),3)
  r12[i]<-round(RMSE(newSRS$BP2,completed2$BP2),3)
  
  # #Output Results into data frame
  final=cbind(b1,b3,b5,b7, b9, b11,b2,b4,b6,b8,b10,b12,r1,r3,r5,r7,r9,r11,r2,r4,r6,r8,r10,r12)
  final1=na.omit(data.frame(final))
  #final2=colMeans(final1)
  setTxtProgressBar(pb, i)
  
}

close(pb)

attach(final1)
#final1

final2<-data.matrix(round(colMeans(final1),3))

final3<-as.data.frame(t(final2))
#View(final3)

Table1<-data.frame(
  Method=c("Method 1","Method 2","Method 3","Method 4", "Method 5" , "Method 6" ),
  tt1bais=c(final3$X1,final3$X2,final3$X3,final3$X4, final3$X5, final3$X6),
  tt2bais=c(final3$X7,final3$X8,final3$X9,final3$X10,final3$X11,final3$X12),
  tt1RMSE=c(final3$X13,final3$X14,final3$X15,final3$X16,final3$X17,final3$X18),
  tt2RMSE=c(final3$X19,final3$X20,final3$X21,final3$X22,final3$X23,final3$X24)
)

#Output
Table1

