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
library(readr)
nlswork <- read_csv("nlswork.csv")
#View(nlswork)

#Checking for missing data
head(nlswork)
sum(is.na(nlswork))
sum(is.na(nlswork$union))
sum(is.na(nlswork$union))
sum(is.na(nlswork$tenure))
sum(is.na(nlswork$wks_work))
sum(is.na(nlswork$ttl_exp))

summary(nlswork$ttl_exp)
summary(nlswork$tenure)

#Creating a dataset with no missing data
full<-na.omit(nlswork)
#summary(full)
#head(full)

#Removing columns that I don't need
full1<-subset.data.frame(select(full,age,tenure,hours,wks_work,ln_wage,ttl_exp ))
#-rownames,-idcode,-year,-birthy_yr,-race,
#-msp,-nev_mar,-grade,-c_city,-south,-ind_code,-occ_code))
#View(full1)

#Checking for correlation between variables
#plot(full1)

#write.csv(full1, "C:\\Users\\alani\\Downloads\\dtTrail3.csv", row.names = FALSE) 


#"C:\Users\alani\Downloads"

#BIAS
#Blank Matrix to put results in 
b1<-matrix(nrow=10,ncol=1) #bias for method 1
b2<-matrix(nrow=10,ncol=1) #bias for method 1
b3<-matrix(nrow=10,ncol=1) #bias for method 2
b4<-matrix(nrow=10,ncol=1) #bias for method 2
b5<-matrix(nrow=10,ncol=1) #bias for method 3
b6<-matrix(nrow=10,ncol=1) ##bias for method 3
b7<-matrix(nrow=10,ncol=1) ##bias for method 4
b8<-matrix(nrow=10,ncol=1) ##bias for method 4
b9<-matrix(nrow=10,ncol=1) ##bias for method 5
b10<-matrix(nrow=10,ncol=1) ##bias for method 5
b11<-matrix(nrow=10,ncol=1) ##bias for method 6
b12<-matrix(nrow=10,ncol=1) ##bias for method 6
#RMSE
r1<-matrix(nrow=10,ncol=1) #RMSE for Method 1
r2<-matrix(nrow=10,ncol=1) #RMSE for Method 1
r3<-matrix(nrow=10,ncol=1) #RMSE for Method 2
r4<-matrix(nrow=10,ncol=1) #RMSE for Method 2
r5<-matrix(nrow=10,ncol=1) #RMSE for Method 3
r6<-matrix(nrow=10,ncol=1) #RMSE for Method 3
r7<-matrix(nrow=10,ncol=1) #RMSE for Method 4
r8<-matrix(nrow=10,ncol=1) #RMSE for Method 4
r9<-matrix(nrow=10,ncol=1) #RMSE for Method 5
r10<-matrix(nrow=10,ncol=1) #RMSE for Method 5
r11<-matrix(nrow=10,ncol=1) #RMSE for Method 6
r12<-matrix(nrow=10,ncol=1) #RMSE for Method 6



pb <- txtProgressBar(max=10, style=3, width=50)

for (i in 1:10)  {

# Next getting a sample of 100
SRS<-sample(nrow(full1),1000,replace=F) #This is a new sample of 100 random data
newSRS1=round(full1[SRS,],4) #Makes a new data frame named newSRS
#View(newSRS1)


#Step 2 creating missing values for scenario 

t1<-delete_MAR_censoring(newSRS1,0.3,"ttl_exp",cols_ctrl = "age")
t2<-delete_MAR_censoring(newSRS1,0.5,"ttl_exp",cols_ctrl = "age")

tt1<-t1$ttl_exp
tt2<-t2$ttl_exp


newSRS2<-as.data.frame(cbind(newSRS1,tt1,tt2))
newSRS3<-as.data.frame(cbind(newSRS1,tt1,tt2))
#newSRS4<-as.data.frame(cbind(newSRS1,tt2))

newSRS3$ID <- 1:nrow(newSRS3)
#View(newSRS3)
#checking missing data
#sum(is.na(newSRS3))
#sum(is.na(newSRS3$tt1))


#Step3 - Imputing with Last One Carried Forward
LOCF<-as.data.frame(na_locf(tt1))
colnames(LOCF)[1]<-"tt1"
LOCF1<-as.data.frame(na_locf(tt2))
colnames(LOCF1)[1]<-"tt2"

#Step 4 - Imputing with Mice Linear Regression
mlr<-mice(newSRS2, m=1, maxit=0, seed = 1234, print= FALSE,remove.collinear = FALSE)

impute=mice::complete(mlr,1)
#View(impute)

# Step 5- Imputing data with 2lonly.pmm

#Getting Predictor Matrix

micelong01 <- mice(newSRS3, maxit = 0)#remove.collinear = FALSE,print= FALSE
meth_micelong1 <- micelong01$method
pred_micelong1 <- micelong01$predictorMatrix


meth_micelong1[c("tt1", "tt2")] <- "2lonly.pmm"

pred_micelong1[, "ID"] <- -2
pred_micelong1[, "ttl_exp"] <- 2

# Next step
micelong1 <- mice(newSRS3, meth = meth_micelong1, pred = pred_micelong1,
                 maxit = 10, seed = 2019, printFlag = FALSE)

imp=mice::complete(micelong1,1)
#View(imp)


#Step 6 - Bayesian Network (Method 2)

# Need subset with only variables listed in dag1
#dm2<-subset.data.frame(select(s2,-id))
dag1a = model2network("[age][hours][wks_work][tenure][ttl_exp|age:hours:wks_work][ln_wage|ttl_exp:age][tt1|ttl_exp][tt2|ttl_exp]")
dfitted1=bn.fit(dag1a,newSRS2)

incomplete1=newSRS2

missing1=matrix(c(sample(nrow(incomplete1),1000),sample(ncol(incomplete1),1000,replace=TRUE)),ncol=2)

incomplete1[missing1] = NA
#head(incomplete1, n = 100 )

completed0=impute(dfitted1,data=incomplete1,method = "parents")
all(complete.cases(completed0))

#all(complete.cases(completed0))
#View(completed0


#Step 7 - Method 3:Bayesian Network with Imputing with Monte Carlo Posterior Inference
completed1a=impute(dfitted1,data = incomplete1,method = "bayes-lw")

#View(completed1a)
#Step 8 - Method 4: Bayesian Network Imputing with exact inference
completed2a=impute(dfitted1,data = incomplete1,method = "exact")
#View(completed2a)


#Bias for all methods
#Method 1
b1[i]<-round(bias(newSRS2$ttl_exp,LOCF$tt1),3)
b2[i]<-round(bias(newSRS2$ttl_exp,LOCF1$tt2),3)
#Method 2
b3[i]<-round(bias(newSRS2$ttl_exp,impute$tt1),3)
b4[i]<-round(bias(newSRS2$ttl_exp,impute$tt2),3)
#Method 3
b5[i]<-round(bias(newSRS2$ttl_exp,imp$tt1),3)
b6[i]<-round(bias(newSRS2$ttl_exp,imp$tt2),3)
#Method 4
b7[i]<-round(bias(newSRS2$ttl_exp,completed0$tt1),3)
b8[i]<-round(bias(newSRS2$ttl_exp,completed0$tt2),3)
#Method 5
b9[i]<-round(bias(newSRS2$ttl_exp,completed1a$tt1),3)
b10[i]<-round(bias(newSRS2$ttl_exp,completed1a$tt1),3)

#Method 6
b11[i]<-round(bias(newSRS2$ttl_exp,completed2a$tt1),3)
b12[i]<-round(bias(newSRS2$ttl_exp,completed2a$tt2),3)

#RMSE for all methods
#Method 1
r1[i]<-round(RMSE(newSRS2$ttl_exp,LOCF$tt1),3)
r2[i]<-round(RMSE(newSRS2$ttl_exp,LOCF1$tt2),3)
#Method 2
r3[i]<-round(RMSE(newSRS2$ttl_exp,impute$tt1),3)
r4[i]<-round(RMSE(newSRS2$ttl_exp,impute$tt2),3)
#Method 3
r5[i]<-round(RMSE(newSRS2$ttl_exp,imp$tt1),3)
r6[i]<-round(RMSE(newSRS2$ttl_exp,imp$tt2),3)
#Method 4
r7[i]<-round(RMSE(newSRS2$ttl_exp,completed0$tt1),3)
r8[i]<-round(RMSE(newSRS2$ttl_exp,completed0$tt2),3)
#Method 5
r9[i]<-round(RMSE(newSRS2$ttl_exp,completed1a$tt1),3)
r10[i]<-round(RMSE(newSRS2$ttl_exp,completed1a$tt1),3)

#Method 6
r11[i]<-round(RMSE(newSRS2$ttl_exp,completed2a$tt1),3)
r12[i]<-round(RMSE(newSRS2$ttl_exp,completed2a$tt2),3)




# #Output Results into data frame
final=cbind(b1,b3,b5,b7, b9, b11,b2,b4,b6,b8,b10,b12,r1,r3,r5,r7,r9,r11,r2,r4,r6,r8,r10,r12)
final1=na.omit(data.frame(final))
#final2=colMeans(final1)
setTxtProgressBar(pb, i)

}

close(pb)

attach(final1)
#final1

final2<-data.matrix(colMeans(final1))

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