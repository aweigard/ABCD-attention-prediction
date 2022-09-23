rm(list=ls())

# load all necessary packages
library(lavaan)
library(semPlot)
library(plotly)
library(magrittr)
library(mice)
library(testthat)
library(car)
library(olsrr)
library(glmnet)
library(e1071)
library(fastDummies)
library(FactoMineR)

#######################################################################
#######################################################################
##### Pre-Processing  #################################################
#######################################################################
#######################################################################

##################################################################
#### Measurement model for attention problems outcome measure ####
##################################################################

# First, we load and reformat a merged data set with measures of 
# attention problems from the year 1 followup visit and a variety 
# of neurocognitive, psychosocial and demographic measures from the 
# earlier (baseline) visit, which we will use to predict prospective 
# attention problems.

# load data 

y1.attn<-read.csv("year1_full_dat.csv",stringsAsFactors = TRUE)

# re-level and name reference levels for each factor

y1.attn$female <- relevel(y1.attn$female, ref = 'no')
y1.attn$race.4level <- relevel(y1.attn$race.4level, ref = 'White')
y1.attn$hisp <- relevel(y1.attn$hisp, ref = 'no')
y1.attn$high.educ <- relevel(y1.attn$high.educ, ref = 'HS Diploma/GED')
y1.attn$married <- relevel(y1.attn$married, ref = 'yes')
y1.attn$household.income <- relevel(y1.attn$household.income, ref = '[>=50K & <100K]')
y1.attn$cash_choice_task <- relevel(y1.attn$cash_choice_task, ref = 'no')

# name specific columns that will be used for prediction

# neurocognitive tests
COG.cols<-c("sst.v","sst.a","sst.Ter","nback.v.avg","nback.v.diff",
            "nback.a.avg","nback.a.diff","nback.Ter.avg","nback.Ter.diff",
            "pea_wiscv_tss","pea_ravlt_total","nihtbx_picvocab_agecorrected",
            "nihtbx_flanker_agecorrected","nihtbx_list_agecorrected",
            "nihtbx_cardsort_agecorrected", "nihtbx_pattern_agecorrected",
            "nihtbx_picture_agecorrected","nihtbx_reading_agecorrected",
            "lmt_scr_perc_correct","cash_choice_task")

# demographic and objective measures
DEM.cols<-c("age","female","race.4level","hisp","high.educ","married",
            "household.income","taylor_ADI","reshist_addr1_p1tot",
            "reshist_addr1_coi_ed_hsgrad", "reshist_addr1_coi_ed_math",
            "reshist_addr1_coi_ed_reading", "reshist_addr1_coi_ed_schpov",
            "reshist_addr1_coi_ed_prxhqece", "reshist_addr1_leadrisk",
            "anthro_waist_cm", "BMI")

# child self-report
SR.cols<-c("neighborhood_crime_y","sch_env","sch_eng","sch_ali",
           "conflict_youth","parent_monit","wkdy_screen","wknd_screen",
           "bis_y_ss_bis_sum", "bis_y_ss_bas_rr", "bis_y_ss_bas_drive",
           "bis_y_ss_bas_fs", "upps_y_ss_negative_urgency",
           "upps_y_ss_lack_of_planning", "upps_y_ss_sensation_seeking",
           "upps_y_ss_positive_urgency","upps_y_ss_lack_of_perseverance")

# child self-report without personality/impulsivity (as construct is similar)
SR.no.imp.cols<-c("neighborhood_crime_y","sch_env","sch_eng","sch_ali",
                  "conflict_youth","parent_monit","wkdy_screen","wknd_screen")

# all predictors
pred.cols<-c(DEM.cols,SR.cols,COG.cols)

# Next we establish the measurement model for the outcome variable. 
# As parent and teacher measures of the childrens' attention problems
# are each affected by specific biases, an ideal "attention problems" 
# measure would identify the common variance between the two raters 
# while modeling rater-specific biases. For this, we use a bifactor 
# model in the structural equation modeling framework. The model 
# has a "General Attention" factor on which both parent and teacher
# rating items load, which is the common factor we are interested in
# predicting, as well as orthogonal parent- and teacher-specific 
# factors, which we estimate only to account for biases unique to the
# specific raters.  

# to prevent family nesting from affecting standard errors, randomly
# sample cases from each family

source("ABCD_family_subsample.R")

y1.SEM<-abcd.fam.subsample(y1.attn)

# columns for measurement model

Attn.cols.final<-c("bpmt_q1", "bpmt_q3", "bpmt_q4", "bpmt_q5",
                   "bpmt_q9", "bpmt_q13", 
                   "cbcl_q04_p","cbcl_q08_p","cbcl_q10_p","cbcl_q41_p",
                   "cbcl_q78_p","cbcl_q93_p","cbcl_q104_p", "cbcl_q01_p",
                   "cbcl_q17_p","cbcl_q61_p" )
cbcl.cols.final<-c("cbcl_q04_p","cbcl_q08_p","cbcl_q10_p","cbcl_q41_p",
                  "cbcl_q78_p","cbcl_q93_p","cbcl_q104_p","cbcl_q13_p",
                  "cbcl_q17_p","cbcl_q45_p","cbcl_q62_p","cbcl_q80_p" )
bpmt.Attention.cols<-c("bpmt_q1","bpmt_q3","bpmt_q4",
                       "bpmt_q5","bpmt_q9","bpmt_q13")

# SEM bifactor model syntax

bf.final <- paste(" Attn =~ ",paste(Attn.cols.final,collapse=" + "),
                 " Parent =~ ",paste(cbcl.cols.final,collapse=" + "),
                 " Teacher =~ ",paste(bpmt.Attention.cols,collapse=" + "),
                 " Attn ~~ 0*Parent 
                      Attn ~~ 0*Teacher
                      Parent ~~ 0*Teacher
                      ",sep="\n")

# fit model

y1.full.model = cfa(bf.final,data=y1.SEM,  
                    std.lv=TRUE,
                    orthogonal=TRUE,
                    estimator = "WLSMV") 

# look at results

summary(y1.full.model,fit=TRUE,standardized=TRUE)

# view schematic of model
semPaths(y1.full.model, "std", edge.label.cex = 0.5, curvePivot = TRUE)

# As shown above, the model displays adequate fit and the standardized 
# factor loadings indicate that the general factor is closely linked to 
# both parent and teacher ratings of attention problems. 

# To avoid overfitting in test data sets (each left out site), we need to 
# fit this outcome measurement model separately in each training fold and 
# then obtain factor scores for the outcome (the general attention factor) 
# for each independent test set using the original model estimated on the 
# training data.


# loop through each training fold

sites<-unique(y1.attn$site_id_l)

training.models<-list()

t<-1

for (s in sites){
  
  t.mod = cfa(bf.final,y1.attn[y1.attn$site_id_l!=s,],  
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  training.models[[t]]<-t.mod 
  
  t<-(t+1)
  
}

names(training.models)<-paste0("lo.",sites)

# add predicted general factor scores from each training set to the data frame

for (s in sites){
  y1.attn[,paste0("lo.",s,".attn")]<-lavPredict(training.models[[paste0("lo.",s)]],
                                                newdata = y1.attn)[,"Attn"]
}

# general attention factors in each fold are close to identical  
fac.cors<-cor(y1.attn[,paste0("lo.",sites,".attn")])

plot_ly(x=sites, y= sites, z = ~fac.cors, type = "heatmap") %>%
  layout(title="Correlations between general attention factor scores 
         estimated in each training fold")

#####################################################################
## Functions for within-fold standardization and imputation #########
#####################################################################

# To avoid overfitting in the test data, we also need to create functions
# that allow A) a composite measurement model, based on z-scores, for one
# of the socio-economic status variables (taylor_ADI), and B) a multiple 
# imputation model to address missing data, to be fit in each training 
# fold and then applied to each independent test site. 

# function for making scaled ADI composite within each fold

gen.ADI.abcd<-function(dat){
  
  #name columns used in composite
  adi.cols<-c("reshist_addr1_adi_edu_h","reshist_addr1_adi_income",
              "reshist_addr1_adi_home_o","reshist_addr1_adi_in_dis",
              "reshist_addr1_adi_unemp","reshist_addr1_adi_pov",
              "reshist_addr1_adi_b138","reshist_addr1_adi_sp",
              "reshist_addr1_adi_ncar")
  
  # z-score
  adi<-dat[,adi.cols]
  means<-apply(adi, 2, function(x) mean(x,na.rm=TRUE))
  sds<-apply(adi, 2, function(x) sd(x,na.rm=TRUE))
  for (c in adi.cols){
    adi[,c]<-(adi[,c]-means[c])/sds[c]
  }
  
  # flip the first 3 reverse-coded items
  adi[,1:3]<-adi[,1:3]*-1
  
  #get summary variable
  
  dat$taylor_ADI<-rowMeans(adi)
  
  attr(dat,"means")<-means
  attr(dat,"sds")<-sds
  
  dat
  
}

#function for applying the training ADI standardization to the test set
test.ADI.abcd<-function(train,test){
  
  #name columns used in composite
  adi.cols<-c("reshist_addr1_adi_edu_h","reshist_addr1_adi_income",
              "reshist_addr1_adi_home_o","reshist_addr1_adi_in_dis",
              "reshist_addr1_adi_unemp","reshist_addr1_adi_pov",
              "reshist_addr1_adi_b138","reshist_addr1_adi_sp",
              "reshist_addr1_adi_ncar")
  
  # z-score
  adi<-test[,adi.cols]
  means<-attr(train,"means")
  sds<-attr(train,"sds")
  for (c in adi.cols){
    adi[,c]<-(adi[,c]-means[c])/sds[c]
  }
  
  # flip the first 3 reverse-coded items
  adi[,1:3]<-adi[,1:3]*-1
  
  #get summary variable
  
  test$taylor_ADI<-rowMeans(adi)
  
  attr(test,"means")<-means
  attr(test,"sds")<-sds
  
  test
  
}

# we will use the mice() function for imputation
# also load add-on mice.reuse() function for imputing new observations 
# in the test data based on the training data imputation model
source("https://raw.githubusercontent.com/prockenschaub/Misc/master/R/mice.reuse/mice.reuse.R")

# make custom function for inputing in each fold:

impute.folds <- function(data,pred.cols,sites){
  
  out<-list()
  
  t<-1
  for (s in sites){
    
    tmp.train<-data[data$site_id_l!=s,]
    tmp.test<-data[data$site_id_l==s,]
    
    tmp.train<-gen.ADI.abcd(dat=tmp.train)
    tmp.test<-test.ADI.abcd(train=tmp.train,test=tmp.test)
    
    tmp.train<-tmp.train[,pred.cols]
    tmp.test<-tmp.test[,pred.cols]
    
    tmp.train.imp <- mice(data = tmp.train,m = 1,print = FALSE)
    tmp.train<-complete(tmp.train.imp)
    
    tmp.test.imp <- suppressWarnings(mice.reuse(tmp.train.imp, 
                                      tmp.test, 
                                      maxit = 1,print = FALSE))
    tmp.test<-tmp.test.imp$`1`
    
    out[[t]]<-list(train=NA,test=NA,imp=NA)
    out[[t]]$train<-tmp.train
    out[[t]]$test<-tmp.test
    out[[t]]$imp<-tmp.train.imp
    
    t<-(t+1)
    
  }
  
  names(out)<-paste0("lo.",sites)
  
  out
  
}
  

# general function to prepare data in each fold
# for modeling:

setup.predicton <- function(fold.data,site,
                            out.data){
  
  #training data
  tmp.dat<-fold.data[[paste0("lo.",site)]]$train
  
  #scale continuous
  s.means<-sapply(tmp.dat[,sapply(tmp.dat,is.numeric)],mean)
  s.sds<-sapply(tmp.dat[,sapply(tmp.dat,is.numeric)],sd)
  for (v in names(s.means)){
    tmp.dat[,v]<-(tmp.dat[,v]-s.means[v])/s.sds[v]}
  
  #dummy variables 
  if(mean(sapply(tmp.dat,is.numeric))<1){
  tmp.dat<-dummy_cols(tmp.dat,remove_first_dummy = TRUE,
                      remove_selected_columns = TRUE)}
  
  # add attention
  tmp.dat$attn<-out.data[out.data$site_id_l!=site,paste0("lo.",site,".attn")]
  
  # test data
  tmp.dat.test<-fold.data[[paste0("lo.",site)]]$test
  
  #scale continuous (same as training)
  for (v in names(s.means)){
    tmp.dat.test[,v]<-(tmp.dat.test[,v]-s.means[v])/s.sds[v]}
  
  # dummy variables
  if(mean(sapply(tmp.dat.test,is.numeric))<1){
  tmp.dat.test<-dummy_cols(tmp.dat.test,remove_first_dummy = TRUE,
                           remove_selected_columns = TRUE)}
  
  out<-list(train=tmp.dat,
            test=tmp.dat.test,
            con.vars=names(s.means))
  
  out
}


# Here we spot check the taylor_ADI composite score functions 
# to make sure they are working properly:

tmp.train<-y1.attn[y1.attn$site_id_l!="site11",]

tmp.test<-y1.attn[y1.attn$site_id_l=="site11",]

tmp.train<-gen.ADI.abcd(dat=tmp.train)

attr(tmp.train,"means")
attr(tmp.train,"sds")
tmp.train$taylor_ADI[1:10]

tmp.train2<-test.ADI.abcd(train=tmp.train,test=tmp.train)

attr(tmp.train2,"means")
attr(tmp.train2,"sds")
tmp.train2$taylor_ADI[1:10]

tmp.test<-test.ADI.abcd(train=tmp.train,test=tmp.test)

attr(tmp.test,"means")
attr(tmp.test,"sds")
tmp.test$taylor_ADI[1:10]

#Here we test the multiple imputation function in an example 
# training fold and assess diagnostics for the success of imputation:

tmp.train<-tmp.train[,pred.cols]

tmp.test<-tmp.test[,pred.cols]

# look at most missing variables
sort(apply(tmp.train,2,function(x) mean(is.na(x))),decreasing = TRUE)[1:20]

# sample imputation model with 3 iterations

tmp.train.imp <- mice(data = tmp.train,
                      m = 3,
                      print = FALSE)

# convergence of chains for key missing variables

plot(tmp.train.imp, c("reshist_addr1_coi_ed_math", "nback.v.avg", "sst.v"))

# density plots to assess how well the imputed data mimics the distributions of real data:

densityplot(tmp.train.imp, ~  reshist_addr1_coi_ed_math)

densityplot(tmp.train.imp, ~  nback.v.avg)

densityplot(tmp.train.imp, ~  sst.v)

densityplot(tmp.train.imp, ~  household.income)

#bivariate plots to asses whether the imputed data replicates bivariate relationships in the observed data:

xyplot(tmp.train.imp,  sst.v ~ pea_wiscv_tss | .imp)

xyplot(tmp.train.imp,  nback.v.avg ~ pea_wiscv_tss  | .imp)

xyplot(tmp.train.imp,   reshist_addr1_coi_ed_math ~ taylor_ADI  | .imp)

xyplot(tmp.train.imp,   household.income ~ taylor_ADI  | .imp)

xyplot(tmp.train.imp,  reshist_addr1_p1tot ~ taylor_ADI  | .imp)


# distribution of data by missingness/propensity score

reshist_addr1_coi_ed_math.na <- is.na(tmp.train$reshist_addr1_coi_ed_math)
sst.v.na <- is.na(tmp.train$sst.v)
nback.v.na <- is.na(tmp.train$nback.v.avg)

train.reshist_addr1_coi_ed_math.na<-complete(tmp.train.imp)
train.reshist_addr1_coi_ed_math.na$reshist_addr1_coi_ed_math<-reshist_addr1_coi_ed_math.na
train.sst.v.na<-complete(tmp.train.imp)
train.sst.v.na$sst.v<-sst.v.na
train.nback.v.na<-complete(tmp.train.imp)
train.nback.v.na$nback.v.avg<-sst.v.na


fit.reshist_addr1_coi_ed_math <- glm(reshist_addr1_coi_ed_math ~., 
                                    train.reshist_addr1_coi_ed_math.na,family = binomial)

fit.sst.v <- glm(sst.v ~., train.sst.v.na,family = binomial)

fit.nback.v <- glm(nback.v.avg ~., train.nback.v.na,family = binomial)

ps.reshist_addr1_coi_ed_math <- rep(predict(fit.reshist_addr1_coi_ed_math, type = "response"),4)
ps.sst.v <- rep(predict(fit.sst.v, type = "response"),4)
ps.nback.v <- rep(predict(fit.nback.v, type = "response"),4)



xyplot(tmp.train.imp, reshist_addr1_coi_ed_math ~ ps.reshist_addr1_coi_ed_math | .imp, pch = c(1, 20), cex = c(0.8, 1.2), 
       xlab = "Probability that reshist_addr1_coi_ed_math is missing", ylab = "reshist_addr1_coi_ed_math")

xyplot(tmp.train.imp, sst.v ~ ps.sst.v  | .imp, pch = c(1, 20), cex = c(0.8, 1.2), 
       xlab = "Probability that sst.v is missing", ylab = "sst.v")

xyplot(tmp.train.imp, sst.v ~ ps.nback.v  | .imp, pch = c(1, 20), cex = c(0.8, 1.2), 
       xlab = "Probability that nback.v.avg is missing", ylab = "nback.v.avg")

# residuals of the relation between propensity score and imputed values:

reshist_addr1_coi_ed_math.comp <- complete(tmp.train.imp, "long", TRUE)$reshist_addr1_coi_ed_math
fit.resid.reshist_addr1_coi_ed_math<-lm(reshist_addr1_coi_ed_math.comp ~ ps.reshist_addr1_coi_ed_math)

sst.v.comp <- complete(tmp.train.imp, "long", TRUE)$sst.v
fit.resid.sst.v<-lm(sst.v.comp ~ ps.sst.v)

nback.v.comp <- complete(tmp.train.imp, "long", TRUE)$nback.v
fit.resid.nback.v<-lm(sst.v.comp ~ ps.nback.v)


densityplot(~residuals(fit.resid.reshist_addr1_coi_ed_math), 
            group = reshist_addr1_coi_ed_math.na, 
            plot.points = FALSE,ref = TRUE, 
            scales = list(y = list(draw = FALSE)), 
            xlab = "Residuals of regression of reshist_addr1_coi_ed_math on propensity score",
            lwd = 2)

densityplot(~residuals(fit.resid.sst.v), group = sst.v.na, plot.points = FALSE,ref = TRUE, scales = list(y = list(draw = FALSE)), xlab = "Residuals of regression of sst.v on propensity score", lwd = 2)

densityplot(~residuals(fit.resid.nback.v), group = nback.v.na, plot.points = FALSE,ref = TRUE, scales = list(y = list(draw = FALSE)), xlab = "Residuals of regression of nback.v on propensity score", lwd = 2)


#apply to test set

tmp.train.imp <- mice(data = tmp.train,m = 1,print = FALSE)

tmp.test.imp <- mice.reuse(tmp.train.imp, tmp.test, maxit = 1,print = FALSE)
tmp.test.imp<-tmp.test.imp$`1`

# These diagnostics indicate that the imputation algorithm converged and was 
# able to replicate several key distributional properties and bivariate 
# relationships present in the observed data. Further, we show that the imputed
# data replicated the relation between the values of several of the most 
# frequently missing variables and their missingness propensity scores. 
# Therefore, we can be confident that the imputed data provide a reasonable 
# approximation of the missing values. 

#######################################################################
## Simple correlations and linear model assumption checks #############
#######################################################################

#correlation matrix and heat map in full data set

y1.attn.fullpred<-gen.ADI.abcd(dat=y1.attn)
y1.attn.fullpred<-y1.attn.fullpred[,c(pred.cols,"Gen.Attn")]

y1.attn.fullpred <- mice(data = y1.attn.fullpred,m = 1,print = FALSE)

y1.attn.fullpred <- complete(y1.attn.fullpred)

all.cors<-cor(model.matrix(~0+., data=y1.attn.fullpred))

plot_ly(x=colnames(all.cors), y= colnames(all.cors), z = ~all.cors, type = "heatmap") %>%
  layout(title="Correlations between all features")

# check variance inflation factor

library(car)

full.model<-lm(Gen.Attn~.,data = y1.attn.fullpred)
summary(full.model)

vif(full.model)

# check normality of residuals 

ols_plot_resid_qq(full.model)

ols_plot_resid_fit(full.model)

#######################################################################
#######################################################################
##### Predictive modeling: full sample ################################
#######################################################################
#######################################################################

#######################################################################
## Create fully imputed training and test data sets for each fold #####
#######################################################################

# Prior to model fitting, we create list objects that contain fully-imputed 
# training data sets and fully-imputed test data sets for 18 iterations of 
# leave-one-site-out cross validation process (3 of the 21 ABCD sites are 
# left out of these analyses because they were pseudo-randomly drawn to be 
# used for final validations of the models prior to publication).

full.imputed<-impute.folds(y1.attn,
                           pred.cols = pred.cols,
                           sites = sites)
save(full.imputed,file="data_imputations_full.RData")


#######################################################################
### Sparse prediction models using individual variables ###############
#######################################################################

# loop to conduct internal cross-validation with ridge, lasso and 
# SVR (both linear and redial basis function kernels) and then 
# test model on left-out site

load("data_imputations_full.RData")

ridge.fits<-list()
ridge.preds<-list()
ridge.train.r<-c()
ridge.test.r<-c()

lasso.fits<-list()
lasso.preds<-list()
lasso.train.r<-c()
lasso.test.r<-c()

svr.linear.fits<-list()
svr.linear.preds<-list()
svr.linear.train.r<-c()
svr.linear.test.r<-c()

svr.rbf.fits<-list()
svr.rbf.preds<-list()
svr.rbf.train.r<-c()
svr.rbf.test.r<-c()

for (s in sites){
  
  tmp.dat<-setup.predicton(full.imputed,s,y1.attn)
  
  ridge.fits[[s]]<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 0)
  ridge.preds[[s]]$train<-predict(ridge.fits[[s]],data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = ridge.fits[[s]]$lambda.1se)
  ridge.preds[[s]]$test<-predict(ridge.fits[[s]],data.matrix(tmp.dat$test), s = ridge.fits[[s]]$lambda.1se)
  ridge.train.r[[s]]<-cor(ridge.preds[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  ridge.test.r[[s]]<-cor(ridge.preds[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  lasso.fits[[s]]<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  lasso.preds[[s]]$train<-predict(lasso.fits[[s]],data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = lasso.fits[[s]]$lambda.1se)
  lasso.preds[[s]]$test<-predict(lasso.fits[[s]],data.matrix(tmp.dat$test), s = lasso.fits[[s]]$lambda.1se)
  lasso.train.r[[s]]<-cor(lasso.preds[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  lasso.test.r[[s]]<-cor(lasso.preds[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  svr.linear.fits[[s]]<-svm(attn~.,
                        data=tmp.dat$train,
                        kernel="linear")
  svr.linear.preds[[s]]$train <- predict(svr.linear.fits[[s]],tmp.dat$train)
  svr.linear.preds[[s]]$test <- predict(svr.linear.fits[[s]],tmp.dat$test)
  svr.linear.train.r[[s]]<-cor(svr.linear.preds[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  svr.linear.test.r[[s]]<-cor(svr.linear.preds[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])

  svr.rbf.fits[[s]]<-svm(attn~.,
                            data=tmp.dat$train,
                            kernel="radial")
  svr.rbf.preds[[s]]$train <- predict(svr.rbf.fits[[s]],tmp.dat$train)
  svr.rbf.preds[[s]]$test <- predict(svr.rbf.fits[[s]],tmp.dat$test)
  svr.rbf.train.r[[s]]<-cor(svr.rbf.preds[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  svr.rbf.test.r[[s]]<-cor(svr.rbf.preds[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
}


ridge.train.r<-as.numeric(ridge.train.r)
ridge.test.r<-as.numeric(ridge.test.r)
lasso.train.r<-as.numeric(lasso.train.r)
lasso.test.r<-as.numeric(lasso.test.r)
svr.linear.train.r<-as.numeric(svr.linear.train.r)
svr.linear.test.r<-as.numeric(svr.linear.test.r)
svr.rbf.train.r<-as.numeric(svr.rbf.train.r)
svr.rbf.test.r<-as.numeric(svr.rbf.test.r)

 save.image(file="sparse_models_6_17_22.RData")
 load("sparse_models_6_17_22.RData")

# little loss in the test set, except for SVR with RBF kernal

mean(ridge.train.r)
mean(ridge.test.r)
mean(ridge.test.r)^2

mean(lasso.train.r)
mean(lasso.test.r)
mean(lasso.test.r)^2

mean(svr.linear.train.r)
mean(svr.linear.test.r)
mean(svr.linear.test.r)^2

mean(svr.rbf.train.r)
mean(svr.rbf.test.r)
mean(svr.rbf.test.r)^2

#plot out

plot.pred<-function(sites,pred){
  
  col=rgb(.75,0,0,.6)
  par(mfrow=c(2,3))
  
  for (s in sites){
    
    actual<-y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")]
    predicted<-pred[[s]]$test
    
    plot(actual,predicted,
         xlab="actual",ylab="predicted",
         main=s,
         col=col,cex.axis=1.25,cex.lab=1.5,cex.main=1.7,pch=16,cex=1.2)
    abline(lm(predicted~actual),lwd=2)
    
  }
  
}

site.ns<-c()
for (s in sites){site.ns[s]<-length(y1.attn[y1.attn$site_id_l==s,"subjectkey"])}

high.n<-names(sort(site.ns,decreasing=TRUE))[1:6]

plot.pred(high.n,ridge.preds)
plot.pred(high.n,lasso.preds)
plot.pred(high.n,svr.rbf.preds)


# look at importance of variables

betaHatRIDGE<-as.data.frame(matrix(NA,nrow=length(sites),ncol = length(rownames(coef(ridge.fits[["site11"]], s = ridge.fits[["site11"]]$lambda.1se)))))
colnames(betaHatRIDGE)<-rownames(coef(ridge.fits[["site11"]], s = ridge.fits[["site11"]]$lambda.1se))
betaHatRIDGE$sites<-sites

betaHatLASSO<-betaHatRIDGE


for (s in sites){
  coefs<-coef(ridge.fits[[s]], s = ridge.fits[[s]]$lambda.1se)
  betaHatRIDGE[betaHatRIDGE$sites==s,1:(length(betaHatRIDGE)-1)]<-c(coefs[,1])
  coefs<-coef(lasso.fits[[s]], s = lasso.fits[[s]]$lambda.1se)
  betaHatLASSO[betaHatLASSO$sites==s,1:(length(betaHatLASSO)-1)]<-c(coefs[,1])
}


betaHatRIDGE.avg<-colMeans(betaHatRIDGE[,colnames(betaHatRIDGE)!="sites"])
betaHatLASSO.avg<-colMeans(betaHatLASSO[,colnames(betaHatLASSO)!="sites"])

sort(abs(betaHatRIDGE.avg),decreasing = TRUE)[1:20]
sort(abs(betaHatLASSO.avg),decreasing = TRUE)[1:20]

#which variables are included by LASSO in all folds?

betaHatLASSO.bin<-betaHatLASSO[,colnames(betaHatLASSO)!="sites"]!=0
betaHatLASSO.bin<-colMeans(betaHatLASSO.bin)

included<-names(betaHatLASSO.bin[betaHatLASSO.bin>.50])

write.csv(data.frame(beta=betaHatLASSO.avg[included]),file="lasso_betas.csv")

#######################################################################
## Principal components analysis and regression #######################
#######################################################################

# conduct and interpret PCA in full sample

# start by imputing in full sample (minus outcome variable) 

y1.attn.PCA<-gen.ADI.abcd(dat=y1.attn)
y1.attn.PCA<-y1.attn.PCA[,pred.cols]

y1.attn.PCA <- mice(data = y1.attn.PCA,m = 1,print = FALSE)

y1.attn.PCA <- complete(y1.attn.PCA)

# make list of non-categorical variables for PCA entry

num.preds<-colnames(y1.attn.PCA)[sapply(y1.attn.PCA,is.numeric)]
cat.preds<-colnames(y1.attn.PCA)[!sapply(y1.attn.PCA,is.numeric)]

# complete PCA

general.pca<-PCA(scale(y1.attn.PCA[,num.preds]),ncp=(length(num.preds)-1))

# interpret general pca

summary(general.pca)
# Eigenvalues
# Dim.1   Dim.2   Dim.3   Dim.4   Dim.5   Dim.6
# Variance               6.051   3.020   2.881   2.412   2.107   1.784
# % of var.             13.154   6.565   6.264   5.243   4.581   3.878
# Cumulative % of var.  13.154  19.719  25.983  31.226  35.807  39.685
# Dim.7   Dim.8   Dim.9  Dim.10  Dim.11  Dim.12
# Variance               1.490   1.317   1.267   1.210   1.168   1.105
# % of var.              3.239   2.864   2.753   2.630   2.539   2.403
# Cumulative % of var.  42.924  45.788  48.541  51.171  53.710  56.113
# Dim.13  Dim.14  Dim.15  Dim.16  Dim.17  Dim.18
# Variance               1.064   1.050   0.961   0.907   0.889   0.879
# % of var.              2.313   2.282   2.088   1.971   1.934   1.910
# Cumulative % of var.  58.426  60.708  62.796  64.767  66.700  68.610
# 

head(general.pca$ind$coord)

# function for PCR regression and tuning

PCA.reg<-function(outcome, PCA.cols, data, folds=10,num.PCs=NA){
  
  if(folds>length(data[,1])){
    stop("number of folds greater than number of observations")}
  
  out<-list()
  
  pred.names<- colnames(data)[colnames(data)!=outcome]
  num.preds<-colnames(data)[colnames(data)%in%PCA.cols & 
                              colnames(data)%in%pred.names]
  cat.preds<-colnames(data)[!colnames(data)%in%PCA.cols& 
                              colnames(data)%in%pred.names]
  
  out$pca<-PCA(scale(data[,num.preds]),
               ncp=(length(num.preds)-1),graph=FALSE)
  
  scaled.comps<-scale(out$pca$ind$coord)
  
  if (is.na(num.PCs)){
    
  cv.vec<-sample(rep(1:folds,(length(data[,1])/folds)+1)[1:length(data[,1])])
  
  nc.range<-1:(length(num.preds)-1)
  
  cv.dat<-expand.grid(nc=nc.range,test.set=1:folds)
  cv.dat$test.r<-NA
  
  for (f in 1:length(cv.dat$nc)){
    train.dat<-data[cv.vec!=cv.dat$test.set[f],c(outcome,cat.preds),drop=FALSE]
    train.dat<-cbind(train.dat,
                     scaled.comps[cv.vec!=cv.dat$test.set[f],
                                  1:cv.dat$nc[f],drop=FALSE])
    test.dat<-data[cv.vec==cv.dat$test.set[f],c(outcome,cat.preds),drop=FALSE]
    test.dat<-cbind(test.dat,
                    scaled.comps[cv.vec==cv.dat$test.set[f],
                                 1:cv.dat$nc[f],drop=FALSE])
    tmp.mod<-lm(paste0(outcome,"~."),train.dat)
    cv.dat$test.r[f]<-cor(predict(tmp.mod,newdata = test.dat),
                          test.dat[,outcome])
  }
  
  out$cv.dat<-cv.dat
  out$cv.means<-tapply(cv.dat$test.r,cv.dat$nc,mean)
  
  max.round<-out$cv.means[round(out$cv.means,2)==max(round(out$cv.means,2))]
  out$winning.nc<-min(as.numeric(names(max.round)))
  
  } else {out$winning.nc<-num.PCs}

  #estimate winnings model
  out$winning.model<-lm(paste0(outcome,"~."),
                        cbind(data[,c(outcome,cat.preds),drop=FALSE],
                              scaled.comps[,1:out$winning.nc,drop=FALSE]))
  
  # project feature weights
  beta.mat<-out$winning.model$coefficients
  beta.mat<-as.matrix(beta.mat[paste0("Dim.",1:out$winning.nc)])
  
  out$fw<-out$pca$var$coord[,1:out$winning.nc]%*%beta.mat
  
  out
  
}

# use PCA regression model to predict in new data

PCA.predict <-function(pcr,newdata){
  
  vars<-row.names(pcr$pca$var$coord)
  
  # component scores in new data
  new.comps<-matrix(dim(as.matrix(pcr$pca$var$coord[,1:pcr$winning.nc])))
  
  loadings<-as.matrix(pcr$pca$var$coord[,1:pcr$winning.nc,drop=FALSE])
  
  for(l in 1:pcr$winning.nc){
    loadings[,l]<-loadings[,l]/sqrt(pcr$pca$eig[l,1])
  }
  
  scores<-as.data.frame(t(t(loadings)%*%t(as.matrix(newdata[,vars]))))
  
  # standardize with reference to training data components
  s.means<-sapply(data.frame(pcr$pca$ind$coord[,1:pcr$winning.nc,drop=FALSE]),mean)
  s.sds<-sapply(data.frame(pcr$pca$ind$coord[,1:pcr$winning.nc,drop=FALSE]),sd)
  for (v in names(s.means)){
    scores[,v]<-(scores[,v]-s.means[v])/s.sds[v]}
  
  
  #predict with new data and new component scores
  new<-cbind(scores,newdata[,!colnames(newdata)%in%vars,drop=FALSE])
  out<-predict(pcr$winning.model,new)
  
}

# run PCA regression

pcareg.fits<-list()
pcareg.preds<-list()
pcareg.train.r<-c()
pcareg.test.r<-c()

for (s in sites){

  tmp.dat<-setup.predicton(full.imputed,s,y1.attn)

  pcareg.fits[[s]]<-PCA.reg(data=tmp.dat$train, outcome = "attn", 
                            PCA.cols = tmp.dat$con.vars)
  pcareg.preds[[s]]$train<-predict(pcareg.fits[[s]]$winning.model)
  pcareg.preds[[s]]$test<-PCA.predict(pcareg.fits[[s]],tmp.dat$test)
  pcareg.train.r[[s]]<-cor(pcareg.preds[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  pcareg.test.r[[s]]<-cor(pcareg.preds[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
}

pca.train.r<-as.numeric(pcareg.train.r)
pca.test.r<-as.numeric(pcareg.test.r)

# feature weights

fw.pca<-data.frame(matrix(NA,length(sites),length(pcareg.fits[[s]]$fw)))
colnames(fw.pca)<-rownames(pcareg.fits[[s]]$fw)
rownames(fw.pca)<-sites

cpw.pca<-data.frame(matrix(NA,length(sites),13))
colnames(cpw.pca)<-names(pcareg.fits[[s]]$winning.model$coefficients[2:14])
rownames(cpw.pca)<-sites

for (s in sites){
  fw.pca[s,]<-pcareg.fits[[s]]$fw
  cpw.pca[s,]<- pcareg.fits[[s]]$winning.model$coefficients[2:14]
}

fw.means<-colMeans(fw.pca)
fw.means[order(abs(fw.means),decreasing = TRUE)]

cpw.means<-colMeans(cpw.pca)
cpw.means[order(abs(cpw.means),decreasing = TRUE)]

save.image(file="pca_regression_models_6_20_22.RData")

# interpret predictive components

pca.interp<-data.frame(lo=sites)

pca.interp$npc<-NA

for (r in 1:length(pca.interp$lo)){
  pca.interp$npc[r]<-pcareg.fits[[sites[r]]]$winning.nc
}

mean(pca.interp$npc)
#[1] 9.5
max(pca.interp$npc)
#[1] 26

pca.interp[,paste("comp",1:26)]<-0

for (r in 1:length(pca.interp$lo)){
  tmp.coefs<-pcareg.fits[[sites[r]]]$winning.model$coefficients
  tmp.coefs<-tmp.coefs[paste0("Dim.",1:pca.interp$npc[r])]
  pca.interp[r,paste("comp",1:pca.interp$npc[r])]<-tmp.coefs
}

plot(colMeans(pca.interp[,paste("comp",1:26)]),
     pch=16,col="darkblue",cex=1.2,ylim=c(c(-.35,.10)),
     ylab="standardized beta",xlab="component")
# lines(colMeans(pca.interp[,paste("comp",1:26)]),
#       col="dodgerblue",lwd=4)
lines(rep(0,26),col="gray",lwd=2,lty=3)
for(r in 1:length(pca.interp$lo)){
  lines(unlist(pca.interp[r,paste("comp",1:26)]),
        col="dodgerblue")
}

########################################################
############## summarize prediction performance ########
########################################################

# Figure 1A:
jpeg("attn_prediction_by_method.jpg",width = 7,height = 5.5,
     units = "in",res = 500)
par(mfrow=c(1,1))
plot(site.ns,lasso.test.r,pch=16,col=rgb(.25,.51,.81,1),cex=1.5,
     xlab="site sample size",ylab="test correlation",
     ylim=c(0,.6),xlim=c(0,760))
points(site.ns,pca.test.r,pch=17,col=rgb(.57,.82,.31,1),cex=1.3)
rect(720,-1,1000,1,col = "lightgray")
ci.lasso<-qt(0.975,df=17)*sd(lasso.test.r)/sqrt(18)
ci.pca<-qt(0.975,df=17)*sd(pca.test.r)/sqrt(18)
arrows(x0=740, y0=mean(lasso.test.r)-ci.lasso, 
       x1=740, y1=mean(lasso.test.r)+ci.lasso, code=3, 
       angle=90, length=0.06, lwd=1)
arrows(x0=770, y0=mean(pca.test.r)-ci.pca, 
       x1=770, y1=mean(pca.test.r)+ci.pca, code=3, 
       angle=90, length=0.06, lwd=1)
points(740,mean(lasso.test.r),pch=16,col=rgb(.25,.51,.81,1),cex=1.5)
points(770,mean(pca.test.r),pch=17,col=rgb(.57,.82,.31,1),cex=1.3)
dev.off()

jpeg("prediction_top_6_sites.jpg",width = 7.7,height = 4.5,
     units = "in",res = 500)
col=rgb(.75,0,0,.2)
par(mfrow=c(2,3),mar=c(2,4,2,2))
for (s in high.n){
  actual<-y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")]
  predicted<-pcareg.preds[[s]]$test
  plot(actual,predicted,
       xlab="",ylab="",
       main="",xlim=c(-1.5,3.6),ylim=c(-1.4,1.4),
       col=col,cex.axis=1.25,cex.lab=1.5,cex.main=1.7,pch=16,cex=1.5)
  abline(lm(predicted~actual),lwd=1.5)
}
dev.off()

mean(pca.test.r)
mean(pca.test.r)^2
mean(pca.test.r)-ci.pca
mean(pca.test.r)+ci.pca

mean(lasso.test.r)
mean(lasso.test.r)^2
mean(lasso.test.r)-ci.lasso
mean(lasso.test.r)+ci.lasso


mean(pca.train.r)
mean(pca.train.r)^2

mean(lasso.train.r)
mean(lasso.train.r)^2

# see if predictions and prediction errors are correlated 
# across PCR and lasso

jpeg("attn_prediction_by_method.jpg",width = 8,height = 12,
     units = "in",res = 300)
par(mfrow=c(5,4))
for (s in sites){
plot(pcareg.preds[[s]]$test,lasso.preds[[s]]$test,
     xlab=paste0(s," PCR pred."),ylab=paste0(s," LASSO pred."),
     main = paste0("r=",
                   round(cor(pcareg.preds[[s]]$test,lasso.preds[[s]]$test),2)))
}
dev.off()

jpeg("attn_prediction_error_by_method.jpg",width = 8,height = 12,
     units = "in",res = 300)
par(mfrow=c(5,4))
for (s in sites){
  PCR.err<-pcareg.preds[[s]]$test-y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")]
  LASSO.err<-lasso.preds[[s]]$test-y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")]
  plot(PCR.err,LASSO.err,
       xlab=paste0(s," PCR error"),ylab=paste0(s," LASSO error"),
       main = paste0("r=",
                     round(cor(PCR.err,LASSO.err),2)))
}
dev.off()
  
#########################################################
#########################################################
##### Comparison of variable domains ####################
#########################################################
#########################################################

#######################################################################
## Create fully imputed training and test data sets for each fold #####
#######################################################################

# impute training/test data with and without variables 
# from each of the three domains

COG.imputed<-impute.folds(y1.attn,
                           pred.cols = COG.cols,
                           sites = sites)

SR.imputed<-impute.folds(y1.attn,
                          pred.cols = SR.cols,
                          sites = sites)

DEM.imputed<-impute.folds(y1.attn,
                         pred.cols = DEM.cols,
                         sites = sites)

COG_SR.imputed<-impute.folds(y1.attn,
                             pred.cols = c(COG.cols,SR.cols),
                             sites = sites)

COG_DEM.imputed<-impute.folds(y1.attn,
                          pred.cols = c(COG.cols,DEM.cols),
                          sites = sites)

SR_DEM.imputed<-impute.folds(y1.attn,
                              pred.cols = c(SR.cols,DEM.cols),
                              sites = sites)

 
 save(COG.imputed,SR.imputed,DEM.imputed,
      COG_SR.imputed,COG_DEM.imputed,SR_DEM.imputed,
      file="data_imputations_domain_analysis.RData")
#load("data_imputations_cog_analysis.RData")

SR.no.imp.imputed<-impute.folds(y1.attn,
                         pred.cols = SR.no.imp.cols,
                         sites = sites)
save(SR.no.imp.imputed,
     file="data_imputations_SR_analysis.RData")


# run LASSO and PCA regression

COG.lasso<-list()
SR.lasso<-list()
DEM.lasso<-list()
COG_SR.lasso<-list()
COG_DEM.lasso<-list()
SR_DEM.lasso<-list()
SR.no.imp.lasso<-list()

COG.pcareg<-list()
SR.pcareg<-list()
DEM.pcareg<-list()
COG_SR.pcareg<-list()
COG_DEM.pcareg<-list()
SR_DEM.pcareg<-list()
SR.no.imp.pcareg<-list()

for (s in sites){

  tmp.dat<-setup.predicton(COG.imputed,s,y1.attn)
  COG.lasso[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  COG.lasso[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  COG.lasso[[s]]$train<-predict(COG.lasso[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = COG.lasso[[s]]$fits$lambda.1se)
  COG.lasso[[s]]$test<-predict(COG.lasso[[s]]$fits,data.matrix(tmp.dat$test), s = COG.lasso[[s]]$fits$lambda.1se)
  COG.lasso[[s]]$train.r<-cor(COG.lasso[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG.lasso[[s]]$test.r<-cor(COG.lasso[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  COG.pcareg[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  COG.pcareg[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  COG.pcareg[[s]]$train<-predict(COG.pcareg[[s]]$fits$winning.model)
  COG.pcareg[[s]]$test<-PCA.predict(COG.pcareg[[s]]$fits,tmp.dat$test)
  COG.pcareg[[s]]$train.r<-cor(COG.pcareg[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG.pcareg[[s]]$test.r<-cor(COG.pcareg[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(SR.imputed,s,y1.attn)
  SR.lasso[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  SR.lasso[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  SR.lasso[[s]]$train<-predict(SR.lasso[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = SR.lasso[[s]]$fits$lambda.1se)
  SR.lasso[[s]]$test<-predict(SR.lasso[[s]]$fits,data.matrix(tmp.dat$test), s = SR.lasso[[s]]$fits$lambda.1se)
  SR.lasso[[s]]$train.r<-cor(SR.lasso[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  SR.lasso[[s]]$test.r<-cor(SR.lasso[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  SR.pcareg[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  SR.pcareg[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  SR.pcareg[[s]]$train<-predict(SR.pcareg[[s]]$fits$winning.model)
  SR.pcareg[[s]]$test<-PCA.predict(SR.pcareg[[s]]$fits,tmp.dat$test)
  SR.pcareg[[s]]$train.r<-cor(SR.pcareg[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  SR.pcareg[[s]]$test.r<-cor(SR.pcareg[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
 
  tmp.dat<-setup.predicton(DEM.imputed,s,y1.attn)
  DEM.lasso[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  DEM.lasso[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  DEM.lasso[[s]]$train<-predict(DEM.lasso[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = DEM.lasso[[s]]$fits$lambda.1se)
  DEM.lasso[[s]]$test<-predict(DEM.lasso[[s]]$fits,data.matrix(tmp.dat$test), s = DEM.lasso[[s]]$fits$lambda.1se)
  DEM.lasso[[s]]$train.r<-cor(DEM.lasso[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  DEM.lasso[[s]]$test.r<-cor(DEM.lasso[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  DEM.pcareg[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  DEM.pcareg[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  DEM.pcareg[[s]]$train<-predict(DEM.pcareg[[s]]$fits$winning.model)
  DEM.pcareg[[s]]$test<-PCA.predict(DEM.pcareg[[s]]$fits,tmp.dat$test)
  DEM.pcareg[[s]]$train.r<-cor(DEM.pcareg[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  DEM.pcareg[[s]]$test.r<-cor(DEM.pcareg[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(COG_SR.imputed,s,y1.attn)
  COG_SR.lasso[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  COG_SR.lasso[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  COG_SR.lasso[[s]]$train<-predict(COG_SR.lasso[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = COG_SR.lasso[[s]]$fits$lambda.1se)
  COG_SR.lasso[[s]]$test<-predict(COG_SR.lasso[[s]]$fits,data.matrix(tmp.dat$test), s = COG_SR.lasso[[s]]$fits$lambda.1se)
  COG_SR.lasso[[s]]$train.r<-cor(COG_SR.lasso[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_SR.lasso[[s]]$test.r<-cor(COG_SR.lasso[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  COG_SR.pcareg[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  COG_SR.pcareg[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  COG_SR.pcareg[[s]]$train<-predict(COG_SR.pcareg[[s]]$fits$winning.model)
  COG_SR.pcareg[[s]]$test<-PCA.predict(COG_SR.pcareg[[s]]$fits,tmp.dat$test)
  COG_SR.pcareg[[s]]$train.r<-cor(COG_SR.pcareg[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_SR.pcareg[[s]]$test.r<-cor(COG_SR.pcareg[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(COG_DEM.imputed,s,y1.attn)
  COG_DEM.lasso[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  COG_DEM.lasso[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  COG_DEM.lasso[[s]]$train<-predict(COG_DEM.lasso[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = COG_DEM.lasso[[s]]$fits$lambda.1se)
  COG_DEM.lasso[[s]]$test<-predict(COG_DEM.lasso[[s]]$fits,data.matrix(tmp.dat$test), s = COG_DEM.lasso[[s]]$fits$lambda.1se)
  COG_DEM.lasso[[s]]$train.r<-cor(COG_DEM.lasso[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_DEM.lasso[[s]]$test.r<-cor(COG_DEM.lasso[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  COG_DEM.pcareg[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  COG_DEM.pcareg[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  COG_DEM.pcareg[[s]]$train<-predict(COG_DEM.pcareg[[s]]$fits$winning.model)
  COG_DEM.pcareg[[s]]$test<-PCA.predict(COG_DEM.pcareg[[s]]$fits,tmp.dat$test)
  COG_DEM.pcareg[[s]]$train.r<-cor(COG_DEM.pcareg[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_DEM.pcareg[[s]]$test.r<-cor(COG_DEM.pcareg[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(SR_DEM.imputed,s,y1.attn)
  SR_DEM.lasso[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  SR_DEM.lasso[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  SR_DEM.lasso[[s]]$train<-predict(SR_DEM.lasso[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = SR_DEM.lasso[[s]]$fits$lambda.1se)
  SR_DEM.lasso[[s]]$test<-predict(SR_DEM.lasso[[s]]$fits,data.matrix(tmp.dat$test), s = SR_DEM.lasso[[s]]$fits$lambda.1se)
  SR_DEM.lasso[[s]]$train.r<-cor(SR_DEM.lasso[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  SR_DEM.lasso[[s]]$test.r<-cor(SR_DEM.lasso[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  SR_DEM.pcareg[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  SR_DEM.pcareg[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  SR_DEM.pcareg[[s]]$train<-predict(SR_DEM.pcareg[[s]]$fits$winning.model)
  SR_DEM.pcareg[[s]]$test<-PCA.predict(SR_DEM.pcareg[[s]]$fits,tmp.dat$test)
  SR_DEM.pcareg[[s]]$train.r<-cor(SR_DEM.pcareg[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  SR_DEM.pcareg[[s]]$test.r<-cor(SR_DEM.pcareg[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(SR.no.imp.imputed,s,y1.attn)
  SR.no.imp.lasso[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  SR.no.imp.lasso[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  SR.no.imp.lasso[[s]]$train<-predict(SR.no.imp.lasso[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = SR.no.imp.lasso[[s]]$fits$lambda.1se)
  SR.no.imp.lasso[[s]]$test<-predict(SR.no.imp.lasso[[s]]$fits,data.matrix(tmp.dat$test), s = SR.no.imp.lasso[[s]]$fits$lambda.1se)
  SR.no.imp.lasso[[s]]$train.r<-cor(SR.no.imp.lasso[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  SR.no.imp.lasso[[s]]$test.r<-cor(SR.no.imp.lasso[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  SR.no.imp.pcareg[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  SR.no.imp.pcareg[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  SR.no.imp.pcareg[[s]]$train<-predict(SR.no.imp.pcareg[[s]]$fits$winning.model)
  SR.no.imp.pcareg[[s]]$test<-PCA.predict(SR.no.imp.pcareg[[s]]$fits,tmp.dat$test)
  SR.no.imp.pcareg[[s]]$train.r<-cor(SR.no.imp.pcareg[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  SR.no.imp.pcareg[[s]]$test.r<-cor(SR.no.imp.pcareg[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
}

save.image("domain_analyses_6_21.RData")

domain.r<-data.frame(COG.train=unlist(lapply(COG.pcareg,FUN=function(x) x$train.r)),
           SR.train=unlist(lapply(SR.pcareg,FUN=function(x) x$train.r)),
           DEM.train=unlist(lapply(DEM.pcareg,FUN=function(x) x$train.r)),
           COG_SR.train=unlist(lapply(COG_SR.pcareg,FUN=function(x) x$train.r)),
           COG_DEM.train=unlist(lapply(COG_DEM.pcareg,FUN=function(x) x$train.r)),
           SR_DEM.train=unlist(lapply(SR_DEM.pcareg,FUN=function(x) x$train.r)),
           COG.test=unlist(lapply(COG.pcareg,FUN=function(x) x$test.r)),
           SR.test=unlist(lapply(SR.pcareg,FUN=function(x) x$test.r)),
           DEM.test=unlist(lapply(DEM.pcareg,FUN=function(x) x$test.r)),
           COG_SR.test=unlist(lapply(COG_SR.pcareg,FUN=function(x) x$test.r)),
           COG_DEM.test=unlist(lapply(COG_DEM.pcareg,FUN=function(x) x$test.r)),
           SR_DEM.test=unlist(lapply(SR_DEM.pcareg,FUN=function(x) x$test.r)) )

colMeans(domain.r)
# COG.train      SR.train     DEM.train  COG_SR.train COG_DEM.train 
# 0.3174459     0.3621923     0.2948020     0.4232799     0.3946140 
# SR_DEM.train      COG.test       SR.test      DEM.test   COG_SR.test 
# 0.4110501     0.3056797     0.3575823     0.2748801     0.4190918 
# COG_DEM.test   SR_DEM.test 
# 0.3849605     0.3919832 

domain.ci<- apply(domain.r,2,
      FUN= function(x) qt(0.975,df=17)*sd(x)/sqrt(18))

colMeans(domain.r)+domain.ci
colMeans(domain.r)-domain.ci

mean(unlist(lapply(SR.no.imp.pcareg,FUN=function(x) x$train.r)))
#[1] 0.2781509
mean(unlist(lapply(SR.no.imp.pcareg,FUN=function(x) x$test.r)))
#[1] 0.2844319


domain.r.lasso<-data.frame(COG.train=unlist(lapply(COG.lasso,FUN=function(x) x$train.r)),
                     SR.train=unlist(lapply(SR.lasso,FUN=function(x) x$train.r)),
                     DEM.train=unlist(lapply(DEM.lasso,FUN=function(x) x$train.r)),
                     COG_SR.train=unlist(lapply(COG_SR.lasso,FUN=function(x) x$train.r)),
                     COG_DEM.train=unlist(lapply(COG_DEM.lasso,FUN=function(x) x$train.r)),
                     SR_DEM.train=unlist(lapply(SR_DEM.lasso,FUN=function(x) x$train.r)),
                     COG.test=unlist(lapply(COG.lasso,FUN=function(x) x$test.r)),
                     SR.test=unlist(lapply(SR.lasso,FUN=function(x) x$test.r)),
                     DEM.test=unlist(lapply(DEM.lasso,FUN=function(x) x$test.r)),
                     COG_SR.test=unlist(lapply(COG_SR.lasso,FUN=function(x) x$test.r)),
                     COG_DEM.test=unlist(lapply(COG_DEM.lasso,FUN=function(x) x$test.r)),
                     SR_DEM.test=unlist(lapply(SR_DEM.lasso,FUN=function(x) x$test.r)) )

colMeans(domain.r.lasso)

# COG.train      SR.train     DEM.train  COG_SR.train COG_DEM.train 
# 0.3117391     0.3541150     0.2804229     0.4204766     0.3855721 
# SR_DEM.train      COG.test       SR.test      DEM.test   COG_SR.test 
# 0.4006495     0.2960681     0.3493571     0.2720220     0.4104044 
# COG_DEM.test   SR_DEM.test 
# 0.3787966     0.3972994 

domain.ci.lasso<- apply(domain.r.lasso,2,
                  FUN= function(x) qt(0.975,df=17)*sd(x)/sqrt(18))

colMeans(domain.r.lasso)+domain.ci.lasso
colMeans(domain.r.lasso)-domain.ci.lasso

mean(unlist(lapply(SR.no.imp.lasso,FUN=function(x) x$train.r)))
#[1] 0.2784817
mean(unlist(lapply(SR.no.imp.lasso,FUN=function(x) x$test.r)))
#[1] 0.2833253



save.image(file="domain_analyses_6_2_22.RData")


#########################################################
#########################################################
##### LASSO predictor stability and efficacy ############
#########################################################
#########################################################

# split-half to identify whether predictors are stable

A.sites<-sample(sites,9,FALSE)
B.sites<-sites[!sites%in%A.sites]
#save(A.sites,B.sites,file="split_half_sites.RData")

y1.attn.A<-gen.ADI.abcd(y1.attn[y1.attn$site_id_l%in%A.sites,])
y1.attn.B<-gen.ADI.abcd(y1.attn[y1.attn$site_id_l%in%B.sites,])

y1.attn.A<-y1.attn.A[,pred.cols]
y1.attn.B<-y1.attn.B[,pred.cols]

# impute each

y1.attn.A <- mice(data = y1.attn.A,m = 1,print = FALSE)
y1.attn.A <- complete(y1.attn.A)

y1.attn.B <- mice(data = y1.attn.B,m = 1,print = FALSE)
y1.attn.B <- complete(y1.attn.B)

# set up for LASSO

s.means<-sapply(y1.attn.A[,sapply(y1.attn.A,is.numeric)],mean)
s.sds<-sapply(y1.attn.A[,sapply(y1.attn.A,is.numeric)],sd)
for (v in names(s.means)){
  y1.attn.A[,v]<-(y1.attn.A[,v]-s.means[v])/s.sds[v]}
y1.attn.A<-dummy_cols(y1.attn.A,remove_first_dummy = TRUE,
                      remove_selected_columns = TRUE)

s.means<-sapply(y1.attn.B[,sapply(y1.attn.B,is.numeric)],mean)
s.sds<-sapply(y1.attn.B[,sapply(y1.attn.B,is.numeric)],sd)
for (v in names(s.means)){
  y1.attn.B[,v]<-(y1.attn.B[,v]-s.means[v])/s.sds[v]}
y1.attn.B<-dummy_cols(y1.attn.B,remove_first_dummy = TRUE,
                      remove_selected_columns = TRUE)

# generate attention measure in each half

t.mod = cfa(bf.final,y1.attn[y1.attn$site_id_l%in%A.sites,],  
            std.lv=TRUE,
            orthogonal=TRUE,
            estimator = "WLSMV")
y1.attn$attn.A<-lavPredict(t.mod,newdata = y1.attn)[,"Attn"]

t.mod = cfa(bf.final,y1.attn[y1.attn$site_id_l%in%B.sites,],  
            std.lv=TRUE,
            orthogonal=TRUE,
            estimator = "WLSMV")
y1.attn$attn.B<-lavPredict(t.mod,newdata = y1.attn)[,"Attn"]


# run LASSO 

lasso.A<-cv.glmnet(data.matrix(y1.attn.A), 
                   y1.attn[y1.attn$site_id_l%in%A.sites,"attn.A"],
                   alpha = 1)

lasso.B<-cv.glmnet(data.matrix(y1.attn.B), 
                   y1.attn[y1.attn$site_id_l%in%B.sites,"attn.B"],
                   alpha = 1)

# features
betaHat.A<-coef(lasso.A, s = lasso.A$lambda.1se)
betaHat.A<-c(betaHat.A[,1])

betaHat.B<-coef(lasso.B, s = lasso.B$lambda.1se)
betaHat.B<-c(betaHat.B[,1])

names(betaHat.A[abs(betaHat.A)>0])
names(betaHat.B[abs(betaHat.B)>0])

both.halves<-intersect(names(betaHat.A[abs(betaHat.A)>0]),
            names(betaHat.B[abs(betaHat.B)>0]))
# [2] "wknd_screen"                   
# [3] "bis_y_ss_bas_drive"            
# [4] "upps_y_ss_lack_of_perseverance"
# [5] "sst.v"                         
# [6] "nback.v.avg"                   
# [7] "pea_wiscv_tss"                 
# [8] "pea_ravlt_total"               
# [9] "nihtbx_list_agecorrected"      
# [10] "nihtbx_reading_agecorrected"   
# [11] "lmt_scr_perc_correct"          
# [12] "female_yes" 

one.half<-unique(c(names(betaHat.A[abs(betaHat.A)>0]),
           names(betaHat.B[abs(betaHat.B)>0])))
# [2] "sch_eng"                       
# [3] "wkdy_screen"                   
# [4] "wknd_screen"                   
# [5] "bis_y_ss_bas_drive"            
# [6] "upps_y_ss_lack_of_planning"    
# [7] "upps_y_ss_lack_of_perseverance"
# [8] "sst.v"                         
# [9] "nback.v.avg"                   
# [10] "pea_wiscv_tss"                 
# [11] "pea_ravlt_total"               
# [12] "nihtbx_list_agecorrected"      
# [13] "nihtbx_picture_agecorrected"   
# [14] "nihtbx_reading_agecorrected"   
# [15] "lmt_scr_perc_correct"          
# [16] "female_yes"                    
# [17] "sch_env"                       
# [18] "parent_monit"                  
# [19] "bis_y_ss_bas_fs"               
# [20] "upps_y_ss_positive_urgency"    
# [21] "nihtbx_cardsort_agecorrected"  
# [22] "nihtbx_pattern_agecorrected"   
# [23] "married_no"   

included.100<-names(betaHatLASSO.bin[betaHatLASSO.bin==1])
included.50<-names(betaHatLASSO.bin[betaHatLASSO.bin>.50])
included.0<-names(betaHatLASSO.bin[betaHatLASSO.bin>0])

sparse <- intersect(both.halves,
          included.100)
sparse <- sparse[sparse!="(Intercept)"]
sparse[sparse=="female_yes"]<-"female"
sparse[sparse=="married_no"]<-"married"
# all included in both halves except "bis_y_ss_bas_drive"

# try prediction with only these sparse variables

sparse.imputed<-impute.folds(y1.attn,
                            pred.cols = sparse,
                            sites = sites)


sparse.lasso<-list()

for (s in sites){
  tmp.dat<-setup.predicton(sparse.imputed,s,y1.attn)
  sparse.lasso[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  sparse.lasso[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  sparse.lasso[[s]]$train<-predict(sparse.lasso[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = sparse.lasso[[s]]$fits$lambda.1se)
  sparse.lasso[[s]]$test<-predict(sparse.lasso[[s]]$fits,data.matrix(tmp.dat$test), s = sparse.lasso[[s]]$fits$lambda.1se)
  sparse.lasso[[s]]$train.r<-cor(sparse.lasso[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  sparse.lasso[[s]]$test.r<-cor(sparse.lasso[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
}



# prediction just as good as full set
mean(unlist(lapply(sparse.lasso,FUN=function(x) x$train.r)))
mean(unlist(lapply(sparse.lasso,FUN=function(x) x$test.r)))


ci.sparse<-qt(0.975,df=17)*sd(unlist(lapply(sparse.lasso,FUN=function(x) x$test.r)))/sqrt(18)

mean(unlist(lapply(sparse.lasso,FUN=function(x) x$test.r)))+ci.sparse
mean(unlist(lapply(sparse.lasso,FUN=function(x) x$test.r)))-ci.sparse

betaHatLASSO.sparse<-as.data.frame(matrix(NA,nrow=length(sites),
                                         ncol = length(rownames(coef(sparse.lasso[["site11"]]$fits, 
                                                                     s = sparse.lasso[["site11"]]$fits$lambda.1se)))))
colnames(betaHatLASSO.sparse)<-rownames(coef(sparse.lasso[["site11"]]$fits, s = sparse.lasso[["site11"]]$fits$lambda.1se))
betaHatLASSO.sparse$sites<-sites

for (s in sites){
  coefs<-coef(sparse.lasso[[s]]$fits, s = sparse.lasso[[s]]$fits$lambda.1se)
  betaHatLASSO.sparse[betaHatLASSO.sparse$sites==s,
                       1:(length(betaHatLASSO.sparse)-1)]<-c(coefs[,1])
}

betaHatLASSO.sparse.avg<-colMeans(betaHatLASSO.sparse[,colnames(betaHatLASSO.sparse)!="sites"])
betaHatLASSO.sparse.avg<-betaHatLASSO.sparse.avg[names(betaHatLASSO.sparse.avg)!="(Intercept)"]

#######################################################################
#######################################################################
##### Sensitivity analysis: medication exclusions #####################
#######################################################################
#######################################################################

## Create fully imputed training and test data sets for each fold #####


nomeds.imputed<-impute.folds(y1.attn[y1.attn$CAS.ADHD==FALSE & 
                                       !is.na(y1.attn$CAS.ADHD),],
                            pred.cols = pred.cols,
                            sites = sites)

#save(nomeds.imputed,file="data_imputations_nomeds.RData")
load("data_imputations_nomeds.RData")

# predictive modeling

y1.attn.nomeds<-y1.attn[y1.attn$CAS.ADHD==FALSE & 
                          !is.na(y1.attn$CAS.ADHD),]

nomeds.lasso<-list()
nomeds.pcareg<-list()


for (s in sites){
  
  tmp.dat<-setup.predicton(nomeds.imputed,s,y1.attn.nomeds)
  nomeds.lasso[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  nomeds.lasso[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), 
                                    y1.attn.nomeds[y1.attn.nomeds$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  nomeds.lasso[[s]]$train<-predict(nomeds.lasso[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = nomeds.lasso[[s]]$fits$lambda.1se)
  nomeds.lasso[[s]]$test<-predict(nomeds.lasso[[s]]$fits,data.matrix(tmp.dat$test), s = nomeds.lasso[[s]]$fits$lambda.1se)
  nomeds.lasso[[s]]$train.r<-cor(nomeds.lasso[[s]]$train,y1.attn.nomeds[y1.attn.nomeds$site_id_l!=s,paste0("lo.",s,".attn")])
  nomeds.lasso[[s]]$test.r<-cor(nomeds.lasso[[s]]$test,y1.attn.nomeds[y1.attn.nomeds$site_id_l==s,paste0("lo.",s,".attn")])
  nomeds.pcareg[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  nomeds.pcareg[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  nomeds.pcareg[[s]]$train<-predict(nomeds.pcareg[[s]]$fits$winning.model)
  nomeds.pcareg[[s]]$test<-PCA.predict(nomeds.pcareg[[s]]$fits,tmp.dat$test)
  nomeds.pcareg[[s]]$train.r<-cor(nomeds.pcareg[[s]]$train,
                                  y1.attn.nomeds[y1.attn.nomeds$site_id_l!=s,paste0("lo.",s,".attn")])
  nomeds.pcareg[[s]]$test.r<-cor(nomeds.pcareg[[s]]$test,
                                 y1.attn.nomeds[y1.attn.nomeds$site_id_l==s,paste0("lo.",s,".attn")])
}


#save.image("medication_sensitivity_analyses.RData")


meds.r<-data.frame(nomeds.train.pca=unlist(lapply(nomeds.pcareg,FUN=function(x) x$train.r)),
                   nomeds.test.pca=unlist(lapply(nomeds.pcareg,FUN=function(x) x$test.r)),
                   nomeds.train.lasso=unlist(lapply(nomeds.lasso,FUN=function(x) x$train.r)),
                   nomeds.test.lasso=unlist(lapply(nomeds.lasso,FUN=function(x) x$test.r)) )
colMeans(meds.r)

ci.lasso.nomeds<-qt(0.975,df=17)*sd(meds.r$nomeds.test.lasso)/sqrt(18)
ci.pca.nomeds<-qt(0.975,df=17)*sd(meds.r$nomeds.test.pca)/sqrt(18)

colMeans(meds.r)["nomeds.test.lasso"]+ci.lasso.nomeds
colMeans(meds.r)["nomeds.test.lasso"]-ci.lasso.nomeds

colMeans(meds.r)["nomeds.test.pca"]+ci.pca.nomeds
colMeans(meds.r)["nomeds.test.pca"]-ci.pca.nomeds


#######################################################################
#######################################################################
##### Plot importance/weights of features #############################
#######################################################################
#######################################################################

load("pca_regression_models_6_20_22.RData")


fw.names<-c("Age","Area Deprivation Index","Area total crimes",
            "Area HS graduation rate","Area 3rd grade math proficiency",
            "Area 3rd grade reading proficiency", "Area school poverty",
            "Area high-quality preschools","Area estimated lead risk",
            "Waist circumference", "BMI", "Neighborhood crime (youth report)",
            "School environment"," School involvement", "School disengagement",
            "Family conflict", "Parental monitoring", "Weekday screen time",
            "Weekend screen time", "BIS summary", "BAS reward responsiveness",
            "BAS drive", "BAS fun seeking", "UPPS negative urgency",
            "UPPS lack of planning", "UPPS sensation seeking", 
            "UPPS positive urgency", "UPPS lack of perseverance", "SST drift rate",
            "SST boundary separation", "SST nondecision time",
            "N-back average drift rate", "N-back drift load effect",
            "N-back average boundary separation","N-back boundary load effect",
            "N-back average nondecision time", "N-back nondecision load effect",
            "WISC matrix reasoning", "RAVLT total correct recalls",
            "NIHTB picture vocabulary", "NIHTB flanker", "NIHTB list sorting",
            "NIHTB card sort", "NIHTB processing speed", 
            "NIHTB picture sequence memory","NIHTB oral reading",
            "Little man task accuracy")
names(fw.names)<-names(fw.means)
fw.cols=c(rep("cornflowerblue",11),
          rep("red2",17),rep("orange",19))
names(fw.cols)<-names(fw.means)

jpeg("con_feature_weights.jpg",width = 8,height = 10,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fw.means[names(sort(abs(fw.means),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fw.means),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fw.means),decreasing = FALSE))],
        border="white",xlim = c(-0.18,0.18),width=0.5)
abline(0,10e10)
dev.off()



cpw.names<-c("Female sex", "Asian race", "Black race", "Other/mixed race",
             "Hispanic ethnicity", "Parent education <HS", 
             "Parent undergraduate degree", "Parent post-graduate degree",
             "Parent some college education", "Parents not married",
             "Parental income <$50k", "Parental income >=$100k",
             "Cash choice 'yes'")
names(cpw.names)<-names(cpw.means)
cpw.cols=c(rep("cornflowerblue",12),
          rep("orange",1))
names(cpw.cols)<-names(cpw.means)

jpeg("cat_feature_weights.jpg",width = 8,height = 4,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(cpw.means[names(sort(abs(cpw.means),decreasing = FALSE))],horiz = TRUE,
        names.arg = cpw.names[names(sort(abs(cpw.means),decreasing = FALSE))],
        col=cpw.cols[names(sort(abs(cpw.means),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()

cpwl<-betaHatLASSO.sparse.avg[names(betaHatLASSO.sparse.avg)%in%names(cpw.means)]
fwl<-betaHatLASSO.sparse.avg[!names(betaHatLASSO.sparse.avg)%in%names(cpw.means)]


jpeg("con_lasso_feature_weights.jpg",width = 8,height = 3.5,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fwl[names(sort(abs(fwl),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fwl),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fwl),decreasing = FALSE))],
        border="white",xlim = c(-0.15,0.15),width=0.5)
abline(0,10e10)
dev.off()

jpeg("cat_lasso_feature_weights.jpg",width = 8,height = 2.1,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(c(0,cpwl[names(sort(abs(cpwl),decreasing = FALSE))]),horiz = TRUE,
        names.arg = c("",cpw.names[names(sort(abs(cpwl),decreasing = FALSE))]),
        col=cpw.cols[names(sort(abs(cpwl),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()

#######################################################################
#######################################################################
##### Look at prediction in lockbox sites #############################
#######################################################################
#######################################################################

# load in and bind with existing data

y1.validation<-read.csv("year1_validation_dat.csv",stringsAsFactors = TRUE)

y1.validation<-rbind(y1.validation,
                     read.csv("year1_full_dat.csv",stringsAsFactors = TRUE))


v.sites<-c("site10","site06","site20")

# re-level and name reference levels for each factor
y1.validation$female <- relevel(y1.validation$female, ref = 'no')
y1.validation$race.4level <- relevel(y1.validation$race.4level, ref = 'White')
y1.validation$hisp <- relevel(y1.validation$hisp, ref = 'no')
y1.validation$high.educ <- relevel(y1.validation$high.educ, ref = 'HS Diploma/GED')
y1.validation$married <- relevel(y1.validation$married, ref = 'yes')
y1.validation$household.income <- relevel(y1.validation$household.income, ref = '[>=50K & <100K]')
y1.validation$cash_choice_task <- relevel(y1.validation$cash_choice_task, ref = 'no')

# generate attention problems scores for the lockbox sites

training.models.v<-list()

t<-1

for (s in v.sites){
  
  t.mod = cfa(bf.final,y1.validation[y1.validation$site_id_l!=s,],  
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  training.models.v[[t]]<-t.mod 
  
  t<-(t+1)
  
}

names(training.models.v)<-paste0("lo.",v.sites)

for (s in v.sites){
  y1.validation[,paste0("lo.",s,".attn")]<-lavPredict(training.models.v[[paste0("lo.",s)]],                                              newdata = y1.validation)[,"Attn"]
}


# prediction

validation.imputed<-impute.folds(y1.validation,
                             pred.cols = pred.cols,
                             sites = v.sites)

#save(validation.imputed,file="data_imputations_validation.RData")
load("data_imputations_nomeds.RData")

validation.lasso<-list()
validation.pcareg<-list()

for (s in v.sites){
  
  tmp.dat<-setup.predicton(validation.imputed,s,y1.validation)
  validation.lasso[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  validation.lasso[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), 
                                        y1.validation[y1.validation$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  validation.lasso[[s]]$train<-predict(validation.lasso[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = validation.lasso[[s]]$fits$lambda.1se)
  validation.lasso[[s]]$test<-predict(validation.lasso[[s]]$fits,data.matrix(tmp.dat$test), s = validation.lasso[[s]]$fits$lambda.1se)
  validation.lasso[[s]]$train.r<-cor(validation.lasso[[s]]$train,y1.validation[y1.validation$site_id_l!=s,paste0("lo.",s,".attn")])
  validation.lasso[[s]]$test.r<-cor(validation.lasso[[s]]$test,y1.validation[y1.validation$site_id_l==s,paste0("lo.",s,".attn")])
  validation.pcareg[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  validation.pcareg[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  validation.pcareg[[s]]$train<-predict(validation.pcareg[[s]]$fits$winning.model)
  validation.pcareg[[s]]$test<-PCA.predict(validation.pcareg[[s]]$fits,tmp.dat$test)
  validation.pcareg[[s]]$train.r<-cor(validation.pcareg[[s]]$train,
                                      y1.validation[y1.validation$site_id_l!=s,paste0("lo.",s,".attn")])
  validation.pcareg[[s]]$test.r<-cor(validation.pcareg[[s]]$test,
                                     y1.validation[y1.validation$site_id_l==s,paste0("lo.",s,".attn")])
}



validation.r<-data.frame(validation.train.pca=unlist(lapply(validation.pcareg,FUN=function(x) x$train.r)),
                         validation.test.pca=unlist(lapply(validation.pcareg,FUN=function(x) x$test.r)),
                         validation.train.lasso=unlist(lapply(validation.lasso,FUN=function(x) x$train.r)),
                         validation.test.lasso=unlist(lapply(validation.lasso,FUN=function(x) x$test.r)) )
validation.r

save.image("validation_analyses.RData")


# also try sparse model

sparse.imputed.v<-impute.folds(y1.validation,
                             pred.cols = sparse,
                             sites = v.sites)


sparse.lasso.v<-list()

for (s in v.sites){
  tmp.dat<-setup.predicton(sparse.imputed.v,s,y1.validation)
  sparse.lasso.v[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  sparse.lasso.v[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.validation[y1.validation$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  sparse.lasso.v[[s]]$train<-predict(sparse.lasso.v[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = sparse.lasso.v[[s]]$fits$lambda.1se)
  sparse.lasso.v[[s]]$test<-predict(sparse.lasso.v[[s]]$fits,data.matrix(tmp.dat$test), s = sparse.lasso.v[[s]]$fits$lambda.1se)
  sparse.lasso.v[[s]]$train.r<-cor(sparse.lasso.v[[s]]$train,(y1.validation[y1.validation$site_id_l!=s,paste0("lo.",s,".attn")]))
  sparse.lasso.v[[s]]$test.r<-cor(sparse.lasso.v[[s]]$test,(y1.validation[y1.validation$site_id_l==s,paste0("lo.",s,".attn")]))
}

lasso.v.r<-data.frame(sparse.lasso.train.lasso=unlist(lapply(sparse.lasso.v,FUN=function(x) x$train.r)),
                      sparse.lasso.test.lasso=unlist(lapply(sparse.lasso.v,FUN=function(x) x$test.r)) )
lasso.v.r

#############################################################################