rm(list=ls())

# load all necessary packages
library(lavaan)
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
#### Read in and format data #####################################
##################################################################

# First, we load and reformat a merged data set with measures of 
# attention problems from the year 1 followup visit and a variety 
# of neurocognitive, psychosocial and demographic measures from the 
# earlier (baseline) visit, which we will use to predict prospective 
# attention problems.

# load data 

b.attn<-read.csv("base_full_dat.csv",stringsAsFactors = TRUE)
y1.attn<-read.csv("year1_full_dat.csv",stringsAsFactors = TRUE)
y2.attn<-read.csv("year2_full_dat.csv",stringsAsFactors = TRUE)

b.full<-read.csv("base_full_full_dat.csv",stringsAsFactors = TRUE)
y1.full<-read.csv("year1_full_full_dat.csv",stringsAsFactors = TRUE)
y2.full<-read.csv("year2_full_full_dat.csv",stringsAsFactors = TRUE)

# re-level and name reference levels for each factor

re_level<-function(dat){
  dat$female <- relevel(dat$female, ref = 'no')
  dat$race.4level <- relevel(dat$race.4level, ref = 'White')
  dat$hisp <- relevel(dat$hisp, ref = 'no')
  dat$high.educ <- relevel(dat$high.educ, ref = 'HS Diploma/GED')
  dat$married <- relevel(dat$married, ref = 'yes')
  dat$household.income <- relevel(dat$household.income, ref = '[>=50K & <100K]')
  dat$cash_choice_task <- relevel(dat$cash_choice_task, ref = 'no')
  dat
}

b.attn<-re_level(b.attn)
y1.attn<-re_level(y1.attn)
y2.attn<-re_level(y2.attn)
b.full<-re_level(b.full)
y1.full<-re_level(y1.full)
y2.full<-re_level(y2.full)

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

##################################################################
#### Measurement model for attention problems outcome measure ####
##################################################################

# To avoid overfitting in test data sets (each left out site), we need to 
# fit this outcome measurement model separately in each training fold and 
# then obtain factor scores for the outcome (the general attention factor) 
# for each independent test set using the original model estimated on the 
# training data.

# define model

cbcl.Attention.cols<-c("cbcl_q01_p","cbcl_q08_p","cbcl_q10_p",
                       "cbcl_q13_p","cbcl_q17_p","cbcl_q41_p",
                       "cbcl_q45_p","cbcl_q61_p","cbcl_q62_p",
                       "cbcl_q80_p")
cbcl.ADHD.cols<-c("cbcl_q04_p","cbcl_q08_p","cbcl_q10_p","cbcl_q41_p",
                  "cbcl_q78_p","cbcl_q93_p","cbcl_q104_p")
bpmt.Attention.cols<-c("bpmt_q1","bpmt_q3","bpmt_q4",
                       "bpmt_q5","bpmt_q9","bpmt_q13")

comp.cbcl.cols<-unique(c(cbcl.ADHD.cols,cbcl.Attention.cols))

cbcl.cols.final<-comp.cbcl.cols[!comp.cbcl.cols%in%c("cbcl_q01_p",
                                                     "cbcl_q61_p")]

Attn.cols.final<-c(bpmt.Attention.cols,comp.cbcl.cols)
Attn.cols.final<-Attn.cols.final[!Attn.cols.final%in%c("cbcl_q13_p",
                                                       "cbcl_q45_p",
                                                       "cbcl_q62_p",
                                                       "cbcl_q80_p")]

all.c<-unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))

bf.final <- paste(" Attn =~ ",paste(Attn.cols.final,collapse=" + "),
                  " Parent =~ ",paste(cbcl.cols.final,collapse=" + "),
                  " Teacher =~ ",paste(bpmt.Attention.cols,collapse=" + "),
                  " Attn ~~ 0*Parent 
                      Attn ~~ 0*Teacher
                      Parent ~~ 0*Teacher
                      ",sep="\n")


# loop through each training fold

sites<-unique(y1.attn$site_id_l)

training.models.b<-list()
training.models.y1<-list()
training.models.y2<-list()

t<-1
for (s in sites){
  
  t.mod = cfa(bf.final,b.attn[b.attn$site_id_l!=s,],  
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  training.models.b[[t]]<-t.mod 
  
  t.mod = cfa(bf.final,y1.attn[y1.attn$site_id_l!=s,],  
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  training.models.y1[[t]]<-t.mod 
  
  t.mod = cfa(bf.final,y2.attn[y2.attn$site_id_l!=s,],  
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  training.models.y2[[t]]<-t.mod 
  
  t<-(t+1)
  
}

names(training.models.b)<-paste0("lo.",sites)
names(training.models.y1)<-paste0("lo.",sites)
names(training.models.y2)<-paste0("lo.",sites)

# add predicted general factor scores from each training set to the data frame

for (s in sites){
  b.attn[,paste0("lo.",s,".attn")]<-lavPredict(training.models.b[[paste0("lo.",s)]],
                                                newdata = b.attn)[,"Attn"]
  y1.attn[,paste0("lo.",s,".attn")]<-lavPredict(training.models.y1[[paste0("lo.",s)]],
                                                newdata = y1.attn)[,"Attn"]
  y2.attn[,paste0("lo.",s,".attn")]<-lavPredict(training.models.y2[[paste0("lo.",s)]],
                                                newdata = y2.attn)[,"Attn"]
}

# general attention factors in each fold are close to identical  
fac.cors.b<-cor(b.attn[,paste0("lo.",sites,".attn")])
fac.cors.y1<-cor(y1.attn[,paste0("lo.",sites,".attn")])
fac.cors.y2<-cor(y2.attn[,paste0("lo.",sites,".attn")])

plot_ly(x=sites, y= sites, z = ~fac.cors.y1, type = "heatmap") %>%
  layout(title="Correlations between general attention factor scores 
         estimated in each training fold")

# for the data sets with missing outcome data, we must impute missing values at
# each fold. 

# we will use the mice() function for imputation
# also load add-on mice.reuse() function for imputing new observations 
# in the test data based on the training data imputation model
source("https://raw.githubusercontent.com/prockenschaub/Misc/master/R/mice.reuse/mice.reuse.R")

# generate imputed data sets and 

trainfull.imp.b<-list()
trainfull.imp.y1<-list()
trainfull.imp.y2<-list()

trainfull.models.y1<-list()
trainfull.models.y2<-list()

t<-1
for (s in sites){
  
  tmp.train.imp <- mice(data = b.full[b.full$site_id_l!=s,all.c],m = 1,print = FALSE)
  trainfull.imp.b[[t]] <- suppressWarnings(mice.reuse(tmp.train.imp, 
                                              b.full[,all.c], 
                                              maxit = 1,print = FALSE))$`1`
  
  tmp.train.imp <- mice(data = y1.full[y1.full$site_id_l!=s,all.c],m = 1,print = FALSE)
  trainfull.imp.y1[[t]] <- suppressWarnings(mice.reuse(tmp.train.imp, 
                                                      y1.full[,all.c], 
                                                      maxit = 1,print = FALSE))$`1`
  
  tmp.train.imp <- mice(data = y2.full[y2.full$site_id_l!=s,all.c],m = 1,print = FALSE)
  trainfull.imp.y2[[t]] <- suppressWarnings(mice.reuse(tmp.train.imp, 
                                                       y2.full[,all.c], 
                                                       maxit = 1,print = FALSE))$`1`
  
  t.mod = cfa(bf.final,trainfull.imp.y1[[t]][y1.full$site_id_l!=s,],
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  trainfull.models.y1[[t]]<-t.mod
  
  t.mod = cfa(bf.final,trainfull.imp.y2[[t]][y2.full$site_id_l!=s,],
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  trainfull.models.y2[[t]]<-t.mod
  

  t<-(t+1)
  
}

names(trainfull.imp.b)<-paste0("lo.",sites)
names(trainfull.imp.y1)<-paste0("lo.",sites)
names(trainfull.imp.y2)<-paste0("lo.",sites)
names(trainfull.models.y1)<-paste0("lo.",sites)
names(trainfull.models.y2)<-paste0("lo.",sites)

# add predicted general factor scores from each training set to the data frame
# include outcome model predictions about baseline site


for (s in sites){
  y1.full[,paste0("lo.",s,".attn")]<-lavPredict(trainfull.models.y1[[paste0("lo.",s)]],
                                                newdata = trainfull.imp.y1[[paste0("lo.",s)]])[,"Attn"]
  y2.full[,paste0("lo.",s,".attn")]<-lavPredict(trainfull.models.y2[[paste0("lo.",s)]],
                                                newdata = trainfull.imp.y2[[paste0("lo.",s)]])[,"Attn"]
  b.full[,paste0("lo.",s,".attn.b.y1m")]<-lavPredict(trainfull.models.y1[[paste0("lo.",s)]],
                                                newdata = trainfull.imp.b[[paste0("lo.",s)]])[,"Attn"]
  b.full[,paste0("lo.",s,".attn.b.y2m")]<-lavPredict(trainfull.models.y2[[paste0("lo.",s)]],
                                                  newdata = trainfull.imp.b[[paste0("lo.",s)]])[,"Attn"]
}

# generate data sets with shared participants across waves 

y1.diff<-y1.full[y1.full$subjectkey%in%b.full$subjectkey,]
y2.diff<-y2.full[y2.full$subjectkey%in%b.full$subjectkey,]

y1.diff<-merge(y1.diff,b.full[,c("subjectkey",colnames(b.full)[grep(".attn.b.y1m",colnames(b.full))])])
y2.diff<-merge(y2.diff,b.full[,c("subjectkey",colnames(b.full)[grep(".attn.b.y2m",colnames(b.full))])])

#generate residuals

for (s in sites){

  reg<-lm(paste0("lo.",s,".attn~lo.",s,".attn.b.y1m"),data=y1.diff[y1.diff$site_id_l!=s,])
  y1.diff[,paste0("lo.",s,".res")]<-(y1.diff[,paste0("lo.",s,".attn")]-predict(reg,newdata=y1.diff))
  
  reg<-lm(paste0("lo.",s,".attn~lo.",s,".attn.b.y2m"),data=y2.diff[y2.diff$site_id_l!=s,])
  y2.diff[,paste0("lo.",s,".res")]<-(y2.diff[,paste0("lo.",s,".attn")]-predict(reg,newdata=y2.diff))
   
}

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
                            out.data,diff=FALSE){
  
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
  
  # add attention or difference score
  if(diff==FALSE){
  tmp.dat$attn<-out.data[out.data$site_id_l!=site,paste0("lo.",site,".attn")]}
  else {
    tmp.dat$attn<-out.data[out.data$site_id_l!=site,paste0("lo.",site,".res")]
  }
  
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
y1.attn.fullpred<-y1.attn.fullpred[,c(pred.cols,"lo.site17.attn")]

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


#correlation matrix and heat map in imputed data set

y1.full.fullpred<-gen.ADI.abcd(dat=y1.full)
y1.full.fullpred<-y1.full.fullpred[,c(pred.cols,"lo.site17.attn")]

y1.full.fullpred <- mice(data = y1.full.fullpred,m = 1,print = FALSE)

y1.full.fullpred <- complete(y1.full.fullpred)

all.cors<-cor(model.matrix(~0+., data=y1.full.fullpred))

plot_ly(x=colnames(all.cors), y= colnames(all.cors), z = ~all.cors, type = "heatmap") %>%
  layout(title="Correlations between all features")

# check variance inflation factor

library(car)

full.model<-lm(lo.site17.attn~.,data = y1.full.fullpred)
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

y1.attn.imputed<-impute.folds(y1.attn,
                           pred.cols = pred.cols,
                           sites = sites)
y1.full.imputed<-impute.folds(y1.full,
                              pred.cols = pred.cols,
                              sites = sites)
y2.full.imputed<-impute.folds(y2.full,
                              pred.cols = pred.cols,
                              sites = sites)
y1.diff.imputed<-impute.folds(y1.diff,
                              pred.cols = pred.cols,
                              sites = sites)
y2.diff.imputed<-impute.folds(y2.diff,
                              pred.cols = pred.cols,
                              sites = sites)
save(y1.attn.imputed,y1.full.imputed,y2.full.imputed,
     y1.diff.imputed,y2.diff.imputed,
     file="data_imputations_all_preds.RData")


#######################################################################
### Sparse prediction models using individual variables ###############
#######################################################################

# loop to conduct internal cross-validation with ridge, lasso and 
# SVR (both linear and radial basis function kernels) and then 
# test model on left-out site

load("data_imputations_all_preds.RData")

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
  
  tmp.dat<-setup.predicton(y1.full.imputed,s,y1.full)
  
  ridge.fits[[s]]<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 0)
  ridge.preds[[s]]$train<-predict(ridge.fits[[s]],data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = ridge.fits[[s]]$lambda.1se)
  ridge.preds[[s]]$test<-predict(ridge.fits[[s]],data.matrix(tmp.dat$test), s = ridge.fits[[s]]$lambda.1se)
  ridge.train.r[[s]]<-cor(ridge.preds[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  ridge.test.r[[s]]<-cor(ridge.preds[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  
  lasso.fits[[s]]<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  lasso.preds[[s]]$train<-predict(lasso.fits[[s]],data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = lasso.fits[[s]]$lambda.1se)
  lasso.preds[[s]]$test<-predict(lasso.fits[[s]],data.matrix(tmp.dat$test), s = lasso.fits[[s]]$lambda.1se)
  lasso.train.r[[s]]<-cor(lasso.preds[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  lasso.test.r[[s]]<-cor(lasso.preds[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  
  svr.linear.fits[[s]]<-svm(attn~.,
                        data=tmp.dat$train,
                        kernel="linear")
  svr.linear.preds[[s]]$train <- predict(svr.linear.fits[[s]],tmp.dat$train)
  svr.linear.preds[[s]]$test <- predict(svr.linear.fits[[s]],tmp.dat$test)
  svr.linear.train.r[[s]]<-cor(svr.linear.preds[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  svr.linear.test.r[[s]]<-cor(svr.linear.preds[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])

  svr.rbf.fits[[s]]<-svm(attn~.,
                            data=tmp.dat$train,
                            kernel="radial")
  svr.rbf.preds[[s]]$train <- predict(svr.rbf.fits[[s]],tmp.dat$train)
  svr.rbf.preds[[s]]$test <- predict(svr.rbf.fits[[s]],tmp.dat$test)
  svr.rbf.train.r[[s]]<-cor(svr.rbf.preds[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  svr.rbf.test.r[[s]]<-cor(svr.rbf.preds[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
}


ridge.train.r<-as.numeric(ridge.train.r)
ridge.test.r<-as.numeric(ridge.test.r)
lasso.train.r<-as.numeric(lasso.train.r)
lasso.test.r<-as.numeric(lasso.test.r)
svr.linear.train.r<-as.numeric(svr.linear.train.r)
svr.linear.test.r<-as.numeric(svr.linear.test.r)
svr.rbf.train.r<-as.numeric(svr.rbf.train.r)
svr.rbf.test.r<-as.numeric(svr.rbf.test.r)

save.image(file="sparse_models.RData")
load("sparse_models.RData")

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

# CI for lasso

ci.lasso<-qt(0.975,df=17)*sd(lasso.test.r)/sqrt(18)
mean(lasso.test.r)+ci.lasso
mean(lasso.test.r)-ci.lasso


#plot out

plot.pred<-function(sites,pred){
  
  col=rgb(.75,0,0,.6)
  par(mfrow=c(2,3))
  
  for (s in sites){
    
    actual<-y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")]
    predicted<-pred[[s]]$test
    
    plot(actual,predicted,
         xlab="actual",ylab="predicted",
         main=s,
         col=col,cex.axis=1.25,cex.lab=1.5,cex.main=1.7,pch=16,cex=1.2)
    abline(lm(predicted~actual),lwd=2)
    
  }
  
}

site.ns<-c()
for (s in sites){site.ns[s]<-length(y1.full[y1.full$site_id_l==s,"subjectkey"])}

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

included<-names(betaHatLASSO.bin[betaHatLASSO.bin==1])

write.csv(data.frame(beta=betaHatLASSO.avg[included]),file="lasso_betas.csv")


# LASSO for the complete y1 sample, imputed y2 sample, and difference scores


lasso.fits.comp<-list()
lasso.preds.comp<-list()
lasso.train.r.comp<-c()
lasso.test.r.comp<-c()

lasso.fits.y2<-list()
lasso.preds.y2<-list()
lasso.train.r.y2<-c()
lasso.test.r.y2<-c()

lasso.fits.y1diff<-list()
lasso.preds.y1diff<-list()
lasso.train.r.y1diff<-c()
lasso.test.r.y1diff<-c()

lasso.fits.y2diff<-list()
lasso.preds.y2diff<-list()
lasso.train.r.y2diff<-c()
lasso.test.r.y2diff<-c()


for (s in sites){
  
  tmp.dat<-setup.predicton(y1.attn.imputed,s,y1.attn)

  lasso.fits.comp[[s]]<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  lasso.preds.comp[[s]]$train<-predict(lasso.fits.comp[[s]],data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = lasso.fits.comp[[s]]$lambda.1se)
  lasso.preds.comp[[s]]$test<-predict(lasso.fits.comp[[s]],data.matrix(tmp.dat$test), s = lasso.fits.comp[[s]]$lambda.1se)
  lasso.train.r.comp[[s]]<-cor(lasso.preds.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  lasso.test.r.comp[[s]]<-cor(lasso.preds.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(y2.full.imputed,s,y2.full)
  
  lasso.fits.y2[[s]]<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y2.full[y2.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  lasso.preds.y2[[s]]$train<-predict(lasso.fits.y2[[s]],data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = lasso.fits.y2[[s]]$lambda.1se)
  lasso.preds.y2[[s]]$test<-predict(lasso.fits.y2[[s]],data.matrix(tmp.dat$test), s = lasso.fits.y2[[s]]$lambda.1se)
  lasso.train.r.y2[[s]]<-cor(lasso.preds.y2[[s]]$train,y2.full[y2.full$site_id_l!=s,paste0("lo.",s,".attn")])
  lasso.test.r.y2[[s]]<-cor(lasso.preds.y2[[s]]$test,y2.full[y2.full$site_id_l==s,paste0("lo.",s,".attn")])

  tmp.dat<-setup.predicton(y1.diff.imputed,s,y1.diff,diff=TRUE)
  
  lasso.fits.y1diff[[s]]<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.diff[y1.diff$site_id_l!=s,paste0("lo.",s,".res")], alpha = 1)
  lasso.preds.y1diff[[s]]$train<-predict(lasso.fits.y1diff[[s]],data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = lasso.fits.y1diff[[s]]$lambda.min)
  lasso.preds.y1diff[[s]]$test<-predict(lasso.fits.y1diff[[s]],data.matrix(tmp.dat$test), s = lasso.fits.y1diff[[s]]$lambda.min)
  lasso.train.r.y1diff[[s]]<-cor(lasso.preds.y1diff[[s]]$train,y1.diff[y1.diff$site_id_l!=s,paste0("lo.",s,".res")])
  lasso.test.r.y1diff[[s]]<-cor(lasso.preds.y1diff[[s]]$test,y1.diff[y1.diff$site_id_l==s,paste0("lo.",s,".res")])
    
  tmp.dat<-setup.predicton(y2.diff.imputed,s,y2.diff,diff=TRUE)
  
  lasso.fits.y2diff[[s]]<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y2.diff[y2.diff$site_id_l!=s,paste0("lo.",s,".res")], alpha = 1)
  lasso.preds.y2diff[[s]]$train<-predict(lasso.fits.y2diff[[s]],data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = lasso.fits.y2diff[[s]]$lambda.min)
  lasso.preds.y2diff[[s]]$test<-predict(lasso.fits.y2diff[[s]],data.matrix(tmp.dat$test), s = lasso.fits.y2diff[[s]]$lambda.min)
  lasso.train.r.y2diff[[s]]<-cor(lasso.preds.y2diff[[s]]$train,y2.diff[y2.diff$site_id_l!=s,paste0("lo.",s,".res")])
  lasso.test.r.y2diff[[s]]<-cor(lasso.preds.y2diff[[s]]$test,y2.diff[y2.diff$site_id_l==s,paste0("lo.",s,".res")])
  
}


lasso.train.r.comp<-as.numeric(lasso.train.r.comp)
lasso.test.r.comp<-as.numeric(lasso.test.r.comp)
mean(lasso.train.r.comp)
mean(lasso.test.r.comp)
mean(lasso.test.r.comp)^2

ci.lasso.comp<-qt(0.975,df=17)*sd(lasso.test.r.comp)/sqrt(18)
mean(lasso.test.r.comp)+ci.lasso.comp
mean(lasso.test.r.comp)-ci.lasso.comp


lasso.train.r.y2<-as.numeric(lasso.train.r.y2)
lasso.test.r.y2<-as.numeric(lasso.test.r.y2)
mean(lasso.train.r.y2)
mean(lasso.test.r.y2)
mean(lasso.test.r.y2)^2

ci.lasso.y2<-qt(0.975,df=17)*sd(lasso.test.r.y2)/sqrt(18)
mean(lasso.test.r.y2)+ci.lasso.y2
mean(lasso.test.r.y2)-ci.lasso.y2


lasso.train.r.y1diff<-as.numeric(lasso.train.r.y1diff)
lasso.test.r.y1diff<-as.numeric(lasso.test.r.y1diff)
mean(lasso.train.r.y1diff)
mean(lasso.test.r.y1diff)
mean(lasso.test.r.y1diff)^2

ci.lasso.y1diff<-qt(0.975,df=17)*sd(lasso.test.r.y1diff)/sqrt(18)
mean(lasso.test.r.y1diff)+ci.lasso.y1diff
mean(lasso.test.r.y1diff)-ci.lasso.y1diff

lasso.train.r.y2diff<-as.numeric(lasso.train.r.y2diff)
lasso.test.r.y2diff<-as.numeric(lasso.test.r.y2diff)
mean(lasso.train.r.y2diff)
mean(lasso.test.r.y2diff)
mean(lasso.test.r.y2diff)^2

ci.lasso.y2diff<-qt(0.975,df=17)*sd(lasso.test.r.y2diff)/sqrt(18)
mean(lasso.test.r.y2diff)+ci.lasso.y2diff
mean(lasso.test.r.y2diff)-ci.lasso.y2diff

# save out for figures
save(lasso.test.r,lasso.test.r.full,
     lasso.test.r.y2,lasso.test.r.y1diff,lasso.test.r.y2diff,
     file="LASSO_test_rs.RData")

save.image(file="sparse_models.RData")
load("sparse_models.RData")

# look at importance of variables

betaHat.comp<-as.data.frame(matrix(NA,nrow=length(sites),ncol = length(rownames(coef(lasso.fits.comp[["site11"]], s = lasso.fits.full[["site11"]]$lambda.1se)))))
colnames(betaHat.comp)<-rownames(coef(lasso.fits.comp[["site11"]], s = lasso.fits.comp[["site11"]]$lambda.1se))
betaHat.comp$sites<-sites

betaHat.y2<-betaHat.comp
betaHat.y1diff<-betaHat.comp
betaHat.y2diff<-betaHat.comp


for (s in sites){
  coefs<-coef(lasso.fits.comp[[s]], s = lasso.fits.comp[[s]]$lambda.1se)
  betaHat.comp[betaHat.comp$sites==s,1:(length(betaHat.comp)-1)]<-c(coefs[,1])
  
  coefs<-coef(lasso.fits.y2[[s]], s = lasso.fits.y2[[s]]$lambda.1se)
  betaHat.y2[betaHat.y2$sites==s,1:(length(betaHat.y2)-1)]<-c(coefs[,1])
  
  coefs<-coef(lasso.fits.y1diff[[s]], s = lasso.fits.y1diff[[s]]$lambda.min)
  betaHat.y1diff[betaHat.y1diff$sites==s,1:(length(betaHat.y1diff)-1)]<-c(coefs[,1])
  
  coefs<-coef(lasso.fits.y2diff[[s]], s = lasso.fits.y2diff[[s]]$lambda.min)
  betaHat.y2diff[betaHat.y2diff$sites==s,1:(length(betaHat.y2diff)-1)]<-c(coefs[,1])
  
}


#which variables are included by LASSO in all folds?

betaHat.comp.bin<-colMeans(betaHat.comp[,colnames(betaHat.comp)!="sites"]!=0)
betaHat.y2.bin<-colMeans(betaHat.y2[,colnames(betaHat.y2)!="sites"]!=0)
betaHat.y1diff.bin<-colMeans(betaHat.y1diff[,colnames(betaHat.y1diff)!="sites"]!=0)
betaHat.y2diff.bin<-colMeans(betaHat.y2diff[,colnames(betaHat.y2diff)!="sites"]!=0)

included.comp<-names(betaHat.comp.bin[betaHat.comp.bin==1])
included.y2<-names(betaHat.y2.bin[betaHat.y2.bin==1])
included.y1diff<-names(betaHat.y1diff.bin[betaHat.y1diff.bin==1])
included.y2diff<-names(betaHat.y2diff.bin[betaHat.y2diff.bin==1])

write.csv(data.frame(beta=betaHat.comp[included.comp]),file="lasso_betas_comp.csv")
write.csv(data.frame(beta=betaHat.y2[included.y2]),file="lasso_betas_y2.csv")
write.csv(data.frame(beta=betaHat.y1diff[included.y1diff]),file="lasso_betas_y1diff.csv")
write.csv(data.frame(beta=betaHat.y2diff[included.y2diff]),file="lasso_betas_y2diff.csv")

save.image(file="sparse_models.RData")
load("sparse_models.RData")

#########################################################
##### LASSO predictor stability and efficacy ############
#########################################################

# split-half to identify whether predictors are stable

A.sites<-sample(sites,9,FALSE)
B.sites<-sites[!sites%in%A.sites]
#save(A.sites,B.sites,file="split_half_sites.RData")

# function to fit LASSO for each half
half.LASSO<-function(dat,sites,imp.out=TRUE,base=NA){
  
  lasso<-gen.ADI.abcd(dat[dat$site_id_l%in%sites,])
  lasso<-lasso[,pred.cols]
  lasso <- mice(data = lasso,m = 1,print = FALSE)
  lasso <- complete(lasso)
  
  s.means<-sapply(lasso[,sapply(lasso,is.numeric)],mean)
  s.sds<-sapply(lasso[,sapply(lasso,is.numeric)],sd)
  for (v in names(s.means)){
    lasso[,v]<-(lasso[,v]-s.means[v])/s.sds[v]}
  lasso<-dummy_cols(lasso,remove_first_dummy = TRUE,
                        remove_selected_columns = TRUE)
  if(is.na(base)){
    if(imp.out==TRUE){
      
        sem.dat <- complete(mice(data = dat[dat$site_id_l%in%sites,all.c],m = 1,print = FALSE))
        
        } else{ sem.dat<-dat[dat$site_id_l%in%sites,all.c] }
  
  t.mod = cfa(bf.final,sem.dat,  
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  dat<-dat[dat$site_id_l%in%sites,]
  dat$Attn<-lavPredict(t.mod,newdata = sem.dat)[,"Attn"]
  
  } else {
    sem.dat <- complete(mice(data = dat[dat$site_id_l%in%sites,all.c],m = 1,print = FALSE))
    
    t.mod = cfa(bf.final,sem.dat,  
                std.lv=TRUE,
                orthogonal=TRUE,
                estimator = "WLSMV")
    
    base.dat <- complete(mice(data = base[base$subjectkey%in%dat$subjectkey & 
                                            base$site_id_l%in%sites,all.c],m = 1,print = FALSE))
    
    t1<-lavPredict(t.mod,newdata = base.dat)[,"Attn"]
    t2<-lavPredict(t.mod,newdata = sem.dat)[,"Attn"]
    
    reg<-lm( t2 ~ t1)
    dat<-dat[dat$site_id_l%in%sites,]
    dat$Attn<-(t2-predict(reg))
    
  }
  
  lasso<-cv.glmnet(data.matrix(lasso), 
                     dat$Attn,
                     alpha = 1)
  lasso
  
}

y1.full.lasso.A<-half.LASSO(dat=y1.full,sites=A.sites)
betaHat.y1.A<-c(coef(y1.full.lasso.A, s = y1.full.lasso.A$lambda.1se)[,1])
y1.full.lasso.B<-half.LASSO(dat=y1.full,sites=B.sites)
betaHat.y1.B<-c(coef(y1.full.lasso.B, s = y1.full.lasso.B$lambda.1se)[,1])

both.halves.y1<-intersect(names(betaHat.y1.A[abs(betaHat.y1.A)>0]),
                          names(betaHat.y1.B[abs(betaHat.y1.B)>0]))

sparsevars.y1<-intersect(both.halves.y1,included)

comp.lasso.A<-half.LASSO(dat=y1.attn,sites=A.sites,imp.out=FALSE)
betaHat.comp.A<-c(coef(comp.lasso.A, s = y1.comp.lasso.A$lambda.1se)[,1])
comp.lasso.B<-half.LASSO(dat=y1.attn,sites=B.sites,imp.out=FALSE)
betaHat.comp.B<-c(coef(comp.lasso.B, s = y1.comp.lasso.B$lambda.1se)[,1])

both.halves.comp<-intersect(names(betaHat.comp.A[abs(betaHat.comp.A)>0]),
                          names(betaHat.comp.B[abs(betaHat.comp.B)>0]))

sparsevars.comp<-intersect(both.halves.comp,included.comp)

y2.full.lasso.A<-half.LASSO(dat=y2.full,sites=A.sites)
betaHat.y2.A<-c(coef(y2.full.lasso.A, s = y2.full.lasso.A$lambda.1se)[,1])
y2.full.lasso.B<-half.LASSO(dat=y2.full,sites=B.sites)
betaHat.y2.B<-c(coef(y2.full.lasso.B, s = y2.full.lasso.B$lambda.1se)[,1])

both.halves.y2<-intersect(names(betaHat.y2.A[abs(betaHat.y2.A)>0]),
                          names(betaHat.y2.B[abs(betaHat.y2.B)>0]))

sparsevars.y2<-intersect(both.halves.y2,included.y2)

y1.diff.lasso.A<-half.LASSO(dat=y1.diff,sites=A.sites,base = b.full)
betaHat.y1diff.A<-c(coef(y1.diff.lasso.A, s = y1.diff.lasso.A$lambda.min)[,1])
y1.diff.lasso.B<-half.LASSO(dat=y1.diff,sites=B.sites,base = b.full)
betaHat.y1diff.B<-c(coef(y1.diff.lasso.B, s = y1.diff.lasso.B$lambda.min)[,1])

both.halves.y1diff<-intersect(names(betaHat.y1diff.A[abs(betaHat.y1diff.A)>0]),
                          names(betaHat.y1diff.B[abs(betaHat.y1diff.B)>0]))

sparsevars.y1diff<-intersect(both.halves.y1diff,included.y1diff)

y2.diff.lasso.A<-half.LASSO(dat=y2.diff,sites=A.sites,base = b.full)
betaHat.y2diff.A<-c(coef(y2.diff.lasso.A, s = y2.diff.lasso.A$lambda.min)[,1])
y2.diff.lasso.B<-half.LASSO(dat=y2.diff,sites=B.sites,base = b.full)
betaHat.y2diff.B<-c(coef(y2.diff.lasso.B, s = y2.diff.lasso.B$lambda.min)[,1])

both.halves.y2diff<-intersect(names(betaHat.y2diff.A[abs(betaHat.y2diff.A)>0]),
                              names(betaHat.y2diff.B[abs(betaHat.y2diff.B)>0]))

sparsevars.y2diff<-intersect(both.halves.y2diff,included.y2diff)


save.image(file="sparse_model_selection.RData")
load("sparse_model_selection.RData")


sparse.y1 <- sparsevars.y1[sparsevars.y1!="(Intercept)"]
sparse.y1[sparse.y1=="female_yes"]<-"female"
sparse.y1[sparse.y1=="married_no"]<-"married"

sparse.comp <- sparsevars.comp[sparsevars.comp!="(Intercept)"]
sparse.comp[sparse.comp=="female_yes"]<-"female"
sparse.comp[sparse.comp=="married_no"]<-"married"

sparse.y2 <- sparsevars.y2[sparsevars.y2!="(Intercept)"]
sparse.y2[sparse.y2=="female_yes"]<-"female"
sparse.y2[sparse.y2=="married_no"]<-"married"

sparse.y1diff <- sparsevars.y1diff [sparsevars.y1diff !="(Intercept)"]
sparse.y1diff [sparse.y1diff =="female_yes"]<-"female"
sparse.y1diff [sparse.y1diff =="married_no"]<-"married"

sparse.y2diff <- sparsevars.y2diff [sparsevars.y2diff !="(Intercept)"]
sparse.y2diff [sparse.y2diff =="female_yes"]<-"female"
sparse.y2diff [sparse.y2diff =="hisp_yes"]<-"hisp"

# try prediction with only these sparse variables

sparse.imp.y1<-impute.folds(y1.full,
                             pred.cols = sparse.y1,
                             sites = sites)

sparse.imp.comp<-impute.folds(y1.attn,
                            pred.cols = sparse.comp,
                            sites = sites)

sparse.imp.y2<-impute.folds(y2.full,
                            pred.cols = sparse.y2,
                            sites = sites)

sparse.imp.y1diff<-impute.folds(y1.diff,
                            pred.cols = sparse.y1diff,
                            sites = sites)

sparse.imp.y2diff<-impute.folds(y2.diff,
                                pred.cols = sparse.y2diff,
                                sites = sites)


sparse.lasso.y1<-list()
sparse.lasso.comp<-list()
sparse.lasso.y2<-list()
sparse.lasso.y1diff<-list()
sparse.lasso.y2diff<-list()

for (s in sites){
  
  tmp.dat<-setup.predicton(sparse.imp.y1,s,y1.full)
  sparse.lasso.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  sparse.lasso.y1[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  sparse.lasso.y1[[s]]$train<-predict(sparse.lasso.y1[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = sparse.lasso.y1[[s]]$fits$lambda.1se)
  sparse.lasso.y1[[s]]$test<-predict(sparse.lasso.y1[[s]]$fits,data.matrix(tmp.dat$test), s = sparse.lasso.y1[[s]]$fits$lambda.1se)
  sparse.lasso.y1[[s]]$train.r<-cor(sparse.lasso.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  sparse.lasso.y1[[s]]$test.r<-cor(sparse.lasso.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(sparse.imp.comp,s,y1.attn)
  sparse.lasso.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  sparse.lasso.comp[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  sparse.lasso.comp[[s]]$train<-predict(sparse.lasso.comp[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = sparse.lasso.comp[[s]]$fits$lambda.1se)
  sparse.lasso.comp[[s]]$test<-predict(sparse.lasso.comp[[s]]$fits,data.matrix(tmp.dat$test), s = sparse.lasso.comp[[s]]$fits$lambda.1se)
  sparse.lasso.compl[[s]]$train.r<-cor(sparse.lasso.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  sparse.lasso.comp[[s]]$test.r<-cor(sparse.lasso.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(sparse.imp.y2,s,y2.full)
  sparse.lasso.y2[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  sparse.lasso.y2[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y2.full[y2.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  sparse.lasso.y2[[s]]$train<-predict(sparse.lasso.y2[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = sparse.lasso.y2[[s]]$fits$lambda.1se)
  sparse.lasso.y2[[s]]$test<-predict(sparse.lasso.y2[[s]]$fits,data.matrix(tmp.dat$test), s = sparse.lasso.y2[[s]]$fits$lambda.1se)
  sparse.lasso.y2[[s]]$train.r<-cor(sparse.lasso.y2[[s]]$train,y2.full[y2.full$site_id_l!=s,paste0("lo.",s,".attn")])
  sparse.lasso.y2[[s]]$test.r<-cor(sparse.lasso.y2[[s]]$test,y2.full[y2.full$site_id_l==s,paste0("lo.",s,".attn")])
 
  tmp.dat<-setup.predicton(sparse.imp.y1diff,s,y1.diff,diff=TRUE)
  sparse.lasso.y1diff[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  sparse.lasso.y1diff[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.diff[y1.diff$site_id_l!=s,paste0("lo.",s,".res")], alpha = 1)
  sparse.lasso.y1diff[[s]]$train<-predict(sparse.lasso.y1diff[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = sparse.lasso.y1diff[[s]]$fits$lambda.min)
  sparse.lasso.y1diff[[s]]$test<-predict(sparse.lasso.y1diff[[s]]$fits,data.matrix(tmp.dat$test), s = sparse.lasso.y1diff[[s]]$fits$lambda.min)
  sparse.lasso.y1diff[[s]]$train.r<-cor(sparse.lasso.y1diff[[s]]$train,y1.diff[y1.diff$site_id_l!=s,paste0("lo.",s,".res")])
  sparse.lasso.y1diff[[s]]$test.r<-cor(sparse.lasso.y1diff[[s]]$test,y1.diff[y1.diff$site_id_l==s,paste0("lo.",s,".res")])
  
  tmp.dat<-setup.predicton(sparse.imp.y2diff,s,y2.diff,diff=TRUE)
  sparse.lasso.y2diff[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  sparse.lasso.y2diff[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y2.diff[y2.diff$site_id_l!=s,paste0("lo.",s,".res")], alpha = 1)
  sparse.lasso.y2diff[[s]]$train<-predict(sparse.lasso.y2diff[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = sparse.lasso.y2diff[[s]]$fits$lambda.min)
  sparse.lasso.y2diff[[s]]$test<-predict(sparse.lasso.y2diff[[s]]$fits,data.matrix(tmp.dat$test), s = sparse.lasso.y2diff[[s]]$fits$lambda.min)
  sparse.lasso.y2diff[[s]]$train.r<-cor(sparse.lasso.y2diff[[s]]$train,y2.diff[y2.diff$site_id_l!=s,paste0("lo.",s,".res")])
  sparse.lasso.y2diff[[s]]$test.r<-cor(sparse.lasso.y2diff[[s]]$test,y2.diff[y2.diff$site_id_l==s,paste0("lo.",s,".res")])
   
}

unlist(lapply(sparse.lasso.y1diff,FUN=function(x) x$test.r))
unlist(lapply(sparse.lasso.y2diff,FUN=function(x) x$test.r))

unlist(lapply(sparse.lasso.y2,FUN=function(x) x$test.r))
unlist(lapply(sparse.lasso.y1,FUN=function(x) x$test.r))

# prediction just as good as full set of predictors

mean(unlist(lapply(sparse.lasso.y1,FUN=function(x) x$train.r)))
mean(unlist(lapply(sparse.lasso.y1,FUN=function(x) x$test.r)))
ci.sparse.lasso.y1<-qt(0.975,df=17)*sd(unlist(lapply(sparse.lasso.y1,FUN=function(x) x$test.r)))/sqrt(18)
mean(unlist(lapply(sparse.lasso.y1,FUN=function(x) x$test.r)))+ci.sparse.lasso.y1
mean(unlist(lapply(sparse.lasso.y1,FUN=function(x) x$test.r)))-ci.sparse.lasso.y1



mean(unlist(lapply(sparse.lasso.comp,FUN=function(x) x$train.r)))
mean(unlist(lapply(sparse.lasso.comp,FUN=function(x) x$test.r)))
ci.sparse.lasso.comp<-qt(0.975,df=17)*sd(unlist(lapply(sparse.lasso.comp,FUN=function(x) x$test.r)))/sqrt(18)
mean(unlist(lapply(sparse.lasso.comp,FUN=function(x) x$test.r)))+ci.sparse.lasso.comp
mean(unlist(lapply(sparse.lasso.comp,FUN=function(x) x$test.r)))-ci.sparse.lasso.comp



mean(unlist(lapply(sparse.lasso.y2,FUN=function(x) x$train.r)))
mean(unlist(lapply(sparse.lasso.y2,FUN=function(x) x$test.r)))

mean(unlist(lapply(sparse.lasso.y1diff,FUN=function(x) x$train.r)))
mean(unlist(lapply(sparse.lasso.y1diff,FUN=function(x) x$test.r)))

mean(unlist(lapply(sparse.lasso.y2diff,FUN=function(x) x$train.r)))
mean(unlist(lapply(sparse.lasso.y2diff,FUN=function(x) x$test.r)))



betaSpar.comp<-as.data.frame(matrix(NA,nrow=length(sites),ncol = length(rownames(coef(sparse.lasso.comp[["site11"]]$fits, s = sparse.lasso.full[["site11"]]$fits$lambda.min)))))
colnames(betaSpar.comp)<-rownames(coef(sparse.lasso.comp[["site11"]]$fits, s =sparse.lasso.comp[["site11"]]$fits$lambda.1se))
betaSpar.comp$sites<-sites

betaSpar.y1<-as.data.frame(matrix(NA,nrow=length(sites),ncol = length(rownames(coef(sparse.lasso.y1[["site11"]]$fits, s = sparse.lasso.y1[["site11"]]$fits$lambda.min)))))
colnames(betaSpar.y1)<-rownames(coef(sparse.lasso.y1[["site11"]]$fits, s =sparse.lasso.y1[["site11"]]$fits$lambda.1se))
betaSpar.y1$sites<-sites

betaSpar.y2<-as.data.frame(matrix(NA,nrow=length(sites),ncol = length(rownames(coef(sparse.lasso.y2[["site11"]]$fits, s = sparse.lasso.y2[["site11"]]$fits$lambda.min)))))
colnames(betaSpar.y2)<-rownames(coef(sparse.lasso.y2[["site11"]]$fits, s =sparse.lasso.y2[["site11"]]$fits$lambda.1se))
betaSpar.y2$sites<-sites

betaSpar.y1diff<-as.data.frame(matrix(NA,nrow=length(sites),ncol = length(rownames(coef(sparse.lasso.y1diff[["site11"]]$fits, s = sparse.lasso.y1diff[["site11"]]$fits$lambda.min)))))
colnames(betaSpar.y1diff)<-rownames(coef(sparse.lasso.y1diff[["site11"]]$fits, s =sparse.lasso.y1diff[["site11"]]$fits$lambda.1se))
betaSpar.y1diff$sites<-sites

betaSpar.y2diff<-as.data.frame(matrix(NA,nrow=length(sites),ncol = length(rownames(coef(sparse.lasso.y2diff[["site11"]]$fits, s = sparse.lasso.y2diff[["site11"]]$fits$lambda.min)))))
colnames(betaSpar.y2diff)<-rownames(coef(sparse.lasso.y2diff[["site11"]]$fits, s =sparse.lasso.y2diff[["site11"]]$fits$lambda.1se))
betaSpar.y2diff$sites<-sites


for (s in sites){

  coefs<-coef(sparse.lasso.y1[[s]]$fits, s = sparse.lasso.y1[[s]]$fits$lambda.1se)
  betaSpar.y1[betaSpar.y1$sites==s,1:(length(betaSpar.y1)-1)]<-c(coefs[,1])
    
  coefs<-coef(sparse.lasso.comp[[s]]$fits, s = sparse.lasso.comp[[s]]$fits$lambda.1se)
  betaSpar.comp[betaSpar.comp$sites==s,1:(length(betaSpar.comp)-1)]<-c(coefs[,1])
  
  coefs<-coef(sparse.lasso.y2[[s]]$fits, s = sparse.lasso.y2[[s]]$fits$lambda.1se)
  betaSpar.y2[betaSpar.y2$sites==s,1:(length(betaSpar.y2)-1)]<-c(coefs[,1])
  
  coefs<-coef(sparse.lasso.y1diff[[s]]$fits, s = sparse.lasso.y1diff[[s]]$fits$lambda.min)
  betaSpar.y1diff[betaSpar.y1diff$sites==s,1:(length(betaSpar.y1diff)-1)]<-c(coefs[,1])
  
  coefs<-coef(sparse.lasso.y2diff[[s]]$fits, s = sparse.lasso.y2diff[[s]]$fits$lambda.min)
  betaSpar.y2diff[betaSpar.y2diff$sites==s,1:(length(betaSpar.y2diff)-1)]<-c(coefs[,1])
  
}



save.image(file="sparse_model_selection.RData")
load("sparse_model_selection.RData")

#######################################################################
## Principal components analysis and regression #######################
#######################################################################

# conduct and interpret PCA in full sample

# start by imputing in complete data sample (minus outcome variable) 

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

pcareg.fits.y1<-list()
pcareg.preds.y1<-list()
pcareg.train.r.y1<-c()
pcareg.test.r.y1<-c()

pcareg.fits.comp<-list()
pcareg.preds.comp<-list()
pcareg.train.r.comp<-c()
pcareg.test.r.comp<-c()

pcareg.fits.y2<-list()
pcareg.preds.y2<-list()
pcareg.train.r.y2<-c()
pcareg.test.r.y2<-c()

pcareg.fits.y1diff<-list()
pcareg.preds.y1diff<-list()
pcareg.train.r.y1diff<-c()
pcareg.test.r.y1diff<-c()

pcareg.fits.y2diff<-list()
pcareg.preds.y2diff<-list()
pcareg.train.r.y2diff<-c()
pcareg.test.r.y2diff<-c()

for (s in sites){

  tmp.dat<-setup.predicton(y1.full.imputed,s,y1.full)
  pcareg.fits.y1[[s]]<-PCA.reg(data=tmp.dat$train, outcome = "attn", 
                            PCA.cols = tmp.dat$con.vars)
  pcareg.preds.y1[[s]]$train<-predict(pcareg.fits.y1[[s]]$winning.model)
  pcareg.preds.y1[[s]]$test<-PCA.predict(pcareg.fits.y1[[s]],tmp.dat$test)
  pcareg.train.r.y1[[s]]<-cor(pcareg.preds.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  pcareg.test.r.y1[[s]]<-cor(pcareg.preds.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(y1.attn.imputed,s,y1.attn)
  pcareg.fits.comp[[s]]<-PCA.reg(data=tmp.dat$train, outcome = "attn", 
                               PCA.cols = tmp.dat$con.vars)
  pcareg.preds.comp[[s]]$train<-predict(pcareg.fits.comp[[s]]$winning.model)
  pcareg.preds.comp[[s]]$test<-PCA.predict(pcareg.fits.comp[[s]],tmp.dat$test)
  pcareg.train.r.comp[[s]]<-cor(pcareg.preds.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  pcareg.test.r.comp[[s]]<-cor(pcareg.preds.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(y2.full.imputed,s,y2.full)
  pcareg.fits.y2[[s]]<-PCA.reg(data=tmp.dat$train, outcome = "attn", 
                               PCA.cols = tmp.dat$con.vars)
  pcareg.preds.y2[[s]]$train<-predict(pcareg.fits.y2[[s]]$winning.model)
  pcareg.preds.y2[[s]]$test<-PCA.predict(pcareg.fits.y2[[s]],tmp.dat$test)
  pcareg.train.r.y2[[s]]<-cor(pcareg.preds.y2[[s]]$train,y2.full[y2.full$site_id_l!=s,paste0("lo.",s,".attn")])
  pcareg.test.r.y2[[s]]<-cor(pcareg.preds.y2[[s]]$test,y2.full[y2.full$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(y1.diff.imputed,s,y1.diff,diff=TRUE)
  pcareg.fits.y1diff[[s]]<-PCA.reg(data=tmp.dat$train, outcome = "attn", 
                               PCA.cols = tmp.dat$con.vars)
  pcareg.preds.y1diff[[s]]$train<-predict(pcareg.fits.y1diff[[s]]$winning.model)
  pcareg.preds.y1diff[[s]]$test<-PCA.predict(pcareg.fits.y1diff[[s]],tmp.dat$test)
  pcareg.train.r.y1diff[[s]]<-cor(pcareg.preds.y1diff[[s]]$train,y1.diff[y1.diff$site_id_l!=s,paste0("lo.",s,".res")])
  pcareg.test.r.y1diff[[s]]<-cor(pcareg.preds.y1diff[[s]]$test,y1.diff[y1.diff$site_id_l==s,paste0("lo.",s,".res")])
  
  tmp.dat<-setup.predicton(y2.diff.imputed,s,y2.diff,diff=TRUE)
  pcareg.fits.y2diff[[s]]<-PCA.reg(data=tmp.dat$train, outcome = "attn", 
                                   PCA.cols = tmp.dat$con.vars)
  pcareg.preds.y2diff[[s]]$train<-predict(pcareg.fits.y2diff[[s]]$winning.model)
  pcareg.preds.y2diff[[s]]$test<-PCA.predict(pcareg.fits.y2diff[[s]],tmp.dat$test)
  pcareg.train.r.y2diff[[s]]<-cor(pcareg.preds.y2diff[[s]]$train,y2.diff[y2.diff$site_id_l!=s,paste0("lo.",s,".res")])
  pcareg.test.r.y2diff[[s]]<-cor(pcareg.preds.y2diff[[s]]$test,y2.diff[y2.diff$site_id_l==s,paste0("lo.",s,".res")])
  
}

pca.train.r.y1<-as.numeric(pcareg.train.r.y1)
pca.test.r.y1<-as.numeric(pcareg.test.r.y1)
pca.train.r.comp<-as.numeric(pcareg.train.r.comp)
pca.test.r.comp<-as.numeric(pcareg.test.r.comp)
pca.train.r.y2<-as.numeric(pcareg.train.r.y2)
pca.test.r.y2<-as.numeric(pcareg.test.r.y2)
pca.train.r.y1diff<-as.numeric(pcareg.train.r.y1diff)
pca.test.r.y1diff<-as.numeric(pcareg.test.r.y1diff)
pca.train.r.y2diff<-as.numeric(pcareg.train.r.y2diff)
pca.test.r.y2diff<-as.numeric(pcareg.test.r.y2diff)

mean(pca.train.r.y1); mean(pca.test.r.y1); mean(pca.test.r.y1)^2
ci.pca.y1<-qt(0.975,df=17)*sd(pca.test.r.y1)/sqrt(18)
mean(pca.test.r.y1)+ci.pca.y1
mean(pca.test.r.y1)-ci.pca.y1

mean(pca.train.r.comp); mean(pca.test.r.comp); mean(pca.test.r.comp)^2
ci.pca.comp<-qt(0.975,df=17)*sd(pca.test.r.comp)/sqrt(18)
mean(pca.test.r.comp)+ci.pca.comp
mean(pca.test.r.comp)-ci.pca.comp

mean(pca.train.r.y2); mean(pca.test.r.y2); mean(pca.test.r.y2)^2
ci.pca.y2<-qt(0.975,df=17)*sd(pca.test.r.y2)/sqrt(18)
mean(pca.test.r.y2)+ci.pca.y2
mean(pca.test.r.y2)-ci.pca.y2


mean(pca.train.r.y1diff); mean(pca.test.r.y1diff); mean(pca.test.r.y1diff)^2
ci.pca.y1diff<-qt(0.975,df=17)*sd(pca.test.r.y1diff)/sqrt(18)
mean(pca.test.r.y1diff)+ci.pca.y1diff
mean(pca.test.r.y1diff)-ci.pca.y1diff


mean(pca.train.r.y2diff); mean(pca.test.r.y2diff); mean(pca.test.r.y2diff)^2
ci.pca.y2diff<-qt(0.975,df=17)*sd(pca.test.r.y2diff)/sqrt(18)
mean(pca.test.r.y2diff)+ci.pca.y2diff
mean(pca.test.r.y2diff)-ci.pca.y2diff

# save out for figures
save(pca.test.r.y1,pca.test.r.comp,
     pca.test.r.y2,pca.test.r.y1diff,pca.test.r.y2diff,
     file="PCA_test_rs.RData")


# feature weights

fw.pca.y1<-data.frame(matrix(NA,length(sites),length(pcareg.fits.y1[[s]]$fw)))
colnames(fw.pca.y1)<-rownames(pcareg.fits.y1[[s]]$fw)
rownames(fw.pca.y1)<-sites
cpw.pca.y1<-data.frame(matrix(NA,length(sites),13))
colnames(cpw.pca.y1)<-names(pcareg.fits.y1[[s]]$winning.model$coefficients[2:14])
rownames(cpw.pca.y1)<-sites

for (s in sites){
  fw.pca.y1[s,]<-pcareg.fits.y1[[s]]$fw
  cpw.pca.y1[s,]<- pcareg.fits.y1[[s]]$winning.model$coefficients[2:14]
}

fw.means.y1<-colMeans(fw.pca.y1)
cpw.means.y1<-colMeans(cpw.pca.y1)

fw.pca.comp<-data.frame(matrix(NA,length(sites),length(pcareg.fits.comp[[s]]$fw)))
colnames(fw.pca.comp)<-rownames(pcareg.fits.comp[[s]]$fw)
rownames(fw.pca.comp)<-sites
cpw.pca.comp<-data.frame(matrix(NA,length(sites),13))
colnames(cpw.pca.comp)<-names(pcareg.fits.comp[[s]]$winning.model$coefficients[2:14])
rownames(cpw.pca.comp)<-sites

for (s in sites){
  fw.pca.comp[s,]<-pcareg.fits.comp[[s]]$fw
  cpw.pca.comp[s,]<- pcareg.fits.comp[[s]]$winning.model$coefficients[2:14]
}

fw.means.comp<-colMeans(fw.pca.comp)
cpw.means.comp<-colMeans(cpw.pca.comp)



fw.pca.y2<-data.frame(matrix(NA,length(sites),length(pcareg.fits.y2[[s]]$fw)))
colnames(fw.pca.y2)<-rownames(pcareg.fits.y2[[s]]$fw)
rownames(fw.pca.y2)<-sites
cpw.pca.y2<-data.frame(matrix(NA,length(sites),13))
colnames(cpw.pca.y2)<-names(pcareg.fits.y2[[s]]$winning.model$coefficients[2:14])
rownames(cpw.pca.y2)<-sites

for (s in sites){
  fw.pca.y2[s,]<-pcareg.fits.y2[[s]]$fw
  cpw.pca.y2[s,]<- pcareg.fits.y2[[s]]$winning.model$coefficients[2:14]
}

fw.means.y2<-colMeans(fw.pca.y2)
cpw.means.y2<-colMeans(cpw.pca.y2)


fw.pca.y1diff<-data.frame(matrix(NA,length(sites),length(pcareg.fits.y1diff[[s]]$fw)))
colnames(fw.pca.y1diff)<-rownames(pcareg.fits.y1diff[[s]]$fw)
rownames(fw.pca.y1diff)<-sites
cpw.pca.y1diff<-data.frame(matrix(NA,length(sites),13))
colnames(cpw.pca.y1diff)<-names(pcareg.fits.y1diff[[s]]$winning.model$coefficients[2:14])
rownames(cpw.pca.y1diff)<-sites

for (s in sites){
  fw.pca.y1diff[s,]<-pcareg.fits.y1diff[[s]]$fw
  cpw.pca.y1diff[s,]<- pcareg.fits.y1diff[[s]]$winning.model$coefficients[2:14]
}

fw.means.y1diff<-colMeans(fw.pca.y1diff)
cpw.means.y1diff<-colMeans(cpw.pca.y1diff)


fw.pca.y2diff<-data.frame(matrix(NA,length(sites),length(pcareg.fits.y2diff[[s]]$fw)))
colnames(fw.pca.y2diff)<-rownames(pcareg.fits.y2diff[[s]]$fw)
rownames(fw.pca.y2diff)<-sites
cpw.pca.y2diff<-data.frame(matrix(NA,length(sites),13))
colnames(cpw.pca.y2diff)<-names(pcareg.fits.y2diff[[s]]$winning.model$coefficients[2:14])
rownames(cpw.pca.y2diff)<-sites

for (s in sites){
  fw.pca.y2diff[s,]<-pcareg.fits.y2diff[[s]]$fw
  cpw.pca.y2diff[s,]<- pcareg.fits.y2diff[[s]]$winning.model$coefficients[2:14]
}

fw.means.y2diff<-colMeans(fw.pca.y2diff)
cpw.means.y2diff<-colMeans(cpw.pca.y2diff)



save.image(file="pca_regression_models.RData")


save(cpw.means.comp,cpw.means.y1,
     cpw.means.y2,cpw.means.y1diff,cpw.means.y2diff,
     fw.means.comp,fw.means.y1,
     fw.means.y2,fw.means.y1diff,fw.means.y2diff,
     file="PCA_weights.RData")


########################################################
########################################################
############## summarize prediction performance ########
########################################################
########################################################

load("PCA_test_rs.RData")
load("LASSO_test_rs.RData")

# Figure 1:
jpeg("attn_prediction_by_method.jpg",width = 7,height = 9,
     units = "in",res = 400)
par(mfrow=c(3,2))

sns<-c(table(y1.full$site_id_l))[sites]
plot(sns,lasso.test.r,pch=16,col=rgb(.25,.51,.81,1),cex=1.5,
     xlab="site sample size",ylab="predicted vs. actual score correlation",
     ylim=c(0,.6),xlim=c(200,1080),main="1-year full sample",)
points(sns,pca.test.r.y1,pch=17,col=rgb(.57,.82,.31,1),cex=1.3)
rect(1000,-1,10000,1,col = "lightgray")
ci.lasso<-qt(0.975,df=17)*sd(lasso.test.r)/sqrt(18)
ci.pca<-qt(0.975,df=17)*sd(pca.test.r.y1)/sqrt(18)
arrows(x0=1030, y0=mean(lasso.test.r)-ci.lasso, 
       x1=1030, y1=mean(lasso.test.r)+ci.lasso, code=3, 
       angle=90, length=0.06, lwd=1)
arrows(x0=1085, y0=mean(pca.test.r.y1)-ci.pca, 
       x1=1085, y1=mean(pca.test.r.y1)+ci.pca, code=3, 
       angle=90, length=0.06, lwd=1)
points(1030,mean(lasso.test.r),pch=16,col=rgb(.25,.51,.81,1),cex=1.5)
points(1085,mean(pca.test.r.y1),pch=17,col=rgb(.57,.82,.31,1),cex=1.3)


sns<-c(table(y1.attn$site_id_l))[sites]
plot(sns,lasso.test.r.comp,pch=16,col=rgb(.25,.51,.81,1),cex=1.5,
     xlab="site sample size",ylab="predicted vs. actual score correlation",
     ylim=c(0,.6),xlim=c(0,820),main="1-year complete data only")
points(sns,pca.test.r.comp,pch=17,col=rgb(.57,.82,.31,1),cex=1.3)
rect(740,-1,10000,1,col = "lightgray")
ci.lasso<-qt(0.975,df=17)*sd(lasso.test.r.comp)/sqrt(18)
ci.pca<-qt(0.975,df=17)*sd(pca.test.r.comp)/sqrt(18)
arrows(x0=770, y0=mean(lasso.test.r.comp)-ci.lasso, 
       x1=770, y1=mean(lasso.test.r.comp)+ci.lasso, code=3, 
       angle=90, length=0.06, lwd=1)
arrows(x0=825, y0=mean(pca.test.r.comp)-ci.pca, 
       x1=825, y1=mean(pca.test.r.comp)+ci.pca, code=3, 
       angle=90, length=0.06, lwd=1)
points(770,mean(lasso.test.r.comp),pch=16,col=rgb(.25,.51,.81,1),cex=1.5)
points(825,mean(pca.test.r.comp),pch=17,col=rgb(.57,.82,.31,1),cex=1.3)


sns<-c(table(y2.full$site_id_l))[sites]
plot(sns,lasso.test.r.y2,pch=16,col=rgb(.25,.51,.81,1),cex=1.5,
     xlab="site sample size",ylab="predicted vs. actual score correlation",
     ylim=c(0,.6),xlim=c(200,1080),main="2-year full sample")
points(sns,pca.test.r.y2,pch=17,col=rgb(.57,.82,.31,1),cex=1.3)
rect(1000,-1,10000,1,col = "lightgray")
ci.lasso<-qt(0.975,df=17)*sd(lasso.test.r.y2)/sqrt(18)
ci.pca<-qt(0.975,df=17)*sd(pca.test.r.y2)/sqrt(18)
arrows(x0=1030, y0=mean(lasso.test.r.y2)-ci.lasso, 
       x1=1030, y1=mean(lasso.test.r.y2)+ci.lasso, code=3, 
       angle=90, length=0.06, lwd=1)
arrows(x0=1085, y0=mean(pca.test.r.y2)-ci.pca, 
       x1=1085, y1=mean(pca.test.r.y2)+ci.pca, code=3, 
       angle=90, length=0.06, lwd=1)
points(1030,mean(lasso.test.r.y2),pch=16,col=rgb(.25,.51,.81,1),cex=1.5)
points(1085,mean(pca.test.r.y2),pch=17,col=rgb(.57,.82,.31,1),cex=1.3)

plot.new()



sns<-c(table(y1.diff$site_id_l))[sites]
plot(sns,lasso.test.r.y1diff,pch=16,col=rgb(.25,.51,.81,1),cex=1.5,
     xlab="site sample size",ylab="predicted vs. actual score correlation",
     ylim=c(0,.6),xlim=c(200,1080),main="1-year change from baseline")
points(sns,pca.test.r.y1diff,pch=17,col=rgb(.57,.82,.31,1),cex=1.3)
rect(1000,-1,10000,1,col = "lightgray")
ci.lasso<-qt(0.975,df=17)*sd(lasso.test.r.y1diff)/sqrt(18)
ci.pca<-qt(0.975,df=17)*sd(pca.test.r.y1diff)/sqrt(18)
arrows(x0=1030, y0=mean(lasso.test.r.y1diff)-ci.lasso, 
       x1=1030, y1=mean(lasso.test.r.y1diff)+ci.lasso, code=3, 
       angle=90, length=0.06, lwd=1)
arrows(x0=1085, y0=mean(pca.test.r.y1diff)-ci.pca, 
       x1=1085, y1=mean(pca.test.r.y1diff)+ci.pca, code=3, 
       angle=90, length=0.06, lwd=1)
points(1030,mean(lasso.test.r.y1diff),pch=16,col=rgb(.25,.51,.81,1),cex=1.5)
points(1085,mean(pca.test.r.y1diff),pch=17,col=rgb(.57,.82,.31,1),cex=1.3)



sns<-c(table(y2.diff$site_id_l))[sites]
plot(sns,lasso.test.r.y2diff,pch=16,col=rgb(.25,.51,.81,1),cex=1.5,
     xlab="site sample size",ylab="predicted vs. actual score correlation",
     ylim=c(0,.6),xlim=c(200,1080),main="2-year change from baseline")
points(sns,pca.test.r.y2diff,pch=17,col=rgb(.57,.82,.31,1),cex=1.3)
rect(1000,-1,10000,1,col = "lightgray")
ci.lasso<-qt(0.975,df=17)*sd(lasso.test.r.y2diff)/sqrt(18)
ci.pca<-qt(0.975,df=17)*sd(pca.test.r.y2diff)/sqrt(18)
arrows(x0=1030, y0=mean(lasso.test.r.y2diff)-ci.lasso, 
       x1=1030, y1=mean(lasso.test.r.y2diff)+ci.lasso, code=3, 
       angle=90, length=0.06, lwd=1)
arrows(x0=1085, y0=mean(pca.test.r.y2diff)-ci.pca, 
       x1=1085, y1=mean(pca.test.r.y2diff)+ci.pca, code=3, 
       angle=90, length=0.06, lwd=1)
points(1030,mean(lasso.test.r.y2diff),pch=16,col=rgb(.25,.51,.81,1),cex=1.5)
points(1085,mean(pca.test.r.y2diff),pch=17,col=rgb(.57,.82,.31,1),cex=1.3)


dev.off()




jpeg("prediction_top_3_sites_y1.jpg",width = 8,height = 15,
     units = "in",res = 500)
col=c(rgb(.75,0,0,.2),rgb(.11,0.56,1,.2),rgb(.2,.8,.2,.2))
names(col)<-high.n[1:3]
par(mfrow=c(5,3))

for (s in high.n[1:3]){
  actual<-y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")]
  predicted<-pcareg.preds.y1[[s]]$test
  plot(actual,predicted,
       xlab="",ylab="",
       main="",xlim=c(-2,3.6),ylim=c(-1.4,1.4),
       col=col[s],cex.axis=1.25,cex.lab=1.5,cex.main=1.7,pch=16,cex=1.5)
  abline(lm(predicted~actual),lwd=1.5)
  
}

for (s in high.n[1:3]){
  actual<-y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")]
  predicted<-pcareg.preds.comp[[s]]$test
  plot(actual,predicted,
       xlab="",ylab="",
       main="",xlim=c(-2,3.6),ylim=c(-1.4,1.4),
       col=col[s],cex.axis=1.25,cex.lab=1.5,cex.main=1.7,pch=16,cex=1.5)
  abline(lm(predicted~actual),lwd=1.5)
  
}


for (s in high.n[1:3]){
  actual<-y2.full[y2.full$site_id_l==s,paste0("lo.",s,".attn")]
  predicted<-pcareg.preds.y2[[s]]$test
  plot(actual,predicted,
       xlab="",ylab="",
       main="",xlim=c(-2,3.6),ylim=c(-1.4,1.4),
       col=col[s],cex.axis=1.25,cex.lab=1.5,cex.main=1.7,pch=16,cex=1.5)
  abline(lm(predicted~actual),lwd=1.5)
  
}

for (s in high.n[1:3]){
  actual<-y1.diff[y1.diff$site_id_l==s,paste0("lo.",s,".res")]
  predicted<-pcareg.preds.y1diff[[s]]$test
  plot(actual,predicted,
       xlab="",ylab="",
       main="",xlim=c(-3,3),ylim=c(-0.4,0.4),
       col=col[s],cex.axis=1.25,cex.lab=1.5,cex.main=1.7,pch=16,cex=1.5)
  abline(lm(predicted~actual),lwd=1.5)
  
}

for (s in high.n[1:3]){
  actual<-y2.diff[y2.diff$site_id_l==s,paste0("lo.",s,".res")]
  predicted<-pcareg.preds.y2diff[[s]]$test
  plot(actual,predicted,
       xlab="",ylab="",
       main="",xlim=c(-3,3),ylim=c(-0.4,0.4),
       col=col[s],cex.axis=1.25,cex.lab=1.5,cex.main=1.7,pch=16,cex=1.5)
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

COG.imputed.comp<-impute.folds(y1.attn,
                           pred.cols = COG.cols,
                           sites = sites)

SR.imputed.comp<-impute.folds(y1.attn,
                          pred.cols = SR.cols,
                          sites = sites)

DEM.imputed.comp<-impute.folds(y1.attn,
                         pred.cols = DEM.cols,
                         sites = sites)

COG_SR.imputed.comp<-impute.folds(y1.attn,
                             pred.cols = c(COG.cols,SR.cols),
                             sites = sites)

COG_DEM.imputed.comp<-impute.folds(y1.attn,
                          pred.cols = c(COG.cols,DEM.cols),
                          sites = sites)

SR_DEM.imputed.comp<-impute.folds(y1.attn,
                              pred.cols = c(SR.cols,DEM.cols),
                              sites = sites)

 
 save(COG.imputed.comp,SR.imputed.comp,DEM.imputed.comp,
      COG_SR.imputed.comp,COG_DEM.imputed.comp,SR_DEM.imputed.comp,
      file="data_imputations_domain_analysis_comp.RData")

 
 COG.imputed.y1<-impute.folds(y1.full,
                                pred.cols = COG.cols,
                                sites = sites)
 
 SR.imputed.y1<-impute.folds(y1.full,
                               pred.cols = SR.cols,
                               sites = sites)
 
 DEM.imputed.y1<-impute.folds(y1.full,
                                pred.cols = DEM.cols,
                                sites = sites)
 
 COG_SR.imputed.y1<-impute.folds(y1.full,
                                   pred.cols = c(COG.cols,SR.cols),
                                   sites = sites)
 
 COG_DEM.imputed.y1<-impute.folds(y1.full,
                                    pred.cols = c(COG.cols,DEM.cols),
                                    sites = sites)
 
 SR_DEM.imputed.y1<-impute.folds(y1.full,
                                   pred.cols = c(SR.cols,DEM.cols),
                                   sites = sites)
 
 
 save(COG.imputed.y1,SR.imputed.y1,DEM.imputed.y1,
      COG_SR.imputed.y1,COG_DEM.imputed.y1,SR_DEM.imputed.y1,
      file="data_imputations_domain_analysis_y1.RData")



# run LASSO and PCA regression, complete data sample

COG.lasso.comp<-list()
SR.lasso.comp<-list()
DEM.lasso.comp<-list()
COG_SR.lasso.comp<-list()
COG_DEM.lasso.comp<-list()
SR_DEM.lasso.comp<-list()
SR.no.imp.lasso.comp<-list()

COG.pcareg.comp<-list()
SR.pcareg.comp<-list()
DEM.pcareg.comp<-list()
COG_SR.pcareg.comp<-list()
COG_DEM.pcareg.comp<-list()
SR_DEM.pcareg.comp<-list()
SR.no.imp.pcareg.comp<-list()

COG.lasso.comp<-list()
SR.lasso.comp<-list()
DEM.lasso.comp<-list()
COG_SR.lasso.comp<-list()
COG_DEM.lasso.comp<-list()
SR_DEM.lasso.comp<-list()
SR.no.imp.lasso.comp<-list()

COG.pcareg.comp<-list()
SR.pcareg.comp<-list()
DEM.pcareg.comp<-list()
COG_SR.pcareg.comp<-list()
COG_DEM.pcareg.comp<-list()
SR_DEM.pcareg.comp<-list()
SR.no.imp.pcareg.comp<-list()


for (s in sites){

  tmp.dat<-setup.predicton(COG.imputed.comp,s,y1.attn)
  COG.lasso.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  COG.lasso.comp[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  COG.lasso.comp[[s]]$train<-predict(COG.lasso.comp[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = COG.lasso.comp[[s]]$fits$lambda.1se)
  COG.lasso.comp[[s]]$test<-predict(COG.lasso.comp[[s]]$fits,data.matrix(tmp.dat$test), s = COG.lasso.comp[[s]]$fits$lambda.1se)
  COG.lasso.comp[[s]]$train.r<-cor(COG.lasso.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG.lasso.comp[[s]]$test.r<-cor(COG.lasso.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  COG.pcareg.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  COG.pcareg.comp[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  COG.pcareg.comp[[s]]$train<-predict(COG.pcareg.comp[[s]]$fits$winning.model)
  COG.pcareg.comp[[s]]$test<-PCA.predict(COG.pcareg.comp[[s]]$fits,tmp.dat$test)
  COG.pcareg.comp[[s]]$train.r<-cor(COG.pcareg.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG.pcareg.comp[[s]]$test.r<-cor(COG.pcareg.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(SR.imputed.comp,s,y1.attn)
  SR.lasso.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  SR.lasso.comp[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  SR.lasso.comp[[s]]$train<-predict(SR.lasso.comp[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = SR.lasso.comp[[s]]$fits$lambda.1se)
  SR.lasso.comp[[s]]$test<-predict(SR.lasso.comp[[s]]$fits,data.matrix(tmp.dat$test), s = SR.lasso.comp[[s]]$fits$lambda.1se)
  SR.lasso.comp[[s]]$train.r<-cor(SR.lasso.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  SR.lasso.comp[[s]]$test.r<-cor(SR.lasso.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  SR.pcareg.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  SR.pcareg.comp[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  SR.pcareg.comp[[s]]$train<-predict(SR.pcareg.comp[[s]]$fits$winning.model)
  SR.pcareg.comp[[s]]$test<-PCA.predict(SR.pcareg.comp[[s]]$fits,tmp.dat$test)
  SR.pcareg.comp[[s]]$train.r<-cor(SR.pcareg.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  SR.pcareg.comp[[s]]$test.r<-cor(SR.pcareg.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
 
  tmp.dat<-setup.predicton(DEM.imputed.comp,s,y1.attn)
  DEM.lasso.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  DEM.lasso.comp[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  DEM.lasso.comp[[s]]$train<-predict(DEM.lasso.comp[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = DEM.lasso.comp[[s]]$fits$lambda.1se)
  DEM.lasso.comp[[s]]$test<-predict(DEM.lasso.comp[[s]]$fits,data.matrix(tmp.dat$test), s = DEM.lasso.comp[[s]]$fits$lambda.1se)
  DEM.lasso.comp[[s]]$train.r<-cor(DEM.lasso.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  DEM.lasso.comp[[s]]$test.r<-cor(DEM.lasso.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  DEM.pcareg.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  DEM.pcareg.comp[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  DEM.pcareg.comp[[s]]$train<-predict(DEM.pcareg.comp[[s]]$fits$winning.model)
  DEM.pcareg.comp[[s]]$test<-PCA.predict(DEM.pcareg.comp[[s]]$fits,tmp.dat$test)
  DEM.pcareg.comp[[s]]$train.r<-cor(DEM.pcareg.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  DEM.pcareg.comp[[s]]$test.r<-cor(DEM.pcareg.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(COG_SR.imputed.comp,s,y1.attn)
  COG_SR.lasso.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  COG_SR.lasso.comp[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  COG_SR.lasso.comp[[s]]$train<-predict(COG_SR.lasso.comp[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = COG_SR.lasso.comp[[s]]$fits$lambda.1se)
  COG_SR.lasso.comp[[s]]$test<-predict(COG_SR.lasso.comp[[s]]$fits,data.matrix(tmp.dat$test), s = COG_SR.lasso.comp[[s]]$fits$lambda.1se)
  COG_SR.lasso.comp[[s]]$train.r<-cor(COG_SR.lasso.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_SR.lasso.comp[[s]]$test.r<-cor(COG_SR.lasso.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  COG_SR.pcareg.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  COG_SR.pcareg.comp[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  COG_SR.pcareg.comp[[s]]$train<-predict(COG_SR.pcareg.comp[[s]]$fits$winning.model)
  COG_SR.pcareg.comp[[s]]$test<-PCA.predict(COG_SR.pcareg.comp[[s]]$fits,tmp.dat$test)
  COG_SR.pcareg.comp[[s]]$train.r<-cor(COG_SR.pcareg.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_SR.pcareg.comp[[s]]$test.r<-cor(COG_SR.pcareg.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(COG_DEM.imputed.comp,s,y1.attn)
  COG_DEM.lasso.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  COG_DEM.lasso.comp[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  COG_DEM.lasso.comp[[s]]$train<-predict(COG_DEM.lasso.comp[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = COG_DEM.lasso.comp[[s]]$fits$lambda.1se)
  COG_DEM.lasso.comp[[s]]$test<-predict(COG_DEM.lasso.comp[[s]]$fits,data.matrix(tmp.dat$test), s = COG_DEM.lasso.comp[[s]]$fits$lambda.1se)
  COG_DEM.lasso.comp[[s]]$train.r<-cor(COG_DEM.lasso.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_DEM.lasso.comp[[s]]$test.r<-cor(COG_DEM.lasso.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  COG_DEM.pcareg.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  COG_DEM.pcareg.comp[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  COG_DEM.pcareg.comp[[s]]$train<-predict(COG_DEM.pcareg.comp[[s]]$fits$winning.model)
  COG_DEM.pcareg.comp[[s]]$test<-PCA.predict(COG_DEM.pcareg.comp[[s]]$fits,tmp.dat$test)
  COG_DEM.pcareg.comp[[s]]$train.r<-cor(COG_DEM.pcareg.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_DEM.pcareg.comp[[s]]$test.r<-cor(COG_DEM.pcareg.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(SR_DEM.imputed.comp,s,y1.attn)
  SR_DEM.lasso.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  SR_DEM.lasso.comp[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  SR_DEM.lasso.comp[[s]]$train<-predict(SR_DEM.lasso.comp[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = SR_DEM.lasso.comp[[s]]$fits$lambda.1se)
  SR_DEM.lasso.comp[[s]]$test<-predict(SR_DEM.lasso.comp[[s]]$fits,data.matrix(tmp.dat$test), s = SR_DEM.lasso.comp[[s]]$fits$lambda.1se)
  SR_DEM.lasso.comp[[s]]$train.r<-cor(SR_DEM.lasso.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  SR_DEM.lasso.comp[[s]]$test.r<-cor(SR_DEM.lasso.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  SR_DEM.pcareg.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  SR_DEM.pcareg.comp[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  SR_DEM.pcareg.comp[[s]]$train<-predict(SR_DEM.pcareg.comp[[s]]$fits$winning.model)
  SR_DEM.pcareg.comp[[s]]$test<-PCA.predict(SR_DEM.pcareg.comp[[s]]$fits,tmp.dat$test)
  SR_DEM.pcareg.comp[[s]]$train.r<-cor(SR_DEM.pcareg.comp[[s]]$train,y1.attn[y1.attn$site_id_l!=s,paste0("lo.",s,".attn")])
  SR_DEM.pcareg.comp[[s]]$test.r<-cor(SR_DEM.pcareg.comp[[s]]$test,y1.attn[y1.attn$site_id_l==s,paste0("lo.",s,".attn")])
  
}


# run LASSO and PCA regression, imputed y1 sample

COG.lasso.y1<-list()
SR.lasso.y1<-list()
DEM.lasso.y1<-list()
COG_SR.lasso.y1<-list()
COG_DEM.lasso.y1<-list()
SR_DEM.lasso.y1<-list()
SR.no.imp.lasso.y1<-list()

COG.pcareg.y1<-list()
SR.pcareg.y1<-list()
DEM.pcareg.y1<-list()
COG_SR.pcareg.y1<-list()
COG_DEM.pcareg.y1<-list()
SR_DEM.pcareg.y1<-list()
SR.no.imp.pcareg.y1<-list()

COG.lasso.y1<-list()
SR.lasso.y1<-list()
DEM.lasso.y1<-list()
COG_SR.lasso.y1<-list()
COG_DEM.lasso.y1<-list()
SR_DEM.lasso.y1<-list()
SR.no.imp.lasso.y1<-list()

COG.pcareg.y1<-list()
SR.pcareg.y1<-list()
DEM.pcareg.y1<-list()
COG_SR.pcareg.y1<-list()
COG_DEM.pcareg.y1<-list()
SR_DEM.pcareg.y1<-list()
SR.no.imp.pcareg.y1<-list()


for (s in sites){
  
  tmp.dat<-setup.predicton(COG.imputed.y1,s,y1.full)
  COG.lasso.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  COG.lasso.y1[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  COG.lasso.y1[[s]]$train<-predict(COG.lasso.y1[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = COG.lasso.y1[[s]]$fits$lambda.1se)
  COG.lasso.y1[[s]]$test<-predict(COG.lasso.y1[[s]]$fits,data.matrix(tmp.dat$test), s = COG.lasso.y1[[s]]$fits$lambda.1se)
  COG.lasso.y1[[s]]$train.r<-cor(COG.lasso.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  COG.lasso.y1[[s]]$test.r<-cor(COG.lasso.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  COG.pcareg.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  COG.pcareg.y1[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  COG.pcareg.y1[[s]]$train<-predict(COG.pcareg.y1[[s]]$fits$winning.model)
  COG.pcareg.y1[[s]]$test<-PCA.predict(COG.pcareg.y1[[s]]$fits,tmp.dat$test)
  COG.pcareg.y1[[s]]$train.r<-cor(COG.pcareg.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  COG.pcareg.y1[[s]]$test.r<-cor(COG.pcareg.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(SR.imputed.y1,s,y1.full)
  SR.lasso.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  SR.lasso.y1[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  SR.lasso.y1[[s]]$train<-predict(SR.lasso.y1[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = SR.lasso.y1[[s]]$fits$lambda.1se)
  SR.lasso.y1[[s]]$test<-predict(SR.lasso.y1[[s]]$fits,data.matrix(tmp.dat$test), s = SR.lasso.y1[[s]]$fits$lambda.1se)
  SR.lasso.y1[[s]]$train.r<-cor(SR.lasso.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  SR.lasso.y1[[s]]$test.r<-cor(SR.lasso.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  SR.pcareg.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  SR.pcareg.y1[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  SR.pcareg.y1[[s]]$train<-predict(SR.pcareg.y1[[s]]$fits$winning.model)
  SR.pcareg.y1[[s]]$test<-PCA.predict(SR.pcareg.y1[[s]]$fits,tmp.dat$test)
  SR.pcareg.y1[[s]]$train.r<-cor(SR.pcareg.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  SR.pcareg.y1[[s]]$test.r<-cor(SR.pcareg.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(DEM.imputed.y1,s,y1.full)
  DEM.lasso.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  DEM.lasso.y1[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  DEM.lasso.y1[[s]]$train<-predict(DEM.lasso.y1[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = DEM.lasso.y1[[s]]$fits$lambda.1se)
  DEM.lasso.y1[[s]]$test<-predict(DEM.lasso.y1[[s]]$fits,data.matrix(tmp.dat$test), s = DEM.lasso.y1[[s]]$fits$lambda.1se)
  DEM.lasso.y1[[s]]$train.r<-cor(DEM.lasso.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  DEM.lasso.y1[[s]]$test.r<-cor(DEM.lasso.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  DEM.pcareg.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  DEM.pcareg.y1[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  DEM.pcareg.y1[[s]]$train<-predict(DEM.pcareg.y1[[s]]$fits$winning.model)
  DEM.pcareg.y1[[s]]$test<-PCA.predict(DEM.pcareg.y1[[s]]$fits,tmp.dat$test)
  DEM.pcareg.y1[[s]]$train.r<-cor(DEM.pcareg.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  DEM.pcareg.y1[[s]]$test.r<-cor(DEM.pcareg.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(COG_SR.imputed.y1,s,y1.full)
  COG_SR.lasso.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  COG_SR.lasso.y1[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  COG_SR.lasso.y1[[s]]$train<-predict(COG_SR.lasso.y1[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = COG_SR.lasso.y1[[s]]$fits$lambda.1se)
  COG_SR.lasso.y1[[s]]$test<-predict(COG_SR.lasso.y1[[s]]$fits,data.matrix(tmp.dat$test), s = COG_SR.lasso.y1[[s]]$fits$lambda.1se)
  COG_SR.lasso.y1[[s]]$train.r<-cor(COG_SR.lasso.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_SR.lasso.y1[[s]]$test.r<-cor(COG_SR.lasso.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  COG_SR.pcareg.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  COG_SR.pcareg.y1[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  COG_SR.pcareg.y1[[s]]$train<-predict(COG_SR.pcareg.y1[[s]]$fits$winning.model)
  COG_SR.pcareg.y1[[s]]$test<-PCA.predict(COG_SR.pcareg.y1[[s]]$fits,tmp.dat$test)
  COG_SR.pcareg.y1[[s]]$train.r<-cor(COG_SR.pcareg.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_SR.pcareg.y1[[s]]$test.r<-cor(COG_SR.pcareg.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(COG_DEM.imputed.y1,s,y1.full)
  COG_DEM.lasso.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  COG_DEM.lasso.y1[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  COG_DEM.lasso.y1[[s]]$train<-predict(COG_DEM.lasso.y1[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = COG_DEM.lasso.y1[[s]]$fits$lambda.1se)
  COG_DEM.lasso.y1[[s]]$test<-predict(COG_DEM.lasso.y1[[s]]$fits,data.matrix(tmp.dat$test), s = COG_DEM.lasso.y1[[s]]$fits$lambda.1se)
  COG_DEM.lasso.y1[[s]]$train.r<-cor(COG_DEM.lasso.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_DEM.lasso.y1[[s]]$test.r<-cor(COG_DEM.lasso.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  COG_DEM.pcareg.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  COG_DEM.pcareg.y1[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  COG_DEM.pcareg.y1[[s]]$train<-predict(COG_DEM.pcareg.y1[[s]]$fits$winning.model)
  COG_DEM.pcareg.y1[[s]]$test<-PCA.predict(COG_DEM.pcareg.y1[[s]]$fits,tmp.dat$test)
  COG_DEM.pcareg.y1[[s]]$train.r<-cor(COG_DEM.pcareg.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  COG_DEM.pcareg.y1[[s]]$test.r<-cor(COG_DEM.pcareg.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(SR_DEM.imputed.y1,s,y1.full)
  SR_DEM.lasso.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  SR_DEM.lasso.y1[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  SR_DEM.lasso.y1[[s]]$train<-predict(SR_DEM.lasso.y1[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = SR_DEM.lasso.y1[[s]]$fits$lambda.1se)
  SR_DEM.lasso.y1[[s]]$test<-predict(SR_DEM.lasso.y1[[s]]$fits,data.matrix(tmp.dat$test), s = SR_DEM.lasso.y1[[s]]$fits$lambda.1se)
  SR_DEM.lasso.y1[[s]]$train.r<-cor(SR_DEM.lasso.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  SR_DEM.lasso.y1[[s]]$test.r<-cor(SR_DEM.lasso.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  SR_DEM.pcareg.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  SR_DEM.pcareg.y1[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  SR_DEM.pcareg.y1[[s]]$train<-predict(SR_DEM.pcareg.y1[[s]]$fits$winning.model)
  SR_DEM.pcareg.y1[[s]]$test<-PCA.predict(SR_DEM.pcareg.y1[[s]]$fits,tmp.dat$test)
  SR_DEM.pcareg.y1[[s]]$train.r<-cor(SR_DEM.pcareg.y1[[s]]$train,y1.full[y1.full$site_id_l!=s,paste0("lo.",s,".attn")])
  SR_DEM.pcareg.y1[[s]]$test.r<-cor(SR_DEM.pcareg.y1[[s]]$test,y1.full[y1.full$site_id_l==s,paste0("lo.",s,".attn")])
  
}

save.image("domain_analyses.RData")

domain.r.comp<-data.frame(COG.train=unlist(lapply(COG.pcareg.comp,FUN=function(x) x$train.r)),
           SR.train=unlist(lapply(SR.pcareg.comp,FUN=function(x) x$train.r)),
           DEM.train=unlist(lapply(DEM.pcareg.comp,FUN=function(x) x$train.r)),
           COG_SR.train=unlist(lapply(COG_SR.pcareg.comp,FUN=function(x) x$train.r)),
           COG_DEM.train=unlist(lapply(COG_DEM.pcareg.comp,FUN=function(x) x$train.r)),
           SR_DEM.train=unlist(lapply(SR_DEM.pcareg.comp,FUN=function(x) x$train.r)),
           COG.test=unlist(lapply(COG.pcareg.comp,FUN=function(x) x$test.r)),
           SR.test=unlist(lapply(SR.pcareg.comp,FUN=function(x) x$test.r)),
           DEM.test=unlist(lapply(DEM.pcareg.comp,FUN=function(x) x$test.r)),
           COG_SR.test=unlist(lapply(COG_SR.pcareg.comp,FUN=function(x) x$test.r)),
           COG_DEM.test=unlist(lapply(COG_DEM.pcareg.comp,FUN=function(x) x$test.r)),
           SR_DEM.test=unlist(lapply(SR_DEM.pcareg.comp,FUN=function(x) x$test.r)) )

colMeans(domain.r.comp)

domain.ci<- apply(domain.r.comp,2,
      FUN= function(x) qt(0.975,df=17)*sd(x)/sqrt(18))

colMeans(domain.r.comp)+domain.ci
colMeans(domain.r.comp)-domain.ci


domain.r.lasso.comp<-data.frame(COG.train=unlist(lapply(COG.lasso.comp,FUN=function(x) x$train.r)),
                     SR.train=unlist(lapply(SR.lasso.comp,FUN=function(x) x$train.r)),
                     DEM.train=unlist(lapply(DEM.lasso.comp,FUN=function(x) x$train.r)),
                     COG_SR.train=unlist(lapply(COG_SR.lasso.comp,FUN=function(x) x$train.r)),
                     COG_DEM.train=unlist(lapply(COG_DEM.lasso.comp,FUN=function(x) x$train.r)),
                     SR_DEM.train=unlist(lapply(SR_DEM.lasso.comp,FUN=function(x) x$train.r)),
                     COG.test=unlist(lapply(COG.lasso.comp,FUN=function(x) x$test.r)),
                     SR.test=unlist(lapply(SR.lasso.comp,FUN=function(x) x$test.r)),
                     DEM.test=unlist(lapply(DEM.lasso.comp,FUN=function(x) x$test.r)),
                     COG_SR.test=unlist(lapply(COG_SR.lasso.comp,FUN=function(x) x$test.r)),
                     COG_DEM.test=unlist(lapply(COG_DEM.lasso.comp,FUN=function(x) x$test.r)),
                     SR_DEM.test=unlist(lapply(SR_DEM.lasso.comp,FUN=function(x) x$test.r)) )

colMeans(domain.r.lasso.comp)


domain.ci<- apply(domain.r.lasso.comp,2,
                  FUN= function(x) qt(0.975,df=17)*sd(x)/sqrt(18))

colMeans(domain.r.lasso.comp)+domain.ci
colMeans(domain.r.lasso.comp)-domain.ci



domain.r.y1<-data.frame(COG.train=unlist(lapply(COG.pcareg.y1,FUN=function(x) x$train.r)),
                          SR.train=unlist(lapply(SR.pcareg.y1,FUN=function(x) x$train.r)),
                          DEM.train=unlist(lapply(DEM.pcareg.y1,FUN=function(x) x$train.r)),
                          COG_SR.train=unlist(lapply(COG_SR.pcareg.y1,FUN=function(x) x$train.r)),
                          COG_DEM.train=unlist(lapply(COG_DEM.pcareg.y1,FUN=function(x) x$train.r)),
                          SR_DEM.train=unlist(lapply(SR_DEM.pcareg.y1,FUN=function(x) x$train.r)),
                          COG.test=unlist(lapply(COG.pcareg.y1,FUN=function(x) x$test.r)),
                          SR.test=unlist(lapply(SR.pcareg.y1,FUN=function(x) x$test.r)),
                          DEM.test=unlist(lapply(DEM.pcareg.y1,FUN=function(x) x$test.r)),
                          COG_SR.test=unlist(lapply(COG_SR.pcareg.y1,FUN=function(x) x$test.r)),
                          COG_DEM.test=unlist(lapply(COG_DEM.pcareg.y1,FUN=function(x) x$test.r)),
                          SR_DEM.test=unlist(lapply(SR_DEM.pcareg.y1,FUN=function(x) x$test.r)) )

colMeans(domain.r.y1)

domain.ci<- apply(domain.r.y1,2,
                  FUN= function(x) qt(0.975,df=17)*sd(x)/sqrt(18))

colMeans(domain.r.y1)+domain.ci
colMeans(domain.r.y1)-domain.ci


domain.r.lasso.y1<-data.frame(COG.train=unlist(lapply(COG.lasso.y1,FUN=function(x) x$train.r)),
                                SR.train=unlist(lapply(SR.lasso.y1,FUN=function(x) x$train.r)),
                                DEM.train=unlist(lapply(DEM.lasso.y1,FUN=function(x) x$train.r)),
                                COG_SR.train=unlist(lapply(COG_SR.lasso.y1,FUN=function(x) x$train.r)),
                                COG_DEM.train=unlist(lapply(COG_DEM.lasso.y1,FUN=function(x) x$train.r)),
                                SR_DEM.train=unlist(lapply(SR_DEM.lasso.y1,FUN=function(x) x$train.r)),
                                COG.test=unlist(lapply(COG.lasso.y1,FUN=function(x) x$test.r)),
                                SR.test=unlist(lapply(SR.lasso.y1,FUN=function(x) x$test.r)),
                                DEM.test=unlist(lapply(DEM.lasso.y1,FUN=function(x) x$test.r)),
                                COG_SR.test=unlist(lapply(COG_SR.lasso.y1,FUN=function(x) x$test.r)),
                                COG_DEM.test=unlist(lapply(COG_DEM.lasso.y1,FUN=function(x) x$test.r)),
                                SR_DEM.test=unlist(lapply(SR_DEM.lasso.y1,FUN=function(x) x$test.r)) )

colMeans(domain.r.lasso.y1)


domain.ci<- apply(domain.r.lasso.y1,2,
                  FUN= function(x) qt(0.975,df=17)*sd(x)/sqrt(18))

colMeans(domain.r.lasso.y1)+domain.ci
colMeans(domain.r.lasso.y1)-domain.ci

#######################################################################
#######################################################################
##### Sensitivity analysis: medication exclusions #####################
#######################################################################
#######################################################################

## Create fully imputed training and test data sets for each fold #####


nomeds.imputed.comp<-impute.folds(y1.attn[y1.attn$CAS.ADHD==FALSE & 
                                       !is.na(y1.attn$CAS.ADHD),],
                            pred.cols = pred.cols,
                            sites = sites)

nomeds.imputed.y1<-impute.folds(y1.full[y1.full$CAS.ADHD==FALSE & 
                                            !is.na(y1.full$CAS.ADHD),],
                                  pred.cols = pred.cols,
                                  sites = sites)

save(nomeds.imputed.full,nomeds.imputed.y1,file="data_imputations_nomeds_sens.RData")

# predictive modeling

y1.attn.nomeds<-y1.attn[y1.attn$CAS.ADHD==FALSE & 
                          !is.na(y1.attn$CAS.ADHD),]
y1.full.nomeds<-y1.full[y1.full$CAS.ADHD==FALSE & 
                          !is.na(y1.full$CAS.ADHD),]

nomeds.lasso.comp<-list()
nomeds.pcareg.comp<-list()

nomeds.lasso.y1<-list()
nomeds.pcareg.y1<-list()


for (s in sites){
  
  tmp.dat<-setup.predicton(nomeds.imputed.comp,s,y1.attn.nomeds)
  nomeds.lasso.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  nomeds.lasso.comp[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), 
                                    y1.attn.nomeds[y1.attn.nomeds$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  nomeds.lasso.comp[[s]]$train<-predict(nomeds.lasso.comp[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = nomeds.lasso.comp[[s]]$fits$lambda.1se)
  nomeds.lasso.comp[[s]]$test<-predict(nomeds.lasso.comp[[s]]$fits,data.matrix(tmp.dat$test), s = nomeds.lasso.comp[[s]]$fits$lambda.1se)
  nomeds.lasso.comp[[s]]$train.r<-cor(nomeds.lasso.comp[[s]]$train,y1.attn.nomeds[y1.attn.nomeds$site_id_l!=s,paste0("lo.",s,".attn")])
  nomeds.lasso.comp[[s]]$test.r<-cor(nomeds.lasso.comp[[s]]$test,y1.attn.nomeds[y1.attn.nomeds$site_id_l==s,paste0("lo.",s,".attn")])
  nomeds.pcareg.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  nomeds.pcareg.comp[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  nomeds.pcareg.comp[[s]]$train<-predict(nomeds.pcareg.comp[[s]]$fits$winning.model)
  nomeds.pcareg.comp[[s]]$test<-PCA.predict(nomeds.pcareg.comp[[s]]$fits,tmp.dat$test)
  nomeds.pcareg.comp[[s]]$train.r<-cor(nomeds.pcareg.comp[[s]]$train,
                                  y1.attn.nomeds[y1.attn.nomeds$site_id_l!=s,paste0("lo.",s,".attn")])
  nomeds.pcareg.comp[[s]]$test.r<-cor(nomeds.pcareg.comp[[s]]$test,
                                 y1.attn.nomeds[y1.attn.nomeds$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(nomeds.imputed.y1,s,y1.full.nomeds)
  nomeds.lasso.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  nomeds.lasso.y1[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), 
                                         y1.full.nomeds[y1.full.nomeds$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  nomeds.lasso.y1[[s]]$train<-predict(nomeds.lasso.y1[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = nomeds.lasso.y1[[s]]$fits$lambda.1se)
  nomeds.lasso.y1[[s]]$test<-predict(nomeds.lasso.y1[[s]]$fits,data.matrix(tmp.dat$test), s = nomeds.lasso.y1[[s]]$fits$lambda.1se)
  nomeds.lasso.y1[[s]]$train.r<-cor(nomeds.lasso.y1[[s]]$train,y1.full.nomeds[y1.full.nomeds$site_id_l!=s,paste0("lo.",s,".attn")])
  nomeds.lasso.y1[[s]]$test.r<-cor(nomeds.lasso.y1[[s]]$test,y1.full.nomeds[y1.full.nomeds$site_id_l==s,paste0("lo.",s,".attn")])
  nomeds.pcareg.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  nomeds.pcareg.y1[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  nomeds.pcareg.y1[[s]]$train<-predict(nomeds.pcareg.y1[[s]]$fits$winning.model)
  nomeds.pcareg.y1[[s]]$test<-PCA.predict(nomeds.pcareg.y1[[s]]$fits,tmp.dat$test)
  nomeds.pcareg.y1[[s]]$train.r<-cor(nomeds.pcareg.y1[[s]]$train,
                                       y1.full.nomeds[y1.full.nomeds$site_id_l!=s,paste0("lo.",s,".attn")])
  nomeds.pcareg.y1[[s]]$test.r<-cor(nomeds.pcareg.y1[[s]]$test,
                                      y1.full.nomeds[y1.full.nomeds$site_id_l==s,paste0("lo.",s,".attn")])
}


save(nomeds.lasso.comp,nomeds.lasso.y1,
     nomeds.pcareg.comp,nomeds.pcareg.y1,file="medication_sensitivity_analyses.RData")


meds.r.comp<-data.frame(nomeds.train.pca=unlist(lapply(nomeds.pcareg.comp,FUN=function(x) x$train.r)),
                   nomeds.test.pca=unlist(lapply(nomeds.pcareg.comp,FUN=function(x) x$test.r)),
                   nomeds.train.lasso=unlist(lapply(nomeds.lasso.comp,FUN=function(x) x$train.r)),
                   nomeds.test.lasso=unlist(lapply(nomeds.lasso.comp,FUN=function(x) x$test.r)) )

meds.r.y1<-data.frame(nomeds.train.pca=unlist(lapply(nomeds.pcareg.y1,FUN=function(x) x$train.r)),
                        nomeds.test.pca=unlist(lapply(nomeds.pcareg.y1,FUN=function(x) x$test.r)),
                        nomeds.train.lasso=unlist(lapply(nomeds.lasso.y1,FUN=function(x) x$train.r)),
                        nomeds.test.lasso=unlist(lapply(nomeds.lasso.y1,FUN=function(x) x$test.r)) )

colMeans(meds.r.y1)

ci.lasso.nomeds.y1<-qt(0.975,df=17)*sd(meds.r.y1$nomeds.test.lasso)/sqrt(18)
ci.pca.nomeds.y1<-qt(0.975,df=17)*sd(meds.r.y1$nomeds.test.pca)/sqrt(18)

colMeans(meds.r.y1)["nomeds.test.lasso"]+ci.lasso.nomeds.y1
colMeans(meds.r.y1)["nomeds.test.lasso"]-ci.lasso.nomeds.y1

colMeans(meds.r.y1)["nomeds.test.pca"]+ci.pca.nomeds.y1
colMeans(meds.r.y1)["nomeds.test.pca"]-ci.pca.nomeds.y1


colMeans(meds.r.comp)

ci.lasso.nomeds.comp<-qt(0.975,df=17)*sd(meds.r.comp$nomeds.test.lasso)/sqrt(18)
ci.pca.nomeds.comp<-qt(0.975,df=17)*sd(meds.r.comp$nomeds.test.pca)/sqrt(18)

colMeans(meds.r.comp)["nomeds.test.lasso"]+ci.lasso.nomeds.comp
colMeans(meds.r.comp)["nomeds.test.lasso"]-ci.lasso.nomeds.comp

colMeans(meds.r.comp)["nomeds.test.pca"]+ci.pca.nomeds.comp
colMeans(meds.r.comp)["nomeds.test.pca"]-ci.pca.nomeds.comp



#######################################################################
#######################################################################
##### Plot importance/weights of features #############################
#######################################################################
#######################################################################

load("sparse_model_selection.RData")
load("PCA_weights.RData")

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
names(fw.names)<-names(fw.means.y1)
fw.cols=c(rep("cornflowerblue",11),
          rep("red2",17),rep("orange",19))
names(fw.cols)<-names(fw.means.y1)

jpeg("con_feature_weights_y1.jpg",width = 8,height = 10,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fw.means.y1[names(sort(abs(fw.means.y1),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fw.means.y1),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fw.means.y1),decreasing = FALSE))],
        border="white",xlim = c(-0.18,0.18),width=0.5)
abline(0,10e10)
dev.off()

jpeg("con_feature_weights_y1_complete.jpg",width = 8,height = 10,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fw.means.comp[names(sort(abs(fw.means.comp),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fw.means.comp),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fw.means.comp),decreasing = FALSE))],
        border="white",xlim = c(-0.18,0.18),width=0.5)
abline(0,10e10)
dev.off()

jpeg("con_feature_weights_y2.jpg",width = 8,height = 10,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fw.means.y2[names(sort(abs(fw.means.y2),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fw.means.y2),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fw.means.y2),decreasing = FALSE))],
        border="white",xlim = c(-0.18,0.18),width=0.5)
abline(0,10e10)
dev.off()

jpeg("con_feature_weights_y1diff.jpg",width = 8,height = 10,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fw.means.y1diff[names(sort(abs(fw.means.y1diff),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fw.means.y1diff),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fw.means.y1diff),decreasing = FALSE))],
        border="white",xlim = c(-0.18,0.18),width=0.5)
abline(0,10e10)
dev.off()

jpeg("con_feature_weights_y2diff.jpg",width = 8,height = 10,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fw.means.y2diff[names(sort(abs(fw.means.y2diff),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fw.means.y2diff),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fw.means.y2diff),decreasing = FALSE))],
        border="white",xlim = c(-0.18,0.18),width=0.5)
abline(0,10e10)
dev.off()


cpw.names<-c("Female sex", "Asian race", "Black race", "Other/mixed race",
             "Hispanic ethnicity", "Parent education <HS", 
             "Parent undergraduate degree", "Parent post-graduate degree",
             "Parent some college education", "Parents not married",
             "Parental income <$50k", "Parental income >=$100k",
             "Cash choice - delayed")
names(cpw.names)<-names(cpw.means.y1)
cpw.cols=c(rep("cornflowerblue",12),
          rep("orange",1))
names(cpw.cols)<-names(cpw.means.y1)

jpeg("cat_feature_weights_y1.jpg",width = 8,height = 4,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(cpw.means.y1[names(sort(abs(cpw.means.y1),decreasing = FALSE))],horiz = TRUE,
        names.arg = cpw.names[names(sort(abs(cpw.means.y1),decreasing = FALSE))],
        col=cpw.cols[names(sort(abs(cpw.means.y1),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()

jpeg("cat_feature_weights_y1_comp.jpg",width = 8,height = 4,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(cpw.means.comp[names(sort(abs(cpw.means.comp),decreasing = FALSE))],horiz = TRUE,
        names.arg = cpw.names[names(sort(abs(cpw.means.comp),decreasing = FALSE))],
        col=cpw.cols[names(sort(abs(cpw.means.comp),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()


jpeg("cat_feature_weights_y2.jpg",width = 8,height = 4,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(cpw.means.y2[names(sort(abs(cpw.means.y2),decreasing = FALSE))],horiz = TRUE,
        names.arg = cpw.names[names(sort(abs(cpw.means.y2),decreasing = FALSE))],
        col=cpw.cols[names(sort(abs(cpw.means.y2),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()


jpeg("cat_feature_weights_y1diff.jpg",width = 8,height = 4,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(cpw.means.y1diff[names(sort(abs(cpw.means.y1diff),decreasing = FALSE))],horiz = TRUE,
        names.arg = cpw.names[names(sort(abs(cpw.means.y1diff),decreasing = FALSE))],
        col=cpw.cols[names(sort(abs(cpw.means.y1diff),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()



jpeg("cat_feature_weights_y2diff.jpg",width = 8,height = 4,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(cpw.means.y2diff[names(sort(abs(cpw.means.y2diff),decreasing = FALSE))],horiz = TRUE,
        names.arg = cpw.names[names(sort(abs(cpw.means.y2diff),decreasing = FALSE))],
        col=cpw.cols[names(sort(abs(cpw.means.y2diff),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()



cpwl.y1<-colMeans(betaSpar.y1[names(betaSpar.y1)%in%names(cpw.means.y1)])
fwl.y1<-colMeans(betaSpar.y1[!names(betaSpar.y1)%in%c(names(cpw.means.y1),"sites","(Intercept)")])

jpeg("con_lasso_feature_weights_y1.jpg",width = 8,height = 4.2,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fwl.y1[names(sort(abs(fwl.y1),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fwl.y1),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fwl.y1),decreasing = FALSE))],
        border="white",xlim = c(-0.15,0.15),width=0.5)
abline(0,10e10)
dev.off()

cpwl.comp<-colMeans(betaSpar.comp[names(betaSpar.comp)%in%names(cpw.means.y1)])
fwl.comp<-colMeans(betaSpar.comp[!names(betaSpar.comp)%in%c(names(cpw.means.y1),"sites","(Intercept)")])

jpeg("con_lasso_feature_weights_y1_comp.jpg",width = 8,height = 4.2,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fwl.comp[names(sort(abs(fwl.comp),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fwl.comp),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fwl.comp),decreasing = FALSE))],
        border="white",xlim = c(-0.15,0.15),width=0.5)
abline(0,10e10)
dev.off()

cpwl.y2<-colMeans(betaSpar.y2[names(betaSpar.y2)%in%names(cpw.means.y1)])
fwl.y2<-colMeans(betaSpar.y2[!names(betaSpar.y2)%in%c(names(cpw.means.y1),"sites","(Intercept)")])

jpeg("con_lasso_feature_weights_y2.jpg",width = 8,height = 4.2,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fwl.y2[names(sort(abs(fwl.y2),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fwl.y2),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fwl.y2),decreasing = FALSE))],
        border="white",xlim = c(-0.15,0.15),width=0.5)
abline(0,10e10)
dev.off()



cpwl.y1diff<-colMeans(betaSpar.y1diff[names(betaSpar.y1diff)%in%names(cpw.means.y1diff)])
fwl.y1diff<-colMeans(betaSpar.y1diff[!names(betaSpar.y1diff)%in%c(names(cpw.means.y1diff),"sites","(Intercept)")])

jpeg("con_lasso_feature_weights_y1diff.jpg",width = 8,height = 3,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fwl.y1diff[names(sort(abs(fwl.y1diff),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fwl.y1diff),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fwl.y1diff),decreasing = FALSE))],
        border="white",xlim = c(-0.15,0.15),width=0.5)
abline(0,10e10)
dev.off()


cpwl.y2diff<-colMeans(betaSpar.y2diff[names(betaSpar.y2diff)%in%names(cpw.means.y2diff)])
fwl.y2diff<-colMeans(betaSpar.y2diff[!names(betaSpar.y2diff)%in%c(names(cpw.means.y2diff),"sites","(Intercept)")])

jpeg("con_lasso_feature_weights_y2diff.jpg",width = 8,height = 2.7,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(fwl.y2diff[names(sort(abs(fwl.y2diff),decreasing = FALSE))],horiz = TRUE,
        names.arg = fw.names[names(sort(abs(fwl.y2diff),decreasing = FALSE))],
        col=fw.cols[names(sort(abs(fwl.y2diff),decreasing = FALSE))],
        border="white",xlim = c(-0.15,0.15),width=0.5)
abline(0,10e10)
dev.off()




jpeg("cat_lasso_feature_weights_y1.jpg",width = 8,height = 2.1,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(c(0,cpwl.y1[names(sort(abs(cpwl.y1),decreasing = FALSE))]),horiz = TRUE,
        names.arg = c("",cpw.names[names(sort(abs(cpwl.y1),decreasing = FALSE))]),
        col=cpw.cols[names(sort(abs(cpwl.y1),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()


jpeg("cat_lasso_feature_weights_y1_comp.jpg",width = 8,height = 2.3,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(c(0,cpwl.comp[names(sort(abs(cpwl.comp),decreasing = FALSE))]),horiz = TRUE,
        names.arg = c("",cpw.names[names(sort(abs(cpwl.comp),decreasing = FALSE))]),
        col=cpw.cols[names(sort(abs(cpwl.comp),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()


jpeg("cat_lasso_feature_weights_y2.jpg",width = 8,height = 2.1,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(c(0,cpwl.y2[names(sort(abs(cpwl.y2),decreasing = FALSE))]),horiz = TRUE,
        names.arg = c("",cpw.names[names(sort(abs(cpwl.y2),decreasing = FALSE))]),
        col=cpw.cols[names(sort(abs(cpwl.y2),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()


jpeg("cat_lasso_feature_weights_y1diff.jpg",width = 8,height = 2.1,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(c(0,cpwl.y1diff[names(sort(abs(cpwl.y1diff),decreasing = FALSE))]),horiz = TRUE,
        names.arg = c("",cpw.names[names(sort(abs(cpwl.y1diff),decreasing = FALSE))]),
        col=cpw.cols[names(sort(abs(cpwl.y1diff),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()

jpeg("cat_lasso_feature_weights_y2diff.jpg",width = 8,height = 2.5,
     units = "in",res = 500)
par(las=2,mar=c(5,15,4,2))
barplot(c(0,cpwl.y2diff[names(sort(abs(cpwl.y2diff),decreasing = FALSE))]),horiz = TRUE,
        names.arg = c("",cpw.names[names(sort(abs(cpwl.y2diff),decreasing = FALSE))]),
        col=cpw.cols[names(sort(abs(cpwl.y2diff),decreasing = FALSE))],
        border="white",xlim = c(-0.30,0.30),width=0.5)
abline(0,10e10)
dev.off()



#######################################################################
#######################################################################
##### Look at prediction in lockbox sites #############################
#######################################################################
#######################################################################

# load in and bind with existing data

comp.validation<-read.csv("year1_validation_dat.csv",stringsAsFactors = TRUE)

comp.validation<-rbind(comp.validation,
                     read.csv("year1_comp_dat.csv",stringsAsFactors = TRUE))

y1.validation<-read.csv("year1_full_validation_dat.csv",stringsAsFactors = TRUE)

y1.validation<-rbind(y1.validation,
                     read.csv("year1_full_dat.csv",stringsAsFactors = TRUE))

y2.validation<-read.csv("year2_full_validation_dat.csv",stringsAsFactors = TRUE)

y2.validation<-rbind(y2.validation,
                     read.csv("year2_full_dat.csv",stringsAsFactors = TRUE))


v.sites<-c("site10","site06","site20")

# re-level and name reference levels for each factor
comp.validation$female <- relevel(comp.validation$female, ref = 'no')
comp.validation$race.4level <- relevel(comp.validation$race.4level, ref = 'White')
comp.validation$hisp <- relevel(comp.validation$hisp, ref = 'no')
comp.validation$high.educ <- relevel(comp.validation$high.educ, ref = 'HS Diploma/GED')
comp.validation$married <- relevel(comp.validation$married, ref = 'yes')
comp.validation$household.income <- relevel(comp.validation$household.income, ref = '[>=50K & <100K]')
comp.validation$cash_choice_task <- relevel(comp.validation$cash_choice_task, ref = 'no')

y1.validation$female <- relevel(y1.validation$female, ref = 'no')
y1.validation$race.4level <- relevel(y1.validation$race.4level, ref = 'White')
y1.validation$hisp <- relevel(y1.validation$hisp, ref = 'no')
y1.validation$high.educ <- relevel(y1.validation$high.educ, ref = 'HS Diploma/GED')
y1.validation$married <- relevel(y1.validation$married, ref = 'yes')
y1.validation$household.income <- relevel(y1.validation$household.income, ref = '[>=50K & <100K]')
y1.validation$cash_choice_task <- relevel(y1.validation$cash_choice_task, ref = 'no')


y2.validation$female <- relevel(y2.validation$female, ref = 'no')
y2.validation$race.4level <- relevel(y2.validation$race.4level, ref = 'White')
y2.validation$hisp <- relevel(y2.validation$hisp, ref = 'no')
y2.validation$high.educ <- relevel(y2.validation$high.educ, ref = 'HS Diploma/GED')
y2.validation$married <- relevel(y2.validation$married, ref = 'yes')
y2.validation$household.income <- relevel(y2.validation$household.income, ref = '[>=50K & <100K]')
y2.validation$cash_choice_task <- relevel(y2.validation$cash_choice_task, ref = 'no')


# generate attention problems scores for the lockbox sites

training.models.v.comp<-list()

t<-1

for (s in v.sites){
  
  t.mod = cfa(bf.final,comp.validation[comp.validation$site_id_l!=s,],  
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  training.models.v.comp[[t]]<-t.mod 
  
  t<-(t+1)
  
}

names(training.models.v.comp)<-paste0("lo.",v.sites)

for (s in v.sites){
  comp.validation[,paste0("lo.",s,".attn")]<-lavPredict(training.models.v.comp[[paste0("lo.",s)]],                                              
                                                      newdata = comp.validation)[,"Attn"]
}



training.models.v.y1<-list()
traincomp.imp.y1<-list()

t<-1

for (s in v.sites){
  
  tmp.train.imp <- mice(data = y1.validation[y1.validation$site_id_l!=s,all.c],m = 1,print = FALSE)
  traincomp.imp.y1[[t]] <- suppressWarnings(mice.reuse(tmp.train.imp, 
                                                       y1.validation[,all.c], 
                                                       maxit = 1,print = FALSE))$`1`
  
  t.mod = cfa(bf.final,y1.validation[y1.validation$site_id_l!=s,],  
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  training.models.v.y1[[t]]<-t.mod 
  
  y1.validation[,paste0("lo.",s,".attn")]<-lavPredict(training.models.v.y1[[t]],                                              
                                                        newdata =traincomp.imp.y1[[t]])[,"Attn"]
  
  t<-(t+1)
  
}



training.models.v.y2<-list()
traincomp.imp.y2<-list()

t<-1

for (s in v.sites){
  
  tmp.train.imp <- mice(data = y2.validation[y2.validation$site_id_l!=s,all.c],m = 1,print = FALSE)
  traincomp.imp.y2[[t]] <- suppressWarnings(mice.reuse(tmp.train.imp, 
                                                       y2.validation[,all.c], 
                                                       maxit = 1,print = FALSE))$`1`
  
  t.mod = cfa(bf.final,y2.validation[y1.validation$site_id_l!=s,],  
              std.lv=TRUE,
              orthogonal=TRUE,
              estimator = "WLSMV")
  training.models.v.y2[[t]]<-t.mod 
  
  y2.validation[,paste0("lo.",s,".attn")]<-lavPredict(training.models.v.y2[[t]],                                              
                                                      newdata =traincomp.imp.y2[[t]])[,"Attn"]
  
  t<-(t+1)
  
}





# prediction

validation.imputed.comp<-impute.folds(comp.validation,
                             pred.cols = pred.cols,
                             sites = v.sites)

validation.imputed.y1<-impute.folds(y1.validation,
                                      pred.cols = pred.cols,
                                      sites = v.sites)


validation.imputed.y2<-impute.folds(y2.validation,
                                    pred.cols = pred.cols,
                                    sites = v.sites)

 save(validation.imputed.comp,validation.imputed.y1,validation.imputed.y2,
 file="data_imputations_validation.RData")

validation.lasso.comp<-list()
validation.pcareg.comp<-list()
validation.lasso.y1<-list()
validation.pcareg.y1<-list()
validation.lasso.y2<-list()
validation.pcareg.y2<-list()

for (s in v.sites){
  
  tmp.dat<-setup.predicton(validation.imputed.comp,s,comp.validation)
  validation.lasso.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  validation.lasso.comp[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), 
                                        comp.validation[comp.validation$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  validation.lasso.comp[[s]]$train<-predict(validation.lasso.comp[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = validation.lasso.comp[[s]]$fits$lambda.1se)
  validation.lasso.comp[[s]]$test<-predict(validation.lasso.comp[[s]]$fits,data.matrix(tmp.dat$test), s = validation.lasso.comp[[s]]$fits$lambda.1se)
  validation.lasso.comp[[s]]$train.r<-cor(validation.lasso.comp[[s]]$train,comp.validation[comp.validation$site_id_l!=s,paste0("lo.",s,".attn")])
  validation.lasso.comp[[s]]$test.r<-cor(validation.lasso.comp[[s]]$test,comp.validation[comp.validation$site_id_l==s,paste0("lo.",s,".attn")])
  validation.pcareg.comp[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  validation.pcareg.comp[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  validation.pcareg.comp[[s]]$train<-predict(validation.pcareg.comp[[s]]$fits$winning.model)
  validation.pcareg.comp[[s]]$test<-PCA.predict(validation.pcareg.comp[[s]]$fits,tmp.dat$test)
  validation.pcareg.comp[[s]]$train.r<-cor(validation.pcareg.comp[[s]]$train,
                                      comp.validation[comp.validation$site_id_l!=s,paste0("lo.",s,".attn")])
  validation.pcareg.comp[[s]]$test.r<-cor(validation.pcareg.comp[[s]]$test,
                                      comp.validation[comp.validation$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(validation.imputed.y1,s,y1.validation)
  validation.lasso.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  validation.lasso.y1[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), 
                                             y1.validation[y1.validation$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  validation.lasso.y1[[s]]$train<-predict(validation.lasso.y1[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = validation.lasso.y1[[s]]$fits$lambda.1se)
  validation.lasso.y1[[s]]$test<-predict(validation.lasso.y1[[s]]$fits,data.matrix(tmp.dat$test), s = validation.lasso.y1[[s]]$fits$lambda.1se)
  validation.lasso.y1[[s]]$train.r<-cor(validation.lasso.y1[[s]]$train,y1.validation[y1.validation$site_id_l!=s,paste0("lo.",s,".attn")])
  validation.lasso.y1[[s]]$test.r<-cor(validation.lasso.y1[[s]]$test,y1.validation[y1.validation$site_id_l==s,paste0("lo.",s,".attn")])
  validation.pcareg.y1[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  validation.pcareg.y1[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  validation.pcareg.y1[[s]]$train<-predict(validation.pcareg.y1[[s]]$fits$winning.model)
  validation.pcareg.y1[[s]]$test<-PCA.predict(validation.pcareg.y1[[s]]$fits,tmp.dat$test)
  validation.pcareg.y1[[s]]$train.r<-cor(validation.pcareg.y1[[s]]$train,
                                           y1.validation[y1.validation$site_id_l!=s,paste0("lo.",s,".attn")])
  validation.pcareg.y1[[s]]$test.r<-cor(validation.pcareg.y1[[s]]$test,
                                          y1.validation[y1.validation$site_id_l==s,paste0("lo.",s,".attn")])
  
  tmp.dat<-setup.predicton(validation.imputed.y2,s,y2.validation)
  validation.lasso.y2[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA)    
  validation.lasso.y2[[s]]$fits<-cv.glmnet(data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), 
                                           y2.validation[y2.validation$site_id_l!=s,paste0("lo.",s,".attn")], alpha = 1)
  validation.lasso.y2[[s]]$train<-predict(validation.lasso.y2[[s]]$fits,data.matrix(tmp.dat$train[,colnames(tmp.dat$train)!="attn"]), s = validation.lasso.y2[[s]]$fits$lambda.1se)
  validation.lasso.y2[[s]]$test<-predict(validation.lasso.y2[[s]]$fits,data.matrix(tmp.dat$test), s = validation.lasso.y2[[s]]$fits$lambda.1se)
  validation.lasso.y2[[s]]$train.r<-cor(validation.lasso.y2[[s]]$train,y2.validation[y2.validation$site_id_l!=s,paste0("lo.",s,".attn")])
  validation.lasso.y2[[s]]$test.r<-cor(validation.lasso.y2[[s]]$test,y2.validation[y2.validation$site_id_l==s,paste0("lo.",s,".attn")])
  validation.pcareg.y2[[s]]<-list(fits=NA,train=NA, test=NA,train.r=NA,test.r=NA) 
  validation.pcareg.y2[[s]]$fits<-PCA.reg(data=tmp.dat$train, outcome = "attn", PCA.cols = tmp.dat$con.vars)
  validation.pcareg.y2[[s]]$train<-predict(validation.pcareg.y2[[s]]$fits$winning.model)
  validation.pcareg.y2[[s]]$test<-PCA.predict(validation.pcareg.y2[[s]]$fits,tmp.dat$test)
  validation.pcareg.y2[[s]]$train.r<-cor(validation.pcareg.y2[[s]]$train,
                                         y2.validation[y2.validation$site_id_l!=s,paste0("lo.",s,".attn")])
  validation.pcareg.y2[[s]]$test.r<-cor(validation.pcareg.y2[[s]]$test,
                                        y2.validation[y2.validation$site_id_l==s,paste0("lo.",s,".attn")])
  
  
}



validation.r<-data.frame(validation.train.pca.comp=unlist(lapply(validation.pcareg.comp,FUN=function(x) x$train.r)),
                         validation.test.pca.comp=unlist(lapply(validation.pcareg.comp,FUN=function(x) x$test.r)),
                         validation.train.lasso.comp=unlist(lapply(validation.lasso.comp,FUN=function(x) x$train.r)),
                         validation.test.lasso.comp=unlist(lapply(validation.lasso.comp,FUN=function(x) x$test.r)), 
                         validation.train.pca.y1=unlist(lapply(validation.pcareg.y1,FUN=function(x) x$train.r)),
                         validation.test.pca.y1=unlist(lapply(validation.pcareg.y1,FUN=function(x) x$test.r)),
                         validation.train.lasso.y1=unlist(lapply(validation.lasso.y1,FUN=function(x) x$train.r)),
                         validation.test.lasso.y1=unlist(lapply(validation.lasso.y1,FUN=function(x) x$test.r)), 
                         validation.train.pca.y2=unlist(lapply(validation.pcareg.y2,FUN=function(x) x$train.r)),
                         validation.test.pca.y2=unlist(lapply(validation.pcareg.y2,FUN=function(x) x$test.r)),
                         validation.train.lasso.y2=unlist(lapply(validation.lasso.y2,FUN=function(x) x$train.r)),
                         validation.test.lasso.y2=unlist(lapply(validation.lasso.y2,FUN=function(x) x$test.r)) )
validation.r

# also try sparse model

sparse.imputed.v<-impute.folds(y1.validation,
                             pred.cols = sparse.y1,
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

save.image("validation_analyses.RData")


#############################################################################
