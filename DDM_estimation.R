rm(list=ls())
source ("dmc/dmc.R")
source("read_abcd_r1.R")

#####################################################################
#####################################################################
############## format data and select subsamples ####################
#####################################################################
#####################################################################


################################################################
####### format n-back data and apply exclusion criteria ########
################################################################

# load trial-level n-back data, merged across all 
# baseline participants (once downloaded from aws)
abcd.nback<-read.csv("nback_revised.csv")

# make DMC-based columns

abcd.nback$S<-as.character(abcd.nback$enback_targettype)
abcd.nback[abcd.nback$S=="target",]$S<-"tar"
abcd.nback[abcd.nback$S=="lure",]$S<-"lur"
abcd.nback[abcd.nback$S=="nonlure",]$S<-"non"
abcd.nback$S<-factor(abcd.nback$S)

abcd.nback$R<-NA
abcd.nback[abcd.nback$S%in%c("tar") & abcd.nback$enback_stim_acc==1 & 
                !is.na(abcd.nback$enback_stim_resp),]$R<-"YES"
abcd.nback[abcd.nback$S%in%c("tar") & abcd.nback$enback_stim_acc==0 & 
                !is.na(abcd.nback$enback_stim_resp),]$R<-"NO"
abcd.nback[abcd.nback$S%in%c("lur","non") & abcd.nback$enback_stim_acc==1 & 
                !is.na(abcd.nback$enback_stim_resp),]$R<-"NO"
abcd.nback[abcd.nback$S%in%c("lur","non") & abcd.nback$enback_stim_acc==0 & 
                !is.na(abcd.nback$enback_stim_resp),]$R<-"YES"
abcd.nback$R<-factor(abcd.nback$R)

abcd.nback$RT<-NA
abcd.nback[!is.na(abcd.nback$R),]$RT<-abcd.nback[!is.na(abcd.nback$R),]$enback_stim_rt/1000

abcd.nback$acc<-NA
abcd.nback[!is.na(abcd.nback$R),]$acc<-abcd.nback[!is.na(abcd.nback$R),]$enback_stim_acc

# separate into baseline and year 2

nback.base<-abcd.nback[abcd.nback$eventname=="baseline_year_1_arm_1",]
nback.y2<-abcd.nback[abcd.nback$eventname=="2_year_follow_up_y_arm_1",]

#compute summary stats

sum.base<-data.frame(unique(nback.base$subject));colnames(sum.base)<-"s"
sum.y2<-data.frame(unique(nback.y2$subject));colnames(sum.y2)<-"s"

sum.base$mrt.0.tar<-NA
sum.base$mrt.0.lur<-NA
sum.base$mrt.0.non<-NA
sum.base$sdrt.0.tar<-NA
sum.base$sdrt.0.lur<-NA
sum.base$sdrt.0.non<-NA
sum.base$acc.0.tar<-NA
sum.base$acc.0.lur<-NA
sum.base$acc.0.non<-NA

sum.base$acc.0.total<-NA
sum.base$p.omit.0<-NA

sum.base$mrt.2.tar<-NA
sum.base$mrt.2.lur<-NA
sum.base$mrt.2.non<-NA
sum.base$sdrt.2.tar<-NA
sum.base$sdrt.2.lur<-NA
sum.base$sdrt.2.non<-NA
sum.base$acc.2.tar<-NA
sum.base$acc.2.lur<-NA
sum.base$acc.2.non<-NA

sum.base$acc.2.total<-NA
sum.base$p.omit.2<-NA

for (r in 1:length(sum.base$s)){
  tmp0<-nback.base[nback.base$subject==sum.base$s[r] & nback.base$enback_loadcon=="0-Back",]
  tmp2<-nback.base[nback.base$subject==sum.base$s[r] & nback.base$enback_loadcon=="2-Back",]
  
  sum.base$mrt.0.tar[r]<-mean(tmp0[tmp0$S=="tar",]$RT,na.rm=TRUE)
  sum.base$mrt.0.lur[r]<-mean(tmp0[tmp0$S=="lur",]$RT,na.rm=TRUE)
  sum.base$mrt.0.non[r]<-mean(tmp0[tmp0$S=="non",]$RT,na.rm=TRUE)
  sum.base$sdrt.0.tar[r]<-sd(tmp0[tmp0$S=="tar",]$RT,na.rm=TRUE)
  sum.base$sdrt.0.lur[r]<-sd(tmp0[tmp0$S=="lur",]$RT,na.rm=TRUE)
  sum.base$sdrt.0.non[r]<-sd(tmp0[tmp0$S=="non",]$RT,na.rm=TRUE)
  sum.base$acc.0.tar[r]<-mean(tmp0[tmp0$S=="tar",]$acc,na.rm=TRUE)
  sum.base$acc.0.lur[r]<-mean(tmp0[tmp0$S=="lur",]$acc,na.rm=TRUE)
  sum.base$acc.0.non[r]<-mean(tmp0[tmp0$S=="non",]$acc,na.rm=TRUE)
  
  sum.base$acc.0.total[r]<-mean(tmp0$acc,na.rm=TRUE)
  sum.base$p.omit.0[r]<-mean(is.na(tmp0$acc))
  
  sum.base$mrt.2.tar[r]<-mean(tmp2[tmp2$S=="tar",]$RT,na.rm=TRUE)
  sum.base$mrt.2.lur[r]<-mean(tmp2[tmp2$S=="lur",]$RT,na.rm=TRUE)
  sum.base$mrt.2.non[r]<-mean(tmp2[tmp2$S=="non",]$RT,na.rm=TRUE)
  sum.base$sdrt.2.tar[r]<-sd(tmp2[tmp2$S=="tar",]$RT,na.rm=TRUE)
  sum.base$sdrt.2.lur[r]<-sd(tmp2[tmp2$S=="lur",]$RT,na.rm=TRUE)
  sum.base$sdrt.2.non[r]<-sd(tmp2[tmp2$S=="non",]$RT,na.rm=TRUE)
  sum.base$acc.2.tar[r]<-mean(tmp2[tmp2$S=="tar",]$acc,na.rm=TRUE)
  sum.base$acc.2.lur[r]<-mean(tmp2[tmp2$S=="lur",]$acc,na.rm=TRUE)
  sum.base$acc.2.non[r]<-mean(tmp2[tmp2$S=="non",]$acc,na.rm=TRUE)
  
  sum.base$acc.2.total[r]<-mean(tmp2$acc,na.rm=TRUE)
  sum.base$p.omit.2[r]<-mean(is.na(tmp2$acc))
  
}

save(sum.base,file="sum_base_fits.RData")
load("sum_base_fits.RData")

# inclusion for each task

inc.0b<-sum.base[sum.base$acc.0.total>=.55 & sum.base$p.omit.0<=.25,"s"]
inc.2b<-sum.base[sum.base$acc.2.total>=.55 & sum.base$p.omit.2<=.25,"s"]

# final groups

base.0.dmc<-nback.base[nback.base$enback_loadcon=="0-Back" & nback.base$subject%in%inc.0b,]
base.0.dmc<-base.0.dmc[,c("subject","S","R","RT")]
colnames(base.0.dmc)<-c("s","S","R","RT")
base.0.dmc$s<-factor(base.0.dmc$s)

base.2.dmc<-nback.base[nback.base$enback_loadcon=="2-Back" & nback.base$subject%in%inc.2b,]
base.2.dmc<-base.2.dmc[,c("subject","S","R","RT")]
colnames(base.2.dmc)<-c("s","S","R","RT")
base.2.dmc$s<-factor(base.2.dmc$s)

# identify RT bin at a which responding rises above chance in baseline

base0.RTs<-base.0.dmc[!is.na(base.0.dmc$RT),]
base0.RTs$acc<-0
base0.RTs[base0.RTs$S=="tar" & base0.RTs$R=="YES",]$acc<-1
base0.RTs[base0.RTs$S=="lur" & base0.RTs$R=="NO",]$acc<-1
base0.RTs[base0.RTs$S=="non" & base0.RTs$R=="NO",]$acc<-1

base2.RTs<-base.2.dmc[!is.na(base.2.dmc$RT),]
base2.RTs$acc<-0
base2.RTs[base2.RTs$S=="tar" & base2.RTs$R=="YES",]$acc<-1
base2.RTs[base2.RTs$S=="lur" & base2.RTs$R=="NO",]$acc<-1
base2.RTs[base2.RTs$S=="non" & base2.RTs$R=="NO",]$acc<-1

acc.by.bin<-data.frame(seq(0,.480,.02),seq(.02,.500,.02))
colnames(acc.by.bin)<-c("lo","hi")
acc.by.bin$n.0<-NA
acc.by.bin$acc.0<-NA
acc.by.bin$n.2<-NA
acc.by.bin$acc.2<-NA

for (r in 1:length(acc.by.bin$lo)){
  RT.acc.0<-base0.RTs[base0.RTs$RT>=acc.by.bin[r,"lo"] & 
                        base0.RTs$RT<acc.by.bin[r,"hi"]  ,]$acc
  RT.acc.2<-base0.RTs[base2.RTs$RT>=acc.by.bin[r,"lo"] & 
                        base2.RTs$RT<acc.by.bin[r,"hi"]  ,]$acc
  acc.by.bin[r,"n.0"]<-length(RT.acc.0)
  acc.by.bin[r,"acc.0"]<-mean(RT.acc.0)
  acc.by.bin[r,"n.2"]<-length(RT.acc.2)
  acc.by.bin[r,"acc.2"]<-mean(RT.acc.2)
}

plot(acc.by.bin$lo,acc.by.bin$acc.0)
plot(acc.by.bin$lo,acc.by.bin$n.0)
plot(acc.by.bin$lo,acc.by.bin$acc.2)
plot(acc.by.bin$lo,acc.by.bin$n.2)


# exclude fast guesses

base.0.dmc<-base.0.dmc[base.0.dmc$RT>=0.200 | is.na(base.0.dmc$RT),]
base.2.dmc<-base.2.dmc[base.2.dmc$RT>=0.200 | is.na(base.2.dmc$RT),]

# save out

save(base.0.dmc,file="base_0_dat_DDM.RData")
save(base.2.dmc,file="base_2_dat_DDM.RData")

################################################################
####### format stop-signal and apply exclusion criteria ########
################################################################

# load trial-level data merged across participants
# (once downloaded from aws and merged) 

abcd.sst.raw<-read.csv("sst_revised.csv",stringsAsFactors = FALSE)

#### format relevant columns for checks and model fitting ####

# make DMC-friendly columns

# subject
abcd.sst.raw$s<-as.factor(abcd.sst.raw$subject) 

# stimulus
abcd.sst.raw$S<-NA
abcd.sst.raw[abcd.sst.raw$sst_stim=="left_arrow",]$S<-"left"
abcd.sst.raw[abcd.sst.raw$sst_stim=="right_arrow",]$S<-"right"
abcd.sst.raw$S<-as.factor(abcd.sst.raw$S)

# response
abcd.sst.raw$R<-NA
abcd.sst.raw[abcd.sst.raw$sst_primaryresp=="left_arrow" & 
               !is.na(abcd.sst.raw$sst_primaryresp),]$R<-"LEFT"
abcd.sst.raw[abcd.sst.raw$sst_primaryresp=="right_arrow" & 
               !is.na(abcd.sst.raw$sst_primaryresp),]$R<-"RIGHT"
abcd.sst.raw[abcd.sst.raw$sst_primaryresp=="" & 
               !is.na(abcd.sst.raw$sst_primaryresp),]$R<-NA
abcd.sst.raw$R<-as.factor(abcd.sst.raw$R)

# stop vs go condition
abcd.sst.raw$SS<-NA
abcd.sst.raw[is.na(abcd.sst.raw$sst_ssd_dur),]$SS<-"GO"
abcd.sst.raw[!is.na(abcd.sst.raw$sst_ssd_dur),]$SS<-"SS"
abcd.sst.raw$SS<-as.factor(abcd.sst.raw$SS)

# SSD (convert to seconds)
abcd.sst.raw$SSD<-abcd.sst.raw$sst_ssd_dur/1000
abcd.sst.raw[is.na(abcd.sst.raw$SSD),]$SSD<-Inf

# RT (convert to seconds)
abcd.sst.raw$RT<-abcd.sst.raw$sst_primaryrt/1000

# Exclude "stop" trials for DDM fits

abcd.sst.raw<-abcd.sst.raw[abcd.sst.raw$SS=="GO",]

#### checks and exclusion criteria #######

# specify summary data frame variables

# subject number and wave
base.subs<-unique(abcd.sst.raw[abcd.sst.raw$eventname=="baseline_year_1_arm_1",]$s)
y2.subs<-unique(abcd.sst.raw[abcd.sst.raw$eventname=="2_year_follow_up_y_arm_1",]$s)

sum.initial<-data.frame(base.subs,"base")
colnames(sum.initial)<-c("s","wave")
sum.tmp<-data.frame(y2.subs,"y2")
colnames(sum.tmp)<-c("s","wave")
sum.initial<-rbind(sum.initial,sum.tmp);rm(sum.tmp)

sum.initial$go.acc<-NA # go trial choice accuracy
sum.initial$p.go.omit<-NA # prob. of omission (non-response) on "go" trials
sum.initial$go.mrt<-NA # mean RT on go trials (correct AND incorrect choices)
sum.initial$go.sdrt<-NA # SD of RT on go trials (correct AND incorrect choices)

# calculate all summary variables for each person
for (s in unique(abcd.sst.raw$s)){
  for (w in unique(abcd.sst.raw$eventname)){
    if (length(abcd.sst.raw[abcd.sst.raw$s==s & abcd.sst.raw$eventname==w,]$s)>0){
      
      tmp<-abcd.sst.raw[abcd.sst.raw$s==s & abcd.sst.raw$eventname==w,]
      wave<-"base";if(w=="2_year_follow_up_y_arm_1"){wave<-"y2"}
      r<-(sum.initial$s==s & sum.initial$wave==wave)
      
      sum.initial[r,]$go.acc<-length(tmp[!is.na(tmp$R) & tmp$S==tolower(tmp$R) & tmp$SS=="GO",]$s)/length(tmp[!is.na(tmp$R) & tmp$SS=="GO",]$s)
      sum.initial[r,]$p.go.omit<-length(tmp[is.na(tmp$RT) & tmp$SS=="GO",]$s)/length(tmp[tmp$SS=="GO",]$s)
      sum.initial[r,]$go.mrt<-mean(tmp[tmp$SS=="GO",]$RT,na.rm = TRUE)
      sum.initial[r,]$go.sdrt<-sd(tmp[tmp$SS=="GO",]$RT,na.rm = TRUE)
    }
  }}

# basic inclusion criteria for DDM

inc.sst.base<-sum.initial[sum.initial$go.acc>=.55 & 
                            sum.initial$p.go.omit<=.25 &
                            sum.initial$wave=="base","s"]

inc.sst.y2<-sum.initial[sum.initial$go.acc>=.55 & 
                          sum.initial$p.go.omit<=.25 &
                          sum.initial$wave=="y2","s"]

# break up into final base and y2 groups

base.sst.dmc<-abcd.sst.raw[abcd.sst.raw$eventname=="baseline_year_1_arm_1"  & 
                             abcd.sst.raw$subject%in%inc.sst.base,]
base.sst.dmc<-base.sst.dmc[,c("subject","S","R","RT")]
colnames(base.sst.dmc)<-c("s","S","R","RT")
base.sst.dmc$s<-factor(base.sst.dmc$s)

y2.sst.dmc<-abcd.sst.raw[abcd.sst.raw$eventname=="2_year_follow_up_y_arm_1"  & 
                           abcd.sst.raw$subject%in%inc.sst.y2,]
y2.sst.dmc<-y2.sst.dmc[,c("subject","S","R","RT")]
colnames(y2.sst.dmc)<-c("s","S","R","RT")
y2.sst.dmc$s<-factor(y2.sst.dmc$s)

# identify RT bin at a which responding rises above chance in baseline

base.RTs<-base.sst.dmc[!is.na(base.sst.dmc$RT),]
base.RTs$acc<-(base.RTs$S==tolower(base.RTs$R))

acc.by.bin<-data.frame(seq(0,.480,.02),seq(.02,.500,.02))
colnames(acc.by.bin)<-c("lo","hi")
acc.by.bin$n<-NA
acc.by.bin$acc<-NA


for (r in 1:length(acc.by.bin$lo)){
  RT.acc<-base.RTs[base.RTs$RT>=acc.by.bin[r,"lo"] & 
                     base.RTs$RT<acc.by.bin[r,"hi"]  ,]$acc
  acc.by.bin[r,"n"]<-length(RT.acc)
  acc.by.bin[r,"acc"]<-mean(RT.acc)
}

plot(acc.by.bin$lo,acc.by.bin$acc)

# exclude fast guesses

base.sst.dmc<-base.sst.dmc[base.sst.dmc$RT>=0.200 | is.na(base.sst.dmc$RT),]
y2.sst.dmc<-y2.sst.dmc[y2.sst.dmc$RT>=0.200 | is.na(y2.sst.dmc$RT),]

# treat RTs outside of the shortest possible response 
# window (>1.700s) as omissions for the sake of the omission model

base.sst.dmc[base.sst.dmc$RT>1.700 & !is.na(base.sst.dmc$RT),]$RT<-NA

y2.sst.dmc[y2.sst.dmc$RT>1.700 & !is.na(y2.sst.dmc$RT),]$RT<-NA

# save out
save(base.sst.dmc,file="base_sst_dat_DDM200.RData")
save(y2.sst.dmc,file="y2_sst_dat_DDM200.RData")


######################################################
####### select hierarchical subsamples (for priors) ##
######################################################

# exclude data from individuals with accuracy <50% and
# omissions > 25% on either task

sum.base[sum.base$acc.0.overall<0.55 | 
           sum.base$acc.0.omit>0.25 | is.na(sum.base$acc.0.overall),
         colnames(sum.base)[grepl(".0.",colnames(sum.base))] ]<-NA

sum.base[sum.base$acc.2.overall<0.55 | 
           sum.base$acc.2.omit>0.25 | is.na(sum.base$acc.2.overall),
         colnames(sum.base)[grepl(".2.",colnames(sum.base))] ]<-NA

sum.base.sst<-sum.initial[sum.initial$go.acc>=.55 & 
                            sum.initial$p.go.omit<=.25 &
                            sum.initial$wave=="base",]

# identify subjects with acceptable data from both n-back tasks and the sst

sum.elig<-sum.base[(is.na(sum.base$acc.0.omit)+is.na(sum.base$acc.2.omit))<1,]
sum.elig<-sum.elig[sum.elig$s%in%sum.base.sst$s,]

# load eligible participant list

ADHD_elig<-read.csv("hier_DDM_eligible_ADHD_proj.csv")

# only look at eligible individuals 

sum.elig<-sum.elig[sum.elig$s%in%ADHD_elig$x,]

# select subsample IDs at random

sub.ids<-sample(sum.elig$s,size = 300,replace = FALSE)

write.csv(sub.ids,
          file="subsample_ADHD_project.csv",
          row.names = FALSE)

#####################################################################
#####################################################################
############## N-back task DDM ######################################
#####################################################################
#####################################################################

#####################################################################
###### individual-level estimation: broad priors ####################
#####################################################################

### load data ###

load("base_0_dat_DDM.RData")

load("base_2_dat_DDM.RData")

#### Design and contaminant (gf) ----
load_model ("DDM","ddm_omit.R")

model <- model.dmc(
  p.map = list(a="1",v=c("S"),z="1",d="1",sz="1",
               sv="1",t0="1",st0="1",
               censor="1",gf="1"), 
  responses = c("NO","YES"),
  match.map = list(M=list(tar="YES",lur="NO",non="NO")),
  factors=list(S=c("lur","non","tar")),
  constants = c(sv=0,sz=0,d=0,censor=2.00), 
  type="rd")

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "a"     "v.lur" "v.non" "v.tar" "z"     "t0"   
# [7] "st0"   "gf"   
# 
# Constants are (see attr(,"constants") ):
#   sz     sv      d censor 
# 0      0      0      2 


p1 <- c(a=1,v.lur=3, v.non = 3, v.tar= 3,
        z=0.5,t0=0.3, st0=0.1,gf=0)

p.prior <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p1,                           
  p2=c(.5,1,1,1,.1,.1,.05,1),
  lower=c(0,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior)) plot.prior(i,p.prior)


zero.base.mdi <- data.model.dmc(base.0.dmc, model)

two.base.mdi <- data.model.dmc(base.2.dmc, model)

sZero.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = zero.base.mdi, thin = 5)
sZero.base.first<-sZero.base[1:4500]
save(sZero.base.first,file="zero_ddm_broad_first.RData")
sZero.base.second<-sZero.base[4501:length(sZero.base)]
save(sZero.base.second,file="zero_ddm_broad_second.RData")

sTwo.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = two.base.mdi, thin = 5)
sTwo.base.first<-sTwo.base[1:4500]
save(sTwo.base.first,file="two_ddm_broad_first.RData")
sTwo.base.second<-sTwo.base[4501:length(sTwo.base)]
save(sTwo.base.second,file="two_ddm_broad_second.RData")

# sampling can be conducted efficiently by using a high-performance 
# computer to estimate parameters for each participant in parallel
# using the RUN.dmc() function (with default parameters) as described
# in the Dynamic Models of Choice tutorial: https://osf.io/pbwx8/


####################################################################
### hierarchical model fits in independent subsample (for priors) ###
#####################################################################

### load data ###

load("base_0_dat_DDM.RData")

load("base_2_dat_DDM.RData")

# independent subsample IDs only
subsample<-read.csv("subsample_ADHD_project.csv")

base.0.dmc<-base.0.dmc[base.0.dmc$s%in%subsample$x,]
base.0.dmc$s<-as.factor(as.character(base.0.dmc$s))

base.2.dmc<-base.2.dmc[base.2.dmc$s%in%subsample$x,]
base.2.dmc$s<-as.factor(as.character(base.2.dmc$s))

#### Design and contaminant (gf) ----
load_model ("DDM","ddm_omit.R")

model <- model.dmc(
  p.map = list(a="1",v=c("S"),z="1",d="1",sz="1",
               sv="1",t0="1",st0="1",
               censor="1",gf="1"), 
  responses = c("NO","YES"),
  match.map = list(M=list(tar="YES",lur="NO",non="NO")),
  factors=list(S=c("lur","non","tar")),
  constants = c(sz = 0, sv = 0, d=0,censor=2.00), 
  type="rd")

p1 <- c(a=1,v.lur=3, v.non = 3, v.tar= 3,
        z=0.5,t0=0.3, st0=.1, gf=0)

p.prior <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p1,                           
  p2=c(.5,1,1,1,.1,.1,.05,1),
  lower=c(0,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

# same prior for mu
p.prior.mu<-p.prior 

# gamma priors for sigma
p1[1:length(p1)] <- rep(1,length(p1))
p.prior.sigma<-prior.p.dmc(
  dists=rep("gamma",length(p1)),
  p1=p1,p2=p1)
# par(mfcol=c(1,3)); for (i in names(p.prior.sigma)) plot.prior(i,p.prior.sigma)

# combine priors for hyperparameters
pp.prior<-list(p.prior.mu,p.prior.sigma)


# make mdi

zero.hier.mdi <- data.model.dmc(base.0.dmc, model)

two.hier.mdi <- data.model.dmc(base.2.dmc, model)


#start points
load("nback_ddm_broad_samples.RData")

hstart.zero <- make.hstart(sZero.base[subsample$x])
theta1.zero <- make.theta1(sZero.base[subsample$x])

hstart.two <- make.hstart(sTwo.base[subsample$x])
theta1.two <- make.theta1(sTwo.base[subsample$x])


# make sampling objects
sZero.hier <- h.samples.dmc(nmc = 40, p.prior=p.prior, 
                            data=zero.hier.mdi, 
                            pp.prior=pp.prior,
                            thin=10,
                            hstart.prior=hstart.zero,
                            theta1=theta1.zero)

save(sZero.hier ,file="zero_hier_ADHD_project.RData")


sTwo.hier <- h.samples.dmc(nmc = 40, p.prior=p.prior, 
                           data=two.hier.mdi, 
                           pp.prior=pp.prior,
                           thin=10,
                           hstart.prior=hstart.two,
                           theta1=theta1.two)

save(sTwo.hier ,file="two_hier_ADHD_project.RData")



# run the following code on a high-performance computer 
# to efficiently generate samples with parallel processing:

rm(list=ls())
source ("dmc/dmc.R")
load_model ("DDM","ddm_omit.R") 

load("zero_hier_ADHD_project.RData")

cores=32

sZero.hier <- h.run.unstuck.dmc(sZero.hier,cores=cores,report=10,p.migrate=0.025,h.p.migrate=0.025)

save.image("zero_hier_ADHD_project.RData")

sZero.hier1 <- h.run.converge.dmc(samples=h.samples.dmc(samples=sZero.hier,nmc=40,thin=10),
                                  nmc=40,cores=cores,report=10,verbose=TRUE)

save.image("zero_hier_ADHD_project.RData")

sZero.hier2 <- h.run.converge.dmc(samples=h.samples.dmc(samples=sZero.hier1,nmc=120,thin=10),
                                  nmc=40,cores=cores,report=10,verbose=TRUE)

save.image("zero_hier_ADHD_project.RData")

rm(list=ls())
source ("dmc/dmc.R")
load_model ("DDM","ddm_omit.R") 

load("two_hier_ADHD_project.RData")

cores=32

sTwo.hier <- h.run.unstuck.dmc(sTwo.hier,cores=cores,report=10,p.migrate=0.025,h.p.migrate=0.025)

save.image("two_hier_ADHD_project.RData")

sTwo.hier1 <- h.run.converge.dmc(samples=h.samples.dmc(samples=sTwo.hier,nmc=40,thin=10),
                                 nmc=40,cores=cores,report=10,verbose=TRUE)

save.image("two_hier_ADHD_project.RData")

sTwo.hier2 <- h.run.converge.dmc(samples=h.samples.dmc(samples=sTwo.hier1,nmc=120,thin=10),
                                 nmc=40,cores=cores,report=10,verbose=TRUE)

save.image("two_hier_ADHD_project.RData")


###### results ######

# read in data

load("zero_hier_ADHD_project.RData")
load("two_hier_ADHD_project.RData")

plot.dmc(sZero.hier2,hyper=TRUE,layout=c(2,2))
gelman.diag.dmc(sZero.hier2,hyper=TRUE)

plot.dmc(sTwo.hier2,hyper=TRUE,layout=c(2,2))
gelman.diag.dmc(sTwo.hier2,hyper=TRUE)

#### model fit

sim.Zero.hier <- h.post.predict.dmc(sZero.hier2,cores=3)
sim.Two.hier <- h.post.predict.dmc(sTwo.hier2,cores=3)

save(sim.Zero.hier,sim.Two.hier,file="nback_hier_fmri_project_pp.RData")

plot.pp.dmc(sim.Zero.hier,layout = c(1,3))

plot.pp.dmc(sim.Two.hier,layout = c(1,3))


##### make informed priors

library(MASS)

# combine all individuals' samples into a single vector for each parameter

zero.samps<-vector("list", length = length(attr(model,"p.vector"))) 
names(zero.samps)<-names(attr(model,"p.vector"))

for (s in 1:length(sZero.hier2)){
  tmp<-sZero.hier2[[s]]
  for (par in names(attr(model,"p.vector"))){
    zero.samps[[par]]<-c(zero.samps[[par]],as.vector(tmp$theta[,par,]))
  }  
}


two.samps<-vector("list", length = length(attr(model,"p.vector"))) 
names(two.samps)<-names(attr(model,"p.vector"))

for (s in 1:length(sTwo.hier2)){
  tmp<-sTwo.hier2[[s]]
  for (par in names(attr(model,"p.vector"))){
    two.samps[[par]]<-c(two.samps[[par]],as.vector(tmp$theta[,par,]))
  }  
}


# upper and lower bounds

lower<-c(0,-Inf,-Inf,-Inf,0,0,0,-Inf)
upper<-c(Inf,Inf,Inf,Inf,1,2,2,Inf) 
names(lower)<-names(p.prior)
names(upper)<-names(p.prior)


h1.gauss.fits.zero<-data.frame(names(attr(model,"p.vector")),NA,NA)
h1.gauss.fits.two<-data.frame(names(attr(model,"p.vector")),NA,NA)
colnames(h1.gauss.fits.zero)<-c("par","m","sd")
colnames(h1.gauss.fits.two)<-c("par","m","sd")

for (par in names(attr(model,"p.vector"))){
  dtnorm0 <- function(X, mean, sd, log = FALSE,low,up) {
    dtnorm(X, mean, sd, low, up,log)}
  
  G.tmp<-fitdistr(zero.samps[[par]], dtnorm0, 
                  start=list(mean=mean(zero.samps[[par]]), 
                             sd=sd(zero.samps[[par]])),
                  low=lower[par],up=upper[par])
  h1.gauss.fits.zero[h1.gauss.fits.zero$par==par,2:3]<-G.tmp$estimate
  
  G.tmp<-fitdistr(two.samps[[par]], dtnorm0, 
                  start=list(mean=mean(two.samps[[par]]), 
                             sd=sd(two.samps[[par]])),
                  low=lower[par],up=upper[par])
  h1.gauss.fits.two[h1.gauss.fits.two$par==par,2:3]<-G.tmp$estimate
}

h1.gauss.fits.zero$m<-as.numeric(h1.gauss.fits.zero$m)
h1.gauss.fits.zero$sd<-as.numeric(h1.gauss.fits.zero$sd)
h1.gauss.fits.two$m<-as.numeric(h1.gauss.fits.two$m)
h1.gauss.fits.two$sd<-as.numeric(h1.gauss.fits.two$sd)


# write out
write.csv(h1.gauss.fits.zero,file = "zero_informed_priors_ADHD_project.csv",row.names = FALSE)
write.csv(h1.gauss.fits.two,file = "two_informed_priors_ADHD_project.csv",row.names = FALSE)


#plot fits

jpeg(filename = "zero_fits_dists_ADHD_project.jpeg",
     units = "in", res = 300,
     width = 12,height = 9)
par(mfrow=c(3,4))
for (par in names(attr(model,"p.vector"))){
  msd<-as.numeric(h1.gauss.fits.zero[h1.gauss.fits.zero$par==par,2:3])
  hist(zero.samps[[par]], prob = TRUE,main = par,xlab= par)
  curve(dtnorm0(x, msd[1], msd[2],low=lower[par],up=upper[par]), 
        col = "red", add = TRUE)
}
dev.off()

jpeg(filename = "two_fits_dists_ADHD_project.jpeg",
     units = "in", res = 300,
     width = 12,height = 9)
par(mfrow=c(3,4))
for (par in names(attr(model,"p.vector"))){
  msd<-as.numeric(h1.gauss.fits.two[h1.gauss.fits.two$par==par,2:3])
  hist(two.samps[[par]], prob = TRUE,main = par,xlab= par)
  curve(dtnorm0(x, msd[1], msd[2],low=lower[par],up=upper[par]), 
        col = "red", add = TRUE)
}
dev.off()


#####################################################################
###### individual-level estimation: informed priors #################
#####################################################################

### load data ###

load("base_0_dat_DDM.RData")

load("base_2_dat_DDM.RData")

subsample<-read.csv("nback_subsample_ADHD_project.csv")

base.0.dmc<-base.0.dmc[!base.0.dmc$s%in%subsample$x,]
base.0.dmc$s<-as.factor(as.character(base.0.dmc$s))

base.2.dmc<-base.2.dmc[!base.2.dmc$s%in%subsample$x,]
base.2.dmc$s<-as.factor(as.character(base.2.dmc$s))

#### Design and contaminant (gf) ----
load_model ("DDM","ddm_omit.R")

model <- model.dmc(
  p.map = list(a="1",v=c("S"),z="1",d="1",sz="1",
               sv="1",t0="1",st0="1",
               censor="1",gf="1"), 
  responses = c("NO","YES"),
  match.map = list(M=list(tar="YES",lur="NO",non="NO")),
  factors=list(S=c("lur","non","tar")),
  constants = c(sz = 0, sv = 0, d=0,censor=2.00), 
  type="rd")

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "a"     "v.lur" "v.non" "v.tar" "z"     "t0"   
# [7] "st0"   "gf"   
# 
# Constants are (see attr(,"constants") ):
#   sz     sv      d censor 
# 0      0      0      2 


hp0<-read.csv("zero_informed_priors_ADHD_project.csv")
rownames(hp0)<-hp0$par

p1.zero <- c(a=hp0["a","m"],v.lur=hp0["v.lur","m"],v.non=hp0["v.non","m"],
             v.tar=hp0["v.tar","m"],z=hp0["z","m"],t0=hp0["t0","m"],
             st0=hp0["st0","m"],gf=hp0["gf","m"])
p2.zero <- c(a=hp0["a","sd"],v.lur=hp0["v.lur","sd"],v.non=hp0["v.non","sd"],
             v.tar=hp0["v.tar","sd"],z=hp0["z","sd"],t0=hp0["t0","sd"],
             st0=hp0["st0","sd"],gf=hp0["gf","sd"])

p.prior.zero <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p1.zero,                           
  p2=p2.zero,
  lower=c(0,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior.zero)) plot.prior(i,p.prior.zero)


hp2<-read.csv("two_informed_priors_ADHD_project.csv")
rownames(hp2)<-hp2$par

p1.two <- c(a=hp2["a","m"],v.lur=hp2["v.lur","m"],v.non=hp2["v.non","m"],
            v.tar=hp2["v.tar","m"],z=hp2["z","m"],t0=hp2["t0","m"],
            st0=hp2["st0","m"],gf=hp2["gf","m"])
p2.two <- c(a=hp2["a","sd"],v.lur=hp2["v.lur","sd"],v.non=hp2["v.non","sd"],
            v.tar=hp2["v.tar","sd"],z=hp2["z","sd"],t0=hp2["t0","sd"],
            st0=hp2["st0","sd"],gf=hp2["gf","sd"])

p.prior.two <- prior.p.dmc(
  dists = rep("tnorm",8),
  p1=p1.two,                           
  p2=p2.two,
  lower=c(0,NA,NA,NA,0,0,0,NA),upper=c(NA,NA,NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior.two)) plot.prior(i,p.prior.two)



zero.base.mdi <- data.model.dmc(base.0.dmc, model)

two.base.mdi <- data.model.dmc(base.2.dmc, model)

sZero.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior.zero,data = zero.base.mdi, thin = 5)
sTwo.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior.two,data = two.base.mdi, thin = 5)

# save(sZero.base,file="zero_DDM_inf_issues.RData")
# save(sTwo.base,file="two_DDM_inf_issues.RData")


sZero.base.first<-sZero.base[1:4500]
save(sZero.base.first,file="zero_ddm_inf_ADHD_first.RData")
sZero.base.second<-sZero.base[4501:length(sZero.base)]
save(sZero.base.second,file="zero_ddm_inf_ADHD_second.RData")

sTwo.base.first<-sTwo.base[1:4500]
save(sTwo.base.first,file="two_ddm_inf_ADHD_first.RData")
sTwo.base.second<-sTwo.base[4501:length(sTwo.base)]
save(sTwo.base.second,file="two_ddm_inf_ADHD_second.RData")

# sampling can be conducted efficiently by using a high-performance 
# computer to estimate parameters for each participant in parallel
# using the RUN.dmc() function (with default parameters) as described
# in the Dynamic Models of Choice tutorial: https://osf.io/pbwx8/

##############################################################
##### generate point estimates (posterior medians)############
##############################################################

#### load in parts, merge them, save out

load("zero_ddm_inf_ADHD_first.RData")
load("zero_ddm_inf_ADHD_second.RData")

sZero.base<-sZero.base.first
sZero.base[4501:(4500+length(sZero.base.second))]<-sZero.base.second
names(sZero.base)<-c(names(sZero.base.first),names(sZero.base.second))

load("two_ddm_inf_ADHD_first.RData")
load("two_ddm_inf_ADHD_second.RData")

sTwo.base<-sTwo.base.first
sTwo.base[4501:(4500+length(sTwo.base.second))]<-sTwo.base.second
names(sTwo.base)<-c(names(sTwo.base.first),names(sTwo.base.second))


save(sZero.base,sTwo.base,file="nback_ddm_inf_ADHD_samples.RData")

#### check convergence, remove non-converged if needed

# load("nback_ddm_inf_ADHD_samples.RData")

rhat.Zero<-gelman.diag.dmc(sZero.base, hyper = FALSE)
rhat.Two<-gelman.diag.dmc(sTwo.base, hyper = FALSE)
save(rhat.Zero,rhat.Two,file="nback_ddm_inf_ABCD_rhat.RData")

#### model fit

# select three sites at random, generate predictions

y1.dat<-read.csv("year1_full_dat.csv")
sample(unique(y1.dat$site_id_l),3,FALSE)
# [1] "site17" "site01" "site05"

abcd_lt01<-read.abcd("abcd_lt01.txt")

s17.ids<-abcd_lt01[abcd_lt01$site_id_l=="site17" & abcd_lt01$eventname=="baseline_year_1_arm_1",]$subjectkey
s01.ids<-abcd_lt01[abcd_lt01$site_id_l=="site01" & abcd_lt01$eventname=="baseline_year_1_arm_1",]$subjectkey
s05.ids<-abcd_lt01[abcd_lt01$site_id_l=="site05" & abcd_lt01$eventname=="baseline_year_1_arm_1",]$subjectkey

save(s17.ids,s01.ids,s05.ids,file="fit_IDs.RData")

load("nback_ddm_inf_ADHD_samples.RData")
load("fit_IDs.RData")

sim.Zero.s17 <- h.post.predict.dmc(sZero.base[names(sZero.base)%in%s17.ids],cores=3)
sim.Zero.s01 <- h.post.predict.dmc(sZero.base[names(sZero.base)%in%s01.ids],cores=3)
sim.Zero.s05 <- h.post.predict.dmc(sZero.base[names(sZero.base)%in%s05.ids],cores=3)

sim.Two.s17 <- h.post.predict.dmc(sTwo.base[names(sTwo.base)%in%s17.ids],cores=3)
sim.Two.s01 <- h.post.predict.dmc(sTwo.base[names(sTwo.base)%in%s01.ids],cores=3)
sim.Two.s05 <- h.post.predict.dmc(sTwo.base[names(sTwo.base)%in%s05.ids],cores=3)

save(sim.Zero.s17,sim.Two.s17,
     sim.Zero.s01,sim.Two.s01,
     sim.Zero.s05,sim.Two.s05,file="nback_ddm_3sites_pp.RData")

plot.pp.dmc(sim.Zero.base,layout = c(1,3))

plot.pp.dmc(sim.Two.base,layout = c(1,3))

#### parameter estimates

# loop to make summary parameters

sZero.base.pars<-sZero.base

for (s in 1:length(sZero.base)){
  t<-sZero.base[[s]]$theta; d<-dim(t)
  Ter<-t[,"t0",]+(t[,"st0",]/2)
  gf.natural<-pnorm(t[,"gf",])
  D <- d + c(0, 2, 0)
  t2<-array(sapply(1:D[3], function(x) cbind(t[,,x], Ter[,x],
                                             gf.natural[,x])), D)
  dimnames(t2)[[2]]<-c(dimnames(t)[[2]],"Ter","gf.natural")
  sZero.base.pars[[s]]$theta<-t2
}

sTwo.base.pars<-sTwo.base

for (s in 1:length(sTwo.base)){
  t<-sTwo.base[[s]]$theta; d<-dim(t)
  Ter<-t[,"t0",]+(t[,"st0",]/2)
  gf.natural<-pnorm(t[,"gf",])
  D <- d + c(0, 2, 0)
  t2<-array(sapply(1:D[3], function(x) cbind(t[,,x], Ter[,x],
                                             gf.natural[,x])), D)
  dimnames(t2)[[2]]<-c(dimnames(t)[[2]],"Ter","gf.natural")
  sTwo.base.pars[[s]]$theta<-t2
}

# save out point estimates 

sZero.base_medians<-lapply(sZero.base.pars,FUN=function(x) apply(x$theta,2,median))
sZero.base_medians<-as.data.frame(t(as.data.frame(sZero.base_medians)))

sTwo.base_medians<-lapply(sTwo.base.pars,FUN=function(x) apply(x$theta,2,median))
sTwo.base_medians<-as.data.frame(t(as.data.frame(sTwo.base_medians)))

write.csv(sZero.base_medians,file="zero_ddm_postmedians_ADHD.csv")

write.csv(sTwo.base_medians,file="two_ddm_postmedians_ADHD.csv")

# merge into one file

tmp0<-read.csv("zero_ddm_postmedians_fmri.csv")
colnames(tmp0)<-c("s",paste0("zero.",colnames(tmp0)[-1]))
tmp2<-read.csv("two_ddm_postmedians_fmri.csv")
colnames(tmp2)<-c("s",paste0("two.",colnames(tmp2)[-1]))

sum.ddm<-merge(tmp0, tmp2,
               by="s",all = TRUE)

write.csv(sum.ddm,file="nback_postmedians_ADHD.csv",row.names = FALSE)


#####################################################################
#####################################################################
############## stop-signal task DDM #################################
#####################################################################
#####################################################################

#####################################################################
###### individual-level estimation: broad priors ####################
#####################################################################

### load data ###

load("base_sst_dat_DDM200.RData")

#### Standard DDM ----
load_model ("DDM","ddm_omit.R")

model <- model.dmc(
  p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1",
                   censor="1",gf="1"),
  match.map = list(M=list(left="LEFT",right="RIGHT")),
  factors=list(S=c("left","right")),
  constants = c(sz = 0, sv = 0, d=0,censor=1.700),
  responses = c("LEFT","RIGHT"),
  type = "rd")


p1 <- c(a=1,v=3, z=0.5,t0=0.3, st0=0.1,gf=0)


p.prior <- prior.p.dmc(
  dists = rep("tnorm",6),
  p1=p1,                           
  p2=c(.5,1,.1,.1,.05,1),
  lower=c(0,NA,0,0,0,NA),upper=c(NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,2)); for (i in names(p.prior)) plot.prior(i,p.prior)


sst.test.mdi <- data.model.dmc(base.sst.dmc, model)

sSST.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = sst.test.mdi, thin = 5)

# break in two and save out

sSST.base.first<-sSST.base[1:4500]
save(sSST.base.first,file="sst_base_ddm2_first.RData")

sSST.base.second<-sSST.base[4501:9342]
save(sSST.base.second,file="sst_base_ddm2_second.RData")


# sampling can be conducted efficiently by using a high-performance 
# computer to estimate parameters for each participant in parallel
# using the RUN.dmc() function (with default parameters) as described
# in the Dynamic Models of Choice tutorial: https://osf.io/pbwx8/


####################################################################
### hierarchical model fits in independent subsample (for priors) ###
#####################################################################


### load data ###

load("base_sst_dat_DDM200.RData")

subsample<-read.csv("subsample_ADHD_project.csv")

base.sst.dmc<-base.sst.dmc[base.sst.dmc$s%in%subsample$x,]
base.sst.dmc$s<-as.factor(as.character(base.sst.dmc$s))

#### Design and contaminant (gf) ----
load_model ("DDM","ddm_omit.R")

model <- model.dmc(
  p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1",
                   censor="1",gf="1"),
  match.map = list(M=list(left="LEFT",right="RIGHT")),
  factors=list(S=c("left","right")),
  constants = c(sz = 0, sv = 0, d=0,censor=1.700),
  responses = c("LEFT","RIGHT"),
  type = "rd")


p1 <- c(a=1,v=3, z=0.5,t0=0.3, st0=0.1,gf=0)


p.prior <- prior.p.dmc(
  dists = rep("tnorm",6),
  p1=p1,                           
  p2=c(.5,1,.1,.1,.05,1),
  lower=c(0,NA,0,0,0,NA),upper=c(NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior)) plot.prior(i,p.prior)

# same prior for mu
p.prior.mu<-p.prior 

# gamma priors for sigma
p1[1:length(p1)] <- rep(1,length(p1))
p.prior.sigma<-prior.p.dmc(
  dists=rep("gamma",length(p1)),
  p1=p1,p2=p1)
# par(mfcol=c(1,3)); for (i in names(p.prior.sigma)) plot.prior(i,p.prior.sigma)

# combine priors for hyperparameters
pp.prior<-list(p.prior.mu,p.prior.sigma)


# make mdi

sst.hier.mdi <- data.model.dmc(base.sst.dmc, model)

#start points
load("sst_ddm_broad_samples.RData")

hstart <- make.hstart(sSST.base[subsample$x])
theta1 <- make.theta1(sSST.base[subsample$x])

# make sampling objects
sSST.hier <- h.samples.dmc(nmc = 40, p.prior=p.prior, 
                           data=sst.hier.mdi, 
                           pp.prior=pp.prior,
                           thin=10,
                           hstart.prior=hstart,
                           theta1=theta1)

save(sSST.hier ,file="sst_hier_ABCD_project.RData")

# run the following code on a high-performance computer 
# to efficiently generate samples with parallel processing:

rm(list=ls())
source ("dmc/dmc.R")
load_model ("DDM","ddm_omit.R") 

load("sst_hier_ABCD_project.RData")

cores=32

sSST.hier <- h.run.unstuck.dmc(sSST.hier,cores=cores,report=10,p.migrate=0.025,h.p.migrate=0.025)

save.image("sst_hier_ABCD_project.RData")

sSST.hier1 <- h.run.converge.dmc(samples=h.samples.dmc(samples=sSST.hier,nmc=40,thin=10),
                                 nmc=40,cores=cores,report=10,verbose=TRUE)

save.image("sst_hier_ABCD_project.RData")

sSST.hier2 <- h.run.converge.dmc(samples=h.samples.dmc(samples=sSST.hier1,nmc=120,thin=10),
                                 nmc=40,cores=cores,report=10,verbose=TRUE)

save.image("sst_hier_ABCD_project.RData")

###### results ######

# read in data

load("sst_hier_ABCD_project.RData")

plot.dmc(sSST.hier2,hyper=TRUE,layout=c(2,2))
gelman.diag.dmc(sZero.hier2,hyper=TRUE)


#### model fit

sim.SST.hier <- h.post.predict.dmc(sSST.hier2,cores=3)

save(sim.SST.hier,file="sst_hier_ADHD_project_pp.RData")
#

plot.pp.dmc(sim.SST.hier,layout = c(1,2))

##### make informed priors

library(MASS)

# combine all individuals' samples into a single vector for each parameter

sst.samps<-vector("list", length = length(attr(model,"p.vector"))) 
names(sst.samps)<-names(attr(model,"p.vector"))

for (s in 1:length(sSST.hier2)){
  tmp<-sSST.hier2[[s]]
  for (par in names(attr(model,"p.vector"))){
    sst.samps[[par]]<-c(sst.samps[[par]],as.vector(tmp$theta[,par,]))
  }  
}

# upper and lower bounds

lower<-c(0,-Inf,0,0,0,-Inf)
upper<-c(Inf,Inf,1,2,2,Inf) 
names(lower)<-names(p.prior)
names(upper)<-names(p.prior)


h1.gauss.fits.sst<-data.frame(names(attr(model,"p.vector")),NA,NA)
colnames(h1.gauss.fits.sst)<-c("par","m","sd")

for (par in names(attr(model,"p.vector"))){
  dtnorm0 <- function(X, mean, sd, log = FALSE,low,up) {
    dtnorm(X, mean, sd, low, up,log)}
  
  G.tmp<-fitdistr(sst.samps[[par]], dtnorm0, 
                  start=list(mean=mean(sst.samps[[par]]), 
                             sd=sd(sst.samps[[par]])),
                  low=lower[par],up=upper[par])
  h1.gauss.fits.sst[h1.gauss.fits.sst$par==par,2:3]<-G.tmp$estimate
  
}

h1.gauss.fits.sst$m<-as.numeric(h1.gauss.fits.sst$m)
h1.gauss.fits.sst$sd<-as.numeric(h1.gauss.fits.sst$sd)


# write out
write.csv(h1.gauss.fits.sst,file = "sst_informed_priors_ADHD_project.csv",row.names = FALSE)


#plot fits

jpeg(filename = "sst_fits_dists_ADHD_project.jpeg",
     units = "in", res = 300,
     width = 12,height = 9)
par(mfrow=c(3,4))
for (par in names(attr(model,"p.vector"))){
  msd<-as.numeric(h1.gauss.fits.sst[h1.gauss.fits.sst$par==par,2:3])
  hist(sst.samps[[par]], prob = TRUE,main = par,xlab= par)
  curve(dtnorm0(x, msd[1], msd[2],low=lower[par],up=upper[par]), 
        col = "red", add = TRUE)
}
dev.off()


#####################################################################
###### individual-level estimation: informed priors #################
#####################################################################

### load data ###

load("base_sst_dat_DDM200.RData")

subsample<-read.csv("subsample_ADHD_project.csv")

base.sst.dmc<-base.sst.dmc[!base.sst.dmc$s%in%subsample$x,]
base.sst.dmc$s<-as.factor(as.character(base.sst.dmc$s))

#### Standard DDM ----
load_model ("DDM","ddm_omit.R")

model <- model.dmc(
  p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1",
                   censor="1",gf="1"),
  match.map = list(M=list(left="LEFT",right="RIGHT")),
  factors=list(S=c("left","right")),
  constants = c(sz = 0, sv = 0, d=0,censor=1.700),
  responses = c("LEFT","RIGHT"),
  type = "rd")


hp0<-read.csv("sst_informed_priors_ADHD_project.csv")
rownames(hp0)<-hp0$par

p1 <- c(a=hp0["a","m"],v=hp0["v","m"],
        z=hp0["z","m"],t0=hp0["t0","m"],
        st0=hp0["st0","m"],gf=hp0["gf","m"])
p2 <- c(a=hp0["a","sd"],v=hp0["v","sd"],
        z=hp0["z","sd"],t0=hp0["t0","sd"],
        st0=hp0["st0","sd"],gf=hp0["gf","sd"])

p.prior <- prior.p.dmc(
  dists = rep("tnorm",6),
  p1=p1,                           
  p2=p2,
  lower=c(0,NA,0,0,0,NA),upper=c(NA,NA,1,2,2,NA)
)
# par(mfcol=c(1,3)); for (i in names(p.prior)) plot.prior(i,p.prior)


sst.base.mdi <- data.model.dmc(base.sst.dmc, model)

sSST.base  <- h.samples.dmc(nmc = 120, p.prior=p.prior,data = sst.base.mdi, thin = 5)

# save(sSST.base,file="sst_DDM_inf_ADHD.RData")

sSST.base.first<-sSST.base[1:4500]
save(sSST.base.first,file="sst_DDM_inf_ADHD_first.RData")
sSST.base.second<-sSST.base[4501:length(sSST.base)]
save(sSST.base.second,file="sst_DDM_inf_ADHD_second.RData")

# sampling can be conducted efficiently by using a high-performance 
# computer to estimate parameters for each participant in parallel
# using the RUN.dmc() function (with default parameters) as described
# in the Dynamic Models of Choice tutorial: https://osf.io/pbwx8/

source("greatlakes.R")
run.greatlakes.dmc("sst_DDM_inf_ADHD_first","DDM","ddm_omit.R","asweigar",
                   "asweigar1",40,wall.hours=10)

source("greatlakes.R")
run.greatlakes.dmc("sst_DDM_inf_ADHD_second","DDM","ddm_omit.R","asweigar",
                   "asweigar1",40,wall.hours=10)


##############################################################
##### generate point estimates (posterior medians)############
##############################################################

#### load in parts, merge them, save out

load("sst_DDM_inf_ADHD_first.RData")
load("sst_DDM_inf_ADHD_second.RData")

sSST.base<-sSST.base.first
sSST.base[4501:(4500+length(sSST.base.second))]<-sSST.base.second
names(sSST.base)<-c(names(sSST.base.first),names(sSST.base.second))

save(sSST.base,file="sst_ddm_inf_ADHD_samples.RData")

#### check convergence, remove non-converged if needed

# load("sst_ddm_inf_ADHD_samples.RData")

rhat.sst<-gelman.diag.dmc(sSST.base, hyper = FALSE)
save(rhat.sst,file="sst_ddm_inf_ADHD_rhat.RData")

# model fit for three sites 

load("sst_ddm_inf_ADHD_samples.RData")
load("fit_IDs.RData")

sim.SST.s17 <- h.post.predict.dmc(sSST.base[names(sSST.base)%in%s17.ids],cores=3)
sim.SST.s01 <- h.post.predict.dmc(sSST.base[names(sSST.base)%in%s01.ids],cores=3)
sim.SST.s05 <- h.post.predict.dmc(sSST.base[names(sSST.base)%in%s05.ids],cores=3)

save(sim.SST.s17,
     sim.SST.s01,
     sim.SST.s05,file="sst_ddm_3sites_pp.RData")

#### summarize parameter estimates

# loop to make summary parameters

sSST.base.pars<-sSST.base

for (s in 1:length(sSST.base)){
  t<-sSST.base[[s]]$theta; d<-dim(t)
  Ter<-t[,"t0",]+(t[,"st0",]/2)
  D <- d + c(0, 1, 0)
  t2<-array(sapply(1:D[3], function(x) cbind(t[,,x], Ter[,x])), D)
  dimnames(t2)[[2]]<-c(dimnames(t)[[2]],"Ter")
  sSST.base.pars[[s]]$theta<-t2
}


# save out point estimates 

sSST.base_medians<-lapply(sSST.base.pars,FUN=function(x) apply(x$theta,2,median))
sSST.base_medians<-as.data.frame(t(as.data.frame(sSST.base_medians)))

write.csv(sSST.base_medians,file="sst_ddm_inf_ADHD_postmedians.csv")

