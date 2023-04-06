rm(list=ls())

# load all necessary packages
library(lavaan)
library(psych)
library(stringi)
library(mice)

#####################################################
######### 1. Select participants ####################
#####################################################

# load function to read ABCD data from a specified directory
# (the path for your local directory containing ABCD data needs
# to be added as the default for the "path" argument in this file)
source("read_abcd.R")

#load in CBCL and BPM-T raw scores
abcd_cbcl01<-read.abcd("abcd_cbcl01.txt")
abcd_bpmt01<-read.abcd("abcd_bpmt01.txt")
acspsw03<-read.abcd("acspsw03.txt")
abcd_lt01<-read.abcd("abcd_lt01.txt")

# relevant family and site data

acspsw03.base<-acspsw03[acspsw03$eventname=="baseline_year_1_arm_1",]
abcd_lt01.base<-abcd_lt01[abcd_lt01$eventname=="baseline_year_1_arm_1",]

family.site<-merge(acspsw03.base[,c("subjectkey","rel_family_id")],
                   abcd_lt01.base[,c("subjectkey","site_id_l")])

# several families have individuals at multiple sites
sort(rowSums(table(family.site[,2:3])==0))[1:10]

# how complete are the columns for the Attention scales?
cbcl.Attention.cols<-c("cbcl_q01_p","cbcl_q08_p","cbcl_q10_p",
                       "cbcl_q13_p","cbcl_q17_p","cbcl_q41_p",
                       "cbcl_q45_p","cbcl_q61_p","cbcl_q62_p",
                       "cbcl_q80_p")
cbcl.ADHD.cols<-c("cbcl_q04_p","cbcl_q08_p","cbcl_q10_p","cbcl_q41_p",
                  "cbcl_q78_p","cbcl_q93_p","cbcl_q104_p")
bpmt.Attention.cols<-c("bpmt_q1","bpmt_q3","bpmt_q4",
                       "bpmt_q5","bpmt_q9","bpmt_q13")

abcd_cbcl01$pc.A1<-rowMeans(!is.na(abcd_cbcl01[,cbcl.Attention.cols]))
abcd_cbcl01$pc.A2<-rowMeans(!is.na(abcd_cbcl01[,cbcl.ADHD.cols]))
abcd_bpmt01$pc.A<-rowMeans(!is.na(abcd_bpmt01[,bpmt.Attention.cols]!=""))

# exclude anyone with missing data
cbcl.inc<-(abcd_cbcl01$pc.A1==1 & abcd_cbcl01$pc.A2==1)
bpmt.inc<-(abcd_bpmt01$pc.A==1)

# break into separate data sets for baseline ("base"),
# year 1 ("y1"), and year 2("y2)
cbcl.base<-abcd_cbcl01$eventname=="baseline_year_1_arm_1"
cbcl.y1<-abcd_cbcl01$eventname=="1_year_follow_up_y_arm_1"
cbcl.y2<-abcd_cbcl01$eventname=="2_year_follow_up_y_arm_1"
bpmt.base<-abcd_bpmt01$eventname=="baseline_year_1_arm_1"
bpmt.y1<-abcd_bpmt01$eventname=="1_year_follow_up_y_arm_1"
bpmt.y2<-abcd_bpmt01$eventname=="2_year_follow_up_y_arm_1"


# included samples with complete data

cbcl.base.subs<-abcd_cbcl01[cbcl.inc & cbcl.base,"subjectkey"]
bpmt.base.subs<-abcd_bpmt01[bpmt.inc & bpmt.base,"subjectkey"]
base.all<-cbcl.base.subs[cbcl.base.subs%in%bpmt.base.subs]

cbcl.y1.subs<-abcd_cbcl01[cbcl.inc & cbcl.y1,"subjectkey"]
bpmt.y1.subs<-abcd_bpmt01[bpmt.inc & bpmt.y1,"subjectkey"]
y1.all<-cbcl.y1.subs[cbcl.y1.subs%in%bpmt.y1.subs]

cbcl.y2.subs<-abcd_cbcl01[cbcl.inc & cbcl.y2,"subjectkey"]
bpmt.y2.subs<-abcd_bpmt01[bpmt.inc & bpmt.y2,"subjectkey"]
y2.all<-cbcl.y2.subs[cbcl.y2.subs%in%bpmt.y2.subs]


# same but regardless of completion

# not bad Ns

length(base.all)
#[1] 3960 

length(y1.all)
# [1] 5915

length(y2.all)
#[1] 4164

# merge data for each

base.dat<-merge(abcd_cbcl01[abcd_cbcl01$subjectkey%in%base.all & cbcl.base,],
            abcd_bpmt01[abcd_bpmt01$subjectkey%in%base.all & bpmt.base,c("subjectkey",bpmt.Attention.cols)])

y1.dat<-merge(abcd_cbcl01[abcd_cbcl01$subjectkey%in%y1.all & cbcl.y1,],
                abcd_bpmt01[abcd_bpmt01$subjectkey%in%y1.all & bpmt.y1,c("subjectkey",bpmt.Attention.cols)])

y2.dat<-merge(abcd_cbcl01[abcd_cbcl01$subjectkey%in%y2.all & cbcl.y2,],
              abcd_bpmt01[abcd_bpmt01$subjectkey%in%y2.all & bpmt.y2,c("subjectkey",bpmt.Attention.cols)])


# merge full sample data sets (regardless of missingness)

base.dat.full<-merge(data.frame(subjectkey=unique(acspsw03$subjectkey)),
                     abcd_cbcl01[abcd_cbcl01$eventname=="baseline_year_1_arm_1",],all.x=TRUE)
base.dat.full<-merge(base.dat.full,
                     abcd_bpmt01[abcd_bpmt01$eventname=="baseline_year_1_arm_1",c("subjectkey",bpmt.Attention.cols)],all.x=TRUE)

y1.dat.full<-merge(data.frame(subjectkey=unique(acspsw03$subjectkey)),
                     abcd_cbcl01[abcd_cbcl01$eventname=="1_year_follow_up_y_arm_1",],all.x=TRUE)
y1.dat.full<-merge(y1.dat.full,
                     abcd_bpmt01[abcd_bpmt01$eventname=="1_year_follow_up_y_arm_1",c("subjectkey",bpmt.Attention.cols)],all.x=TRUE)


y2.dat.full<-merge(data.frame(subjectkey=unique(acspsw03$subjectkey)),
                   abcd_cbcl01[abcd_cbcl01$eventname=="2_year_follow_up_y_arm_1",],all.x=TRUE)
y2.dat.full<-merge(y2.dat.full,
                   abcd_bpmt01[abcd_bpmt01$eventname=="2_year_follow_up_y_arm_1",c("subjectkey",bpmt.Attention.cols)],all.x=TRUE)


# add site and family

base.dat<-merge(base.dat,family.site)
y1.dat<-merge(y1.dat,family.site)
y2.dat<-merge(y2.dat,family.site)

base.dat.full<-merge(base.dat.full,family.site)
y1.dat.full<-merge(y1.dat.full,family.site)
y2.dat.full<-merge(y2.dat.full,family.site)


# remove anyone who has siblings across sites

multi.site.fams<-sort(rowSums(table(family.site[,c("rel_family_id","site_id_l")])==0))[1:7]
multi.site.fams<-names(multi.site.fams)

base.dat<-base.dat[!base.dat$rel_family_id%in%multi.site.fams,]
y1.dat<-y1.dat[!y1.dat$rel_family_id%in%multi.site.fams,]
y2.dat<-y2.dat[!y2.dat$rel_family_id%in%multi.site.fams,]

base.dat.full<-base.dat.full[!base.dat.full$rel_family_id%in%multi.site.fams,]
y1.dat.full<-y1.dat.full[!y1.dat.full$rel_family_id%in%multi.site.fams,]
y2.dat.full<-y2.dat.full[!y2.dat.full$rel_family_id%in%multi.site.fams,]

#remove site 22 (small n)
base.dat<-base.dat[base.dat$site_id_l!="site22",]
y1.dat<-y1.dat[y1.dat$site_id_l!="site22",]
y2.dat<-y2.dat[y2.dat$site_id_l!="site22",]

base.dat.full<-base.dat.full[base.dat.full$site_id_l!="site22",]
y1.dat.full<-y1.dat.full[y1.dat.full$site_id_l!="site22",]
y2.dat.full<-y2.dat.full[y2.dat.full$site_id_l!="site22",]

# remove anyone with completely missing data across both CBCL and BPMT

all.c<-unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))

base.dat.full<-base.dat.full[rowSums(is.na(base.dat.full[,all.c]))<length(all.c),]
y1.dat.full<-y1.dat.full[rowSums(is.na(y1.dat.full[,all.c]))<length(all.c),]
y2.dat.full<-y2.dat.full[rowSums(is.na(y2.dat.full[,all.c]))<length(all.c),]

# determine who can be included in hierarchical
# fitting subgroup for DDM analyses

hier.ddm<-family.site

# include only singletons

singles<-table(hier.ddm$rel_family_id)
singles<-singles[singles==1]
singles<-names(singles)

hier.ddm<-hier.ddm[hier.ddm$rel_family_id%in%singles,]

# remove anyone from held-out sites

hier.ddm<-hier.ddm[!hier.ddm$site_id_l%in%v.sites,]

# remove anyone in current prediction data sets

hier.ddm<-hier.ddm[!hier.ddm$subjectkey%in%base.all,]
hier.ddm<-hier.ddm[!hier.ddm$subjectkey%in%y1.all,]
hier.ddm<-hier.ddm[!hier.ddm$subjectkey%in%y2.all,]

hier.ddm<-hier.ddm$subjectkey

write.csv(hier.ddm,
          file="hier_DDM_eligible_ADHD_proj.csv",
          row.names = FALSE)

# make sure prior sample individuals are removed from all data sets 
prior.ids<-read.csv("subsample_ADHD_project.csv")

base.dat$subjectkey[base.dat$subjectkey%in%prior.ids$x]
y1.dat$subjectkey[y1.dat$subjectkey%in%prior.ids$x]
y2.dat$subjectkey[y2.dat$subjectkey%in%prior.ids$x]

base.dat.full<-base.dat.full[!base.dat.full$subjectkey%in%prior.ids$x,]
y1.dat.full<-y1.dat.full[!y1.dat.full$subjectkey%in%prior.ids$x,]
y2.dat.full<-y2.dat.full[!y2.dat.full$subjectkey%in%prior.ids$x,]

###########################################
### Select training and validation data ###
###########################################

# year 1 has better Ns per site by a lot
base.ns<-table(base.dat$site_id_l)
y1.ns<-table(y1.dat$site_id_l)

# sample one site from each 
qs<-quantile(y1.ns,probs = c(0,.25,.50,.75))
qs
sample(y1.ns[y1.ns>qs[2] & y1.ns<=qs[3]],1)
# site10 
sample(y1.ns[y1.ns>qs[3] & y1.ns<=qs[4]],1)
# site06 
sample(y1.ns[y1.ns>qs[4]],1)
# site20 

# IDs for prediction and validation (i.e., "lockbox") sites 

v.sites<-c("site10","site06","site20")
p.sites<-names(y1.ns)[!names(y1.ns)%in%v.sites]

base.p.dat<-base.dat[base.dat$site_id_l%in%p.sites,]
y1.p.dat<-y1.dat[y1.dat$site_id_l%in%p.sites,]
y2.p.dat<-y2.dat[y2.dat$site_id_l%in%p.sites,]

base.v.dat<-base.dat[base.dat$site_id_l%in%v.sites,]
y1.v.dat<-y1.dat[y1.dat$site_id_l%in%v.sites,]
y2.v.dat<-y2.dat[y2.dat$site_id_l%in%v.sites,]


base.full.p.dat<-base.dat.full[base.dat.full$site_id_l%in%p.sites,]
y1.full.p.dat<-y1.dat.full[y1.dat.full$site_id_l%in%p.sites,]
y2.full.p.dat<-y2.dat.full[y2.dat.full$site_id_l%in%p.sites,]

base.full.v.dat<-base.dat.full[base.dat.full$site_id_l%in%v.sites,]
y1.full.v.dat<-y1.dat.full[y1.dat.full$site_id_l%in%v.sites,]
y2.full.v.dat<-y2.dat.full[y2.dat.full$site_id_l%in%v.sites,]

# convert columns of interest to numeric

for (c in c(bpmt.Attention.cols,cbcl.Attention.cols,cbcl.ADHD.cols)){
  base.p.dat[,c]<-as.numeric(base.p.dat[,c])
  base.v.dat[,c]<-as.numeric(base.v.dat[,c])
  y1.p.dat[,c]<-as.numeric(y1.p.dat[,c])
  y1.v.dat[,c]<-as.numeric(y1.v.dat[,c])
  y2.p.dat[,c]<-as.numeric(y2.p.dat[,c])
  y2.v.dat[,c]<-as.numeric(y2.v.dat[,c])
  base.full.p.dat[,c]<-as.numeric(base.full.p.dat[,c])
  base.full.v.dat[,c]<-as.numeric(base.full.v.dat[,c])
  y1.full.p.dat[,c]<-as.numeric(y1.full.p.dat[,c])
  y1.full.v.dat[,c]<-as.numeric(y1.full.v.dat[,c])
  y2.full.p.dat[,c]<-as.numeric(y2.full.p.dat[,c])
  y2.full.v.dat[,c]<-as.numeric(y2.full.v.dat[,c])}


#save out each
save(base.p.dat,y1.p.dat,y1.p.dat,
     file="prediction_data.RData")
save(base.v.dat,y1.v.dat,y2.v.dat,
     file="validation_data.RData")
save(base.full.p.dat,y1.full.p.dat,y2.full.p.dat,
     file="prediction_data_full.RData")
save(base.full.v.dat,y1.full.v.dat,y2.full.v.dat,
     file="validation_data_full.RData")



###############################################################
######### 2. Compute Demographic Variables ####################
###############################################################

# load in relevant data structures
acspsw03<-read.abcd("acspsw03.txt")
pdem02<-read.abcd("pdem02.txt")
abcd_lt01<-read.abcd("abcd_lt01.txt")

#baseline values to match up with pdem
acspsw03.base<-acspsw03[acspsw03$eventname=="baseline_year_1_arm_1" ,]
abcd_lt01.base<-abcd_lt01[abcd_lt01$eventname=="baseline_year_1_arm_1" ,]

#merge
ABCD.DEAP.demo<-merge(acspsw03.base,pdem02,by.x = "subjectkey",
                      by.y="subjectkey")
ABCD.DEAP.demo<-merge(ABCD.DEAP.demo,
                      abcd_lt01.base[,c("subjectkey","site_id_l" )],
                      by.x = "subjectkey",
                      by.y="subjectkey")


# take the steps outlined in the ABCD scripts below to produce
# the same codings for demographics as are used in DEAP

#https://github.com/ABCD-STUDY/analysis-nda/blob/master/notebooks/general/core_demographics3.0.R
#https://github.com/ABCD-STUDY/analysis-nda/blob/master/notebooks/general/categorical_extension3.0.R

# Age
ABCD.DEAP.demo$age = ABCD.DEAP.demo$interview_age.x

# Female (gender)
ABCD.DEAP.demo$female = ABCD.DEAP.demo$sex.x
ABCD.DEAP.demo$female[which(ABCD.DEAP.demo$female=="")]=NA
ABCD.DEAP.demo$female = factor(as.numeric(ABCD.DEAP.demo$female == "F"), 
                               levels = 0:1, labels = c("no", "yes") ) 

# Household income
household.income = ABCD.DEAP.demo$demo_comb_income_v2
household.income[household.income == "1"] = 1 # "[<50K]"
household.income[household.income == "2"] = 1 # "[<50K]"
household.income[household.income == "3"] = 1 # "[<50K]"
household.income[household.income == "4"] = 1 # "[<50K]"
household.income[household.income == "5"] = 1 # "[<50K]"
household.income[household.income == "6"] = 1 # "[<50K]"
household.income[household.income == "7"] = 2 # "[>=50K & <100K]"
household.income[household.income == "8"] = 2 # "[>=50K & <100K]"
household.income[household.income == "9"] = 3 # "[>=100K]"
household.income[household.income == "10"] = 3 # "[>=100K]"
household.income[household.income == "777"] = NA
household.income[household.income == "999"] = NA
household.income[household.income %in% c(NA, "999", "777")] = NA
ABCD.DEAP.demo$household.income = factor( household.income, levels= 1:3, labels = c("[<50K]", "[>=50K & <100K]", "[>=100K]") )

# highest education
high.educ1 = ABCD.DEAP.demo$demo_prnt_ed_v2
high.educ2 = ABCD.DEAP.demo$demo_prtnr_ed_v2
high.educ1[which(high.educ1 == "999")] = NA
high.educ2[which(high.educ2 == "999")] = NA
high.educ1[which(high.educ1 == "777")] = NA
high.educ2[which(high.educ2 == "777")] = NA
high.educ = pmax(as.numeric(as.character(high.educ1)), as.numeric(as.character(high.educ2)), na.rm=T)
idx <- which(high.educ %in% 0:12, arr.ind = TRUE)
high.educ[idx] = 1 # "< HS Diploma"
idx <- which(high.educ %in% 13:14, arr.ind = TRUE)
high.educ[idx] = 2 # "HS Diploma/GED"
idx <- which(high.educ %in% 15:17, arr.ind = TRUE)
high.educ[idx] = 3 # "Some College"
idx <- which(high.educ == 18, arr.ind = TRUE)
high.educ[idx] = 4 # "Bachelor"
idx <- which(high.educ %in% 19:21, arr.ind = TRUE)
high.educ[idx] = 5 # "Post Graduate Degree"
high.educ[which(high.educ == "999")]=NA
high.educ[which(high.educ == "777")]=NA
ABCD.DEAP.demo$high.educ = factor( high.educ, levels= 1:5, 
                                   labels = c("< HS Diploma","HS Diploma/GED",
                                              "Some College","Bachelor","Post Graduate Degree") )


# marital status
married = rep(NA, length(ABCD.DEAP.demo$demo_prnt_marital_v2))
married[ABCD.DEAP.demo$demo_prnt_marital_v2 == 1] = 1
married[ABCD.DEAP.demo$demo_prnt_marital_v2 %in% 2:6] = 0
ABCD.DEAP.demo$married = factor( married, levels= 0:1, labels = c("no", "yes") )

# Add another variable that also includes couples that just live together. 
married.livingtogether = rep(NA, length(ABCD.DEAP.demo$demo_prnt_marital_v2))
married.livingtogether[ABCD.DEAP.demo$demo_prnt_marital_v2 %in% c(1,6)] = 1
married.livingtogether[ABCD.DEAP.demo$demo_prnt_marital_v2 %in% 2:5] = 0
ABCD.DEAP.demo$married.or.livingtogether = factor( married.livingtogether, levels= 0:1, labels = c("no", "yes") )

# race

ABCD.DEAP.demo$white=as.numeric(ABCD.DEAP.demo$demo_race_a_p___10)

ABCD.DEAP.demo$black=as.numeric(ABCD.DEAP.demo$demo_race_a_p___11)

ABCD.DEAP.demo$asian = 0
ABCD.DEAP.demo[ABCD.DEAP.demo$demo_race_a_p___18==1 | ABCD.DEAP.demo$demo_race_a_p___19==1  |
                 ABCD.DEAP.demo$demo_race_a_p___20==1  |  ABCD.DEAP.demo$demo_race_a_p___21==1 | 
                 ABCD.DEAP.demo$demo_race_a_p___22==1 | ABCD.DEAP.demo$demo_race_a_p___23==1 |
                 ABCD.DEAP.demo$demo_race_a_p___24==1,"asian"]<-1

ABCD.DEAP.demo$aian = 0
ABCD.DEAP.demo[ABCD.DEAP.demo$demo_race_a_p___12==1 | ABCD.DEAP.demo$demo_race_a_p___13==1,"aian"]<-1

ABCD.DEAP.demo$nhpi = 0
ABCD.DEAP.demo[ABCD.DEAP.demo$demo_race_a_p___14==1 |  ABCD.DEAP.demo$demo_race_a_p___15==1 |
                 ABCD.DEAP.demo$demo_race_a_p___16==1 | ABCD.DEAP.demo$demo_race_a_p___17==1,"nhpi"]<-1

ABCD.DEAP.demo$other = as.numeric(ABCD.DEAP.demo$demo_race_a_p___25)

mixed.num = rowSums(ABCD.DEAP.demo[,c("white","black","asian","aian","nhpi","other")])

ABCD.DEAP.demo$mixed=NA
ABCD.DEAP.demo[mixed.num<=1,"mixed"]=0
ABCD.DEAP.demo[mixed.num>1,"mixed"]=1

# Race 4 level

ABCD.DEAP.demo$race.4level=NA
ABCD.DEAP.demo[ABCD.DEAP.demo$white==1,]$race.4level=1
ABCD.DEAP.demo[ABCD.DEAP.demo$black==1,]$race.4level=2
ABCD.DEAP.demo[ABCD.DEAP.demo$asian==1,]$race.4level=3
ABCD.DEAP.demo[ABCD.DEAP.demo$aian==1,]$race.4level=4
ABCD.DEAP.demo[ABCD.DEAP.demo$nhpi==1,]$race.4level=4
ABCD.DEAP.demo[ABCD.DEAP.demo$other==1,]$race.4level=4
ABCD.DEAP.demo[ABCD.DEAP.demo$mixed==1,]$race.4level=4

ABCD.DEAP.demo$race.4level<- factor(ABCD.DEAP.demo$race.4level,
                                    levels = 1:4,
                                    labels = c("White","Black","Asian","Other/Mixed"))


# Race 6 level

ABCD.DEAP.demo$race.6level=NA
ABCD.DEAP.demo[ABCD.DEAP.demo$white==1,]$race.6level=1
ABCD.DEAP.demo[ABCD.DEAP.demo$black==1,]$race.6level=2
ABCD.DEAP.demo[ABCD.DEAP.demo$asian==1,]$race.6level=3
ABCD.DEAP.demo[ABCD.DEAP.demo$aian==1,]$race.6level=4
ABCD.DEAP.demo[ABCD.DEAP.demo$nhpi==1,]$race.6level=4
ABCD.DEAP.demo[ABCD.DEAP.demo$other==1,]$race.6level=5
ABCD.DEAP.demo[ABCD.DEAP.demo$mixed==1,]$race.6level=6

ABCD.DEAP.demo$race.6level<- factor(ABCD.DEAP.demo$race.6level,
                                    levels = 1:6,
                                    labels = c("White","Black","Asian","AIAN/NHPI","Other","Mixed"))



# Hispanic ethnicity

ABCD.DEAP.demo$hisp=NA
ABCD.DEAP.demo[!is.na(ABCD.DEAP.demo$demo_ethn_v2) & ABCD.DEAP.demo$demo_ethn_v2==1,"hisp"]=1
ABCD.DEAP.demo[!is.na(ABCD.DEAP.demo$demo_ethn_v2) & ABCD.DEAP.demo$demo_ethn_v2==2,"hisp"]=0

ABCD.DEAP.demo$hisp<- factor(ABCD.DEAP.demo$hisp,
                             levels=0:1,
                             labels= c("no","yes"))



# write out only relevant variables
ABCD.DEAP.demo.small<-ABCD.DEAP.demo[,c("subjectkey","rel_family_id","site_id_l",
                                        "age","female","household.income","high.educ",
                                        "married","married.or.livingtogether",
                                        "race.4level","race.6level","hisp")]
write.csv(ABCD.DEAP.demo.small,file="DEAP_reconstructed_demo_r4.csv",row.names = FALSE)


###############################################################
######### 3. Compute Stimulant Med Use Variable ###############
###############################################################

# name relevant child medication columns in ABCD data

med.cols.medsy01<-c(paste0("med",seq(1,15,1),"_rxnorm_p"))

med.cols.plus01<-c(paste0("pls2_p_med",seq(1,5,1),"_rxnorm"),
                   paste0("pls3_med",seq(1,2,1),"_rxnorm_p"))

med.cols.abcd_plus01<-c(paste0("pls1_med",seq(1,5,1),"_rxnorm"),
                        paste0("pls2_med",seq(1,6,1),"_rxnorm"),
                        paste0("pls3_med",seq(1,5,1)), 
                        paste0("pls4_med",seq(1,5,1)))

med.cols.all<-c(med.cols.medsy01,med.cols.plus01,med.cols.abcd_plus01)

# load med data sources

medsy01<-read.abcd("medsy01.txt")
medsy01<-medsy01[medsy01$eventname=="baseline_year_1_arm_1",]
plus01<-read.abcd("plus01.txt")
plus01<-plus01[plus01$eventname=="baseline_year_1_arm_1",]
abcd_plus01<-read.abcd("abcd_plus01.txt")
abcd_plus01<-abcd_plus01[abcd_plus01$eventname=="baseline_year_1_arm_1",]

# load drug names of from RxNorm

# all drugs classified in RxNorm as:
# "Centrally acting sympathomimetics"
RxNorm.CAS<-list.files("Med_classes/RXnorm_CAS/")

# add two other drugs FDA approved for ADHD treatment:
# Clonidine and Guanfacine
RxNorm.non.stims<-list.files("Med_classes/RXnorm_approved_nonstims/")

# retired "non-drug" categories present in the ABCD data
RxNorm.retired<-read.csv("Med_classes/RXnorm_likely_meds/likely.csv")
RxNorm.retired.non.stims<-RxNorm.retired[RxNorm.retired$CAS==FALSE,]
RxNorm.retired.stims<-RxNorm.retired[RxNorm.retired$CAS==TRUE,]

# generate lists of all possible formulations of the
# stimulant and non-stimulant medications

CAS.list<-list()

for (g in 1:length(RxNorm.CAS)){
  tmp<-read.csv(paste0("Med_classes/RxNorm_all_meds/",RxNorm.CAS[g]),skip = 1,row.names = NULL)
  colnames(tmp)<-c(colnames(tmp)[-1],"gen");tmp$gen<-RxNorm.CAS[g]
  CAS.list[[g]]<-tmp
}

CAS.list<-do.call(rbind,CAS.list)

non.stim.list<-list()

for (g in 1:length(RxNorm.non.stims)){
  tmp<-read.csv(paste0("Med_classes/RxNorm_all_meds/",RxNorm.non.stims[g]),skip = 1,row.names = NULL)
  colnames(tmp)<-c(colnames(tmp)[-1],"gen");tmp$gen<-RxNorm.non.stims[g]
  non.stim.list[[g]]<-tmp
}

non.stim.list<-do.call(rbind,non.stim.list)

# however, some rxcui numbers are not specific to one medication, 
# including those with "Term Type" DF and DFG

# not present in lists: GPCK, BPCK, PSN, SY, TMSY, ET

# delete those term types from all lists:

CAS.list<-CAS.list[CAS.list$termType!="DF" & CAS.list$termType!="DFG",]
non.stim.list<-non.stim.list[non.stim.list$termType!="DF" & non.stim.list$termType!="DFG",]

# finally, add relevant retired rxcuis the their respective lists

CAS.retired<-data.frame(NA,RxNorm.retired.stims$rxcui,
                        RxNorm.retired.stims$name,NA,NA,NA,
                        RxNorm.retired.stims$active)
colnames(CAS.retired)<-colnames(CAS.list)
CAS.list<-rbind(CAS.list,CAS.retired)


non.stim.retired<-data.frame(NA,RxNorm.retired.non.stims$rxcui,
                             RxNorm.retired.non.stims$name,NA,NA,NA,
                             RxNorm.retired.non.stims$active)
colnames(non.stim.retired)<-colnames(non.stim.list)
non.stim.list<-rbind(non.stim.list,non.stim.retired)


# loop to generate medication use variables

ABCD.med.data<-data.frame(medsy01$subjectkey);colnames(ABCD.med.data)<-"subjectkey"
ABCD.med.data$CAS.ADHD<-FALSE
ABCD.med.data$non.stim.ADHD<-FALSE

for (s in ABCD.med.data$subjectkey){
  meds<-as.numeric(stri_extract_first_regex(medsy01[medsy01$subjectkey==s,med.cols.medsy01],"[0-9]+"))
  meds<-c(meds,as.numeric(stri_extract_first_regex(plus01[plus01$subjectkey==s,med.cols.plus01],"[0-9]+")))
  meds<-c(meds,as.numeric(stri_extract_first_regex(abcd_plus01[abcd_plus01$subjectkey==s,med.cols.abcd_plus01],"[0-9]+")))
  CAS.meds<-meds[meds%in%CAS.list$rxcui]
  non.stim.meds<-meds[meds%in%non.stim.list$rxcui]
  if(length(CAS.meds)>0){ABCD.med.data[ABCD.med.data$subjectkey==s,]$CAS.ADHD<-TRUE}
  if(length(non.stim.meds)>0){ABCD.med.data[ABCD.med.data$subjectkey==s,]$non.stim.ADHD<-TRUE}
}

ABCD.med.data$either.ADHD<-as.logical((ABCD.med.data$CAS.ADHD+ABCD.med.data$non.stim.ADHD)>0)

ABCD.med.data$only.non.stims<-as.logical((ABCD.med.data$CAS.ADHD-ABCD.med.data$non.stim.ADHD)==-1)

write.csv(ABCD.med.data,"ABCD_ADHD_meds.csv",row.names = FALSE)

#####################################################
######### 4. Compile  Features ######################
#####################################################


load("prediction_data.RData")
load("validation_data.RData")
load("prediction_data_full.RData")
load("validation_data_full.RData")

#demographic data, previously converted to DEAP categories

demo<-read.csv("DEAP_reconstructed_demo_r4.csv")
demo.cols<-c("subjectkey","age","female","race.4level","hisp","high.educ",
             "married","household.income")

base<-merge(base.p.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)
year1<-merge(y1.p.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)
year2<-merge(y2.p.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)

base.full<-merge(base.full.p.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)
year1.full<-merge(y1.full.p.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)
year2.full<-merge(y2.full.p.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)

base.v<-merge(base.v.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)
year1.v<-merge(y1.v.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)
year2.v<-merge(y2.v.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)
year1.full.v<-merge(y1.full.v.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)
year2.full.v<-merge(y2.full.v.dat,demo[,demo.cols],by="subjectkey",all.x=TRUE)

# neighborhood poverty

abcd_rhds01<-read.abcd("abcd_rhds01.txt")
adi1.cols<-colnames(abcd_rhds01)[24:40]
abcd_rhds01[,adi1.cols]<-apply(abcd_rhds01[,adi1.cols], 2, function(x) as.numeric(x))
abcd_rhds01<-as.data.frame(abcd_rhds01)
adi1<-abcd_rhds01[abcd_rhds01$eventname=="baseline_year_1_arm_1",c("subjectkey",adi1.cols)]

base<-merge(base,adi1,by="subjectkey",all.x=TRUE)
year1<-merge(year1,adi1,by="subjectkey",all.x=TRUE)
year2<-merge(year2,adi1,by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,adi1,by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,adi1,by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,adi1,by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,adi1,by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,adi1,by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,adi1,by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,adi1,by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,adi1,by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,adi1,by="subjectkey",all.x=TRUE)

# neighborhood crime, education system, and lead exposure data

neighborhood<-abcd_rhds01[abcd_rhds01$eventname=="baseline_year_1_arm_1",
                    c("subjectkey","reshist_addr1_p1tot",
                      "reshist_addr1_coi_ed_hsgrad","reshist_addr1_coi_ed_math",
                      "reshist_addr1_coi_ed_reading","reshist_addr1_coi_ed_schpov",
                      "reshist_addr1_coi_ed_prxhqece","reshist_addr1_leadrisk")]

neighborhood[,2:8]<-apply(neighborhood[,2:8], 2, function(x) as.numeric(x))

base<-merge(base,neighborhood,by="subjectkey",all.x=TRUE)
year1<-merge(year1,neighborhood,by="subjectkey",all.x=TRUE)
year2<-merge(year2,neighborhood,by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,neighborhood,by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,neighborhood,by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,neighborhood,by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,neighborhood,by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,neighborhood,by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,neighborhood,by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,neighborhood,by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,neighborhood,by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,neighborhood,by="subjectkey",all.x=TRUE)

# basic needs un-affordability

pdem02<-read.abcd("pdem02.txt")

BNU.cols<-paste0("demo_fam_exp",1:7,"_v2")

pdem02[,BNU.cols]<-pdem02[,BNU.cols]<-apply(pdem02[,BNU.cols], 2, function(x) gsub("777", NA, x) )
pdem02[,BNU.cols]<-apply(pdem02[,BNU.cols], 2, function(x) as.numeric(x))

pdem02$basic_needs_sum<-rowSums(pdem02[,BNU.cols])

base<-merge(base,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)
year1<-merge(year1,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)
year2<-merge(year2,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,pdem02[,c("subjectkey","basic_needs_sum")],by="subjectkey",all.x=TRUE)


# perceived neighborhood safety, parent and youth

abcd_pnsc01<-read.abcd("abcd_pnsc01.txt")
abcd_pnsc01<-abcd_pnsc01[abcd_pnsc01$eventname=="baseline_year_1_arm_1",]

abcd_nsc01<-read.abcd("abcd_nsc01.txt")
abcd_nsc01<-abcd_nsc01[abcd_nsc01$eventname=="baseline_year_1_arm_1",]

abcd_pnsc01<-merge(abcd_pnsc01,
                   abcd_nsc01[,c("subjectkey","neighborhood_crime_y")])

crime.cols<-c("neighborhood1r_p",
              "neighborhood2r_p",
              "neighborhood3r_p",
              "neighborhood_crime_y")

abcd_pnsc01[,crime.cols]<-apply(abcd_pnsc01[,crime.cols], 2, function(x) as.numeric(x))

cor(abcd_pnsc01[,crime.cols],use="complete")

abcd_pnsc01$neighborhood_p_sum<-rowMeans(abcd_pnsc01[,crime.cols[1:3]])

crime.inc<-c("subjectkey","neighborhood_p_sum","neighborhood_crime_y")

base<-merge(base,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)
year1<-merge(year1,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)
year2<-merge(year2,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,abcd_pnsc01[,crime.inc],by="subjectkey",all.x=TRUE)

# school environment

srpf01<-read.abcd("srpf01.txt")

env.cols<-paste0("school_",c(2:7),"_y")
eng.cols<-paste0("school_",c(8,9,10,12),"_y")
ali.cols<-paste0("school_",c(15,17),"_y")

srpf01[,c(env.cols,eng.cols,ali.cols)]<-apply(
  srpf01[,c(env.cols,eng.cols,ali.cols)], 
  2, function(x) as.numeric(x))

srpf01$sch_env<-rowMeans(srpf01[,env.cols])
srpf01$sch_eng<-rowMeans(srpf01[,eng.cols])
srpf01$sch_ali<-rowMeans(srpf01[,ali.cols])

sch.cols<-c("subjectkey","sch_env",
            "sch_eng","sch_ali")

school<-srpf01[srpf01$visit=="baseline_year_1_arm_1",sch.cols]

base<-merge(base,school,by="subjectkey",all.x=TRUE)
year1<-merge(year1,school,by="subjectkey",all.x=TRUE)
year2<-merge(year2,school,by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,school,by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,school,by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,school,by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,school,by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,school,by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,school,by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,school,by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,school,by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,school,by="subjectkey",all.x=TRUE)


# SST and Nback scanner task behavior 

ABCD_DDM_base<-read.csv("sst_ddm_inf_ADHD_postmedians.csv")
colnames(ABCD_DDM_base)<-c("subjectkey",paste0("sst.",colnames(ABCD_DDM_base)[-1]))

tmp<-read.csv("nback_postmedians_ADHD.csv")
colnames(tmp)<-c("subjectkey",colnames(tmp)[-1])

ABCD_DDM_base<-merge(ABCD_DDM_base,
                     tmp,
                     all = TRUE,by="subjectkey"); rm(tmp)

ABCD_DDM_base$zero.v.avg<-rowMeans(ABCD_DDM_base[,c("zero.v.lur","zero.v.non",
                                                    "zero.v.tar")])

ABCD_DDM_base$two.v.avg<-rowMeans(ABCD_DDM_base[,c("two.v.lur",
                                                   "two.v.non","two.v.tar")])

ABCD_DDM_base$nback.v.avg<-rowMeans(ABCD_DDM_base[,c("zero.v.lur","zero.v.non",
                                                     "zero.v.tar","two.v.lur",
                                                     "two.v.non","two.v.tar")])

ABCD_DDM_base$nback.v.diff<-(ABCD_DDM_base$zero.v.avg-ABCD_DDM_base$two.v.avg)

ABCD_DDM_base$nback.a.avg<-rowMeans(ABCD_DDM_base[,c("zero.a","two.a")])

ABCD_DDM_base$nback.a.diff<-(ABCD_DDM_base$zero.a-ABCD_DDM_base$two.a)

ABCD_DDM_base$nback.Ter.avg<-rowMeans(ABCD_DDM_base[,c("zero.Ter","two.Ter")])

ABCD_DDM_base$nback.Ter.diff<-(ABCD_DDM_base$zero.Ter-ABCD_DDM_base$two.Ter)

ddm.cols<-c("subjectkey","sst.v","sst.a","sst.Ter",
            "nback.v.avg","nback.v.diff",
            "nback.a.avg", "nback.a.diff",
            "nback.Ter.avg", "nback.Ter.diff")

ddm<-ABCD_DDM_base[,ddm.cols]

base<-merge(base,ddm,by="subjectkey",all.x=TRUE)
year1<-merge(year1,ddm,by="subjectkey",all.x=TRUE)
year2<-merge(year2,ddm,by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,ddm,by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,ddm,by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,ddm,by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,ddm,by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,ddm,by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,ddm,by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,ddm,by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,ddm,by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,ddm,by="subjectkey",all.x=TRUE)

# other cognitive tasks

abcd_ps01<-read.abcd("abcd_ps01.txt")
abcd_tbss01<-read.abcd("abcd_tbss01.txt")
lmtp201<-read.abcd("lmtp201.txt")
cct01<-read.abcd("cct01.txt")

pearson.cols<-c("pea_wiscv_tss",
                "pea_ravlt_sd_trial_i_tc" , 
                "pea_ravlt_sd_trial_ii_tc",
                "pea_ravlt_sd_trial_iii_tc",
                "pea_ravlt_sd_trial_iv_tc",
                "pea_ravlt_sd_trial_v_tc",
                "pea_ravlt_sd_listb_tc",
                "pea_ravlt_sd_trial_vi_tc",
                "pea_ravlt_ld_trial_vii_tc")

nihtb.cols<-c("nihtbx_picvocab_agecorrected","nihtbx_flanker_agecorrected",
              "nihtbx_list_agecorrected","nihtbx_cardsort_agecorrected",
              "nihtbx_pattern_agecorrected","nihtbx_picture_agecorrected",
              "nihtbx_reading_agecorrected")

pearson<-abcd_ps01[abcd_ps01$eventname=="baseline_year_1_arm_1",c("subjectkey",pearson.cols)]
pearson[,pearson.cols]<-apply(pearson[,pearson.cols], 2, function(x) as.numeric(x))
pearson$pea_ravlt_total<-rowSums(pearson[,c("pea_ravlt_sd_trial_i_tc" , 
                                                "pea_ravlt_sd_trial_ii_tc",
                                                "pea_ravlt_sd_trial_iii_tc",
                                                "pea_ravlt_sd_trial_iv_tc",
                                                "pea_ravlt_sd_trial_v_tc",
                                                "pea_ravlt_sd_listb_tc",
                                                "pea_ravlt_sd_trial_vi_tc",
                                                "pea_ravlt_ld_trial_vii_tc")])
pearson.cols<-c("pea_wiscv_tss","pea_ravlt_total")

nihtb<-abcd_tbss01[abcd_tbss01$eventname=="baseline_year_1_arm_1",c("subjectkey",nihtb.cols)]
nihtb[,nihtb.cols]<-apply(nihtb[,nihtb.cols], 2, function(x) as.numeric(x))

lmt<-lmtp201[lmtp201$eventname=="baseline_year_1_arm_1",c("subjectkey","lmt_scr_perc_correct")]
lmt[,"lmt_scr_perc_correct"]<-as.numeric(lmt[,"lmt_scr_perc_correct"])

cct<-cct01[cct01$eventname=="baseline_year_1_arm_1",c("subjectkey","cash_choice_task")]
cct[,"cash_choice_task"]<-as.numeric(cct[,"cash_choice_task"])
cct[cct$cash_choice_task==3 & !is.na(cct$cash_choice_task),"cash_choice_task"]<-NA
cct$cash_choice_task<-(cct$cash_choice_task-1)
cct$cash_choice_task<- factor(cct$cash_choice_task,
                             levels=0:1,
                             labels= c("no","yes"))


other.cog<-merge(pearson,nihtb,by="subjectkey",all=TRUE)
other.cog<-merge(other.cog,lmt,by="subjectkey",all=TRUE)
other.cog<-merge(other.cog,cct,by="subjectkey",all=TRUE)


base<-merge(base,other.cog,by="subjectkey",all.x=TRUE)
year1<-merge(year1,other.cog,by="subjectkey",all.x=TRUE)
year2<-merge(year2,other.cog,by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,other.cog,by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,other.cog,by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,other.cog,by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,other.cog,by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,other.cog,by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,other.cog,by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,other.cog,by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,other.cog,by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,other.cog,by="subjectkey",all.x=TRUE)

# Family conflict (parent and youth)

abcd_fes01<-read.abcd("abcd_fes01.txt")
fes02<-read.abcd("fes02.txt")

p.conf.cols<-paste0("fam_enviro",1:9,
                    c("","r","","r","","","r","","r"),
                    "_p")
y.conf.cols<-paste0("fes_youth_q",1:9)

abcd_fes01[,y.conf.cols]<-apply(abcd_fes01[,y.conf.cols], 2, function(x) as.numeric(x))
fes02[,p.conf.cols]<-apply(fes02[,p.conf.cols], 2, function(x) as.numeric(x))


abcd_fes01$conflict_youth<-rowSums(abcd_fes01[,y.conf.cols])
fes02$conflict_parent<-rowSums(fes02[,p.conf.cols])

abcd_fes01<-abcd_fes01[abcd_fes01$eventname=="baseline_year_1_arm_1",]
fes02<-fes02[fes02$eventname=="baseline_year_1_arm_1" & fes02$dataset_id==34291,]

fes<-merge(abcd_fes01,fes02,by="subjectkey",all=TRUE)

base<-merge(base,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)
year1<-merge(year1,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)
year2<-merge(year2,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,fes[,c("subjectkey","conflict_youth","conflict_parent")],by="subjectkey",all.x=TRUE)


# Parental monitoring

pmq01<-read.abcd("pmq01.txt")

pmq.cols<-paste0("parent_monitor_q",1:5,"_y")

pmq01[,pmq.cols]<-apply(pmq01[,pmq.cols], 2, function(x) as.numeric(x))

pmq01$parent_monit<-rowMeans(pmq01[,pmq.cols])

pmq01<-pmq01[pmq01$eventname=="baseline_year_1_arm_1",]

base<-merge(base,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)
year1<-merge(year1,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)
year2<-merge(year2,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,pmq01[,c("subjectkey","parent_monit")],by="subjectkey",all.x=TRUE)

# BIS/BAS and UPPS

abcd_mhy02<-read.abcd("abcd_mhy02.txt")

bis.bas<-c("bis_y_ss_bis_sum","bis_y_ss_bas_rr","bis_y_ss_bas_drive","bis_y_ss_bas_fs")
upps<-c("upps_y_ss_negative_urgency","upps_y_ss_lack_of_planning",
        "upps_y_ss_sensation_seeking","upps_y_ss_positive_urgency",
        "upps_y_ss_lack_of_perseverance")


abcd_mhy02[,c(bis.bas,upps)]<-apply(abcd_mhy02[,c(bis.bas,upps)], 2, function(x) as.numeric(x))

abcd_mhy02<-abcd_mhy02[abcd_mhy02$eventname=="baseline_year_1_arm_1",]

base<-merge(base,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)
year1<-merge(year1,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)
year2<-merge(year2,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,abcd_mhy02[,c("subjectkey",bis.bas,upps)],by="subjectkey",all.x=TRUE)


# BMI and waist circumference

abcd_ant01<-read.abcd("abcd_ant01.txt")

abcd_ant01$anthroheightcalc<-as.numeric(abcd_ant01$anthroheightcalc)*0.0254
abcd_ant01$anthroweightcalc<-as.numeric(abcd_ant01$anthroweight1lb)*0.45359237
abcd_ant01$anthro_waist_cm<-as.numeric(abcd_ant01$anthro_waist_cm)
abcd_ant01$BMI<-abcd_ant01$anthroweightcalc/abcd_ant01$anthroheightcalc^2

#trim implausible BMI values (clearly height/weight measurement errors)
abcd_ant01[(abcd_ant01$BMI>100 | abcd_ant01$BMI<10) & !is.na(abcd_ant01$BMI),]$BMI<-NA

abcd_ant01<-abcd_ant01[abcd_ant01$eventname=="baseline_year_1_arm_1",]

base<-merge(base,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)


base<-merge(base,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)
year1<-merge(year1,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)
year2<-merge(year2,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,abcd_ant01[,c("subjectkey","anthro_waist_cm","BMI")],by="subjectkey",all.x=TRUE)

# screen media

abcd_stq01<-read.abcd("abcd_stq01.txt")

wkdy.cols<-paste0("screen",c(1:5,""),"_wkdy_y")
wknd.cols<-paste0("screen",c(7:12),"_wknd_y")

abcd_stq01[,c(wkdy.cols,wknd.cols)]<-apply(abcd_stq01[,c(wkdy.cols,wknd.cols)], 2, function(x) as.numeric(x))

abcd_stq01$wkdy_screen<-rowSums(abcd_stq01[,wkdy.cols])
abcd_stq01$wknd_screen<-rowSums(abcd_stq01[,wknd.cols])

abcd_stq01<-abcd_stq01[abcd_stq01$eventname=="baseline_year_1_arm_1",]

base<-merge(base,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)
year1<-merge(year1,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)
year2<-merge(year2,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,abcd_stq01[,c("subjectkey","wkdy_screen","wknd_screen")],by="subjectkey",all.x=TRUE)

# Stim meds 

meds<-read.csv("ABCD_ADHD_meds.csv")

base<-merge(base,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)

base<-merge(base,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)
year1<-merge(year1,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)
year2<-merge(year2,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)

base.full<-merge(base.full,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)
year1.full<-merge(year1.full,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)
year2.full<-merge(year2.full,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)

base.v<-merge(base.v,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)
year1.v<-merge(year1.v,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)
year2.v<-merge(year2.v,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)

base.full.v<-merge(base.full.v,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)
year1.full.v<-merge(year1.full.v,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)
year2.full.v<-merge(year2.full.v,meds[,c("subjectkey","CAS.ADHD")],by="subjectkey",all.x=TRUE)

# save out data sets

write.csv(base,file="base_full_dat.csv",row.names = F)
write.csv(year1,file="year1_full_dat.csv",row.names = F)
write.csv(year2,file="year2_full_dat.csv",row.names = F)

write.csv(base.full,file="base_full_full_dat.csv",row.names = F)
write.csv(year1.full,file="year1_full_full_dat.csv",row.names = F)
write.csv(year2.full,file="year2_full_full_dat.csv",row.names = F)

write.csv(base.v,file="base_validation_dat.csv",row.names = F)
write.csv(year1.v,file="year1_validation_dat.csv",row.names = F)
write.csv(year2.v,file="year2_validation_dat.csv",row.names = F)

write.csv(base.full.v,file="base_full_validation_dat.csv",row.names = F)
write.csv(year1.full.v,file="year1_full_validation_dat.csv",row.names = F)
write.csv(year2.full.v,file="year2_full_validation_dat.csv",row.names = F)


#####################################################
######### 5. Bifactor Measurement Model #############
#####################################################

# fit bifactor models

load("prediction_data.RData")

# to prevent family nesting from affecting standard errors, 
# randomly sample cases from each family
source("ABCD_family_subsample.R")

base.SEM<-abcd.fam.subsample(base.p.dat)
y1.SEM<-abcd.fam.subsample(y1.p.dat)
save(base.SEM,y1.SEM,file="SEM_data_final.RData")

full.cbcl.cols<-unique(c(cbcl.ADHD.cols,cbcl.Attention.cols))

bf.comp <- paste(" Attn =~ ",paste(c(bpmt.Attention.cols,comp.cbcl.cols),collapse=" + "),
                 " Parent =~ ",paste(comp.cbcl.cols,collapse=" + "),
                 " Teacher =~ ",paste(bpmt.Attention.cols,collapse=" + "),
                 " Attn ~~ 0*Parent 
                      Attn ~~ 0*Teacher
                      Parent ~~ 0*Teacher
                      ",sep="\n")


base.model.comp = cfa(bf.comp,data=base.SEM,  
                      std.lv=TRUE,
                      orthogonal=TRUE,
                      estimator = "WLSMV") 
summary(base.model.comp,fit=TRUE,standardized=TRUE)


y1.model.comp = cfa(bf.comp,data=y1.SEM,  
                    std.lv=TRUE,
                    orthogonal=TRUE,
                    estimator = "WLSMV") 
summary(y1.model.comp,fit=TRUE,standardized=TRUE)

sum.paper<-summary(y1.model.comp,fit=TRUE,standardized=TRUE)
write.csv(sum.paper$PE,file="loadings.csv")


# remove loadings for items with weak (<0.30) relations
# with each factor in larger (y1) sample

# weak on general attention factor 
# (makes sense, as not core ADHD Sx.): 
# cbcl_q13_p - Confused or seems to be in a fog
# cbcl_q45_p - Nervous, highstrung, or tense
# cbcl_q62_p - Poorly coordinated or clumsy
# cbcl_q80_p - Stares blankly

# weak on parent specific factor 
# (these are core ADHD Sx.): 
# cbcl_q01_p - Acts too young for his/her age 
# cbcl_q61_p - Poor school work

cbcl.cols.final<-comp.cbcl.cols[!comp.cbcl.cols%in%c("cbcl_q01_p",
                                                     "cbcl_q61_p")]

Attn.cols.final<-c(bpmt.Attention.cols,comp.cbcl.cols)
Attn.cols.final<-Attn.cols.final[!Attn.cols.final%in%c("cbcl_q13_p",
                                                       "cbcl_q45_p",
                                                       "cbcl_q62_p",
                                                       "cbcl_q80_p")]

bf.final <- paste(" Attn =~ ",paste(Attn.cols.final,collapse=" + "),
                  " Parent =~ ",paste(cbcl.cols.final,collapse=" + "),
                  " Teacher =~ ",paste(bpmt.Attention.cols,collapse=" + "),
                  " Attn ~~ 0*Parent 
                      Attn ~~ 0*Teacher
                      Parent ~~ 0*Teacher
                      ",sep="\n")


base.model.final = cfa(bf.final,data=base.SEM,  
                       std.lv=TRUE,
                       orthogonal=TRUE,
                       estimator = "WLSMV") 
summary(base.model.final,fit=TRUE,standardized=TRUE)


y1.model.final = cfa(bf.final,data=y1.SEM,  
                     std.lv=TRUE,
                     orthogonal=TRUE,
                     estimator = "WLSMV") 
summary(y1.model.final,fit=TRUE,standardized=TRUE)

sum.paper.constrained<-summary(y1.model.final,fit=TRUE,standardized=TRUE)
write.csv(sum.paper.constrained$PE,file="loadings_constrained.csv")



# factor scores across all models and data sets

base.SEM$Gen.bmod<-lavPredict(base.model.final)[,"Attn"]
base.SEM$Parent.bmod<-lavPredict(base.model.final)[,"Parent"]
base.SEM$Teacher.bmod<-lavPredict(base.model.final)[,"Teacher"]

base.SEM$Gen.y1mod<-lavPredict(y1.model.final,newdata = base.SEM)[,"Attn"]
base.SEM$Parent.y1mod<-lavPredict(y1.model.final,newdata = base.SEM)[,"Parent"]
base.SEM$Teacher.y1mod<-lavPredict(y1.model.final,newdata = base.SEM)[,"Teacher"]

y1.SEM$Gen.y1mod<-lavPredict(y1.model.final)[,"Attn"]
y1.SEM$Parent.y1mod<-lavPredict(y1.model.final)[,"Parent"]
y1.SEM$Teacher.y1mod<-lavPredict(y1.model.final)[,"Teacher"]

y1.SEM$Gen.y1mod.unconstrained<-lavPredict(y1.model.comp)[,"Attn"]
cor(y1.SEM$Gen.y1mod.unconstrained,y1.SEM$Gen.y1mod)

y1.SEM$Gen.bmod<-lavPredict(base.model.final,newdata = y1.SEM)[,"Attn"]
y1.SEM$Parent.bmod<-lavPredict(base.model.final,newdata = y1.SEM)[,"Parent"]
y1.SEM$Teacher.bmod<-lavPredict(base.model.final,newdata = y1.SEM)[,"Teacher"]

# relations between general factors from each of the models 
# (y1, baseline) when applied to each data set

cor(base.SEM$Gen.bmod,base.SEM$Gen.y1mod)
# [1] 0.9984602

cor(y1.SEM$Gen.bmod,y1.SEM$Gen.y1mod)
#[1] 0.9984623

# test-retest reliability

ids.both<-base.SEM$subjectkey[base.SEM$subjectkey %in% y1.SEM$subjectkey]
fs.cols<-c("subjectkey","Gen.bmod","Gen.y1mod")

base.ids.both<-base.SEM[base.SEM$subjectkey%in%ids.both,fs.cols]
y1.ids.both<-y1.SEM[y1.SEM$subjectkey%in%ids.both,fs.cols]

colnames(base.ids.both)<-c("subjectkey","Gen.bmod.base","Gen.y1mod.base")
colnames(y1.ids.both)<-c("subjectkey","Gen.bmod.y1","Gen.y1mod.y1")

test.retest<-merge(base.ids.both,y1.ids.both,by="subjectkey")

#model fit separately to each sample
ICC(test.retest[,c("Gen.bmod.base","Gen.y1mod.y1")])

#base model fit to both samples
ICC(test.retest[,c("Gen.bmod.base","Gen.bmod.y1")])

#y1 model fit to both samples
ICC(test.retest[,c("Gen.y1mod.base","Gen.y1mod.y1")])

# fit to entire baseline and y1 prediction samples

base.model.full = cfa(bf.final,data=base.p.dat,  
                      std.lv=TRUE,
                      orthogonal=TRUE,
                      estimator = "WLSMV") 
summary(base.model.full,fit=TRUE,standardized=TRUE)


y1.model.full = cfa(bf.final,data=y1.p.dat,  
                    std.lv=TRUE,
                    orthogonal=TRUE,
                    estimator = "WLSMV") 
summary(y1.model.full,fit=TRUE,standardized=TRUE)

base.p.dat$Gen.Attn<-lavPredict(base.model.full)[,"Attn"]
base.p.dat$Parent.Attn<-lavPredict(base.model.full)[,"Parent"]
base.p.dat$Teacher.Attn<-lavPredict(base.model.full)[,"Teacher"]

y1.p.dat$Gen.Attn<-lavPredict(y1.model.full)[,"Attn"]
y1.p.dat$Parent.Attn<-lavPredict(y1.model.full)[,"Parent"]
y1.p.dat$Teacher.Attn<-lavPredict(y1.model.full)[,"Teacher"]

write.csv(base.p.dat,file="baseline_attention_fac.csv",row.names = FALSE)
write.csv(y1.p.dat,file="year1_attention_fac.csv",row.names = FALSE)


#####################################################
######### 6. Multiple Imputation Analyses ###########
#####################################################

# read in

base<-read.csv("base_full_dat.csv")
year1<-read.csv("year1_full_dat.csv")
year2<-read.csv("year2_full_dat.csv")

base.full<-read.csv("base_full_full_dat.csv")
year1.full<-read.csv("year1_full_full_dat.csv")
year2.full<-read.csv("year2_full_full_dat.csv")

cbcl.Attention.cols<-c("cbcl_q01_p","cbcl_q08_p","cbcl_q10_p",
                       "cbcl_q13_p","cbcl_q17_p","cbcl_q41_p",
                       "cbcl_q45_p","cbcl_q61_p","cbcl_q62_p",
                       "cbcl_q80_p")
cbcl.ADHD.cols<-c("cbcl_q04_p","cbcl_q08_p","cbcl_q10_p","cbcl_q41_p",
                  "cbcl_q78_p","cbcl_q93_p","cbcl_q104_p")
bpmt.Attention.cols<-c("bpmt_q1","bpmt_q3","bpmt_q4",
                       "bpmt_q5","bpmt_q9","bpmt_q13")

# identify people with missing teacher data

year1.full$t.missing<-is.na(rowMeans(year1.full[,bpmt.Attention.cols]))
year2.full$t.missing<-is.na(rowMeans(year2.full[,bpmt.Attention.cols]))

#randomly remove teacher data from 1000 to test imputation efficacy

year1.full.1000<-year1.full
year1.full.1000$f.miss<-FALSE
y1.fm<-sample(year1.full.1000[year1.full.1000$t.missing==FALSE,]$subjectkey,1000,FALSE)
year1.full.1000[year1.full.1000$subjectkey%in%y1.fm,]$f.miss<-TRUE
year1.full.1000[year1.full.1000$f.miss,bpmt.Attention.cols]<-NA


save(year1.full,year2.full,year1.full.1000,file="missing_data_randomized.RData")
load("missing_data_randomized.RData")

# Impute with MICE and run standard checks

year1.imp <- mice(data = year1.full[,unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))],
                      m = 5,
                      print = FALSE)

year2.imp <- mice(data = year2.full[,unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))],
                  m = 5,
                  print = FALSE)

year1.1000.imp <- mice(data = year1.full.1000[,unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))],
                  m = 5,
                  print = FALSE)


### standard MICE checks on BPMT variables ###

# density plots to assess how well the imputed data mimics the distributions of real data:

densityplot(year1.imp, ~  bpmt_q1)
densityplot(year1.imp, ~  bpmt_q3)
densityplot(year1.imp, ~  bpmt_q4)
densityplot(year1.imp, ~  bpmt_q5)
densityplot(year1.imp, ~  bpmt_q9)
densityplot(year1.imp, ~  bpmt_q13)

densityplot(year2.imp, ~  bpmt_q1)
densityplot(year2.imp, ~  bpmt_q3)
densityplot(year2.imp, ~  bpmt_q4)
densityplot(year2.imp, ~  bpmt_q5)
densityplot(year2.imp, ~  bpmt_q9)
densityplot(year2.imp, ~  bpmt_q13)

# distribution of data by missingness/propensity score

bpmt_q4.na1 <- is.na(year1.full$bpmt_q4)
bpmt_q13.na1 <- is.na(year1.full$bpmt_q4)

train.bpmt_q4.na1<-complete(year1.1000.imp)
train.bpmt_q4.na1$miss<-bpmt_q4.na1

train.bpmt_q13.na1<-fulllete(year1.1000.imp)
train.bpmt_q13.na1$miss<-bpmt_q13.na1


fit.bpmt_q4.na1 <- glm(miss ~., train.bpmt_q4.na1,family = binomial)
ps.bpmt_q4.na1 <- rep(predict(fit.bpmt_q4.na1, type = "response"),4)

fit.bpmt_q13.na1 <- glm(miss ~., train.bpmt_q13.na1,family = binomial)
ps.bpmt_q13.na1 <- rep(predict(fit.bpmt_q13.na1, type = "response"),4)

xyplot(year1.imp , bpmt_q4 ~ ps.bpmt_q4.na1 | .imp, pch = c(1, 20), cex = c(0.8, 1.2), xlab = "Probability that bpmt_q4 is missing", ylab = "bpmt_q4")
xyplot(year1.imp , bpmt_q13 ~ ps.bpmt_q13.na1 | .imp, pch = c(1, 20), cex = c(0.8, 1.2), xlab = "Probability that bpmt_q13 is missing", ylab = "bpmt_q13")

# check correlation matrix for imputed versus non-imputed variables of interest

v.int<-c("bpmt_q4","bpmt_q5","bpmt_q13","cbcl_q08_p","cbcl_q80_p")


# mostly slightly lower in imputed data
miss.cm<-cor(complete(year1.imp,action = 0)[,v.int],use="pairwise.complete")
miss.cm 

miss.cm-cor(complete(year1.imp,action = 1)[,v.int])
cor(complete(year1.imp,action = 1)[,v.int])

miss.cm-cor(complete(year1.imp,action = 2)[,v.int])
cor(complete(year1.imp,action = 2)[,v.int])

miss.cm-cor(complete(year1.imp,action = 3)[,v.int])
cor(complete(year1.imp,action = 3)[,v.int])

miss.cm2<-cor(complete(year2.imp,action = 0)[,v.int],use="pairwise.complete")
miss.cm2 

miss.cm-miss.cm2

miss.cm2-cor(complete(year2.imp,action = 1)[,v.int])
cor(complete(year2.imp,action = 1)[,v.int])

miss.cm2-cor(complete(year2.imp,action = 2)[,v.int])
cor(complete(year2.imp,action = 2)[,v.int])

miss.cm2-cor(complete(year2.imp,action = 3)[,v.int])
cor(complete(year2.imp,action = 3)[,v.int])

### How well does it do with the 1000 who were randomly treated as missing?

# very low correlations between imputed and actual bpmt items

cor(complete(year1.1000.imp,action = 1)[year1.full.1000$f.miss,"bpmt_q4"],
     year1.full[year1.full.1000$f.miss,"bpmt_q4"])

cor(complete(year1.1000.imp,action = 1)[year1.full.1000$f.miss,"bpmt_q13"],
    year1.full[year1.full.1000$f.miss,"bpmt_q13"])

cor(complete(year1.1000.imp,action = 1)[year1.full.1000$f.miss,"bpmt_q5"],
    year1.full[year1.full.1000$f.miss,"bpmt_q5"])

# some systematic reductions in correlations between cbcl and bpmt items

cm.1000.y1<-cor(year1.full[year1.full.1000$f.miss,v.int],use="pairwise.complete")
cm.1000.y1

cor(complete(year1.1000.imp,action = 1)[year1.full.1000$f.miss,v.int])
cm.1000.y1-cor(complete(year1.1000.imp,action = 1)[year1.full.1000$f.miss,v.int])

cor(complete(year1.1000.imp,action = 2)[year1.full.1000$f.miss,v.int])
cm.1000.y1-cor(complete(year1.1000.imp,action = 2)[year1.full.1000$f.miss,v.int])

cor(complete(year1.1000.imp,action = 3)[year1.full.1000$f.miss,v.int])
cm.1000.y1-cor(complete(year1.1000.imp,action = 3)[year1.full.1000$f.miss,v.int])

### how does this affect factor structure and factor scores? #####

# data sets for comparison

y1.sem<-year1.full[,unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))]
inc.sem<-!is.na(rowMeans(y1.sem))
y1.sem<-y1.sem[inc.sem,]

y1.sem.imp<-complete(year1.imp,action = 1)[,unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))]

y1.sem.1000<-year1.full[year1.full.1000$f.miss,unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))]
inc.1000<-!is.na(rowMeans(y1.sem.1000))
y1.sem.1000<-y1.sem.1000[inc.1000,]

y1.sem.1000.imp<-complete(year1.1000.imp,action = 1)[year1.full.1000$f.miss,unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))]
y1.sem.1000.imp<-y1.sem.1000.imp[inc.1000,]

#variables for SEM

comp.cbcl.cols<-unique(c(cbcl.ADHD.cols,cbcl.Attention.cols))

cbcl.cols.final<-comp.cbcl.cols[!comp.cbcl.cols%in%c("cbcl_q01_p",
                                                     "cbcl_q61_p")]

Attn.cols.final<-c(bpmt.Attention.cols,comp.cbcl.cols)
Attn.cols.final<-Attn.cols.final[!Attn.cols.final%in%c("cbcl_q13_p",
                                                       "cbcl_q45_p",
                                                       "cbcl_q62_p",
                                                       "cbcl_q80_p")]

bf.final <- paste(" Attn =~ ",paste(Attn.cols.final,collapse=" + "),
                  " Parent =~ ",paste(cbcl.cols.final,collapse=" + "),
                  " Teacher =~ ",paste(bpmt.Attention.cols,collapse=" + "),
                  " Attn ~~ 0*Parent 
                      Attn ~~ 0*Teacher
                      Parent ~~ 0*Teacher
                      ",sep="\n")

# compare full and imputed data set models


y1.sem.model = cfa(bf.final,data=y1.sem,  
                   std.lv=TRUE,
                   orthogonal=TRUE,
                   estimator = "WLSMV") 
y1.sem.sum<-summary(y1.sem.model,fit=TRUE,standardized=TRUE)

y1.sem.imp.model = cfa(bf.final,data=y1.sem.imp,  
                       std.lv=TRUE,
                       orthogonal=TRUE,
                       estimator = "WLSMV") 
y1.sem.imp.sum<-summary(y1.sem.imp.model,fit=TRUE,standardized=TRUE)

# teacher items contribute systematically less
plot(y1.sem.sum$PE$std.all[1:6],y1.sem.imp.sum$PE$std.all[1:6],
     xlim=c(.2,.8),ylim=c(.2,.8),
     xlab="listwise deletion",
     ylab="imputation",col="blue");abline(0,1)
points(y1.sem.sum$PE$std.all[6:16],y1.sem.imp.sum$PE$std.all[6:16],col="red")

# randomly selected 1000 missing people

y1.sem.1000.model = cfa(bf.final,data=y1.sem.1000,  
                        std.lv=TRUE,
                        orthogonal=TRUE,
                        estimator = "WLSMV") 
y1.sem.1000.sum<-summary(y1.sem.1000.model,fit=TRUE,standardized=TRUE)

y1.sem.1000.imp.model = cfa(bf.final,data=y1.sem.1000.imp,  
                            std.lv=TRUE,
                            orthogonal=TRUE,
                            estimator = "WLSMV") 
y1.sem.1000.imp.sum<-summary(y1.sem.1000.imp.model,fit=TRUE,standardized=TRUE)

# teacher items contribute even less when imputed

plot(y1.sem.1000.sum$PE$std.all[1:6],y1.sem.1000.imp.sum$PE$std.all[1:6],
     xlim=c(.2,.8),ylim=c(.2,.8),
     xlab="listwise deletion",
     ylab="imputation",col="blue");abline(0,1)
points(y1.sem.1000.sum$PE$std.all[6:16],y1.sem.1000.imp.sum$PE$std.all[6:16],col="red")


y1.sem.full1000.imp<-complete(year1.1000.imp,action = 1)[,unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))]

y1.sem.full1000.imp.model = cfa(bf.final,
                                data=y1.sem.full1000.imp,  
                            std.lv=TRUE,
                            orthogonal=TRUE,
                            estimator = "WLSMV") 

y1.sem.full1000.imp.sum<-summary(y1.sem.full1000.imp.model,fit=TRUE,standardized=TRUE)

plot(y1.sem.sum$PE$std.all[1:6],y1.sem.full1000.imp.sum$PE$std.all[1:6],
     xlim=c(.2,.8),ylim=c(.2,.8),
     xlab="listwise deletion",
     ylab="imputation",col="blue");abline(0,1)
points(y1.sem.sum$PE$std.all[6:16],y1.sem.full1000.imp.sum$PE$std.all[6:16],col="red")

# factor scores 

y1.sem$Attn<-lavPredict(y1.sem.model)[,"Attn"]
y1.sem$Parent<-lavPredict(y1.sem.model)[,"Parent"]
y1.sem$Teacher<-lavPredict(y1.sem.model)[,"Teacher"]

y1.sem.imp$Attn<-lavPredict(y1.sem.imp.model)[,"Attn"]
y1.sem.imp$Parent<-lavPredict(y1.sem.imp.model)[,"Parent"]
y1.sem.imp$Teacher<-lavPredict(y1.sem.imp.model)[,"Teacher"]

y1.sem.full1000.imp$Attn<-lavPredict(y1.sem.full1000.imp.model)[,"Attn"]
y1.sem.full1000.imp$Parent<-lavPredict(y1.sem.full1000.imp.model)[,"Parent"]
y1.sem.full1000.imp$Teacher<-lavPredict(y1.sem.full1000.imp.model)[,"Teacher"]


cor(y1.sem[,c("Attn","Parent","Teacher")])
cor(y1.sem.imp[,c("Attn","Parent","Teacher")])
cor(y1.sem.full1000.imp[,c("Attn","Parent","Teacher")])

plot(y1.sem.imp$Attn[inc.sem],y1.sem$Attn);cor(y1.sem.imp$Attn[inc.sem],y1.sem$Attn)
cor(y1.sem.imp$Attn[inc.sem],y1.sem$Attn);cor(y1.sem.imp$Attn[inc.sem],y1.sem$Attn)
#[1] 0.9996569

plot(y1.sem.full1000.imp$Attn[year1.full.1000$f.miss],y1.sem.imp$Attn[year1.full.1000$f.miss])
cor(y1.sem.full1000.imp$Attn[year1.full.1000$f.miss],y1.sem.imp$Attn[year1.full.1000$f.miss])
#[1] 0.9637267

# Figure for supplemental

jpeg(filename = "imputed_complete_comparison.jpg",width = 5,height = 5,units = "in",res=300)

plot(y1.sem.imp$Attn[year1.full.1000$f.miss],y1.sem.full1000.imp$Attn[year1.full.1000$f.miss],
     xlab="ADHD Sx., complete teacher data",
     ylab="ADHD Sx., missing teacher data",pch=16,
    col=rgb(.13,.59,.89,.5))
abline(lm(y1.sem.full1000.imp$Attn[year1.full.1000$f.miss]~y1.sem.imp$Attn[year1.full.1000$f.miss]))

dev.off()

save.image("imputation_sensitivity_ANL.RData")

#############################################################################
######### 7. Factor Score Robustness and Test-Retest Correlations ###########
#############################################################################
                                           
# make sure constrained model is robust when fit to all imputed data sets

source("ABCD_family_subsample.R")
base.SEM.full<-abcd.fam.subsample(base.full)
y1.SEM.full<-year1.full[year1.full$subjectkey%in%base.SEM.full$subjectkey,]
y2.SEM.full<-year2.full[year2.full$subjectkey%in%base.SEM.full$subjectkey,]

# what proportion had teacher data imputed?

sum(y1.SEM.full$subjectkey%in%year1.full[year1.full$t.missing,]$subjectkey)
sum(y2.SEM.full$subjectkey%in%year2.full[year2.full$t.missing,]$subjectkey)

all.c<-unique(c(cbcl.Attention.cols,cbcl.ADHD.cols,bpmt.Attention.cols))

base.SEM.full<-cbind(subjectkey=base.SEM.full$subjectkey,
                     complete(mice(data = base.SEM.full[,all.c],
                             m = 1,print = FALSE,)))
y1.SEM.full<-cbind(subjectkey=y1.SEM.full$subjectkey,
                     complete(mice(data = y1.SEM.full[,all.c],
                                   m = 1,print = FALSE,)))
y2.SEM.full<-cbind(subjectkey=y2.SEM.full$subjectkey,
                   complete(mice(data = y2.SEM.full[,all.c],
                                 m = 1,print = FALSE,)))

save(base.SEM.full,y1.SEM.full,y2.SEM.full,file="SEM_data_full_final.RData")

comp.cbcl.cols<-unique(c(cbcl.ADHD.cols,cbcl.Attention.cols))
cbcl.cols.final<-comp.cbcl.cols[!comp.cbcl.cols%in%c("cbcl_q01_p",
                                                     "cbcl_q61_p")]

Attn.cols.final<-c(bpmt.Attention.cols,comp.cbcl.cols)
Attn.cols.final<-Attn.cols.final[!Attn.cols.final%in%c("cbcl_q13_p",
                                                       "cbcl_q45_p",
                                                       "cbcl_q62_p",
                                                       "cbcl_q80_p")]

bf.final <- paste(" Attn =~ ",paste(Attn.cols.final,collapse=" + "),
                  " Parent =~ ",paste(cbcl.cols.final,collapse=" + "),
                  " Teacher =~ ",paste(bpmt.Attention.cols,collapse=" + "),
                  " Attn ~~ 0*Parent 
                      Attn ~~ 0*Teacher
                      Parent ~~ 0*Teacher
                      ",sep="\n")


base.full.model.final = cfa(bf.final,data=base.SEM.full,  
                       std.lv=TRUE,
                       orthogonal=TRUE,
                       estimator = "WLSMV") 
summary(base.full.model.final,fit=TRUE,standardized=TRUE)


y1.full.model.final = cfa(bf.final,data=y1.SEM.full,  
                     std.lv=TRUE,
                     orthogonal=TRUE,
                     estimator = "WLSMV") 
summary(y1.full.model.final,fit=TRUE,standardized=TRUE)

y2.full.model.final = cfa(bf.final,data=y2.SEM.full,  
                          std.lv=TRUE,
                          orthogonal=TRUE,
                          estimator = "WLSMV") 
summary(y2.full.model.final,fit=TRUE,standardized=TRUE)



write.csv(summary(y1.full.model.final,fit=TRUE,standardized=TRUE)$PE,
          file="loadings_constrained_full_y1.csv")

write.csv(summary(y2.full.model.final,fit=TRUE,standardized=TRUE)$PE,
          file="loadings_constrained_full_y2.csv")


# factor scores across all models and data sets 

base.SEM.full$Gen.bmod<-lavPredict(base.full.model.final)[,"Attn"]

base.SEM.full$Gen.y1mod<-lavPredict(y1.full.model.final,newdata = base.SEM.full)[,"Attn"]

base.SEM.full$Gen.y2mod<-lavPredict(y2.full.model.final,newdata = base.SEM.full)[,"Attn"]

y1.SEM.full$Gen.bmod<-lavPredict(base.full.model.final,newdata = y1.SEM.full)[,"Attn"]

y1.SEM.full$Gen.y1mod<-lavPredict(y1.full.model.final)[,"Attn"]

y1.SEM.full$Gen.y2mod<-lavPredict(y2.full.model.final,newdata = y1.SEM.full)[,"Attn"]

y1.SEM.full$Gen.y1c<-lavPredict(y1.model.final,newdata = y1.SEM.full)[,"Attn"]


y2.SEM.full$Gen.bmod<-lavPredict(base.full.model.final,newdata = y2.SEM.full)[,"Attn"]

y2.SEM.full$Gen.y1mod<-lavPredict(y1.full.model.final,newdata = y2.SEM.full)[,"Attn"]

y2.SEM.full$Gen.y2mod<-lavPredict(y2.full.model.final)[,"Attn"]

y2.SEM.full$Gen.y1c<-lavPredict(y1.model.final,newdata = y2.SEM.full)[,"Attn"]

# relations between general factors from each of the models 
# (y1, baseline) when applied to each data set

cor(base.SEM.full[,c("Gen.bmod","Gen.y1mod","Gen.y2mod")])
                                           
cor(y1.SEM.full[,c("Gen.bmod","Gen.y1mod","Gen.y2mod","Gen.y1c")])

cor(y2.SEM.full[,c("Gen.bmod","Gen.y1mod","Gen.y2mod","Gen.y1c")])

# test-retest correlation

ids.both.y1<-base.SEM.full$subjectkey[base.SEM.full$subjectkey %in% y1.SEM.full$subjectkey]
ids.both.y2<-base.SEM.full$subjectkey[base.SEM.full$subjectkey %in% y2.SEM.full$subjectkey]
fs.cols.y1<-c("subjectkey","Gen.bmod","Gen.y1mod")
fs.cols.y2<-c("subjectkey","Gen.bmod","Gen.y2mod")

base.ids.both.y1<-base.SEM.full[base.SEM.full$subjectkey%in%ids.both.y1,fs.cols.y1]
y1.ids.both<-y1.SEM.full[y1.SEM.full$subjectkey%in%ids.both.y1,fs.cols.y1]

base.ids.both.y2<-base.SEM.full[base.SEM.full$subjectkey%in%ids.both.y2,fs.cols.y2]
y2.ids.both<-y2.SEM.full[y2.SEM.full$subjectkey%in%ids.both.y2,fs.cols.y2]

colnames(base.ids.both.y1)<-c("subjectkey","Gen.bmod.base","Gen.y1mod.base")
colnames(y1.ids.both)<-c("subjectkey","Gen.bmod.y1","Gen.y1mod.y1")

colnames(base.ids.both.y2)<-c("subjectkey","Gen.bmod.base","Gen.y2mod.base")
colnames(y2.ids.both)<-c("subjectkey","Gen.bmod.y2","Gen.y2mod.y2")


test.retest.y1<-merge(base.ids.both.y1,y1.ids.both,by="subjectkey")
test.retest.y2<-merge(base.ids.both.y2,y2.ids.both,by="subjectkey")

#model fit separately to each sample

ICC(test.retest.y1[,c("Gen.bmod.base","Gen.y1mod.y1")])

ICC(test.retest.y2[,c("Gen.bmod.base","Gen.y2mod.y2")])

# outcome model fit to both samples

ICC(test.retest.y1[,c("Gen.y1mod.base","Gen.y1mod.y1")])

ICC(test.retest.y2[,c("Gen.y2mod.base","Gen.y2mod.y2")])

cor(test.retest.y1[-1])

cor(test.retest.y2[-1])

