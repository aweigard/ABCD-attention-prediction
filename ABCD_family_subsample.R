abcd.fam.subsample<-function(data,id="subjectkey",family="rel_family_id"){
  data$select<-FALSE
  fam.ids<-unique(data[,family])
  cnts<-table(as.character(data[,family]))
  singles<-names(cnts[cnts==1])
  multiples<-names(cnts[cnts>1])
  
  data[data[,family]%in%singles,"select"]<-TRUE
  
  for (f in multiples){
    tmp<-sample(c(1,rep(0,length(data[data[,family]==f,family])-1)),
                size = length(data[data[,family]==f,family]))
    c<-1
    for (s in data[data[,family]==f,id]){
      data[data[,id]==s,"select"]<-as.logical(tmp[c])
      c<-(c+1)
    }
  }
  
  out<-data[data$select==TRUE,]
  out
}