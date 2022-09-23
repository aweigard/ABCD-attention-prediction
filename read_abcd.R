# function for loading standard ABCD .txt files

read.abcd<-function(file,
                    path="" # enter path for your local data storage directory 
){
  dat<-read.delim(paste0(path,file),na.strings=c("","NA"))
  labs<-as.character(dat[1,])
  names(labs)<-colnames(dat)
  dat<-data.frame(dat[-1,])
  row.names(dat)<-1:length(dat[,1])
  dat<-as.data.frame(apply(dat, 2, function(x) gsub("^$|^ $", NA, x)))
  comment(dat)<-labs
  dat
}