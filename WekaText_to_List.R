



library(recommenderlab)
setwd("/workdir/Ligon/Thermal/8x8thermal/April2020_CoolStuff/TextOut")

TextFiles<-list.files(,pattern=".txt")
# create progress bar
pb <- txtProgressBar(min = 1, max = length(TextFiles), style = 3)
TextList<-list()
for(t in 1:length(TextFiles)){
  name<-TextFiles[t]
  name<-gsub(".txt","",name)
  x<-read.table(TextFiles[t])
  x<-as.matrix(x)
  a<-dropNA(x)
  b<-name
  
  TextList[[t]]<-list(b,a)
  
  # update progress bar
  setTxtProgressBar(pb, t)
    
}

fullnames<-names(TextList)
fullnames<-gsub(".txt","",fullnames)
names(TextList)<-fullnames

setwd("/workdir/Ligon/Thermal/8x8thermal/April2020_CoolStuff")
save(TextList,file="TextList.Rdata")
