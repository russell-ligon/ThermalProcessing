



library(recommenderlab)


TextFiles<-list.files(,pattern=".txt")

TextList<-list()
for(t in 1:length(TextFiles)){
  name<-TextFiles[t]
  name<-gsub(".txt","",name)
  x<-read.table(TextFiles[t])
  x<-as.matrix(x)
  TextList[[t]]<-dropNA(x)
}

fullnames<-names(TextList)
fullnames<-gsub(".txt","",fullnames)
names(TextList)<-fullnames

setwd("/workdir/Ligon/Thermal/8x8thermal/April2020_CoolStuff")
save(TextList,file="TextList.Rdata")