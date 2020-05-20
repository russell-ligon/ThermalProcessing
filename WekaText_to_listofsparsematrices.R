




library(recommenderlab)

setwd(dir<-"/workdir/Ligon/Thermal/8x8thermal/April2020_CoolStuff/morningframes/wekaText")


textfiles<-list.files(,pattern=".txt")
main.list.dropped<-list()
for(t in 1:length(textfiles)){
  frame<-as.matrix(read.table(textfiles[t],header=FALSE))
  main.list.dropped[[t]]<-dropNA(frame)
}
framenames<-gsub(".txt","",textfiles)
names(main.list.dropped)<-framenames

save(main.list.dropped,file="morningframes.Rdata")

