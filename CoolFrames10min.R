#function to return daily second from military time (xx-xx-xx)
dailysecondfinder<-function(militarytime){
  
  t.hour<-unlist(lapply(militarytime,function(x) as.numeric(as.character(strsplit(x,"-")[[1]][1]))))
  t.min <-unlist(lapply(militarytime,function(x) as.numeric(as.character(strsplit(x,"-")[[1]][2]))))
  t.sec <-unlist(lapply(militarytime,function(x) as.numeric(as.character(strsplit(x,"-")[[1]][3]))))
  
  #calculate daily.sec values
  dailysecond<-t.sec+(t.min*60)+(t.hour*60*60)
  return(dailysecond)
}

#Plot matrix for quick vizualization
library(colorRamps)
colfunc<-colorRampPalette(rev(c("red","orange","yellow","springgreen","aquamarine","mediumblue","darkorchid4","gray40","black")))

vizframe<-function(matrix){
  matrix[is.na(matrix)]<-0
  t3.d<-matrix
  t3.d.r<-t3.d[rev(1:nrow(t3.d)),]
  image(1:ncol(t3.d.r), 1:nrow(t3.d.r), t(t3.d.r), 
        #breaks=c(0,1,20,50,90,300,500),
        col = colfunc(60),
        zlim=c(60,80),
        xlim=c(1,640),ylim=c(1,480),
        axes = FALSE,asp=1, xlab='',ylab='',main="") 
}



#function to rearange vector into matrix, so we can apply over whole dataframe
thermalarrange<-function(arow,nr,nc){
  matrix(arow,nrow=nr,ncol=nc,byrow=TRUE)
}


mainwd<-"D:/MICE/c57validation"
setwd(mainwd)

infoX<-read.csv("infobox.csv",header = TRUE);infoX<-infoX[,-5]
infoX$Filepath<-as.character(infoX$Filepath)

coolframes<-"D:/MICE/c57validation/coolframes"



TrialArenaDays<-unique(infoX$Filepath)

for(q in 1:length(TrialArenaDays)){
  
  trialdayname<-strsplit(TrialArenaDays[q],"\\",fixed=TRUE)[[1]]
  trialdayname<-trialdayname[length(trialdayname)]
  
  batchnames<-list.files(TrialArenaDays[q],pattern='Batch',full.names=TRUE)
  megas<-list.files(TrialArenaDays[q],pattern='Summary',full.names=TRUE)
  
  for(b in 1:length(batchnames)){
    batchnombres<-read.csv(batchnames[b])
    batchnombres<-cbind(batchnombres,b)
    batchnombres$megarow<-c(1:nrow(batchnombres))
    if(b==1){
      fullbatchnames<-batchnombres
    } else {
      fullbatchnames<-rbind(fullbatchnames,batchnombres)
    }
  }
  fullbatchnames[,1]<-gsub("n_","",fullbatchnames[,1])
  fullbatchnames[,1]<-gsub("Record_","",fullbatchnames[,1])
  
  fullbatchnames$time<-sapply(strsplit(as.character(fullbatchnames[,1]),"_"),"[[",2)#takes 2nd element from each stringsplit
  
  #calculate daily.sec values
  fullbatchnames$dailysecond<-dailysecondfinder(fullbatchnames$time)
  
  #keeps only obs between noon and 10pm
  fullbatchnames<-fullbatchnames[which(fullbatchnames$dailysecond>43199 & fullbatchnames$dailysecond<79201),] 
  
  #keep only frames collected every 10 minutes = 600 seconds
  keepbatchnames<-fullbatchnames[which(fullbatchnames$dailysecond %in% seq(fullbatchnames$dailysecond[1],fullbatchnames$dailysecond[nrow(fullbatchnames)],600)),]
  
  
  ###check/make directories for saving plot images
  newDir<-paste(coolframes,trialdayname,sep = "/")
  if (dir.exists(newDir)){
    setwd(newDir)
  } else {
    dir.create(newDir)
    setwd(newDir)
  }
  
  for(ef in 1:nrow(keepbatchnames)){
    frame<-keepbatchnames[ef,]
    
    boogs<-as.matrix(fread(megas[frame$b],sep=",",header=F,nrows=1,skip = (frame$megarow-1)))
    
    tempmat<-thermalarrange(boogs,480,640)
    write.table(tempmat,file=paste0(frame$Name,".csv"),sep=";",row.names = FALSE,col.names =FALSE)
    # png(paste(frame$Name,".png",sep=''),bg="transparent",width=640,height=480)#use if making transparent pngs
    # vizframe(tempmat)
    # dev.off()
  }
  
  setwd(mainwd)
  
}