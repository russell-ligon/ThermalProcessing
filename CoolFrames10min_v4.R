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


library(data.table)
#library(vroom)
#devtools::install_github("collectivemedia/tictoc")
library(tictoc)
library(pryr)


mainwd<-"H:/MICE/c57validation"
setwd(mainwd)

infoX<-read.csv("infobox.csv",header = TRUE);infoX<-infoX[,-5]
infoX$Filepath<-as.character(infoX$Filepath)

coolframes<-"H:/MICE/c57validation/coolframes"



TrialArenaDays<-unique(infoX$Filepath)
TrialArenaDays<-TrialArenaDays[!is.na(TrialArenaDays)]


for(q in 1:length(TrialArenaDays)){
#for(q in c(41)){  
  trialdayname<-strsplit(TrialArenaDays[q],"\\",fixed=TRUE)[[1]]
  trialdayname<-trialdayname[length(trialdayname)]
  
  print(trialdayname)
  
  batchnames<-list.files(TrialArenaDays[q],pattern='Batch',full.names=TRUE)
  megas<-list.files(TrialArenaDays[q],pattern='Summary',full.names=TRUE)
  
  if(length(batchnames)>0){
    which.columns.the.names.in<-NCOL(as.matrix(fread(batchnames[1],sep=",",header=F,skip=1,nrows=3)))
    
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
    
    if(which.columns.the.names.in==2){
      fullbatchnames<-fullbatchnames[,-1]
    }
    
    
    fullbatchnames[,1]<-gsub("n_","",fullbatchnames[,1])
    fullbatchnames[,1]<-gsub("Record_","",fullbatchnames[,1])
    fullbatchnames[,1]<-gsub("new_","",fullbatchnames[,1])
    
    
    date<-sapply(strsplit(as.character(fullbatchnames[,1]),"_"),"[[",1)#takes 1st element from each stringsplit
    dateadder<-ifelse(date==unique(date)[1],0,80000)
    
    
    fullbatchnames$time<-sapply(strsplit(as.character(fullbatchnames[,1]),"_"),"[[",2)#takes 2nd element from each stringsplit
    
    #calculate daily.sec values
    fullbatchnames$dailysecond<-dailysecondfinder(fullbatchnames$time)
    fullbatchnames$dailysecond<-fullbatchnames$dailysecond+dateadder
    
    #keeps only obs between noon and 10pm
    fullbatchnames<-fullbatchnames[which(fullbatchnames$dailysecond>43199 & fullbatchnames$dailysecond<79201),] 
    
    
    
    #keep only frames collected every 10 minutes = 600 seconds
    # BUT if a match isn't found exactly 10 minutes later (saving issues), keep looking for next closest frame in the future
    compseq<-seq(fullbatchnames$dailysecond[1],fullbatchnames$dailysecond[nrow(fullbatchnames)],600)
    
    for(cs in 1:length(compseq)){
      matches<-which(fullbatchnames$dailysecond==compseq[cs])
      if(length(matches)>0){
        keepseq<-matches
      } else {
        j=0
        while(length(matches)==0){
          j<-j+1
          matches<-which(fullbatchnames$dailysecond==(compseq[cs]+j))
        }
        keepseq<-matches
      }
      if(cs==1){
        fullkeep<-keepseq
      } else {
        fullkeep<-c(fullkeep,keepseq)
      }
    }
    
    #keepbatchnames<-fullbatchnames[which(fullbatchnames$dailysecond %in% seq(fullbatchnames$dailysecond[1],fullbatchnames$dailysecond[nrow(fullbatchnames)],600)),]
    
    keepbatchnames<-fullbatchnames[fullkeep,]
    for(oneeach in 1:length(unique(keepbatchnames$dailysecond))){
      tempbatchnames<-keepbatchnames[which(keepbatchnames$dailysecond==(unique(keepbatchnames$dailysecond)[oneeach])),]
      tempbatchnames<-tempbatchnames[1,]
      if(oneeach==1){
        keeper<-tempbatchnames
      } else {
        keeper<-rbind(keeper,tempbatchnames)
      }
    }
    
    keepbatchnames<-keeper
    rm(fullbatchnames)
    
    ###check/make directories for saving plot images
    newDir<-paste(coolframes,trialdayname,sep = "/")
    if (dir.exists(newDir)){
      setwd(newDir)
    } else {
      dir.create(newDir)
      setwd(newDir)
    }
    
    # setwd(TrialArenaDays[q])
    
    for(ms in 1:length(unique(keepbatchnames$b))){
      thismega<-keepbatchnames[which(keepbatchnames$b==unique(keepbatchnames$b)[ms]),] 
      
      megafile<-megas[thismega$b[1]]
      megafile<-gsub("/","\\",megafile,fixed=TRUE)
      
      want_rows<-thismega$megarow  
      
      zNames<-thismega$Name
      
      #Can use the cmd argument in 'fread' to control the read-in file
      #this can be coupled with the 'sed' call, to ONLY read in specific rows
      # from the Summary_XXX.csvs 
      # BUT, this is somehow slower than reading in the whole thing, and then subsetting
      #https://stackoverflow.com/questions/61348204/read-specific-non-consecutive-rows-using-data-tablefread-equivalent-to-the
      # writeLines(paste0(c(1, 1+want_rows), "p"), "commands.sed")
      # 
      # comndname<-paste0(newDir,"\\","commands.sed")
      # 
      # commandline<-paste0("sed -n -f ",comndname," ",megafile)
      # 
      # look<- as.matrix(fread(cmd=commandline,sep=",",header=F,fill=TRUE));print(object.size(look),units="Mb")}
      
      #delete commands sed file
      #unlink("commands.sed")
      
      # benchmark("sed"={
      #   writeLines(paste0(c(1, 1+want_rows), "p"), "commands.sed")
      #   comndname<-paste0(newDir,"\\","commands.sed")
      #   commandline<-paste0("sed -n -f ",comndname," ",megafile)
      #   look<- as.matrix(fread(cmd=commandline,sep=",",header=F,fill=TRUE))
      #   unlink("commands.sed")
      #   
      #   },"fread"={
      #   look2<-as.matrix(fread(megafile,sep=",",header=F,fill=TRUE))
      #   #print(object.size(look2),units="Mb")
      #   look2a<-look2[want_rows,,drop=FALSE]
      #  },
      #   replications=10,
      #   columns=c("test","replications","elapsed","relative","user.self","sys.self"))
      # 
      # 
      
      ################
      # look2<-as.matrix(fread(megafile,sep=",",header=F,fill=TRUE))
      # #print(object.size(look2),units="Mb")
      # look2a<-look2[want_rows,,drop=FALSE]
      # 
      # rm(look2)
      # 
      # library(bigmemory)
      # bm.mainframe <- as.big.matrix(look2,type="double",
      #                               backingfile="mega.bin",descriptorfile="mega.desc")
      # datadesc<-dget("mega.desc")
      
      writeLines(paste0(c(1, 1+want_rows), "p"), "commands.sed")
      comndname<-paste0(newDir,"\\","commands.sed")
      commandline<-paste0("sed -n -f ",comndname," ",megafile)
      look2<- as.matrix(fread(cmd=commandline,sep=",",header=F,fill=TRUE))
      unlink("commands.sed")
      
      print(mem_used())
      
      for(ef in 1:nrow(look2)){
        frame<-look2[ef,]
        fname<-zNames[ef]
        #boogs<-as.matrix(fread(megas[frame$b],sep=",",header=F,nrows=1,skip = (frame$megarow-1)))
        
        tempmat<-thermalarrange(frame,480,640)
        write.table(tempmat,file=paste0(fname,".csv"),sep=";",row.names = FALSE,col.names =FALSE)
        # png(paste(frame$Name,".png",sep=''),bg="transparent",width=640,height=480)#use if making transparent pngs
        # vizframe(tempmat)
        # dev.off()
      }
      
      #############
      if (file.exists("mega.bin"))
        file.remove("mega.bin")
      
      if (file.exists("mega.desc")) 
        file.remove("mega.desc")
      
      if (file.exists("bm.mainframe"))
        rm(bm.mainframe)
      
      
      if (file.exists("look2"))
        rm(look2)
      
      gc()
      ##################
    }
    
  }
  
  
  
  
  
 
  
  setwd(mainwd)
  
}



nowfolders<-list.dirs(,full.names = TRUE)
nowfolders<-nowfolders[-1]

for(nf in 1:length(nowfolders)){
  thesefiles<-list.files(nowfolders[nf],full.names = TRUE)
  
  NAfile<-thesefiles[grep("NA",thesefiles)]
  if(length(NAfile)>0){
    unlink(NAfile)
    print(paste(nf,nowfolders[nf],NAfile))
  }
  
}
