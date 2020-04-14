


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("The directory to be pee-analyzed [1] needs to be identified, 
       as does the directory where output is to go [2],
       as does framerate[3] (1 or 4), as does nightonly[4] (y/n), 
       as does folder search term [5], as does the number of cores to use [6],
       as does the number of seconds to look back for evidence of warming (try 20) [7],
       as does the temp increase required for warming (recommend 1)[8],
       a
       as does whether to keep checking for new folders (check/nocheck) [9],
       as does the distance (in mm) determining whether pee pixels are clustered together (5 recommended) [10],
       as does the year of the trial [11],
       as does whether or not to use coordinates to clip analyses (y or n) [12],
       as does whether to make plots or not (yes or no) [13], 
       as does the # of frames forward to look for cooling [14], 
       as does the temperature change for detecting cooling [15]",
       call.=FALSE)
} 
if (length(args)!=15) {
  stop("The directory to be pee-analyzed [1] needs to be identified, 
       as does the directory where output is to go [2],
       as does framerate[3] (1 or 4), as does nightonly[4] (y/n), 
       as does folder search term [5], as does the number of cores to use [6],
       as does the number of seconds to look back for warming (try 20 for starters) [7],
       as does the temp increase required for warming (recommend 1)[8],
      
       as does whether to keep checking for new folders (check/nocheck) [9],
       as does the distance (in mm) determining whether pee pixels are clustered together (5 recommended) [10],
       as does the year of the trial [11],
       as does whether or not to use coordinates to clip analyses (y or n) [12],
       as does whether to make plots or not (yes or no) [13], 
       as does the # of frames forward to look for cooling [14], 
       as does the temperature change for detecting cooling [15]",call.=FALSE)
}
#########Pee cluster has 3 columns of detections (viewable in gifs)
############### Pee Snapshots have 3 rows of detections (viewable in pdfs)


library(zoo)
library(data.table)
library(snow)
library(doParallel)
library(biganalytics)
library(raster)
library(sp)
library(recommenderlab)
library(stringr)

#FUNCTIONS
source("DetectionFunctions.R", chdir = F)


colfunc.noblack<-colorRampPalette(rev(c("red","orange","yellow","springgreen","aquamarine","mediumblue","darkorchid4","gray85")))





# dirmain<-"/local/workdir/Ligon/Thermal/8x8thermal/C57mega"
# savedir<-"/local/workdir/Ligon/Thermal/Mar19_2020_outdata"
# framerate<-1
# nightonly<-"y"
# folderclassifier<-"_A"
# howmancoresshouldiuse<-10
# 
# windowsize<-20
# 
# temperature.increase<-1
# 
# checkfornewfolders<-"check"
# cluserdistance.mm<-5#Window size for rolling correlation
# yearoftrial<-"2018"
# trimtocoors<-"y"
# makerawpixelplots<-"no"
# 
# warmwindow<-180
# temperature.decrease<-1



dirmain<-args[1]
savedir<-args[2]
framerate<-as.numeric(as.character(args[3]))
nightonly<-as.character(args[4])
folderclassifier<-as.character(args[5])
howmancoresshouldiuse<-as.numeric(as.character(args[6]))

windowsize<-as.numeric(as.character(args[7]))
temperature.increase<-as.numeric(as.character(args[8]))
checkfornewfolders<-as.character(args[9])
cluserdistance.mm<-as.numeric(as.character(args[10]))#Window size for rolling correlation
yearoftrial<-as.character(args[11])
trimtocoors<-as.character(args[12])
makerawpixelplots<-as.character(args[13])

warmwindow<-as.numeric(as.character(args[14]))
temperature.decrease<-as.numeric(as.character(args[15]))




MB1pxcm<-9.1#pix/cm
MB2pxcm<-9.5#pix/cm
a13distance<-2.6472 #pix/com
a24distance<-2.615 #pix/com

setwd(savedir)
simplename<-strsplit(dirmain,"/")[[1]][length(strsplit(dirmain,"/")[[1]])]
zz <- file(paste(simplename,".Rout",sep=""),open="at")
#zz1<-file(paste(simplename,"1.Rout",sep=""),open="wt")
sink(zz, append=TRUE, type="message" )
sink(zz, append=TRUE, type="output" )

verybeginningtime<-Sys.time()

setwd(dirmain)

print(paste("Folders in ",dirmain," to be analyzed",sep=''))
print(paste("Summary data to be saved in ",savedir,sep=''))
print(paste("Frame rate is ",framerate," frames/sec",sep=''))
print(paste("Analyzing night only? ", nightonly,sep=''))
print(paste("Analyze folders with ", folderclassifier," in name", sep=''))
print(paste("Use",howmancoresshouldiuse,"cores"))
#print(paste("Final.Desired.Htz =",FDH,"(relevant for slope analyses"))

# print(paste("Type of correlation looking for=",typeofcorr))
#print(paste("Window size for median smoothing (k parameter in rollmedian)=",median.win.size))
#dirmain<-("C:/Users/Rusty/Amazon Drive/MICE/Thermal")
print(paste(windowsize,"= frames back to compare for *warming*"))
print(paste("Temperature increase necessary count as *warm* =",temperature.increase))
print(paste("Check for new folders with while loop =",checkfornewfolders))
print(paste("Pee pixels within",cluserdistance.mm, "mm will be clustered together as part of the same cluster"))
print(paste("Trial was run in", yearoftrial))
print(paste("Make raw pee detection pngs? ", makerawpixelplots))

print(paste("Seconds in the future to look for cool: ",warmwindow))
print(paste("Temp drop in the future to look for cooling: ",temperature.decrease))

# setwd(dir<-"/local/workdir/Ligon/Thermal/Summary_Thermal_MB1_OFT_NY2-142_072518")
# savedir<-"/local/workdir/Ligon/Thermal/PeeSummaries"
# 
# /workdir/Ligon/Thermal/Summary_Thermal_MB1_OFT_NY2-142_072518
######################################
##########################################################
##########################################################
##########################################################
##########################################################


PastFrame<-windowsize#Number of frames (i.e. rows of chunks) to look into past
FuturoFrame<-warmwindow #Number of frames (i.e. rows of chunks to look into future)
temprise<-temperature.increase#Temperature that categorizes a "warming" (from a mouse, or urination)
tempdrop<-temperature.decrease#Temperature change that consitutes a 'cool future'
minpix<-1 #Minimum number of pixels to be considered a 'cluster'







#########################################
foldersanalyzed<-0
iteration<-1
yearoftrial<-paste0(yearoftrial,"-")




coordinatefiles<-list.files("/local/workdir/Ligon/Thermal/8x8thermal/Arena_coordinates",full.names=TRUE,pattern="D1.")#creates list (2 items) with the megaframe, and its rowname file



if(checkfornewfolders=="check"){
  while(checkfornewfolders=="check"){
    
      #only on the first run do we need to define shortfolders and allfolders,
      # afterwards, it is set by evaluating the difference in folder names
      if(iteration==1){  
        allfolders<-list.dirs(path = dirmain, full.names = TRUE, recursive = TRUE)
        
        for(af in 1:length(allfolders)){
          sf<-unlist(strsplit(allfolders[af],split='/'))
          sf1<-sf[length(sf)]
          if(af==1){
            shortfolders<-sf1
          } else {
            shortfolders<-rbind(shortfolders,sf1)
          }
        }
        
        pos = grep(folderclassifier, shortfolders[,1])
        shortfolders<-shortfolders[pos]
        allfolders<-allfolders[pos]
        print("Folders to be analyzed")
        print(allfolders)
      } else {
        print("Folders to be analyzed")
        print(allfolders)
      }
        
        
      #loop through all folders existing *currently*
      for(afn in 1:length(allfolders)){
        print(paste(afn," out of ",length(allfolders),sep=''))
        setwd(dirmain)
        setwd(dir<-allfolders[afn])
        #foldyname<-shortfolders[afn]
        
        folderidentifier<-strsplit(dir,"/")[[1]][length(strsplit(dir,"/")[[1]])]#Gets name of folder being processed
        print(folderidentifier)
        trialidentifier<-strsplit(folderidentifier,"_")[[1]][1]
        camera<-strsplit(folderidentifier,"_")[[1]][2]
        day<-strsplit(folderidentifier,"_")[[1]][3]
        
        #directories<-list.dirs(dir) #Gets any sub-directories 
        BigCSVs<-list.files(dir,full.names=TRUE,pattern="Summary")# To date, I've named these megaframes "Summary_xxxxxx.csv"
        rownombres<-list.files(dir,full.names=FALSE,pattern="Batch")# To date. I've been creating separate files with the frame IDs (from the original ; csvs)
        # 
        # #read in and save full, single list of framenames
        # for(rn in 1:length(rownombres)){
        #   n1<-read.csv(rownombres[[rn]])
        #   n1n<-n1[,"Name",drop=FALSE]
        #   if(rn==1){
        #     fullframenames<-n1n
        #   } else {
        #     fullframenames<-rbind(fullframenames,n1n)
        #   }
        # }
        # write.csv(fullframenames,file=paste(folderidentifier,"fullframenames.csv",sep=''))
        ###################
        
        if(length(rownombres)==0 | length(BigCSVs)==0){
          print(paste("Folder ",folderidentifier," doesn't have the *right* csvs",sep=''))
        } else {
          
          relevantfilelength<-length(rownombres) # how many megaframes in this folder? #currently -1 to get rid of partial transfer
          #numerics<-sprintf("%05d",seq(0,relevantfilelength-1,1)) #creates 5 value (e.g. 00001) names
          nume.1<-unlist(strsplit(rownombres,".csv"))
          nume.2<-lapply(strsplit(nume.1,"_"),function(x) x[[2]])
          numerics<-unlist(nume.2)
          
          
          megabox<-strsplit(folderidentifier,"_")[[1]][length(strsplit(folderidentifier,"_")[[1]])]
          pix.cm<-ifelse(megabox=="MB1",MB1pxcm,MB2pxcm)
          
          #if the last bit of the folder identifier does NOT have "MB" in it, its length will be 0, 
          # and it is an 8x8 trial
          if(length(grep("MB",megabox))==0){
            #if there is an A13 in the folderidentifier then its length will be 1, and we assign 
            #the pix.cm of the a13distance (else, a24distance)
            pix.cm<-ifelse(length(grep("A13",folderidentifier))==1,a13distance,a24distance)
          }
          
          
          #cluserdistance.mm<-5
          divisor<-10/cluserdistance.mm #gives divisor for calculation below
          specifdist<-round((pix.cm/divisor),digits=3) #replace this 2 if you want to change the clustering distance. 
          #dividing by 2 means that pee pixels beyond 0.5cm will be clustered separately so:
          realspecificdistance<-paste0(cluserdistance.mm,"mm")
          
          folderidentifier<-paste(folderidentifier,
                                  realspecificdistance,sep=".")
          print(folderidentifier)
          
          
          ############
          #Create empty matrix for pee indexes
          #MASTER.pee.matrix<-matrix(nrow=700000,ncol=101)
          
          
          #cockatoo loop
          #loops through all megaframes (i.e. chunked data files)
          #For each loop, pulls megaframe AND adds necessary rows from previous megaframe (so that future temps can be evaluated)
          #THEN it applies SweepPeeDetector (including pee.profiler function) to look for pee profile at every pixel location (columns)
          # and every time (rows)
          
          # Next, pee coordinates are converted into a distance matrix, and those pee-pixels within
          # the threshold distance (specificDistance) are clustered. Once clustered, a sub-loop goes through
          # each cluster and finds the centroid, adds unclustered size (= n pixels) and mean temp and binds to pee.matrix
          
          
          
          startcols<-seq(1,307200,640)
          endcols<-seq(640,307200,640)
          
          
           # startcols<-seq(1,307200,6400)
           # endcols<-seq(6400,307200,6400)
          
          start.time <- Sys.time()
          
          
          
          # end.time <- Sys.time()
          # time.taken <- end.time - start.time
          # time.taken
          
          
          #startfileindex<-8
          startfileindex<-1
          
  
          
          
          startingtimeforthisfolder <- Sys.time()
      
          
            matcher<-paste0(camera,".*",trialidentifier)
            matcher1<-paste0(trialidentifier,".*",camera)
            matchma<-c(matcher,matcher1)
            

            TheseCoords<-coordinatefiles[grepl(paste(matchma,collapse="|"), coordinatefiles)]

            
          #if we are not trimming to coordinates (e.g. for volume calibration/validation), let's define relevant stuff here
          if(trimtocoors=="n"){
            TheseCoords<-list(1,2) #makes TheseCoords a 2 element list, so that the length condition check below is "TRUE"
            L.coords<-data.frame(cbind(c(0,638,638,0),c(0,0,480,480)));colnames(L.coords)<-c("X","Y")
            R.coords<-data.frame(cbind(c(639,640,640,639),c(0,0,480,480)));colnames(R.coords)<-c("X","Y")
          }  
          
            
          #Only run the analysis if the "trial" has matching arena coords for trimming, IF we are using trimtocoors='y'    
          if(length(TheseCoords)>0){
            
            if(trimtocoors=="y"){
              L.coords<-read.csv(TheseCoords[grep("L.c",TheseCoords)])
              R.coords<-read.csv(TheseCoords[grep("R.c",TheseCoords)])              
            }

            
            
            cockatoocount<-1
          for(cockatoo in startfileindex:relevantfilelength){
            print(paste("Folder: ",folderidentifier,sep=''))
            getem<-list.files(dir,full.names=TRUE,pattern=numerics[cockatoo])#creates list (2 items) with the megaframe, and its rowname file
            
  
            #print(getem)
            
            
            which.columns.the.names.in<-NCOL(as.matrix(fread(getem[1],sep=",",header=F,skip=1,nrows=3)))
            
            
            frameIDs<-c(as.matrix(fread(getem[1],sep=",",header=F,skip=1,select=which.columns.the.names.in)))#Alphabetically, the frameIDs are item [1] in the list
              #print("Taking first and only column for Summary Names.a")
            
              if(length(grep("2018",frameIDs[1]))==1){
                yearoftrial<-"2018-"
              }
              if(length(grep("2019",frameIDs[1]))==1){
                yearoftrial<-"2019-"
              }
            
            
              frameinfo<-strsplit(frameIDs,yearoftrial)
              frameinfo2<-lapply(frameinfo,function(x) strsplit(x,"_"))
              f1a<-frameinfo2[[1]]#take the first timestamp for this megacsv
              
              checkthetime<-f1a[[length(f1a)]][length(f1a[[length(f1a)]])]
              
              fixduplicatissue<-0
              
              rowstodrop<-grep("_0",str_sub(frameIDs,-2,))
           
               #if there are 1 or more rows that need to be dropped because of this duplication issue
              if(length(rowstodrop)>0){
                fixduplicatissue<-1
                
                #this will only happen if the first frame is on a duplicated time ending in another "_0"
                if(checkthetime=="0"){
                    checkthetime<-f1a[[length(f1a)]][2] #take the 2nd element of this timestamp, split on "_"
                }
              }
             
              
              
             
              
              
              # print(paste("which.columns.the.names.in =",which.columns.the.names.in ))
              # print(head(frameIDs))
              # print(paste("frameinfo=",frameinfo))
              # 
              # print(paste("length of f1a=", length(f1a)))
              # print(paste("f1a = ", f1a))
              
              print(paste("MegaCSV starting time =", checkthetime))
              
              
              
            # out <- tryCatch(frameIDs<-c(as.matrix(fread(getem[1],sep=",",header=F,skip=1,select=2))),error = function(e) { cat('In         error handler\n'); print(e); e })#function checks if fread returns error
            # 
            # #If any of the values are "error", returns TRUE (and we skip)
            # if(any(class(out) == "error")){
            #   frameIDs<-c(as.matrix(fread(getem[1],sep=",",header=F,skip=1,select=1)))#Alphabetically, the frameIDs are item [1] in the list
            #   #print("Taking first and only column for Summary Names.a")
            #   frameinfo<-strsplit(frameIDs,yearoftrial)
            #   frameinfo2<-lapply(frameinfo,function(x) strsplit(x,"_"))
            #   f1a<-frameinfo2[[1]]
            #   checkthetime<-f1a[[1]][length(f1a[[1]])]
            # } else {
            #   frameIDs<-c(as.matrix(fread(getem[1],sep=",",header=F,skip=1,select=2)))#Alphabetically, the frameIDs are item [1] in the list
            #   frameinfo<-strsplit(frameIDs,yearoftrial)
            #   frameinfo2<-lapply(frameinfo,function(x) strsplit(x,"_"))
            #   f1a<-frameinfo2[[1]]
            #   checkthetime<-f1a[[2]][c(2:3)]
            # }
            
            
            #Alphabetically, the frameIDs are item [1] in the list
            #print(frameIDs)
            
  
            
            startinghourforthismegaframe<-as.numeric(as.character(unlist(strsplit(unlist(strsplit(checkthetime[1],"_")),"-"))[1]))
            startingminforthismegaframe<-as.numeric(as.character(unlist(strsplit(unlist(strsplit(checkthetime[1],"_")),"-"))[1]))
            
            #This if/else statement is meant to skip all megacsvs (reading/processing) that start before noon (i.e. 12)
            #or start after 10pm (i.e. 22)
            if(nightonly=="y" & 
               (startinghourforthismegaframe<12 | 
                startinghourforthismegaframe>21)){
              print("skipping megaframe because the start time is not in the dark")
            } else {
              
              print("bout to fread...")
              ####
              mainframe <- as.matrix(fread(getem[2],sep=",",header=F,fill=TRUE))#Reads the megaframe in
              #mainframe <- as.matrix(fread(getem[2],sep=",",header=F,nrows=50,skip=380))#Reads the megaframe in
              #boogs<-as.matrix(fread(getem[2],sep=",",header=F,nrows=1,skip = 1999))#Reads the megaframe in
              frameIDs.keep<-frameIDs
              print("fread in...")
              
              origrow<-nrow(mainframe)
              
              if(fixduplicatissue==1){
                mainframe<-mainframe[-rowstodrop,]
                frameIDs.keep<-frameIDs[-rowstodrop]
              }
              
              #get orig size to tell how far to analyze on composite matrix
              
              
                #################################################
                #if not the last megaframe, pull rows needed from next megaframe (for future, cool analysis)
                if(cockatoo!=relevantfilelength){
                 
                  Fframes<-list.files(dir,full.names=TRUE,pattern=numerics[cockatoo+1])#Fframes lists filenames for the next frame
                 
                  #if there is a filename/dup issue, with seconds duplicated, then fixduplicate issue will be 1 (set earlier)
                  # and we'll pull in extra rows from the next frame to make sure we have enough
                  if(fixduplicatissue==1){
                    futFRAME.IDs<-c(as.matrix(fread(Fframes[1],sep=",",header=F,
                                                    skip=1,
                                                    select=which.columns.the.names.in,
                                                    nrows=(2*FuturoFrame))))#Alphabetically, the frameIDs are item [1] in the list
                    addfutureframe <- as.matrix(fread(Fframes[2],sep=",",header=F,
                                                      #skip=rowstopull,
                                                      nrows=(2*FuturoFrame)))
                    
                    
                    
                    nextrowstodrop<-grep("_0",str_sub(futFRAME.IDs,-2,))
                    addfutureframe<-addfutureframe[-nextrowstodrop,]
                    futFRAME.IDs<-futFRAME.IDs[-nextrowstodrop]
                    
                  } else {
                    
  
                    futFRAME.IDs<-c(as.matrix(fread(Fframes[1],sep=",",header=F,
                                                  skip=1,
                                                  select=which.columns.the.names.in,
                                                  nrows=FuturoFrame)))#Alphabetically, the frameIDs are item [1] in the list
                    
                    addfutureframe <- as.matrix(fread(Fframes[2],sep=",",header=F,
                                                #skip=rowstopull,
                                                nrows=FuturoFrame))
                    
                  }    
                  
                    mainframe<-rbind(mainframe,addfutureframe)#combines mainframe and addrame
                    row.names(mainframe)<-c(frameIDs.keep,futFRAME.IDs)#adds row.names 
                    frameIDs.keep<-row.names(mainframe)
                  
                  
                  

                } else {
                  #Nothing! No need to modify mainframe
                }
                
              
              

                
                #if not the first mega frame, pull rows needed from previous megaframe (for pee analyses)
                if(cockatoo!=1){
                  prevframes<-list.files(dir,full.names=TRUE,pattern=numerics[cockatoo-1])#prevframes lists filenames for the previous frame
                  
                  #number of previous frames will always be 2000
                  rowstopull<-2000-PastFrame
                  
                  
                    # out2 <- tryCatch( frameIDs2<-c(as.matrix(fread(prevframes[1],sep=",",header=F,skip=rowstopull,nrows=PastFrame))),
                    #                  error = function(e) { cat('In         error handler\n'); print(e); e })#function checks if fread returns error
                    # 
                    #If any of the values are "error", returns TRUE (and we skip)
                    # if(any(class(out2) == "error")){
                    #   frameIDs2<-c(as.matrix(fread(prevframes[1],sep=",",header=F,skip=rowstopull,select=1,nrows=PastFrame)))#only reads in necessary rows from previous megaframe
                    #   #only reads in necessary rows from previous megaframe
                    #   #print("Taking first and only column for SummaryNames.c")
                    # } else {
                    #   frameIDs2<-c(as.matrix(fread(prevframes[1],sep=",",header=F,skip=rowstopull,select=2,nrows=PastFrame)))#only reads in necessary rows from previous megaframe
                    # }
                  
                  if(fixduplicatissue==1){
                    rowstopull2<-2000-PastFrame*2
                    frameIDs2<-c(as.matrix(fread(prevframes[1],sep=",",header=F,skip=(1+rowstopull2),select=which.columns.the.names.in,nrows=PastFrame)))#Alphabetically, the frameIDs are item [1] in the list
                    addframe <- as.matrix(fread(prevframes[2],sep=",",header=F,skip=rowstopull2,nrows=PastFrame))
                    
                    prevrowstodrop<-grep("_0",str_sub(frameIDs2,-2,))
                    addframe<-addframe[-prevrowstodrop,]
                    prestartedframes<-nrows(addframe)
                    frameIDs2<-frameIDs2[-prevrowstodrop]
                    
                  } else {
                    frameIDs2<-c(as.matrix(fread(prevframes[1],sep=",",header=F,skip=(1+rowstopull),select=which.columns.the.names.in,nrows=PastFrame)))#Alphabetically, the frameIDs are item [1] in the list
                    addframe <- as.matrix(fread(prevframes[2],sep=",",header=F,skip=rowstopull,nrows=PastFrame))
                    
                  }
                  
                  mainframe<-rbind(addframe,mainframe)#combines mainframe and addrame
                  row.names(mainframe)<-c(frameIDs2,frameIDs.keep)#adds row.names 

                    
  
                } else {
                  origrow<-nrow(mainframe)
                  row.names(mainframe)<-frameIDs.keep
                }
                
  
                
                thesearethenames<-row.names(mainframe)
                
                
  
                
                ################
                #########################
                ############################################
                ##########################################################################
                #################################################################################################
                ##########################################################################
                ############################################
                ##########################
                ################
                #############
                if (file.exists("mega.bin"))
                  file.remove("mega.bin")
                
                if (file.exists("mega.desc")) 
                  file.remove("mega.desc")
                ##################
                
                library(bigmemory)
                bm.mainframe <- as.big.matrix(mainframe,type="double",
                                              backingfile="mega.bin",descriptorfile="mega.desc")
                datadesc<-dget("mega.desc")
                
                gc()
                #garbage collector to clean up memory
                #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                #$$$$$$$$$$$$$$$$$$$$$$$$$$
                #Generates 307200 rows where the values in column 1 are possible y coordinates (1-480)
                #and the values in column 2 are possible x coordinates (1-640)
                possCombo<-expand.grid(seq(1,480,1),seq(1,640,1));colnames(possCombo)<-c("row","col")
                possCombo[,1]<-as.integer(possCombo[,1])
                possCombo[,2]<-as.integer(possCombo[,2])
                ###########
               
                
                # flagme<-1
                # mainframe2<-attach.big.matrix(datadesc)
                # main.list<-list()
                # for(viva in startframe:fullon){
                #   main.list[[flagme]]<-runitall(submainframe=mainframe2[c((viva-PastFrame),viva,(viva+FuturoFrame)),],Lpolygoncoords=L.coords,Rpolygoncoords=R.coords)
                #   flagme<-flagme+1
                # }
                # 
                # 
                # 
                
                
                
                
                
                
                thismanycores<-howmancoresshouldiuse
                #RB <<- parallel::makeCluster(thismanycores, outfile="")
                RB  <<- parallel::makeCluster(thismanycores, outfile=paste(savedir,"/",folderidentifier,"RunBig.txt",sep=''))
                
                #  Ccl<-makeCluster(6)
                registerDoParallel(RB,cores=thismanycores)
                clusterEvalQ(RB, library(data.table)) ## Discard result
                print("Starting to run foreach FullON processor")
                
                startframe<-1+PastFrame
                fullon<-nrow(mainframe)-FuturoFrame
  
                #startframe<-1995
                #fullon<-startframe+1
                
                # templist<-list()
                # for(viva in startframe:fullon){
                #   print(viva)
                #   mainframe2<-attach.big.matrix(datadesc)
                #   templist[[viva]]<-runitall(submainframe=mainframe2[c((viva-PastFrame),viva,(viva+FuturoFrame)),],Lpolygoncoords=L.coords,Rpolygoncoords=R.coords)
                # 
                # }
                
                main.list<- foreach(viva=startframe:fullon, .packages=c('sp','bigmemory'),
                                    .inorder=TRUE) %dopar% {
                  #print(viva)
                  #require(bigmemory)
                  mainframe2<-attach.big.matrix(datadesc)
                  runitall(submainframe=mainframe2[c((viva-PastFrame),viva,(viva+FuturoFrame)),],Lpolygoncoords=L.coords,Rpolygoncoords=R.coords)
                }
                parallel::stopCluster(RB)
                
                names(main.list)<-thesearethenames[c(startframe:fullon)]
  
                ##############################################################################################
                ##############################################################################################
                #######################
                if(length(main.list)>1){
                  rm(mainframe)
                }
                
                #######################
                ############################################
                #
                urine<-t(dplyr::bind_rows(lapply(main.list,function(x){listpeedetails(x)})))
                colnames(urine)<-colnames(listpeedetails(main.list[[1]]))
                rownames(urine)<-names(main.list)
                
  
                library(recommenderlab)
                main.list.dropped<-lapply(main.list,function(x){dropNA(x)})
                names(main.list.dropped)<-names(main.list)
                
                
                #print(object.size(main.list), units="Mb")
                #print(object.size(main.list.dropped), units="Mb")
                

                rm(main.list)
                
                
                #if need to convert to array
                #########################
                #newarray<-array(as.numeric(unlist(Left.mega.pee.list)), dim=c(480, 640, length(Left.mega.pee.list)))
                
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                ###########################################################
                gc()
                #garbage collector to clean up memory
  
                library(rlist)
                
                
                
                if(cockatoocount<2){ #only for first megacsv
                  #print(paste("cockatoo=1, dim(mega.pee.matrix) ",dim(mega.pee.matrix),"class(mega.pee.matrix) ",class(mega.pee.matrix),sep=''))
  
                  MASTER.pee.matrix<-urine
                   
                  peeframes<-main.list.dropped
  
                  
                } else { #for all remaining megacsvs
                  MASTER.pee.matrix<-rbind(MASTER.pee.matrix,urine)
                  peeframes<-c(peeframes,main.list.dropped)
                }
                
                setwd(savedir)
                cockatoocount<-cockatoocount+1
                
                save(MASTER.pee.matrix,file=paste(folderidentifier,"MASTER.pee.matrix.Rdata",sep="."))
                #MASTER.pee.matrix is a large matrix (nrows = total number of frames analyzed X 101 columns)
                print("ZZ")
                save(peeframes,file=paste(folderidentifier,"Snapshots.Rdata",sep="_")) 
                #Snapshotsummaries are summed pee detections per timeframe of a megacsv
                #this is typically 2000 frames/ 4htz, = 500 seconds, 8 minutes and 20 seconds
                # --- could be useful for plotting dynamically changing urine allocation over time
                
                setwd(dir)
                
                if (file.exists("mainframe"))
                  rm(mainframe)
                
                if (file.exists("bm.mainframe"))
                  rm(bm.mainframe)
                
                if (file.exists("bm.peematrix"))
                  rm(bm.peematrix)
                
                gc()
                
                
                print(paste(cockatoo," megacsv down!",sep=""))
                print(Sys.time())
                
                print("Total size of workspace =")
                print(object.size(x=lapply(ls(), get)), units="Mb")
                
                print("Ten biggest items in workspace =")
                print(tail(sort(sapply(ls(),function(x){format(object.size(get(x)), units="Mb")}))))
                
              
            }#closes up if/else for megacsvs that start before noon or after 10pm
            
            
            if(cockatoo==relevantfilelength){
              folderend.time <- Sys.time()
              print(paste("started at",startingtimeforthisfolder))
              print(paste("finished at ",folderend.time))
            }
            
          }#closes cockatoo loop
          end.time <- Sys.time()
          time.taken <- end.time - start.time
          time.taken
         
          
          setwd(savedir)
          save(MASTER.pee.matrix,file=paste(folderidentifier,"MASTER.pee.matrix.Rdata",sep="."))
          save(peeframes,file=paste(folderidentifier,"Snapshots.Rdata",sep="_")) 
          
          
          if(makerawpixelplots=="yes"){
            onecore<-Sys.time()
            ###check/make directories for saving plot images
            whereplotgoes<-getwd()
            pv<-paste0(trialidentifier,"_",camera,"_",day,"_pixelviz")
            newDir<-paste(whereplotgoes,pv,sep = "/")
            if (dir.exists(newDir)){
              setwd(newDir)
            } else {
              dir.create(newDir)
              setwd(newDir)
            }
            
            
            
            ################
            cumulativepeeplotter<-function(peef,cumframe,framename,mousecolors){
              
              colfunc<-colorRampPalette(rev(c("red","orange","yellow","springgreen","aquamarine","mediumblue","darkorchid4","gray40","black")))
              
              
              ####################################
              #FRAMENAME
              called<-framename
              
              t2.d<-peef
              t2.d.r<-t2.d[rev(1:nrow(t2.d)),]
              
              t3.d<-cumframe
              t3.d.r<-t3.d[rev(1:nrow(t3.d)),]
              
  
              
              png(paste("Cum_LR_",called,".png",sep=''),bg="transparent",width=1200,height=800)#use if making transparent pngs
              image(1:ncol(t3.d.r), 1:nrow(t3.d.r), t(t3.d.r), col = colfunc(60),zlim=c(0,200),
                    axes = FALSE,asp=1, xlab='',ylab='',main="")  
              image(1:ncol(t2.d.r), 1:nrow(t2.d.r), t(t2.d.r), col =  mousecolors, 
                    breaks=c(0,11,31,101,121), axes = FALSE,asp=1, xlab='',ylab='',main="",add=TRUE)   
              text(275,20,called,cex=.6)
              dev.off()
            }
            #############
            mousecolors<-c("goldenrod3", "darkorange", "darkorchid3","blue3")
            if(camera=="A13"){mousecolors<-c("goldenrod3", "darkorange", "darkorchid3","blue3")}
            if(camera=="A24"){mousecolors<-c("goldenrod3", "darkorange", "darkgreen","darkred")}
            
            llist<-length(peeframes)
            
            for(rsrl in 1:llist){
              frnm<-names(peeframes)[rsrl]
              if(rsrl==1){
                cumFRAME<-growl<-dropNA2matrix(peeframes[[rsrl]])
                cumFRAME[cumFRAME>40]<-NA #sets pixels with values greater than 40 (mouse and clusters) to NA
                cumFRAME[cumFRAME<40]<-1 #sets pixels with values less than 40 (pee) to 1
              } else {
                bingo<-growl<-dropNA2matrix(peeframes[[rsrl]])
                bingo[bingo>40]<-NA #sets pixels with values greater than 40 (mouse and clusters) to NA
                bingo[bingo<40]<-1 #sets pixels with values less than 40 (pee) to 1
                bingo[is.na(bingo)]<-0
                
                cumFRAME[is.na(cumFRAME)]<-0
                cumFRAME<-bingo+cumFRAME
                cumFRAME[cumFRAME==0]<-NA
                
              }
              cumulativepeeplotter(peef=growl,cumframe=cumFRAME,framename=frnm,mousecolors)
            }
            
            
            plot.time<-Sys.time()-onecore
            plot.time
            }
            
            #setwd(whereplotgoes)
          
          
            }
          
          
          
          # Shutdown cluster neatly
          # if(exists("RB")){
          #   if(!is.null(RB)) {
          #     parallel::stopCluster(RB)
          #     RB <- c()
          #   }
          # }
  
        } #if batchnames missing
        
        foldersanalyzed<-foldersanalyzed+1
      } #end of afn loop
      
      
      allfolders.2<-list.dirs(path = dirmain, full.names = TRUE, recursive = TRUE)
      
      for(af in 1:length(allfolders.2)){
        sf<-unlist(strsplit(allfolders.2[af],split='/'))
        sf1<-sf[length(sf)]
        if(af==1){
          shortfolders.2<-sf1
        } else {
          shortfolders.2<-rbind(shortfolders.2,sf1)
        }
      }
      
      pos = grep(folderclassifier, shortfolders.2[,1])
      shortfolders.2<-shortfolders.2[pos]
      allfolders.2<-allfolders.2[pos]
      # print("Folders to be analyzed")
      # print(allfolders)
      
      #identifies any folders that exist in allfolders.2 that WERE NOT in allfolders
      newfolders<-setdiff(allfolders.2,allfolders)
      new.short.folders<-setdiff(shortfolders.2,shortfolders)
      
      #if the number of 'new' folders is 0 
      # AND
      # the number of folders found now (length(allfolders.2)) is the same as the number of folders analyzed
      # then checkfornewfolders will be changed to 'nocheck' and the while loop can exit, OTHERWISE, it will run again
      if(length(newfolders)==0 & length(allfolders.2)==foldersanalyzed){
        checkfornewfolders="nocheck"
      } 
      
      allfolders<-newfolders
      shortfolders<-new.short.folders
      iteration<-iteration+1
  
    }#close while loop
  
  ## reset message sink and close the file connection
  sink(type="output")
  sink(type="message")
  close(zz)

  
  
  print(paste("Folders in ",dirmain," to be analyzed",sep=''))
  print(paste("Summary data to be saved in ",savedir,sep=''))
  print(paste("Frame rate is ",framerate," frames/sec",sep=''))
  print(paste("Analyzing night only? ", nightonly,sep=''))
  #print(paste("Downsampling to 1htz? ", downsample,sep=''))
  print(paste("Analyze folders with ", folderclassifier," in name", sep=''))
  print(paste("Use",howmancoresshouldiuse,"cores"))
  #print(paste("Final.Desired.Htz =",FDH,"(relevant for slope analyses"))
  
  # print(paste("Type of correlation looking for=",typeofcorr))
  #print(paste("Window size for median smoothing (k parameter in rollmedian)=",median.win.size))
  #dirmain<-("C:/Users/Rusty/Amazon Drive/MICE/Thermal")
  print(paste(windowsize,"= frames back to compare for *warming*"))
  #print(paste("Proportion of frames that need to have negative slope in",windowsize,"frames=",propneeded))
  #print(paste("Pee detection strategy", pee.detection.strategy,"(c=correlation, n=negativeslope, cn=correlation and negativeslope)"))
  #print(paste("Correlation coefficient threshold=",CCT))
  #print(paste("Width of rolling correlation window (in frames)=",rollwidth))
  print(paste("Check for new folders with while loop =",checkfornewfolders))
  
  
  
  now<-Sys.time()
  print("*******************************************")
  print(paste("Started at", verybeginningtime,"    Finished at",now))
  print(now-verybeginningtime)
} else {
  #
}









# Shutdown cluster neatly
# if(exists("RB")){
#   if(!is.null(RB)) {
#     parallel::stopCluster(RB)
#     RB <- c()
#   }
# }







