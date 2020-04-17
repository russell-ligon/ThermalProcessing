
library("lubridate")
library("trajr")
library("sp")

setwd("C:/Users/Rusty/Amazon Drive/MICE")
source("MouseFunctions.R", chdir = F)


IsSide<-function(p.a.x,p.a.y, p.b.x,p.b.y, p.c.x, p.c.y){
  #(p.a.x,p.a.y), (p.b.x,p.b.y), (p.c.x,p.c.y) each are a points with an x and y coordinate
  #a = line point 1; b = line point 2; c = point to check against
  #https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line
  return ((p.b.x - p.a.x)*(p.c.y - p.a.y) - (p.b.y - p.a.y)*(p.c.x - p.a.x)) > 0;
}

getpositionalvalues<-function(coordinateset){
  
  minx<-min(coordinateset[,1])
  maxx<-max(coordinateset[,1])
  miny<-min(coordinateset[,2])
  maxy<-max(coordinateset[,2])
  
  ###############################
  botleft<-coordinateset[which(coordinateset[,1]<=(minx+30)),]
  botleft<-botleft[which(botleft[,2]==min(botleft[,2])),]
  
  botright<-coordinateset[which(coordinateset[,1]>=(maxx-30)),]
  botright<-botright[which(botright[,2]==min(botright[,2])),]
  
  topleft<-coordinateset[which(coordinateset[,1]<=(minx+30)),]
  topleft<-topleft[which(topleft[,2]==max(topleft[,2])),]
  
  topright<-coordinateset[which(coordinateset[,1]>=(maxx-30)),]
  topright<-topright[which(topright[,2]==max(topright[,2])),]
  
  if(class(botleft)=='matrix'){if(nrow(botleft)>1){botleft<-botleft[1,]}}
  if(class(botright)=='matrix'){if(nrow(botright)>1){botright<-botright[1,]}}
  if(class(topleft)=='matrix'){if(nrow(topleft)>1){topleft<-topleft[1,]}}
  if(class(topright)=='matrix'){if(nrow(topright)>1){topright<-topright[1,]}}
  
  vertz<-rbind(botleft,botright,topleft,topright)
  return(vertz)
}

setwd("C:/Users/Rusty/Amazon Drive/MICE/Thermal/IRdata")
load("IRfileslist.Rdata")

setwd("C:/Users/Rusty/Amazon Drive/MICE/Thermal/Next/8x8PeeData/PeeDetected")

MetaData<-read.csv("1_8x8_SQ_master.csv")
MetaData$quad<-paste0("a",MetaData$X8x8.quad)
MetaData$quad<-as.factor(MetaData$quad)
MetaData<-MetaData[-which(colnames(MetaData)=="notes")]
MetaData$X8x8.quad[which(MetaData$X8x8.quad=="not there")]<-NA
MetaData$X8x8.quad<-factor(paste0("a",MetaData$X8x8.quad))
MetaData$X8x8_quad<-MetaData$X8x8.quad
MetaData$id<-paste(MetaData$strain,MetaData$id,sep='.')


vizdir<-"C:/Users/Rusty/Amazon Drive/MICE/Thermal/Next/8x8PeeData/Visualize"


dir<-"C:/Users/Rusty/Amazon Drive/MICE/Thermal/Next/8x8PeeData/PeeDetected"
setwd(dir)

ROIs<-read.csv("ROI_coordinates.csv")
ROIs$arena<-tolower(ROIs$arena)
ROIs<-ROIs[,c(1:13)]
ROIs<-na.omit(ROIs)

HelperData<-read.csv("HelperData_AG.csv")
HelperData$ID.left<-ifelse(HelperData$ID.left=="A1","a1",
                           ifelse(HelperData$ID.left=="A2","a2",
                                  ifelse(HelperData$ID.left=="A3","a3",
                                         ifelse(HelperData$ID.left=="A4","a4",HelperData$ID.left))))
HelperData$ID.right<-ifelse(HelperData$ID.right=="A1","a1",
                           ifelse(HelperData$ID.right=="A2","a2",
                                  ifelse(HelperData$ID.right=="A3","a3",
                                         ifelse(HelperData$ID.right=="A4","a4",HelperData$ID.right))))
HelperData<-HelperData[,c("Trial", "Arena", "ID.left", "ID.right", "TopScreen", 
                          "botdivide.x", "botdivide.y", "topdivide.x", "topdivide.y", "toplefthoriz.x", 
                          "toplefthoriz.y", "botlefthoriz.x", "botlefthoriz.y", "toprighthoriz.x", 
                          "toprighthoriz.y", "botrighthoriz.x", "botrighthoriz.y")]

folders<-list.dirs(dir,full.names=TRUE,recursive = FALSE)# 
folders<-folders[-1]
AllTrialIDs<-unlist(lapply(folders, function(x) {strsplit(x,"/")[[1]][length(strsplit(x,"/")[[1]])]}))



#COOrdinates for polygons defining 'hard' boundaries for pee detection
coordinatefiles<-list.files("C:/Users/Rusty/Amazon Drive/MICE/Thermal/Arena_coordinates",full.names=TRUE,pattern="D1.")#creates list (2 items) with the megaframe, and its rowname file



#makeplots, not yet re-integrated with current pipeline
#makeplots<-"yes"
makeplots<-"no"

a13distance<-2.6472 #pix/com
a24distance<-2.615 #pix/com

#centeradjust="no"
centeradjust="yes"

#wingplots="stacked"
wingplots="climbinglines"

#Do you want to drop pee clusters that are only a single pixel in size?
#---------------------
#dropsinglepixelclusters=TRUE
dropsinglepixelclusters=FALSE

#Do you want to implement a time/distance threshold? This threshold drops pee detections
# that are too far from from any location
# that the mouse has recently visited
#--------------------------------------
#implementdistancethreshold=TRUE
implementdistancethreshold=FALSE

trialflag<-1
roimismatchflag<-1
#alltrials will be a list with n elements, corresponding to the # of folders analyzed (each corresponding to 1 trial)
# each element of alltrials will be a list, with m elements, corresponding to the # of days analyzed for that trial
alltrials<-list()
for(trial in 1:length(folders)){
#for(trial in 10:10){  
  print(paste("Analyzing trial ",trial," out of ",length(folders),sep=''))
  ClusterData<-list.files((folders[[trial]]),full.names=TRUE,pattern="(MASTER).*\\.Rdata$")# To date, I've named these megaframes "Summary_xxxxxx.csv"
  
  #Pull Trial-specific helper & meta info
  TrialID<-AllTrialIDs[trial]
  print(paste0("(",TrialID,")"))
  HelperTrial<-HelperData[which(HelperData$Trial==TrialID),]
  MetaTrial<-MetaData[which(MetaData$trial==TrialID),];MetaTrial<-MetaTrial[,!grepl("OFT",colnames(MetaTrial))]#Drop OFT data
  #MetaTrial$masschange<-MetaTrial$X8x8.post.mass.g.-MetaTrial$X8x8.mass.g
  
  #Pull Trial-specific tracks
  trackfiles<-list.files(folders[[trial]],full.names = TRUE,pattern="(fixed).*\\.csv$")
  ROIfiles<-trackfiles[grepl("region",trackfiles)]
  trackfiles<-setdiff(trackfiles,ROIfiles)
  
  #IR social data
  trialfiles<-grep(TrialID,names(IRfileslist))
  D1.ir.index<-names(IRfileslist)[trialfiles][grep("D1",names(IRfileslist)[trialfiles])]
  D2.ir.index<-names(IRfileslist)[trialfiles][grep("D2",names(IRfileslist)[trialfiles])]
  
  D1.IRdata<-IRfileslist[[which(names(IRfileslist)==D1.ir.index)]]
  D2.IRdata<-IRfileslist[[which(names(IRfileslist)==D2.ir.index)]]
  
  #Get names of Arena-specific, trial-specific files to read in
  A13s<-ClusterData[grepl("A13",ClusterData)]
  A24s<-ClusterData[grepl("A24",ClusterData)]
  Alist<-list(A13s,A24s)
  
 
  #Get trial, camera, and side-specific coordinates
  TheseCoords<-coordinatefiles[grepl(TrialID, coordinatefiles)]
  CoordsA13<-TheseCoords[grepl("A13", TheseCoords)]
  CoordsA24<-TheseCoords[grepl("A24", TheseCoords)]
  
  skipA13<-skipA24<-"no"
  
  if(length(CoordsA13)>0){
    A13.left<-read.csv(CoordsA13[grep("L.c",CoordsA13)])
    A13.right<-read.csv(CoordsA13[grep("R.c",CoordsA13)])  
  } else {
    skipA13<-"yes"
  }

  if(length(CoordsA24)>0){
    A24.left<-read.csv(CoordsA24[grep("L.c",CoordsA24)])
    A24.right<-read.csv(CoordsA24[grep("R.c",CoordsA24)])
  } else {
    skipA24<-"yes"
  }
  
  dayflag<-1
  daylist<-list()
    for(day in 1:2){ #this length should equal the number of days
      print(paste("Day", day))
      #txtProgressBar(min=1,max=length(Just.Frames),initial=teddy,width = 80,style=3)
      
      
      #Set Day
      p.day<-ifelse(day==1,"D1","D2")
      # p.day<-ifelse(day==1,"D1",
      #           ifelse(day==2,"D2",
      #               ifelse(day==3,"D3","D4")))
      # 
      if(p.day=="D1"){IRdata<-D1.IRdata}
      if(p.day=="D2"){IRdata<-D2.IRdata}
      
      #get actual, trial-specific roi data
      trialspecific.ROI<-ROIs[which(ROIs$trial==TrialID & ROIs$day==p.day),]
      trialspecific.ROI$arena<-factor(trialspecific.ROI$arena)
      
      #########
      rowissueA13<-0
      rowissueA24<-0
      cam13daymissing<-0
      cam24daymissing<-0
      fixafter<-0
      #load cluster pee info
      Camera.13<-Alist[[1]]
      if(any(grepl(p.day,Camera.13))){
        c13<-load(Camera.13[grepl(p.day,Camera.13)])
        Camera.13.day<-get(c13);if(exists("thisday")){rm(thisday)};if(exists("MASTER.pee.matrix")){rm(MASTER.pee.matrix)}
        if(nrow(Camera.13.day)>700000){Camera.13.day<-Camera.13.day[c(700001:nrow(Camera.13.day)),];
        row.names(Camera.13.day)<-seq(43200,(43200-1+(nrow(Camera.13.day))));nameissueA13<-"yes";print("WONKY MASTER PEE A13")}
        Camera.13.day<-data.frame(Camera.13.day)
      } else {
        print(paste0(TrialID,": Data for Camera13, ",p.day," missing"))
        cam13daymissing<-1
      }
      
      
      Camera.24<-Alist[[2]]
      if(any(grepl(p.day,Camera.24))){
        c24<-load(Camera.24[grepl(p.day,Camera.24)])
        Camera.24.day<-get(c24);if(exists("thisday")){rm(thisday)};if(exists("MASTER.pee.matrix")){rm(MASTER.pee.matrix)}
        if(nrow(Camera.24.day)>700000){Camera.24.day<-Camera.24.day[c(700001:nrow(Camera.24.day)),];
        row.names(Camera.24.day)<-seq(43200,(43200-1+nrow(Camera.24.day)));nameissueA24<-"yes";print("WONKY MASTER PEE A24")}
        Camera.24.day<-data.frame(Camera.24.day)
      } else {
        print(paste0(TrialID,": Data for Camera24, ",p.day," missing"))
        cam24daymissing<-1
      }
      
      
      if(skipA13=="yes"){cam13daymissing<-1}
      if(skipA24=="yes"){cam24daymissing<-1}#codes data from this camera as missing, based on earlier 'skipA24' if coordinates missing
      
      
      if(cam13daymissing==1 & cam24daymissing==1){
        print(paste0(TrialID,": Data for BOTH cameras, ",p.day," missing...skipping"))
      } else { 
      
        if(cam13daymissing==0){  
          ##### deal with rownames
          if(length(grep("2018",row.names(Camera.13.day)[1]))==1){
           A13.frames<-data.frame(unlist(lapply(strsplit(row.names(Camera.13.day),"ecord_2018-"),function(x){x[[2]]})))
           colnames(A13.frames)<-"Name";
          }
          if(length(grep("2019",row.names(Camera.13.day)[1]))==1){
            A13.frames<-data.frame(unlist(lapply(strsplit(row.names(Camera.13.day),"ecord_2019-"),function(x){x[[2]]})))
            colnames(A13.frames)<-"Name";
          }
          if(length(grep("43200",row.names(Camera.13.day)[1]))==1){
            A13.frames<-data.frame(unlist(lapply(strsplit(row.names(Camera.13.day),"ecord"),function(x){x[[1]]})))
            colnames(A13.frames)<-"Name";
            rowissueA13<-1
          }
        }
      #
        if(cam24daymissing==0){  
          if(length(grep("2018",row.names(Camera.24.day)[1]))==1){
            A24.frames<-data.frame(unlist(lapply(strsplit(row.names(Camera.24.day),"ecord_2018-"),function(x){x[[2]]})))
            colnames(A24.frames)<-"Name";
          }
          if(length(grep("2019",row.names(Camera.24.day)[1]))==1){
            A24.frames<-data.frame(unlist(lapply(strsplit(row.names(Camera.24.day),"ecord_2019-"),function(x){x[[2]]})))
            colnames(A24.frames)<-"Name";
          }
          if(length(grep("43200",row.names(Camera.24.day)[1]))==1){
            A24.frames<-data.frame(unlist(lapply(strsplit(row.names(Camera.24.day),"ecord"),function(x){x[[1]]})))
            colnames(A24.frames)<-"Name";
            rowissueA24<-1
          }
        }
      #pull track,roi info
      alltracks<-trackfiles[grepl(p.day,trackfiles)]
      allrois<-ROIfiles[grepl(p.day,ROIfiles)]
      
      #Match up track info with full frame name info
      if(cam13daymissing==0){  
        
          fix.a1.track<-0
        
          a1.track<-read.csv(alltracks[grepl("a1-",alltracks,ignore.case = TRUE)]);a1.track$info<-as.character(a1.track$info);a1.track$info[3]<-as.character("a1")
          #if track file doesn't start at first frame of video (position ==0), add in those filler rows
          if(a1.track$position[1]!=0){
            startpos.a1<-a1.track$position[1]
            temp.a1<-cbind(seq(0,(startpos.a1-1)),NA,NA);colnames(temp.a1)<-c("position","x0","y0")
            nw.a1<-rbind(temp.a1,a1.track[,c("position","x0","y0")])
            nw.inf<-c(a1.track$info,rep("",nrow(nw.a1)-nrow(a1.track)))
            a1.track<-cbind(nw.a1,nw.inf);colnames(a1.track)<-c("position","x0","y0","info")
            a1.track$info<-as.character(a1.track$info)            
            fix.a1.track<-1
            
          }
          
          a1.track$info[1]<-gsub("2019","2018",a1.track$info[1])#deal with 2019 issue
            a1.start<-strsplit(a1.track$info[1],"ecord_2018-")[[1]][2]
            if(nchar(a1.start)>14){a1.start<-substr(a1.start,0,14)}
            #deals with issue if rownames are WONKY and replaced with daily second value
            if(rowissueA13==1){
              timestring<-as.numeric(strsplit(strsplit(a1.start,"_")[[1]][2],"-")[[1]])
              a1.start.time<-timestring[1]*60*60+timestring[2]*60+timestring[3]
              a1.track$position<-A13.frames$Name[c(which(A13.frames$Name==a1.start.time):(-1+nrow(a1.track)+which(A13.frames$Name==a1.start.time)))]
            } else {
              namematchlength<-length(which(A13.frames$Name==a1.start)) #check if we can find the time in the track data IN the peedata
              if(namematchlength>0){
                  a1.aligner.coords<-a1.track$position
                  a1.track$position<-A13.frames$Name[c(which(A13.frames$Name==a1.start):(-1+nrow(a1.track)+which(A13.frames$Name==a1.start)))]
                  
                  #found it
              } else {
                #it's missing, so we need to change our input track/roi data to match the pee data
                timestring.track<-as.numeric(strsplit(strsplit(a1.start,"_")[[1]][2],"-")[[1]])
                timestring.pee<-as.numeric(strsplit(strsplit(as.character(A13.frames$Name[1]),"_")[[1]][2],"-")[[1]])
                
                a1.track.time<-timestring.track[1]*60*60+timestring.track[2]*60+timestring.track[3]
                a1.pee.time<-timestring.pee[1]*60*60+timestring.pee[2]*60+timestring.pee[3]
                numtodump<-a1.pee.time-a1.track.time
                
                    a1.info.orig<-a1.track$info[c(1:3)]
                    a1.track<-a1.track[(numtodump+1):nrow(a1.track),]
                    a1.track$info[c(1:3)]<-a1.info.orig
                    a1.aligner.coords<-a1.track$position
                    a1.track$position<-A13.frames$Name[c(1:nrow(a1.track))]
                    
               }
            }
            
            if(length(allrois[grepl("a1-",allrois,ignore.case = TRUE)])==1){
              a1.roi<-read.csv(allrois[grepl("a1-",allrois,ignore.case = TRUE)]);colnames(a1.roi)[2]<-"a1.roi"
              
              if(fix.a1.track==1 | a1.aligner.coords[1]!=0){
                a1.track.2<-cbind(a1.track,a1.aligner.coords)
                a1.hope<-merge(a1.roi,a1.track.2,by.x="position",by.y="a1.aligner.coords",all.y=TRUE)
                a1.roi<-a1.hope
              }
              
              if(nrow(a1.track)==nrow(a1.roi)){
                a1.track$roi<-a1.roi$a1.roi
              } else {
                  print(paste("A1 rois and track info don't match for ",TrialID,sep=''))
                }
              
            } else {
              a1.track$roi<-NA
            }
            
            fix.a3.track<-0
            
            a3.track<-read.csv(alltracks[grepl("a3-",alltracks,ignore.case = TRUE)]);a3.track$info<-as.character(a3.track$info);a3.track$info[3]<-as.character("a3")
            #if track file doesn't start at first frame of video (position ==0), add in those filler rows
            if(a3.track$position[1]!=0){
              startpos.a3<-a3.track$position[1]
              temp.a3<-cbind(seq(0,(startpos.a3-1)),NA,NA);colnames(temp.a3)<-c("position","x0","y0")
              nw.a3<-rbind(temp.a3,a3.track[,c("position","x0","y0")])
              nw.inf<-c(a3.track$info,rep("",nrow(nw.a3)-nrow(a3.track)))
              a3.track<-cbind(nw.a3,nw.inf);colnames(a3.track)<-c("position","x0","y0","info")
              a3.track$info<-as.character(a3.track$info)            
              fix.a3.track<-1
              
            }
            
            a3.track$info[1]<-gsub("2019","2018",a3.track$info[1])#deal with 2019 issue
            a3.start<-strsplit(a3.track$info[1],"ecord_2018-")[[1]][2]
            if(nchar(a3.start)>14){a3.start<-substr(a3.start,0,14)}
            #deals with issue if rownames are WONKY and replaced with daily second value
            if(rowissueA13==1){
              timestring<-as.numeric(strsplit(strsplit(a3.start,"_")[[1]][2],"-")[[1]])
              a3.start.time<-timestring[1]*60*60+timestring[2]*60+timestring[3]
              a3.track$position<-A13.frames$Name[c(which(A13.frames$Name==a3.start.time):(-1+nrow(a3.track)+which(A13.frames$Name==a3.start.time)))]
            } else {
              namematchlength<-length(which(A13.frames$Name==a3.start)) #check if we can find the time in the track data IN the peedata
              if(namematchlength>0){
                a3.aligner.coords<-a3.track$position
                a3.track$position<-A13.frames$Name[c(which(A13.frames$Name==a3.start):(-1+nrow(a3.track)+which(A13.frames$Name==a3.start)))]
                #found it
              } else {
                #it's missing, so we need to change our input track/roi data to match the pee data
                timestring.track<-as.numeric(strsplit(strsplit(a3.start,"_")[[1]][2],"-")[[1]])
                timestring.pee<-as.numeric(strsplit(strsplit(as.character(A13.frames$Name[1]),"_")[[1]][2],"-")[[1]])
                
                a3.track.time<-timestring.track[1]*60*60+timestring.track[2]*60+timestring.track[3]
                a3.pee.time<-timestring.pee[1]*60*60+timestring.pee[2]*60+timestring.pee[3]
                numtodump<-a3.pee.time-a3.track.time
                
                a3.info.orig<-a3.track$info[c(1:3)]
                a3.track<-a3.track[(numtodump+1):nrow(a3.track),]
                a3.track$info[c(1:3)]<-a3.info.orig
                a3.aligner.coords<-a3.track$position
                a3.track$position<-A13.frames$Name[c(1:nrow(a3.track))]
                
              }
            }
            
            if(length(allrois[grepl("a3-",allrois,ignore.case = TRUE)])==1){
              a3.roi<-read.csv(allrois[grepl("a3-",allrois,ignore.case = TRUE)]);colnames(a3.roi)[2]<-"a3.roi"
              
              if(fix.a3.track==1 | a3.aligner.coords[1]!=0){
                a3.track.2<-cbind(a3.track,a3.aligner.coords)
                a3.hope<-merge(a3.roi,a3.track.2,by.x="position",by.y="a3.aligner.coords",all.y=TRUE)
                a3.roi<-a3.hope
              }
              
              if(nrow(a3.track)==nrow(a3.roi)){
                a3.track$roi<-a3.roi$a3.roi
              } else {
                print(paste("A3 rois and track info don't match for ",TrialID,sep=''))
              }
              
            } else {
              a3.track$roi<-NA
            }
            
      }
       ####
       if(cam24daymissing==0){
         
         ###########################
         fix.a2.track<-0
         
         a2.track<-read.csv(alltracks[grepl("a2-",alltracks,ignore.case = TRUE)]);a2.track$info<-as.character(a2.track$info);a2.track$info[3]<-as.character("a2")
         #if track file doesn't start at first frame of video (position ==0), add in those filler rows
         if(a2.track$position[1]!=0){
           startpos.a2<-a2.track$position[1]
           temp.a2<-cbind(seq(0,(startpos.a2-1)),NA,NA);colnames(temp.a2)<-c("position","x0","y0")
           nw.a2<-rbind(temp.a2,a2.track[,c("position","x0","y0")])
           nw.inf<-c(a2.track$info,rep("",nrow(nw.a2)-nrow(a2.track)))
           a2.track<-cbind(nw.a2,nw.inf);colnames(a2.track)<-c("position","x0","y0","info")
           a2.track$info<-as.character(a2.track$info)            
           fix.a2.track<-1
           
         }
         
         a2.track$info[1]<-gsub("2019","2018",a2.track$info[1])#deal with 2019 issue
         a2.start<-strsplit(a2.track$info[1],"ecord_2018-")[[1]][2]
         if(nchar(a2.start)>14){a2.start<-substr(a2.start,0,14)}
         #deals with issue if rownames are WONKY and replaced with daily second value
         if(rowissueA24==1){
           timestring<-as.numeric(strsplit(strsplit(a2.start,"_")[[1]][2],"-")[[1]])
           a2.start.time<-timestring[1]*60*60+timestring[2]*60+timestring[3]
           a2.track$position<-A24.frames$Name[c(which(A24.frames$Name==a2.start.time):(-1+nrow(a2.track)+which(A24.frames$Name==a2.start.time)))]
         } else {
           namematchlength<-length(which(A24.frames$Name==a2.start)) #check if we can find the time in the track data IN the peedata
           if(namematchlength>0){
             a2.aligner.coords<-a2.track$position
             a2.track$position<-A24.frames$Name[c(which(A24.frames$Name==a2.start):(-1+nrow(a2.track)+which(A24.frames$Name==a2.start)))]
             #found it
           } else {
             timestring.track<-as.numeric(strsplit(strsplit(a2.start,"_")[[1]][2],"-")[[1]])
             timestring.pee<-as.numeric(strsplit(strsplit(as.character(A24.frames$Name[1]),"_")[[1]][2],"-")[[1]])
             
             a2.track.time<-timestring.track[1]*60*60+timestring.track[2]*60+timestring.track[3]
             a2.pee.time<-timestring.pee[1]*60*60+timestring.pee[2]*60+timestring.pee[3]
             numtodump<-a2.pee.time-a2.track.time
             
             a2.info.orig<-a2.track$info[c(1:3)]
             a2.track<-a2.track[(numtodump+1):nrow(a2.track),]
             a2.track$info[c(1:3)]<-a2.info.orig
             a2.aligner.coords<-a2.track$position
             a2.track$position<-A24.frames$Name[c(1:nrow(a2.track))]
             
           }
         }
         
         if(length(allrois[grepl("a2-",allrois,ignore.case = TRUE)])==1){
           a2.roi<-read.csv(allrois[grepl("a2-",allrois,ignore.case = TRUE)]);colnames(a2.roi)[2]<-"a2.roi"
           
           if(fix.a2.track==1 | a2.aligner.coords[1]!=0){
             a2.track.2<-cbind(a2.track,a2.aligner.coords)
             a2.hope<-merge(a2.roi,a2.track.2,by.x="position",by.y="a2.aligner.coords",all.y=TRUE)
             a2.roi<-a2.hope
           }
           
           if(nrow(a2.track)==nrow(a2.roi)){
             a2.track$roi<-a2.roi$a2.roi
           } else {
             print(paste("A2 rois and track info don't match for ",TrialID,sep=''))
           }
           
         } else {
           a2.track$roi<-NA
         }
         
         #################################
         fix.a4.track<-0
         
         a4.track<-read.csv(alltracks[grepl("a4-",alltracks,ignore.case = TRUE)]);a4.track$info<-as.character(a4.track$info);a4.track$info[3]<-as.character("a4")
         #if track file doesn't start at first frame of video (position ==0), add in those filler rows
         if(a4.track$position[1]!=0){
           startpos.a4<-a4.track$position[1]
           temp.a4<-cbind(seq(0,(startpos.a4-1)),NA,NA);colnames(temp.a4)<-c("position","x0","y0")
           nw.a4<-rbind(temp.a4,a4.track[,c("position","x0","y0")])
           nw.inf<-c(a4.track$info,rep("",nrow(nw.a4)-nrow(a4.track)))
           a4.track<-cbind(nw.a4,nw.inf);colnames(a4.track)<-c("position","x0","y0","info")
           a4.track$info<-as.character(a4.track$info)            
           fix.a4.track<-1
           
         }
         
         a4.track$info[1]<-gsub("2019","2018",a4.track$info[1])#deal with 2019 issue
         a4.start<-strsplit(a4.track$info[1],"ecord_2018-")[[1]][2]
         if(nchar(a4.start)>14){a4.start<-substr(a4.start,0,14)}
         #deals with issue if rownames are WONKY and replaced with daily second value
         if(rowissueA24==1){
           timestring<-as.numeric(strsplit(strsplit(a4.start,"_")[[1]][2],"-")[[1]])
           a4.start.time<-timestring[1]*60*60+timestring[2]*60+timestring[3]
           a4.track$position<-A13.frames$Name[c(which(A13.frames$Name==a4.start.time):(-1+nrow(a4.track)+which(A13.frames$Name==a4.start.time)))]
         } else {
           namematchlength<-length(which(A13.frames$Name==a4.start)) #check if we can find the time in the track data IN the peedata
           if(namematchlength>0){
             a4.aligner.coords<-a4.track$position
             a4.track$position<-A13.frames$Name[c(which(A13.frames$Name==a4.start):(-1+nrow(a4.track)+which(A13.frames$Name==a4.start)))]
             #found it
           } else {
             timestring.track<-as.numeric(strsplit(strsplit(a4.start,"_")[[1]][2],"-")[[1]])
             timestring.pee<-as.numeric(strsplit(strsplit(as.character(A24.frames$Name[1]),"_")[[1]][2],"-")[[1]])
             
             a4.track.time<-timestring.track[1]*60*60+timestring.track[2]*60+timestring.track[3]
             a4.pee.time<-timestring.pee[1]*60*60+timestring.pee[2]*60+timestring.pee[3]
             numtodump<-a4.pee.time-a4.track.time
             
             a4.info.orig<-a4.track$info[c(1:3)]
             a4.track<-a4.track[(numtodump+1):nrow(a4.track),]
             a4.track$info[c(1:3)]<-a4.info.orig
             a4.aligner.coords<-a4.track$position
             a4.track$position<-A24.frames$Name[c(1:nrow(a4.track))]             
           }
         }
         
         if(length(allrois[grepl("a4-",allrois,ignore.case = TRUE)])==1){
           a4.roi<-read.csv(allrois[grepl("a4-",allrois,ignore.case = TRUE)]);colnames(a4.roi)[2]<-"a4.roi"
           
           if(fix.a4.track==1 | a4.aligner.coords[1]!=0){
             a4.track.2<-cbind(a4.track,a4.aligner.coords)
             a4.hope<-merge(a4.roi,a4.track.2,by.x="position",by.y="a4.aligner.coords",all.y=TRUE)
             a4.roi<-a4.hope
           }
           
           if(nrow(a4.track)==nrow(a4.roi)){
             a4.track$roi<-a4.roi$a4.roi
           } else {
             print(paste("A4 rois and track info don't match for ",TrialID,sep=''))
           }
           
         } else {
           a4.track$roi<-NA
         }
         

       }
      
      
      if(cam13daymissing==0){
        a1.track$roi<-as.character(a1.track$roi);
        a3.track$roi<-as.character(a3.track$roi);
        a1.track$roi[(a1.track$roi=='')]<-'generalmid'
        a3.track$roi[(a3.track$roi=='')]<-'generalmid'
        a1.track$roi<-as.factor(a1.track$roi);
        a3.track$roi<-as.factor(a3.track$roi);
      }   
        
      if(cam24daymissing==0){ 
        a2.track$roi<-as.character(a2.track$roi);
        a4.track$roi<-as.character(a4.track$roi);  
        a2.track$roi[(a2.track$roi=='')]<-'generalmid'
        a4.track$roi[(a4.track$roi=='')]<-'generalmid'
        a2.track$roi<-as.factor(a2.track$roi);
        a4.track$roi<-as.factor(a4.track$roi);
      }


      
      if(cam13daymissing==0 & cam24daymissing==0){tracklist<-list(a1.track,a2.track,a3.track,a4.track)}
      if(cam13daymissing==1 & cam24daymissing==0){tracklist<-list(a2.track,a4.track)}
      if(cam13daymissing==0 & cam24daymissing==1){tracklist<-list(a1.track,a3.track)}
      
      if(cam13daymissing==0 & cam24daymissing==0){
        if(rowissueA24==1 & rowissueA13==1){
          
          for(TL in 1:length(tracklist)){
  
            focus<-tracklist[[TL]]
            trackerguide<-focus$info[c(1:3)]
            quad.ID<-trackerguide[3]
            focus<-focus[,-4]
            colnames(focus)<-c("dailysecond",(paste0(quad.ID,colnames(focus)[c(2:4)])))
            if(TL==1){
              combo_locations<-focus
            } else {
              combo_locations<-merge(combo_locations,focus,by="dailysecond",all=TRUE)
            }
            
          }
          
          
          nameissueA13<-"yes"#if issue detected, this is set to "yes" below
          nameissueA24<-"yes"
        } else {
          #pull track info to attach to pee info
          #creates dataframe called "combo_locations" with quad-specific location data
          for(TL in 1:length(tracklist)){
            focus<-tracklist[[TL]]
            trackerguide<-focus$info[c(1:3)]
            focus<-focus[,-4]
            focus<-focus[,c(2,3,1,4)]
            reftime<-strsplit(trackerguide[1],"ecord_2018-")[[1]][2]#reference time
            Htz<-1#trackerguide[2]
            quad.ID<-trackerguide[3]
            
            
            
            #if reftime has a subsecond,need to get it 'aligned' to full second
            
              focus$t.stamp<-apply(focus[,3,drop=FALSE],1,function(x) {strsplit(x,"_")[[1]][2]})
              focus$t.hour<-apply(focus[,5,drop=FALSE],1,function(x) as.numeric(as.character(strsplit(x,"-")[[1]][1])))
              focus$t.min<-apply(focus[,5,drop=FALSE],1,function(x) as.numeric(as.character(strsplit(x,"-")[[1]][2])))
              focus$t.sec<-apply(focus[,5,drop=FALSE],1,function(x) as.numeric(as.character(strsplit(x,"-")[[1]][3])))
              
              #Finds first instance of a given time (not including the sub-measurements/second)
              focus$unique.time<-!duplicated(focus$t.stamp)      
              ####
              #subsets mainframe to analyze only at 1htz intervals
              streamlined.focus<-focus[which(focus$unique.time),]
              
              #calculate daily.sec values
              streamlined.focus$dailysecond<-streamlined.focus$t.sec+(streamlined.focus$t.min*60)+(streamlined.focus$t.hour*60*60)
    
              #keep only coordinates, roi, and dailysecond
              streamlined.focus<-streamlined.focus[,c(1,2,10,4)]
              colnames(streamlined.focus)[c(1,2,4)]<-paste(quad.ID,colnames(streamlined.focus)[c(1,2,4)],sep='')
            
            
            if(TL==1){
              combo_locations<-streamlined.focus
            } else {
              combo_locations<-merge(combo_locations,streamlined.focus,by="dailysecond",all=TRUE)
            }
          }
        
          nameissueA13<-"no"#if issue detected, this is set to "yes" below
          nameissueA24<-"no"
        }
      } else {
        if(rowissueA24==1 | rowissueA13==1){
          
          for(TL in 1:length(tracklist)){
            
            focus<-tracklist[[TL]]
            trackerguide<-focus$info[c(1:3)]
            quad.ID<-trackerguide[3]
            focus<-focus[,-4]
            colnames(focus)<-c("dailysecond",(paste0(quad.ID,colnames(focus)[c(2:4)])))
            if(TL==1){
              combo_locations<-focus
            } else {
              combo_locations<-merge(combo_locations,focus,by="dailysecond",all=TRUE)
            }
            
          }
          
          
          nameissueA13<-"yes"#if issue detected, this is set to "yes" below
          nameissueA24<-"yes"
        } else {
          #pull track info to attach to pee info
          #creates dataframe called "combo_locations" with quad-specific location data
          for(TL in 1:length(tracklist)){
            focus<-tracklist[[TL]]
            trackerguide<-focus$info[c(1:3)]
            focus<-focus[,-4]
            focus<-focus[,c(2,3,1,4)]
            reftime<-strsplit(trackerguide[1],"ecord_2018-")[[1]][2]#reference time
            Htz<-1#trackerguide[2]
            quad.ID<-trackerguide[3]
            
            
            
            #if reftime has a subsecond,need to get it 'aligned' to full second
            
            focus$t.stamp<-apply(focus[,3,drop=FALSE],1,function(x) {strsplit(x,"_")[[1]][2]})
            focus$t.hour<-apply(focus[,5,drop=FALSE],1,function(x) as.numeric(as.character(strsplit(x,"-")[[1]][1])))
            focus$t.min<-apply(focus[,5,drop=FALSE],1,function(x) as.numeric(as.character(strsplit(x,"-")[[1]][2])))
            focus$t.sec<-apply(focus[,5,drop=FALSE],1,function(x) as.numeric(as.character(strsplit(x,"-")[[1]][3])))
            
            #Finds first instance of a given time (not including the sub-measurements/second)
            focus$unique.time<-!duplicated(focus$t.stamp)      
            ####
            #subsets mainframe to analyze only at 1htz intervals
            streamlined.focus<-focus[which(focus$unique.time),]
            
            #calculate daily.sec values
            streamlined.focus$dailysecond<-streamlined.focus$t.sec+(streamlined.focus$t.min*60)+(streamlined.focus$t.hour*60*60)
            
            #keep only coordinates, roi, and dailysecond
            streamlined.focus<-streamlined.focus[,c(1,2,10,4)]
            colnames(streamlined.focus)[c(1,2,4)]<-paste(quad.ID,colnames(streamlined.focus)[c(1,2,4)],sep='')
            
            
            if(TL==1){
              combo_locations<-streamlined.focus
            } else {
              combo_locations<-merge(combo_locations,streamlined.focus,by="dailysecond",all=TRUE)
            }
          }
          
          nameissueA13<-"no"#if issue detected, this is set to "yes" below
          nameissueA24<-"no"
        }
      }
      
      if(cam13daymissing==0){  
        if(nameissueA13=="no"){  
          #get frame info, camera13
          if(length(grep("2018",rownames(Camera.13.day)[1]))==1){
            inf1<-(lapply(strsplit(rownames(Camera.13.day),"ecord_2018-"),function(x) {x[2]}))
          } else {
            inf1<-(lapply(strsplit(rownames(Camera.13.day),"ecord_2019-"),function(x) {x[2]}))
          }
          
          inf2<-(lapply(inf1,function(x){strsplit(x,"_")}))
          inf3<-(lapply(inf2,function(x) {x[[1]][2]}))
          hour<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][1]}))))
          minutes<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][2]}))))
          second<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][3]}))))
          dailyseconds<-second+(minutes*60)+(hour*60*60)
          frame.info<-paste(p.day,unlist(lapply(inf2,function(x) {x[[1]][2]})),sep=".")
          Camera.13.day$frameinfo<-frame.info
          Camera.13.day$dailysecond<-dailyseconds
        } else {
          Camera.13.day$frameinfo<-"abc"
          Camera.13.day$dailysecond<-as.numeric(as.character(row.names(Camera.13.day)))
        }
      } else {
        fixafter<-1
      }
      
      if(cam24daymissing==0){
        if(nameissueA24=="no"){    
          #get frame info, camera24
          if(length(grep("2018",rownames(Camera.24.day)[1]))==1){
            inf1<-(lapply(strsplit(rownames(Camera.24.day),"ecord_2018-"),function(x) {x[2]}))
          } else {
            inf1<-(lapply(strsplit(rownames(Camera.24.day),"ecord_2019-"),function(x) {x[2]}))
          }
    
          inf2<-(lapply(inf1,function(x){strsplit(x,"_")}))
          inf3<-(lapply(inf2,function(x) {x[[1]][2]}))
          hour<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][1]}))))
          minutes<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][2]}))))
          second<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][3]}))))
          dailyseconds<-second+(minutes*60)+(hour*60*60)
          frame.info<-paste(p.day,unlist(lapply(inf2,function(x) {x[[1]][2]})),sep=".")
          Camera.24.day$frameinfo<-frame.info
          Camera.24.day$dailysecond<-dailyseconds
        } else {
          Camera.24.day$frameinfo<-"abc"
          Camera.24.day$dailysecond<-as.numeric(as.character(row.names(Camera.24.day)))
        }        
      } else {
        Camera.24.day<-Camera.13.day
        Camera.24.day[,c(1:101)]<-NA
      }

      if(fixafter==1){
        Camera.13.day<-Camera.24.day
        Camera.13.day[,c(1:101)]<-NA
      }
      
      
      orig.singlecam.names<-colnames(Camera.24.day)
      
      ###################
      #unique colnames by camera
      colnames(Camera.13.day)<-paste(colnames(Camera.13.day),".A13",sep='');colnames(Camera.13.day)[103]<-"dailysecond"
      colnames(Camera.24.day)<-paste(colnames(Camera.24.day),".A24",sep='');colnames(Camera.24.day)[103]<-"dailysecond"
      
      #Full time (noon to 10pm)
      fulltimerunner<-as.data.frame(seq((12*60*60),(22*60*60),1));colnames(fulltimerunner)[1]<-"dailysecond"
      
      #this merge temporally aligns the data from both cameras
      wookA13<-merge(fulltimerunner,Camera.13.day,by="dailysecond",all.x=TRUE)
      wookA24<-merge(fulltimerunner,Camera.24.day,by="dailysecond",all.x=TRUE)
      
      positpeeA13<-merge(wookA13,combo_locations,by="dailysecond",all.x=TRUE)
      positpeeA24<-merge(wookA24,combo_locations,by="dailysecond",all.x=TRUE)
      
      if(cam13daymissing==1){
        colnames(positpeeA13)[c(104:109)]<-c("a1x0", "a1y0", "a1roi", "a3x0", "a3y0", "a3roi")
        #positpeeA13[,match(c("a1x0", "a1y0", "a1roi", "a3x0", "a3y0", "a3roi"),colnames(positpeeA13))]<-NA
      }
      if(cam24daymissing==1){
        colnames(positpeeA24)[c(104:109)]<-c("a2x0", "a2y0", "a2roi", "a4x0", "a4y0", "a4roi")
        #positpeeA24[,match(c("a2x0", "a2y0", "a2roi", "a4x0", "a4y0", "a4roi"),colnames(positpeeA24))]<-NA
      }

      ####################################
      #positionANDpee<-merge(wook2,combo_locations,by="dailysecond",all.x=TRUE)
      #beepr::beep(4)
      ####################################
      # c("dailysecond", "a1x0", "a1y0", "a1roi", "a3x0", "a3y0", "a3roi",
      #   "frameinfo.A13", "Nclusts.A13", "Clustsize.A13", "mTemp.A13",
      #   "x1.A13", "y1.A13", "Clustsize.1.A13", "mTemp.1.A13", "x1.1.A13",
      #   "y1.1.A13", "Clustsize.2.A13", "mTemp.2.A13", "x1.2.A13", "y1.2.A13",
      #   "Clustsize.3.A13", "mTemp.3.A13", "x1.3.A13", "y1.3.A13", "Clustsize.4.A13",
      #   "mTemp.4.A13", "x1.4.A13", "y1.4.A13", "Clustsize.5.A13", "mTemp.5.A13",
      #   "x1.5.A13", "y1.5.A13", "Clustsize.6.A13", "mTemp.6.A13", "x1.6.A13",
      #   "y1.6.A13", "Clustsize.7.A13", "mTemp.7.A13", "x1.7.A13", "y1.7.A13",
      #   "Clustsize.8.A13", "mTemp.8.A13", "x1.8.A13", "y1.8.A13", "Clustsize.9.A13",
      #   "mTemp.9.A13", "x1.9.A13", "y1.9.A13", "Clustsize.10.A13", "mTemp.10.A13",
      #   "x1.10.A13", "y1.10.A13", "Clustsize.11.A13", "mTemp.11.A13",
      #   "x1.11.A13", "y1.11.A13", "Clustsize.12.A13", "mTemp.12.A13",
      #   "x1.12.A13", "y1.12.A13", "Clustsize.13.A13", "mTemp.13.A13",
      #   "x1.13.A13", "y1.13.A13", "Clustsize.14.A13", "mTemp.14.A13",
      #   "x1.14.A13", "y1.14.A13", "Clustsize.15.A13", "mTemp.15.A13",
      #   "x1.15.A13", "y1.15.A13", "Clustsize.16.A13", "mTemp.16.A13",
      #   "x1.16.A13", "y1.16.A13", "Clustsize.17.A13", "mTemp.17.A13",
      #   "x1.17.A13", "y1.17.A13", "Clustsize.18.A13", "mTemp.18.A13",
      #   "x1.18.A13", "y1.18.A13", "Clustsize.19.A13", "mTemp.19.A13",
      #   "x1.19.A13", "y1.19.A13", "Clustsize.20.A13", "mTemp.20.A13",
      #   "x1.20.A13", "y1.20.A13", "Clustsize.21.A13", "mTemp.21.A13",
      #   "x1.21.A13", "y1.21.A13", "Clustsize.22.A13", "mTemp.22.A13",
      #   "x1.22.A13", "y1.22.A13", "Clustsize.23.A13", "mTemp.23.A13",
      #   "x1.23.A13", "y1.23.A13", "Clustsize.24.A13", "mTemp.24.A13",
      #   "x1.24.A13", "y1.24.A13")
      alignedA13<-positpeeA13[,c("dailysecond", "a1x0", "a1y0", "a1roi", "a3x0", "a3y0", "a3roi",
                                 "frameinfo.A13", "Nclusts.A13", "Clustsize.A13", "mTemp.A13",
                                 "x1.A13", "y1.A13", "Clustsize.1.A13", "mTemp.1.A13", "x1.1.A13",
                                 "y1.1.A13", "Clustsize.2.A13", "mTemp.2.A13", "x1.2.A13", "y1.2.A13",
                                 "Clustsize.3.A13", "mTemp.3.A13", "x1.3.A13", "y1.3.A13", "Clustsize.4.A13",
                                 "mTemp.4.A13", "x1.4.A13", "y1.4.A13", "Clustsize.5.A13", "mTemp.5.A13",
                                 "x1.5.A13", "y1.5.A13", "Clustsize.6.A13", "mTemp.6.A13", "x1.6.A13",
                                 "y1.6.A13", "Clustsize.7.A13", "mTemp.7.A13", "x1.7.A13", "y1.7.A13",
                                 "Clustsize.8.A13", "mTemp.8.A13", "x1.8.A13", "y1.8.A13", "Clustsize.9.A13",
                                 "mTemp.9.A13", "x1.9.A13", "y1.9.A13", "Clustsize.10.A13", "mTemp.10.A13",
                                 "x1.10.A13", "y1.10.A13", "Clustsize.11.A13", "mTemp.11.A13",
                                 "x1.11.A13", "y1.11.A13", "Clustsize.12.A13", "mTemp.12.A13",
                                 "x1.12.A13", "y1.12.A13", "Clustsize.13.A13", "mTemp.13.A13",
                                 "x1.13.A13", "y1.13.A13", "Clustsize.14.A13", "mTemp.14.A13",
                                 "x1.14.A13", "y1.14.A13", "Clustsize.15.A13", "mTemp.15.A13",
                                 "x1.15.A13", "y1.15.A13", "Clustsize.16.A13", "mTemp.16.A13",
                                 "x1.16.A13", "y1.16.A13", "Clustsize.17.A13", "mTemp.17.A13",
                                 "x1.17.A13", "y1.17.A13", "Clustsize.18.A13", "mTemp.18.A13",
                                 "x1.18.A13", "y1.18.A13", "Clustsize.19.A13", "mTemp.19.A13",
                                 "x1.19.A13", "y1.19.A13", "Clustsize.20.A13", "mTemp.20.A13",
                                 "x1.20.A13", "y1.20.A13", "Clustsize.21.A13", "mTemp.21.A13",
                                 "x1.21.A13", "y1.21.A13", "Clustsize.22.A13", "mTemp.22.A13",
                                 "x1.22.A13", "y1.22.A13", "Clustsize.23.A13", "mTemp.23.A13",
                                 "x1.23.A13", "y1.23.A13", "Clustsize.24.A13", "mTemp.24.A13",
                                 "x1.24.A13", "y1.24.A13")]
      #alignedA13<-positionANDpee[,c(1,206:208,212:214,103,2:102)]
      colnames(alignedA13)[c(9:109)]<-orig.singlecam.names[c(1:101)]
      
      
      # c("dailysecond", "a2x0", "a2y0", "a2roi", "a4x0", "a4y0", "a4roi",
      #   "frameinfo.A24", "Nclusts.A24", "Clustsize.A24", "mTemp.A24",
      #   "x1.A24", "y1.A24", "Clustsize.1.A24", "mTemp.1.A24", "x1.1.A24",
      #   "y1.1.A24", "Clustsize.2.A24", "mTemp.2.A24", "x1.2.A24", "y1.2.A24",
      #   "Clustsize.3.A24", "mTemp.3.A24", "x1.3.A24", "y1.3.A24", "Clustsize.4.A24",
      #   "mTemp.4.A24", "x1.4.A24", "y1.4.A24", "Clustsize.5.A24", "mTemp.5.A24",
      #   "x1.5.A24", "y1.5.A24", "Clustsize.6.A24", "mTemp.6.A24", "x1.6.A24",
      #   "y1.6.A24", "Clustsize.7.A24", "mTemp.7.A24", "x1.7.A24", "y1.7.A24",
      #   "Clustsize.8.A24", "mTemp.8.A24", "x1.8.A24", "y1.8.A24", "Clustsize.9.A24",
      #   "mTemp.9.A24", "x1.9.A24", "y1.9.A24", "Clustsize.10.A24", "mTemp.10.A24",
      #   "x1.10.A24", "y1.10.A24", "Clustsize.11.A24", "mTemp.11.A24",
      #   "x1.11.A24", "y1.11.A24", "Clustsize.12.A24", "mTemp.12.A24",
      #   "x1.12.A24", "y1.12.A24", "Clustsize.13.A24", "mTemp.13.A24",
      #   "x1.13.A24", "y1.13.A24", "Clustsize.14.A24", "mTemp.14.A24",
      #   "x1.14.A24", "y1.14.A24", "Clustsize.15.A24", "mTemp.15.A24",
      #   "x1.15.A24", "y1.15.A24", "Clustsize.16.A24", "mTemp.16.A24",
      #   "x1.16.A24", "y1.16.A24", "Clustsize.17.A24", "mTemp.17.A24",
      #   "x1.17.A24", "y1.17.A24", "Clustsize.18.A24", "mTemp.18.A24",
      #   "x1.18.A24", "y1.18.A24", "Clustsize.19.A24", "mTemp.19.A24",
      #   "x1.19.A24", "y1.19.A24", "Clustsize.20.A24", "mTemp.20.A24",
      #   "x1.20.A24", "y1.20.A24", "Clustsize.21.A24", "mTemp.21.A24",
      #   "x1.21.A24", "y1.21.A24", "Clustsize.22.A24", "mTemp.22.A24",
      #   "x1.22.A24", "y1.22.A24", "Clustsize.23.A24", "mTemp.23.A24",
      #   "x1.23.A24", "y1.23.A24", "Clustsize.24.A24", "mTemp.24.A24",
      #   "x1.24.A24", "y1.24.A24")
      alignedA24<-positpeeA24[,c("dailysecond", "a2x0", "a2y0", "a2roi", "a4x0", "a4y0", "a4roi",
                                 "frameinfo.A24", "Nclusts.A24", "Clustsize.A24", "mTemp.A24",
                                 "x1.A24", "y1.A24", "Clustsize.1.A24", "mTemp.1.A24", "x1.1.A24",
                                 "y1.1.A24", "Clustsize.2.A24", "mTemp.2.A24", "x1.2.A24", "y1.2.A24",
                                 "Clustsize.3.A24", "mTemp.3.A24", "x1.3.A24", "y1.3.A24", "Clustsize.4.A24",
                                 "mTemp.4.A24", "x1.4.A24", "y1.4.A24", "Clustsize.5.A24", "mTemp.5.A24",
                                 "x1.5.A24", "y1.5.A24", "Clustsize.6.A24", "mTemp.6.A24", "x1.6.A24",
                                 "y1.6.A24", "Clustsize.7.A24", "mTemp.7.A24", "x1.7.A24", "y1.7.A24",
                                 "Clustsize.8.A24", "mTemp.8.A24", "x1.8.A24", "y1.8.A24", "Clustsize.9.A24",
                                 "mTemp.9.A24", "x1.9.A24", "y1.9.A24", "Clustsize.10.A24", "mTemp.10.A24",
                                 "x1.10.A24", "y1.10.A24", "Clustsize.11.A24", "mTemp.11.A24",
                                 "x1.11.A24", "y1.11.A24", "Clustsize.12.A24", "mTemp.12.A24",
                                 "x1.12.A24", "y1.12.A24", "Clustsize.13.A24", "mTemp.13.A24",
                                 "x1.13.A24", "y1.13.A24", "Clustsize.14.A24", "mTemp.14.A24",
                                 "x1.14.A24", "y1.14.A24", "Clustsize.15.A24", "mTemp.15.A24",
                                 "x1.15.A24", "y1.15.A24", "Clustsize.16.A24", "mTemp.16.A24",
                                 "x1.16.A24", "y1.16.A24", "Clustsize.17.A24", "mTemp.17.A24",
                                 "x1.17.A24", "y1.17.A24", "Clustsize.18.A24", "mTemp.18.A24",
                                 "x1.18.A24", "y1.18.A24", "Clustsize.19.A24", "mTemp.19.A24",
                                 "x1.19.A24", "y1.19.A24", "Clustsize.20.A24", "mTemp.20.A24",
                                 "x1.20.A24", "y1.20.A24", "Clustsize.21.A24", "mTemp.21.A24",
                                 "x1.21.A24", "y1.21.A24", "Clustsize.22.A24", "mTemp.22.A24",
                                 "x1.22.A24", "y1.22.A24", "Clustsize.23.A24", "mTemp.23.A24",
                                 "x1.23.A24", "y1.23.A24", "Clustsize.24.A24", "mTemp.24.A24",
                                 "x1.24.A24", "y1.24.A24")]
      #alignedA24<-positionANDpee[,c(1,209:211,215:217,205,104:204)]
      colnames(alignedA24)[c(9:109)]<-orig.singlecam.names[c(1:101)]
      
      ######################################################
      colnames(alignedA13)[seq(12,109,4)]<-"true.y"
      colnames(alignedA13)[seq(13,109,4)]<-"true.x"
      colnames(alignedA24)[seq(12,109,4)]<-"true.y"
      colnames(alignedA24)[seq(13,109,4)]<-"true.x"
      ##############
      trim.alignedA13<-alignedA13[,-9]#drop Nclusts
      trim.alignedA24<-alignedA24[,-9]#drop Nclusts
      
      #convert to long format, with a column identifying frame
      #OLD VERSION, KEEPS MULTIPLE ENTRIES PER TIME INTERVAL, REGARDLESS OF EMPTINESS
            # bseq<-seq(9,108,4)
            # for(b in 1:25){
            #   print(paste(b," of 25",sep=''))
            #   H<-bseq[b]+3
            #   if(b==1){
            #     longA13<- trim.alignedA13[,c(1:8,bseq[b]:H)]
            #     longA24<- trim.alignedA24[,c(1:8,bseq[b]:H)]
            #     
            #   } else {
            #     ta13<-trim.alignedA13[,c(1:8,bseq[b]:H)]
            #     colnames(ta13)<-colnames(longA13)
            #     longA13<-rbind(longA13,ta13)
            #     
            #     ta24<-trim.alignedA24[,c(1:8,bseq[b]:H)]
            #     colnames(ta24)<-colnames(longA24)
            #     longA24<-rbind(longA24,ta24)
            #   }
            # }
      ####
      bseq<-seq(9,108,4)
      for(b in 1:25){
        #print(paste(b," of 25",sep=''))
        H<-bseq[b]+3
        if(b==1){
          longA13<- trim.alignedA13[,c(1:8,bseq[b]:H)]
          longA13$nNas<-apply(longA13,1,function(x) sum(is.na(x)))
          longA24<- trim.alignedA24[,c(1:8,bseq[b]:H)]
          longA24$nNas<-apply(longA24,1,function(x) sum(is.na(x)))
          
        } else {
          ta13<-trim.alignedA13[,c(1:8,bseq[b]:H)]
          #
          ta13$nNas<-apply(ta13,1,function(x) sum(is.na(x)))
          ta13<-ta13[which(ta13$nNas<4),]
          colnames(ta13)<-colnames(longA13)

          #ta13<-ta13[,-13]
          ##
          longA13<-rbind(longA13,ta13)
          
          ta24<-trim.alignedA24[,c(1:8,bseq[b]:H)]
          #
          ta24$nNas<-apply(ta24,1,function(x) sum(is.na(x)))
          ta24<-ta24[which(ta24$nNas<4),]
          colnames(ta24)<-colnames(longA24)
          #ta24<-ta24[,-13]
          ##
          longA24<-rbind(longA24,ta24)
        }
      }
      ###
      longA13<-longA13[order(longA13$dailysecond),]
      longA24<-longA24[order(longA24$dailysecond),]
      

      
      full.long.A13<-longA13 #keep as 'whole'
      full.long.A24<-longA24
      
      if(dropsinglepixelclusters==TRUE){
        print("dropping single pixel clusters")
        singlepixelsA13<-which(longA13$Clustsize==1)
        singlepixelsA24<-which(longA24$Clustsize==1)
        #Drop single pixel clusters
        for(a13px in singlepixelsA13){
          tofix<-longA13[a13px,,drop=FALSE]
          tofix[,c(9:13)]<-c(NA,NA,NA,NA,10)
          longA13[a13px,]<-tofix
        }
        for(a24px in singlepixelsA24){
          tofix<-longA24[a24px,,drop=FALSE]
          tofix[,c(9:13)]<-c(NA,NA,NA,NA,10)
          longA24[a24px,]<-tofix
        }
      }

      
      
      #droping by NAs in these named columns no longer works??
     # longA13.peeonly<-na.omit(longA13,cols=c("Clustsize","mTemp","true.y","true.x"))
     # longA24.peeonly<-na.omit(longA24,cols=c("Clustsize","mTemp","true.y","true.x"))
      

      
      # longA13.peeonly<-longA13
      # longA13.peeonly$nNas<-apply(longA13,1,function(x) sum(is.na(x)))
      # longA13.peeonly<-longA13.peeonly[which(longA13.peeonly$nNas<4),]
      # 
      # longA24.peeonly<-longA24
      # longA24.peeonly$nNas<-apply(longA24,1,function(x) sum(is.na(x)))
      # longA24.peeonly<-longA24.peeonly[which(longA24.peeonly$nNas<4),]
      # 
      # 
      # longA13<-longA13.peeonly
      # longA24<-longA24.peeonly
      
      # alignedA13[,seq(4,100,4)]<-(-1)*alignedA13[,seq(4,100,4)]#A13, make y-axis negative (0,-480)
      # alignedA24[,seq(4,100,4)]<-abs(alignedA24[,seq(4,100,4)]-480)#A24, flip y-axis, turn 1 to 479, and turn 480 to 0
      
      A13.reference<-HelperTrial[which(HelperTrial$Arena=="13"),]
      A13.reference$sex.left<-MetaTrial$sex[which(A13.reference$ID.left==MetaTrial$X8x8_quad)]
      A13.reference$sex.right<-MetaTrial$sex[which(A13.reference$ID.right==MetaTrial$X8x8_quad)]
      
      
      A24.reference<-HelperTrial[which(HelperTrial$Arena=="24"),]
      A24.reference$sex.left<-MetaTrial$sex[which(A24.reference$ID.left==MetaTrial$X8x8_quad)]
      A24.reference$sex.right<-MetaTrial$sex[which(A24.reference$ID.right==MetaTrial$X8x8_quad)]
      
      
      
      #GEOMETRIC RESTRICTIONS
      ##############################
      ############################################
      ################################################################
      if(A13.reference$TopScreen=="inner"){
        #SPLITS LEFT/RIGHT
        
        longA13$quad<-ifelse(IsSide(A13.reference$topdivide.x,A13.reference$topdivide.y,
                             A13.reference$botdivide.x,A13.reference$botdivide.y,
                             longA13$true.x,longA13$true.y)<0,
                             "Right",
                             "Left")
        #################
        #LEFT SIDE
        #keep below left top
        longA13$keep1<-ifelse(longA13$quad=="Left",
                              ifelse(IsSide(A13.reference$toplefthoriz.x,A13.reference$toplefthoriz.y,
                                            A13.reference$topdivide.x,A13.reference$topdivide.y,
                                            longA13$true.x,longA13$true.y)>0,
                                     "in","out"),NA)
        #keep right of edge
        longA13$keep2<-ifelse(longA13$quad=="Left",
                              ifelse(IsSide(A13.reference$botlefthoriz.x,A13.reference$botlefthoriz.y,
                                            A13.reference$toplefthoriz.x,A13.reference$toplefthoriz.y,
                                            longA13$true.x,longA13$true.y)>0,
                                     "in","out"),NA)
        #keep up of bottom
        longA13$keep3<-ifelse(longA13$quad=="Left",
                              ifelse(IsSide(A13.reference$botdivide.x,A13.reference$botdivide.y,
                                            A13.reference$botlefthoriz.x,A13.reference$toplefthoriz.y,
                                            longA13$true.x,longA13$true.y)>0,
                                     "in","out"),NA)
        #############################
        #################
        #RIGHT SIDE
        #Keep below right top
        longA13$keep1<-ifelse(longA13$quad=="Left",longA13$keep1,
                              ifelse(IsSide(A13.reference$topdivide.x,A13.reference$topdivide.y,
                                            A13.reference$toprighthoriz.x,A13.reference$toprighthoriz.y,
                                            longA13$true.x,longA13$true.y)>0,"in","out"))
        #keep left of edge
        longA13$keep2<-ifelse(longA13$quad=="Left",longA13$keep2,
                              ifelse(IsSide(A13.reference$toprighthoriz.x,A13.reference$toprighthoriz.y,
                                            A13.reference$botrighthoriz.x,A13.reference$botrighthoriz.y,
                                            longA13$true.x,longA13$true.y)>0,"in","out"))
        #keep up of bottom
        longA13$keep3<-ifelse(longA13$quad=="Left",longA13$keep3,
                              ifelse(IsSide(A13.reference$botrighthoriz.x,A13.reference$botrighthoriz.y,
                                            A13.reference$botdivide.x,A13.reference$botdivide.y,
                                            longA13$true.x,longA13$true.y)>0,"in","out"))
        
        #########################################################
        longA13$keep<-ifelse(apply(longA13[,c(15:17)],1,function(x) any(x=="out")),"out","in")
        longA13<-longA13[,-c(15:17)]
        
        keeplongA13<-longA13
        keeplongA13[which(keeplongA13$keep=="out"),c(9:12,14)]<-NA
        
        keeplongA13$quad<-ifelse(keeplongA13$quad=="Left",
                                 as.character(A13.reference$ID.left),
                                 as.character(A13.reference$ID.right))
        
      }
      if(A13.reference$TopScreen=="outer"){
        longA13$quad<-ifelse(IsSide(A13.reference$topdivide.x,A13.reference$topdivide.y,
                                    A13.reference$botdivide.x,A13.reference$botdivide.y,
                                    longA13$true.x,longA13$true.y)<0,
                             "Right",
                             "Left")

        #################
        #LEFT SIDE
        #keep below left top
        longA13$keep1<-ifelse(longA13$quad=="Left",
                             ifelse(IsSide(A13.reference$toplefthoriz.x,A13.reference$toplefthoriz.y,
                                           A13.reference$topdivide.x,A13.reference$topdivide.y,
                                           longA13$true.x,longA13$true.y)>0,
                                    "in","out"),NA)
        #keep right of edge
        longA13$keep2<-ifelse(longA13$quad=="Left",
                             ifelse(IsSide(A13.reference$botlefthoriz.x,A13.reference$botlefthoriz.y,
                                           A13.reference$toplefthoriz.x,A13.reference$toplefthoriz.y,
                                           longA13$true.x,longA13$true.y)>0,
                                    "in","out"),NA)
        #keep up of bottom
        longA13$keep3<-ifelse(longA13$quad=="Left",
                             ifelse(IsSide(A13.reference$botdivide.x,A13.reference$botdivide.y,
                                           A13.reference$botlefthoriz.x,A13.reference$toplefthoriz.y,
                                           longA13$true.x,longA13$true.y)>0,
                                    "in","out"),NA)
        #############################
        #################
        #RIGHT SIDE
        #Keep below right top
        longA13$keep1<-ifelse(longA13$quad=="Left",longA13$keep1,
                             ifelse(IsSide(A13.reference$topdivide.x,A13.reference$topdivide.y,
                                           A13.reference$toprighthoriz.x,A13.reference$toprighthoriz.y,
                                           longA13$true.x,longA13$true.y)>0,"in","out"))
        #keep left of edge
        longA13$keep2<-ifelse(longA13$quad=="Left",longA13$keep2,
                             ifelse(IsSide(A13.reference$toprighthoriz.x,A13.reference$toprighthoriz.y,
                                           A13.reference$botrighthoriz.x,A13.reference$botrighthoriz.y,
                                           longA13$true.x,longA13$true.y)>0,"in","out"))
        #keep up of bottom
        longA13$keep3<-ifelse(longA13$quad=="Left",longA13$keep3,
                             ifelse(IsSide(A13.reference$botrighthoriz.x,A13.reference$botrighthoriz.y,
                                           A13.reference$botdivide.x,A13.reference$botdivide.y,
                                           longA13$true.x,longA13$true.y)>0,"in","out"))
        
        #########################################################
        longA13$keep<-ifelse(apply(longA13[,c(15:17)],1,function(x) any(x=="out")),"out","in")
        longA13<-longA13[,-c(15:17)]
        
        keeplongA13<-longA13
        keeplongA13[which(keeplongA13$keep=="out"),c(9:12,14)]<-NA
        
        keeplongA13$quad<-ifelse(keeplongA13$quad=="Left",
                                 as.character(A13.reference$ID.left),
                                 as.character(A13.reference$ID.right))
        
      }
      
      if(A24.reference$TopScreen=="inner"){
        #SPLITS LEFT/RIGHT
        longA24$quad<-ifelse(IsSide(A24.reference$topdivide.x,A24.reference$topdivide.y,
                                    A24.reference$botdivide.x,A24.reference$botdivide.y,
                                    longA24$true.x,longA24$true.y)<0,
                             "Right",
                             "Left")

        #################
        #LEFT SIDE
        #keep below left top
        longA24$keep1<-ifelse(longA24$quad=="Left",
                             ifelse(IsSide(A24.reference$toplefthoriz.x,A24.reference$toplefthoriz.y,
                                           A24.reference$topdivide.x,A24.reference$topdivide.y,
                                           longA24$true.x,longA24$true.y)>0,
                                    "in","out"),NA)
        #keep right of edge
        longA24$keep2<-ifelse(longA24$quad=="Left",
                             ifelse(IsSide(A24.reference$botlefthoriz.x,A24.reference$botlefthoriz.y,
                                           A24.reference$toplefthoriz.x,A24.reference$toplefthoriz.y,
                                           longA24$true.x,longA24$true.y)>0,
                                    "in","out"),NA)
        #keep up of bottom
        longA24$keep3<-ifelse(longA24$quad=="Left",
                             ifelse(IsSide(A24.reference$botdivide.x,A24.reference$botdivide.y,
                                           A24.reference$botlefthoriz.x,A24.reference$toplefthoriz.y,
                                           longA24$true.x,longA24$true.y)>0,
                                    "in","out"),NA)
        
        #################
        #RIGHT SIDE
        #Keep below right top
        longA24$keep1<-ifelse(longA24$quad=="Left",longA24$keep1,
                             ifelse(IsSide(A24.reference$topdivide.x,A24.reference$topdivide.y,
                                           A24.reference$toprighthoriz.x,A24.reference$toprighthoriz.y,
                                           longA24$true.x,longA24$true.y)>0,"in","out"))
        #keep left of edge
        longA24$keep2<-ifelse(longA24$quad=="Left",longA24$keep2,
                             ifelse(IsSide(A13.reference$toprighthoriz.x,A13.reference$toprighthoriz.y,
                                           A13.reference$botrighthoriz.x,A13.reference$botrighthoriz.y,
                                           longA24$true.x,longA24$true.y)>0,"in","out"))
        #keep up of bottom
        longA24$keep3<-ifelse(longA24$quad=="Left",longA24$keep3,
                             ifelse(IsSide(A13.reference$botrighthoriz.x,A13.reference$botrighthoriz.y,
                                           A13.reference$botdivide.x,A13.reference$botdivide.y,
                                           longA24$true.x,longA24$true.y)>0,"in","out"))
        
        #########################################################
        longA24$keep<-ifelse(apply(longA24[,c(15:17)],1,function(x) any(x=="out")),"out","in")
        longA24<-longA24[,-c(15:17)]
        
        keeplongA24<-longA24
        keeplongA24[which(keeplongA24$keep=="out"),c(9:12,14)]<-NA
        
        keeplongA24$quad<-ifelse(keeplongA24$quad=="Left",
                                 as.character(A24.reference$ID.left),
                                 as.character(A24.reference$ID.right))
        
      }
      if(A24.reference$TopScreen=="outer"){
        longA24$quad<-ifelse(IsSide(A24.reference$topdivide.x,A24.reference$topdivide.y,
                                    A24.reference$botdivide.x,A24.reference$botdivide.y,
                                    longA24$true.x,longA24$true.y)<0,
                             "Right",
                             "Left")
        #################
        #LEFT SIDE
        #keep below left top
        longA24$keep1<-ifelse(longA24$quad=="Left",
                              ifelse(IsSide(A24.reference$toplefthoriz.x,A24.reference$toplefthoriz.y,
                                            A24.reference$topdivide.x,A24.reference$topdivide.y,
                                            longA24$true.x,longA24$true.y)>0,
                                     "in","out"),NA)
        #keep right of edge
        longA24$keep2<-ifelse(longA24$quad=="Left",
                              ifelse(IsSide(A24.reference$botlefthoriz.x,A24.reference$botlefthoriz.y,
                                            A24.reference$toplefthoriz.x,A24.reference$toplefthoriz.y,
                                            longA24$true.x,longA24$true.y)>0,
                                     "in","out"),NA)
        #keep up of bottom
        longA24$keep3<-ifelse(longA24$quad=="Left",
                              ifelse(IsSide(A24.reference$botdivide.x,A24.reference$botdivide.y,
                                            A24.reference$botlefthoriz.x,A24.reference$toplefthoriz.y,
                                            longA24$true.x,longA24$true.y)>0,
                                     "in","out"),NA)
        
        #################
        #RIGHT SIDE
        #Keep below right top
        longA24$keep1<-ifelse(longA24$quad=="Left",longA24$keep1,
                              ifelse(IsSide(A24.reference$topdivide.x,A24.reference$topdivide.y,
                                            A24.reference$toprighthoriz.x,A24.reference$toprighthoriz.y,
                                            longA24$true.x,longA24$true.y)>0,"in","out"))
        #keep left of edge
        longA24$keep2<-ifelse(longA24$quad=="Left",longA24$keep2,
                              ifelse(IsSide(A13.reference$toprighthoriz.x,A13.reference$toprighthoriz.y,
                                            A13.reference$botrighthoriz.x,A13.reference$botrighthoriz.y,
                                            longA24$true.x,longA24$true.y)>0,"in","out"))
        #keep up of bottom
        longA24$keep3<-ifelse(longA24$quad=="Left",longA24$keep3,
                              ifelse(IsSide(A13.reference$botrighthoriz.x,A13.reference$botrighthoriz.y,
                                            A13.reference$botdivide.x,A13.reference$botdivide.y,
                                            longA24$true.x,longA24$true.y)>0,"in","out"))
        
        #########################################################
        longA24$keep<-ifelse(apply(longA24[,c(15:17)],1,function(x) any(x=="out")),"out","in")
        longA24<-longA24[,-c(15:17)]
        
        
        keeplongA24<-longA24
        keeplongA24[which(keeplongA24$keep=="out"),c(9:12,14)]<-NA
        
        keeplongA24$quad<-ifelse(keeplongA24$quad=="Left",
                                 as.character(A24.reference$ID.left),
                                 as.character(A24.reference$ID.right))
        
      }
      
      
      #DROP 'pee'instances if they don't fall in the polygon of their quad
      keeplongA13[which(keeplongA13$quad==A13.reference$ID.left & 
                          !(point.in.polygon(keeplongA13$true.x,
                                            keeplongA13$true.y,
                                            A13.left$X,
                                            A13.left$Y))),c(9:12,14)]<-NA
      keeplongA13[which(keeplongA13$quad==A13.reference$ID.right & 
                          !(point.in.polygon(keeplongA13$true.x,
                                             keeplongA13$true.y,
                                             A13.right$X,
                                             A13.right$Y))),c(9:12,14)]<-NA
      
      
      keeplongA24[which(keeplongA24$quad==A24.reference$ID.left & 
                          !(point.in.polygon(keeplongA24$true.x,
                                             keeplongA24$true.y,
                                             A24.left$X,
                                             A24.left$Y))),c(9:12,14)]<-NA
      keeplongA24[which(keeplongA24$quad==A24.reference$ID.right & 
                          !(point.in.polygon(keeplongA24$true.x,
                                             keeplongA24$true.y,
                                             A24.right$X,
                                             A24.right$Y))),c(9:12,14)]<-NA
      
      
      
      
      keeplongA13<-keeplongA13[,c(1:12,14:15,13)]
      keeplongA24<-keeplongA24[,c(1:12,14:15,13)]
      
      keeplongA13<-keeplongA13[order(keeplongA13$dailysecond),]
      keeplongA24<-keeplongA24[order(keeplongA24$dailysecond),]
      
      #A13
      #FIND MINIMUM DISTANCE WITHIN 5 seconds prior to pee detection
      for(rsrl in 1:nrow(keeplongA13)){

        #if there is a valid pee
        if(!is.na(keeplongA13[rsrl,"quad"])){
          #print(rsrl)
          #for backcalculating locations, ignore if backlocations are unavailable
          if(rsrl>5){
            backlook<-rsrl-5
          } else {
            backlook<-rsrl
          }
          look.13<-keeplongA13[c(backlook:rsrl),]
          reflook13<-keeplongA13[rsrl,]
          
          if(reflook13[,"quad"]=="a1"){
            look.13$distZ<-sqrt((look.13$a1x0-reflook13$true.x)^2 + (look.13$a1y0-reflook13$true.y)^2)
            if(all(is.na(look.13$distZ))){ #if all distances are missing, return distance of 1000
              keeplongA13$distance[rsrl]<-10000
            } else {
              keeplongA13$distance[rsrl]<-min(look.13$distZ,na.rm = TRUE)
            }
          }else{
            look.13$distZ<-sqrt((look.13$a3x0-reflook13$true.x)^2 + (look.13$a3y0-reflook13$true.y)^2)
            if(all(is.na(look.13$distZ))){ #if all distances are missing, return distance of 1000
              keeplongA13$distance[rsrl]<-10000
            } else {
              keeplongA13$distance[rsrl]<-min(look.13$distZ,na.rm = TRUE)
            }
          }
        } else {
          keeplongA13$distance[rsrl]<-NA
        }
        
      }


      
      #A24
      #FIND MINIMUM DISTANCE WITHIN 5 seconds prior to pee detection
      for(rsrl in 1:nrow(keeplongA24)){
        
        #if there is a valid pee
        if(!is.na(keeplongA24[rsrl,"quad"])){
          #print(rsrl)
          #for backcalculating locations, ignore if backlocations are unavailable
          if(rsrl>5){
            backlook<-rsrl-5
          } else {
            backlook<-rsrl
          }
          look.24<-keeplongA24[c(backlook:rsrl),]
          reflook24<-keeplongA24[rsrl,]
          
          if(reflook24[,"quad"]=="a2"){
            look.24$distZ<-sqrt((look.24$a2x0-reflook24$true.x)^2 + (look.24$a2y0-reflook24$true.y)^2)
            if(all(is.na(look.24$distZ))){ #if all distances are missing, return distance of 1000
              keeplongA24$distance[rsrl]<-10000
            } else {
              keeplongA24$distance[rsrl]<-min(look.24$distZ,na.rm = TRUE)
            }
          }else{
            look.24$distZ<-sqrt((look.24$a4x0-reflook24$true.x)^2 + (look.24$a4y0-reflook24$true.y)^2)
            if(all(is.na(look.24$distZ))){ #if all distances are missing, return distance of 1000
              keeplongA24$distance[rsrl]<-10000
            } else {
              keeplongA24$distance[rsrl]<-min(look.24$distZ,na.rm = TRUE)
            }
          }
        } else {
          keeplongA24$distance[rsrl]<-NA
        }
        
      }
      
      
      
      if(implementdistancethreshold==TRUE){
        print(paste0("Pee->Mouse distance threshold enforced"))
        #discard detections that were farther than 10cm
        keeplongA13[which(keeplongA13$distance==10000),c(9:14,16)]<-NA
        keeplongA13$distance<-keeplongA13$distance/a13distance
        keeplongA13[which(keeplongA13$distance>10),c(9:14,16)]<-NA
        
        #discard detections that were farther than 10cm
        keeplongA24[which(keeplongA24$distance==10000),c(9:14,16)]<-NA
        keeplongA24$distance<-keeplongA24$distance/a24distance
        keeplongA24[which(keeplongA24$distance>10),c(9:14,16)]<-NA
      }

      


      
      keeplongA13$a1roi<-as.factor(keeplongA13$a1roi);
      keeplongA13$a3roi<-as.factor(keeplongA13$a3roi);
      keeplongA24$a2roi<-as.factor(keeplongA24$a2roi);
      keeplongA24$a4roi<-as.factor(keeplongA24$a4roi);
      

      
      library(forcats)
      library(dplyr)
      keeplongA13 = keeplongA13 %>% mutate_if(is.factor,
                                                fct_explicit_na,
                                                  na_level = 'missing')
      keeplongA24 = keeplongA24 %>% mutate_if(is.factor,
                                                fct_explicit_na,
                                                  na_level = 'missing')

      
      
      keeplongA13$quad<-as.factor(keeplongA13$quad)
      keeplongA24$quad<-as.factor(keeplongA24$quad)
      
      
      
      kt13<-keeplongA13
      kt24<-keeplongA24
      
      #####clean up pee detections, keeping only the first in a string of consecutive pee detections (but also dropping single instance detections FIRST)
      

      
      for(tt in 1:2){
        if(tt==1){
          peefile<-keeplongA13
          specifdistA<-a13distance
        }
        if(tt==2){
          peefile<-keeplongA24
          specifdistA<-a24distance
        }
        specifdistB<-specifdistA*2
        peefile$ispee<-ifelse(!is.na(peefile$Clustsize),1,0)
        peefile$second<-peefile$dailysecond
        #will go through each unique second and group pee events if they are within a certain distance of each other,
        # so that a single pee event will not be recorded multiply
        # this loop written on July 4, 2019 :)
        
        duplicatedseconds<-peefile[duplicated(peefile$second)|duplicated(peefile$second,fromLast = TRUE),]
        startingrowflag<-1
        
        #only run duplicated fixer if duplicates detected
        if(nrow(duplicatedseconds)>0){
          for(usa in 1:length(unique(duplicatedseconds$second))){
            
            sec.value<-unique(duplicatedseconds$second)[usa]
            where.is.this.sec<-which(peefile$second==sec.value)
            howmanduplicates<-length(where.is.this.sec)
            goto<-where.is.this.sec[1]-1
            newkeepA<-peefile[c(startingrowflag:goto),]
            newkeepC<-peefile[c((goto+howmanduplicates+1):nrow(peefile)),]
            
            sec1<-duplicatedseconds[which(duplicatedseconds$second==sec.value),]
            
            sec1NAs<-sec1[which(is.na(sec1$Clustsize)),]
            
            
            newkeepA$group<-NA
            newkeepC$group<-NA
            sec1$group<-NA
            #if there are multiple rows/second, then there might be inflated detections (i.e. multiple detections for a single pee event)
            #this keeps only a single pee-event PER group (so it can keep multiple pee events/second, but only if they are spatial distinct)
            if(nrow(sec1)>1){
              if(length(which(is.na(sec1$Clustsize))>0)){
                sec1a<-sec1[-which(is.na(sec1$Clustsize)),]
                sec1$group<-1
                if(nrow(sec1a)==0){
                  sec1<-sec1[1,]
                  sec1$group<-NA
                } else {
                  sec1<-sec1a
                }
                
              }
              
              if(nrow(sec1)>1){
                sec.dist.matrix<-dist(sec1[,c("true.x","true.y")])#uses x/y coordinates to create distance matrix
                sec.clust <-fastcluster::hclust(sec.dist.matrix, method="single", members=NULL) #Clusters pee based on distance matrix
                sec.cutted<-cutree(sec.clust,h=specifdistA)#Determines cutoff for what is a cluster (based on 'specificDistance' based on folder title (MB1,MB2))
                sec1$group<-sec.cutted
                
                #keepgroupvec<-sec.cutted
                
                for(independence in 1:length(unique(sec1$group))){
                  this.group<-sec1[which(sec1$group==unique(sec1$group)[independence]),]
                  keep.from.this.group<-this.group[which.max(this.group$Clustsize),]
                  
                  if(independence==1){
                    kg<-keep.from.this.group
                  } else {
                    kg<-rbind(kg,keep.from.this.group)
                  }
                  
                }
                gv<-which(colnames(kg)=="group")
                sec1<-kg
                #sec1<-kg[,-gv]#drop column corresponding to second-specific grouping
              }
              #sec1<-rbind(sec1NAs,sec1)
            }
            
            
            
            peefile<-rbind(newkeepA,sec1,newkeepC)
            
            if(usa==1){
              newSUB<-sec1
            } else {
              newSUB<-rbind(newSUB,sec1)
            }
          }
        } else {
          newSUB<-peefile
        }
        
        #plug in group values back into peefile (for concurrently detected pee spots that are sufficiently separated spatially)
        for(fw in 1:nrow(newSUB)){
          plugin<-newSUB[fw,]
          indx<-which(peefile$dailysecond==plugin$dailysecond & peefile$mTemp==plugin$mTemp & peefile$distance==plugin$distance)
          peefile[indx,]<-plugin
        }
        
        newkeep<-peefile
        
        #find and then uniquely label consecutive groups of pee-detections  
        vs<-rle(newkeep$keep)
        is.na.rle <- rle(is.na(newkeep$keep))
        
        counting<-1  
        for(Am in 1:length(vs$values)){
          if(!is.na(vs$values[Am])){
            vs$values[[Am]]<-paste0(counting,vs$values[[Am]])
            counting<-counting+1
          }
        }
        newkeep$grouping<-inverse.rle(vs)  
        
        counting2<-1
        for(Merica in 1:length(is.na.rle$values)){
          if(is.na.rle$values[Merica]){
            is.na.rle$values[[Merica]]<-paste0(counting2,is.na.rle$values[[Merica]])
            counting2<-counting2+1
          }
        }
        newkeep$grouping2<-inverse.rle(is.na.rle) 
        
        newkeep$maingroup<-ifelse(grepl('TRUE',newkeep$grouping2),newkeep$grouping2,newkeep$grouping)
        
        groupingANDgrouping2<-c(which(colnames(newkeep)=="grouping"),which(colnames(newkeep)=="grouping2"))
        
        newkeep<-newkeep[,-groupingANDgrouping2]#eliminate grouping and grouping2 columns
        
        
        ####  
        for(usa2 in 1:length(unique(newkeep$maingroup))){ #start at 2 because first 'grouping' is NA values
          grpA<-newkeep[which(newkeep$maingroup==unique(newkeep$maingroup)[usa2]),]
          
          #if maingroup does NOT have 'in' in it, it belongs to an NA group, and we 'else' and keep the whole thing
          if(grepl('in',grpA$maingroup)[1]){
              grpA<-grpA[which(!is.na(grpA$Clustsize)),]#now drop clusters with size = NA, assigned if they weren't in the right quad
              nrgoups<-length(unique(grpA$group))
              useme<-grpA
              
              
              
              if(nrow(grpA)<2){
                #if there is ONLY A SINGLE row/group, then that is likely a misdetection (real pee has a time course and is likely to be detected in sequential frames)   
                grpA[1,c("Clustsize","mTemp","true.y","true.x","distance","ispee","keep")]<-NA

              } else {
              
              for(whiskerdoo in 1:nrgoups){
                thisparticulargroup<-unique(useme$group)[whiskerdoo]
                
                #a consecutive string of pees, none of which had been detected at the same second, wont already have 'group' ids
                #and thus will have NA as group, for which we need to evaluate with is.na NOT ==
                if(is.na(thisparticulargroup)){
                  grp<-useme[which(is.na(useme$group)),]
                } else {
                  grp<-useme[which(useme$group==thisparticulargroup),]
                }
                

                
                
              #if there are multiple rows/group, then there might be inflated detections (i.e. multiple detections for a single pee event)
              if(nrow(grp)>1){
                if(length(which(is.na(grp$Clustsize)))>0){
                  grp<-grp[-which(is.na(grp$Clustsize)),] #unlikely to have groups containing NA values, but just in case
                }
                
                grp.dist.matrix<-dist(grp[,c("true.x","true.y")])#uses x/y coordinates to create distance matrix
                grp.clust <-fastcluster::hclust(grp.dist.matrix, method="single", members=NULL) #Clusters pee based on distance matrix
                grp.cutted<-cutree(grp.clust,h=specifdistB)#Determines cutoff for what is a cluster (based on 'specificDistance' based on folder title (MB1,MB2))
                grp$group<-grp.cutted
                
                for(independence in 1:length(unique(grp$group))){
                  this.group<-grp[which(grp$group==unique(grp$group)[independence]),]
                  this.group[-which.max(this.group$Clustsize),c("Clustsize","mTemp","true.y","true.x","distance","ispee","keep")]<-NA
                  
                  if(independence==1){
                    kg<-this.group
                  } else {
                    kg<-rbind(kg,this.group)
                  }
                }
                disgroup<-which(colnames(kg)=="group")
                grp<-kg
                #grp<-kg[,-disgroup]#drop column corresponding to consecutive grouping
              } 
                
                if(whiskerdoo==1){
                  grpA<-grp
                } else {
                  grpA<-rbind(grpA,grp)
                }
            }#close whiskerdoo
              
              }
              
          } else {
            #do nothing, just keep all NA rows
          }  
          if(usa2==1){
            newkeep2<-grpA
          } else {
            newkeep2<-rbind(newkeep2,grpA)
          }
        } 
        
        newkeep3<-newkeep2[order(newkeep2$second),]
        
        
        if(tt==1){
          keeplongA13<-newkeep3
        }
        if(tt==2){
          keeplongA24<-newkeep3
        }
        
      }

      zza13<-keeplongA13
      zza24<-keeplongA24
      
      
      #if there are any rows with NA in the dailysecond column, drop those rows
      if(length(which(is.na(keeplongA13$dailysecond)))>0){
         keeplongA13<-keeplongA13[-which(is.na(keeplongA13$dailysecond)),]
      }
      if(length(which(is.na(keeplongA24$dailysecond)))>0){
        keeplongA24<-keeplongA24[-which(is.na(keeplongA24$dailysecond)),]
      }
      
      
      
      
      
      ######################################################
      library(dplyr)
      Just.Frames2<-as.data.frame(unique(keeplongA24$dailysecond));colnames(Just.Frames2)<-"dailysecond"
      
      
      
      #a2peeroll
      a2.totpee.info<-subset(keeplongA24,(keeplongA24$quad=="a2" |is.na(keeplongA24$quad) ))
      a2.totpee<-a2.totpee.info %>%
                  group_by(dailysecond) %>%
                    summarise_at("Clustsize", sum, na.rm = TRUE)
      a2spots<-a2.totpee.info %>%
                  group_by(dailysecond,a2roi) %>%
                    summarise_at(c("a2x0","a2y0"), mean, na.rm = TRUE)
      a2spots$a2x0[is.nan(a2spots$a2x0)]<-NA;a2spots$a2y0[is.nan(a2spots$a2y0)]<-NA;
      a2info<-merge(a2spots,a2.totpee,by="dailysecond")      
      peeroller<-merge(Just.Frames2,a2info,by="dailysecond",all.x=TRUE)
      colnames(peeroller)[5]<-"a2.pee.vol"
      peeroller$a2.roll.pee<-rollapply(peeroller$a2.pee.vol,width=20,function(x) sum(x,na.rm=TRUE),partial=TRUE)

      #a4peeroll
      a4.totpee.info<-subset(keeplongA24,(keeplongA24$quad=="a4" |is.na(keeplongA24$quad) ))
      a4.totpee<-a4.totpee.info %>%
                  group_by(dailysecond) %>%
                    summarise_at("Clustsize", sum, na.rm = TRUE)
      a4spots<-a4.totpee.info %>%
                    group_by(dailysecond,a4roi) %>%
                      summarise_at(c("a4x0","a4y0"), mean, na.rm = TRUE)
      a4spots$a4x0[is.nan(a4spots$a4x0)]<-NA;a4spots$a4y0[is.nan(a4spots$a4y0)]<-NA;
      a4info<-merge(a4spots,a4.totpee,by="dailysecond")            
      peeroller<-merge(peeroller,a4info,by="dailysecond",all.x=TRUE)
      colnames(peeroller)[10]<-"a4.pee.vol"
      peeroller$a4.roll.pee<-rollapply(peeroller$a4.pee.vol,width=20,function(x) sum(x,na.rm=TRUE),partial=TRUE)
      #a1peeroll
      a1.totpee.info<-subset(keeplongA13,(keeplongA13$quad=="a1" |is.na(keeplongA13$quad) ))
      a1.totpee<-a1.totpee.info %>%
                  group_by(dailysecond) %>%
                    summarise_at("Clustsize", sum, na.rm = TRUE)
      a1spots<-a1.totpee.info %>%
        group_by(dailysecond,a1roi) %>%
        summarise_at(c("a1x0","a1y0"), mean, na.rm = TRUE)
      a1spots$a1x0[is.nan(a1spots$a1x0)]<-NA;a1spots$a1y0[is.nan(a1spots$a1y0)]<-NA;
      a1info<-merge(a1spots,a1.totpee,by="dailysecond")            
      peeroller<-merge(peeroller,a1info,by="dailysecond",all.x=TRUE)
      colnames(peeroller)[15]<-"a1.pee.vol"
      peeroller$a1.roll.pee<-rollapply(peeroller$a1.pee.vol,width=20,function(x) sum(x,na.rm=TRUE),partial=TRUE)
      #a3peeroll
      a3.totpee.info<-subset(keeplongA13,(keeplongA13$quad=="a3" |is.na(keeplongA13$quad) ))
      a3.totpee<-a3.totpee.info %>%
        group_by(dailysecond) %>%
        summarise_at("Clustsize", sum, na.rm = TRUE)
      a3spots<-a3.totpee.info %>%
        group_by(dailysecond,a3roi) %>%
        summarise_at(c("a3x0","a3y0"), mean, na.rm = TRUE)
      a3spots$a3x0[is.nan(a3spots$a3x0)]<-NA;a3spots$a3y0[is.nan(a3spots$a3y0)]<-NA;
      a3info<-merge(a3spots,a3.totpee,by="dailysecond")            
      peeroller<-merge(peeroller,a3info,by="dailysecond",all.x=TRUE)
      colnames(peeroller)[20]<-"a3.pee.vol"
      peeroller$a3.roll.pee<-rollapply(peeroller$a3.pee.vol,width=20,function(x) sum(x,na.rm=TRUE),partial=TRUE)
      ###################################################
      keep.peeroller<-peeroller
      
      
      peeroller[,c("a2.pee.vol", "a2.roll.pee", 
                   "a4.pee.vol", "a4.roll.pee", 
                   "a1.pee.vol", "a1.roll.pee", 
                   "a3.pee.vol", "a3.roll.pee")][is.na(peeroller[,c("a2.pee.vol", "a2.roll.pee", 
                                                                    "a4.pee.vol", "a4.roll.pee", 
                                                                    "a1.pee.vol", "a1.roll.pee", 
                                                                    "a3.pee.vol", "a3.roll.pee")])]<-0
      
      summary(peeroller)
      
      
      peeroller<-peeroller[,c("dailysecond", 
                              "a1roi", "a1x0", "a1y0", "a1.pee.vol", "a1.roll.pee",
                              "a3roi", "a3x0","a3y0", "a3.pee.vol", "a3.roll.pee",
                              "a2roi", "a2x0", "a2y0", "a2.pee.vol", "a2.roll.pee", 
                              "a4roi", "a4x0", "a4y0", "a4.pee.vol", "a4.roll.pee")]
      
      #IMPUTE MISSING X/Y COORDINATES FOR MOUSE LOCATIONS
      peeroller[,c("a1x0", "a1y0", 
                   "a3x0","a3y0", "a2x0", "a2y0", "a4x0", "a4y0")]<-apply(peeroller[,c("a1x0", "a1y0", 
                               "a3x0","a3y0", "a2x0", "a2y0", "a4x0", "a4y0")],2,function(X){approxfun(seq_along(X),X)(seq_along(X))})
      
      

      A1trj <- TrajFromCoords(track=peeroller[,c("dailysecond","a1x0", "a1y0")],xCol="a1x0",yCol="a1y0",timeCol="dailysecond",fps=1,spatial="pixels")
      A2trj <- TrajFromCoords(track=peeroller[,c("dailysecond","a2x0", "a2y0")],xCol="a2x0",yCol="a2y0",timeCol="dailysecond",fps=1,spatial="pixels")
      A3trj <- TrajFromCoords(track=peeroller[,c("dailysecond","a3x0", "a3y0")],xCol="a3x0",yCol="a3y0",timeCol="dailysecond",fps=1,spatial="pixels")
      A4trj <- TrajFromCoords(track=peeroller[,c("dailysecond","a4x0", "a4y0")],xCol="a4x0",yCol="a4y0",timeCol="dailysecond",fps=1,spatial="pixels")
      trajlist<-list(A1trj,A2trj,A3trj,A4trj)
      
      for(Cairo in 1:4){
        
        quadname<-c("a1","a2","a3","a4")[Cairo]
        
        thistrack<-trajlist[[Cairo]]
        derivs<-TrajDerivatives(thistrack)
        ii<-merge(thistrack,cbind(derivs$speed,derivs$speedTimes),by.x="displacementTime",by.y="V2",all.x=TRUE)
        colnames(ii)[ncol(ii)]<-"speed"
        accels<-cbind(derivs$accelerationTimes,derivs$acceleration);colnames(accels)<-c("displacementTime","acceleration")
        iii<-merge(ii,accels,by="displacementTime",all.x=TRUE)
        colnames(iii)[2]<-"dailysecond";iii<-as.data.frame(iii)
        iv<-iii[,c("dailysecond","speed","acceleration")]
        iv$rollspeed<-rollapply(iv$speed,width=20,function(x) mean(x,na.rm=TRUE),partial=TRUE)
        iv$rollacceleration<-rollapply(iv$acceleration,width=20,function(x) mean(x,na.rm=TRUE),partial=TRUE)
        
        iv$displacement<-Mod(iii$displacement)
        iv$displacement[is.na(iv$displacement)]<-0
        
          
        if(quadname=="a1" | quadname=="a3")
          iv$displacement<-iv$displacement/a13distance
        
        if(quadname=="a2" | quadname=="a4")
          iv$displacement<-iv$displacement/a24distance
        
        colnames(iv)[c(2:6)]<-c(paste(quadname,colnames(iv)[c(2:6)],sep='.'))
        
        
        peeroller<-merge(peeroller,iv,by="dailysecond",all.x=TRUE)
        
      }
      
      
      peeroller[,c("a1.displacement", "a2.displacement", 
                   "a3.displacement", "a4.displacement")][is.na(peeroller[,c("a1.displacement", "a2.displacement", 
                                                                             "a3.displacement", "a4.displacement")])]<-0
      
      peeroller$a1cumpee<-cumsum(peeroller$a1.pee.vol)
      peeroller$a2cumpee<-cumsum(peeroller$a2.pee.vol)
      peeroller$a3cumpee<-cumsum(peeroller$a3.pee.vol)
      peeroller$a4cumpee<-cumsum(peeroller$a4.pee.vol)
      
      peeroller$a1cumdist<-cumsum(peeroller$a1.displacement)
      peeroller$a2cumdist<-cumsum(peeroller$a2.displacement)
      peeroller$a3cumdist<-cumsum(peeroller$a3.displacement)
      peeroller$a4cumdist<-cumsum(peeroller$a4.displacement)
      
      
      if(cam13daymissing==1){
        #colnames(positpeeA13)[c(104:109)]<-c("a1x0", "a1y0", "a1roi", "a3x0", "a3y0", "a3roi")
        peeroller[,match(c("a1x0", "a1y0", "a1roi", "a3x0", "a3y0", "a3roi","a1.displacement", "a3.displacement","a1cumdist","a3cumdist"),
                         colnames(peeroller))]<-NA
      }
      if(cam24daymissing==1){
        #colnames(positpeeA24)[c(104:109)]<-c("a2x0", "a2y0", "a2roi", "a4x0", "a4y0", "a4roi")
        peeroller[,match(c("a2x0", "a2y0", "a2roi", "a4x0", "a4y0", "a4roi","a2.displacement", "a4.displacement","a2cumdist","a4cumdist"),
                         colnames(peeroller))]<-NA
      }
      
      
      back.that.pee.up<-peeroller
      
      #PHYSICAL LOCATION COUNTER
      #######################################################
      #A1
        #puts 1s for instances of roi, and 0s for errything else
        back.that.pee.up$a1_a2_barrier<-ifelse(is.na(back.that.pee.up$a1roi),0,
                                               ifelse(back.that.pee.up$a1roi=="a1_a2_barrier",1,0))
        back.that.pee.up$a1_a2_barrier<-cumsum(back.that.pee.up$a1_a2_barrier)
        #
        back.that.pee.up$a1_a3_barrier<-ifelse(is.na(back.that.pee.up$a1roi),0,
                                               ifelse(back.that.pee.up$a1roi=="a1_a3_barrier",1,0))
        back.that.pee.up$a1_a3_barrier<-cumsum(back.that.pee.up$a1_a3_barrier)
        #
        back.that.pee.up$a1_central_corner<-ifelse(is.na(back.that.pee.up$a1roi),0,
                                               ifelse(back.that.pee.up$a1roi=="a1_central_corner",1,0))
        back.that.pee.up$a1_central_corner<-cumsum(back.that.pee.up$a1_central_corner)
        #
        back.that.pee.up$a1_water<-ifelse(is.na(back.that.pee.up$a1roi),0,
                                                   ifelse((back.that.pee.up$a1roi=="a1_a3_water"|back.that.pee.up$a1roi=="a1_a2_water"),1,0))
        back.that.pee.up$a1_water<-cumsum(back.that.pee.up$a1_water)
      
      #A3
        #puts 1s for instances of roi, and 0s for errything else
        back.that.pee.up$a3_a1_barrier<-ifelse(is.na(back.that.pee.up$a3roi),0,
                                               ifelse(back.that.pee.up$a3roi=="a3_a1_barrier",1,0))
        back.that.pee.up$a3_a1_barrier<-cumsum(back.that.pee.up$a3_a1_barrier)
        #
        back.that.pee.up$a3_a4_barrier<-ifelse(is.na(back.that.pee.up$a3roi),0,
                                               ifelse(back.that.pee.up$a3roi=="a3_a4_barrier",1,0))
        back.that.pee.up$a3_a4_barrier<-cumsum(back.that.pee.up$a3_a4_barrier)
        #
        back.that.pee.up$a3_central_corner<-ifelse(is.na(back.that.pee.up$a3roi),0,
                                                   ifelse(back.that.pee.up$a3roi=="a3_central_corner",1,0))
        back.that.pee.up$a3_central_corner<-cumsum(back.that.pee.up$a3_central_corner)
        #
        back.that.pee.up$a3_water<-ifelse(is.na(back.that.pee.up$a3roi),0,
                                          ifelse((back.that.pee.up$a3roi=="a3_a1_water"|back.that.pee.up$a3roi=="a3_a4_water"),1,0))
        back.that.pee.up$a3_water<-cumsum(back.that.pee.up$a3_water)
      
      #A2
        #puts 1s for instances of roi, and 0s for errything else
        back.that.pee.up$a2_a1_barrier<-ifelse(is.na(back.that.pee.up$a2roi),0,
                                               ifelse(back.that.pee.up$a2roi=="a2_a1_barrier",1,0))
        back.that.pee.up$a2_a1_barrier<-cumsum(back.that.pee.up$a2_a1_barrier)
        #
        back.that.pee.up$a2_a4_barrier<-ifelse(is.na(back.that.pee.up$a2roi),0,
                                               ifelse(back.that.pee.up$a2roi=="a2_a4_barrier",1,0))
        back.that.pee.up$a2_a4_barrier<-cumsum(back.that.pee.up$a2_a4_barrier)
        #
        back.that.pee.up$a2_central_corner<-ifelse(is.na(back.that.pee.up$a2roi),0,
                                                   ifelse(back.that.pee.up$a2roi=="a2_central_corner",1,0))
        back.that.pee.up$a2_central_corner<-cumsum(back.that.pee.up$a2_central_corner)
        #
        back.that.pee.up$a2_water<-ifelse(is.na(back.that.pee.up$a2roi),0,
                                          ifelse((back.that.pee.up$a2roi=="a2_a1_water"|back.that.pee.up$a2roi=="a2_a4_water"),1,0))
        back.that.pee.up$a2_water<-cumsum(back.that.pee.up$a2_water)

      #A4
        #puts 1s for instances of roi, and 0s for errything else
        back.that.pee.up$a4_a2_barrier<-ifelse(is.na(back.that.pee.up$a4roi),0,
                                               ifelse(back.that.pee.up$a4roi=="a4_a2_barrier",1,0))
        back.that.pee.up$a4_a2_barrier<-cumsum(back.that.pee.up$a4_a2_barrier)
        #
        back.that.pee.up$a4_a3_barrier<-ifelse(is.na(back.that.pee.up$a4roi),0,
                                               ifelse(back.that.pee.up$a4roi=="a4_a3_barrier",1,0))
        back.that.pee.up$a4_a3_barrier<-cumsum(back.that.pee.up$a4_a3_barrier)
        #
        back.that.pee.up$a4_central_corner<-ifelse(is.na(back.that.pee.up$a4roi),0,
                                                   ifelse(back.that.pee.up$a4roi=="a4_central_corner",1,0))
        back.that.pee.up$a4_central_corner<-cumsum(back.that.pee.up$a4_central_corner)
        #
        back.that.pee.up$a4_water<-ifelse(is.na(back.that.pee.up$a4roi),0,
                                          ifelse((back.that.pee.up$a4roi=="a4_a3_water"|back.that.pee.up$a4roi=="a4_a2_water"),1,0))
        back.that.pee.up$a4_water<-cumsum(back.that.pee.up$a4_water)
      
      ####################
        #PEE LOCATION COUNTER
        #######################################################
        #A1
        #puts 1s for instances of roi and pee, and 0s for errything else
        back.that.pee.up$a1_a2_barrier.pee<-ifelse(is.na(back.that.pee.up$a1roi),0,
                                               ifelse((back.that.pee.up$a1roi=="a1_a2_barrier" & back.that.pee.up$a1.pee.vol>0),1,0))
        back.that.pee.up$a1_a2_barrier.pee<-cumsum(back.that.pee.up$a1_a2_barrier.pee)
        #
        back.that.pee.up$a1_a3_barrier.pee<-ifelse(is.na(back.that.pee.up$a1roi),0,
                                               ifelse((back.that.pee.up$a1roi=="a1_a3_barrier"& back.that.pee.up$a1.pee.vol>0),1,0))
        back.that.pee.up$a1_a3_barrier.pee<-cumsum(back.that.pee.up$a1_a3_barrier.pee)
        #
        back.that.pee.up$a1_central_corner.pee<-ifelse(is.na(back.that.pee.up$a1roi),0,
                                                   ifelse((back.that.pee.up$a1roi=="a1_central_corner" & back.that.pee.up$a1.pee.vol>0),1,0))
        back.that.pee.up$a1_central_corner.pee<-cumsum(back.that.pee.up$a1_central_corner.pee)
        #
        back.that.pee.up$a1_water.pee<-ifelse(is.na(back.that.pee.up$a1roi),0,
                                          ifelse(((back.that.pee.up$a1roi=="a1_a3_water"|back.that.pee.up$a1roi=="a1_a2_water") & back.that.pee.up$a1.pee.vol>0),1,0))
        back.that.pee.up$a1_water.pee<-cumsum(back.that.pee.up$a1_water.pee)
        
        #A3
        #puts 1s for instances of roi, and 0s for errything else
        back.that.pee.up$a3_a1_barrier.pee<-ifelse(is.na(back.that.pee.up$a3roi),0,
                                               ifelse((back.that.pee.up$a3roi=="a3_a1_barrier" & back.that.pee.up$a3.pee.vol>0),1,0))
        back.that.pee.up$a3_a1_barrier.pee<-cumsum(back.that.pee.up$a3_a1_barrier.pee)
        #
        back.that.pee.up$a3_a4_barrier.pee<-ifelse(is.na(back.that.pee.up$a3roi),0,
                                               ifelse((back.that.pee.up$a3roi=="a3_a4_barrier" & back.that.pee.up$a3.pee.vol>0),1,0))
        back.that.pee.up$a3_a4_barrier.pee<-cumsum(back.that.pee.up$a3_a4_barrier.pee)
        #
        back.that.pee.up$a3_central_corner.pee<-ifelse(is.na(back.that.pee.up$a3roi),0,
                                                   ifelse((back.that.pee.up$a3roi=="a3_central_corner" & back.that.pee.up$a3.pee.vol>0),1,0))
        back.that.pee.up$a3_central_corner.pee<-cumsum(back.that.pee.up$a3_central_corner.pee)
        #
        back.that.pee.up$a3_water.pee<-ifelse(is.na(back.that.pee.up$a3roi),0,
                                          ifelse(((back.that.pee.up$a3roi=="a3_a1_water"|back.that.pee.up$a3roi=="a3_a4_water") & back.that.pee.up$a3.pee.vol>0),1,0))
        back.that.pee.up$a3_water.pee<-cumsum(back.that.pee.up$a3_water.pee)
        
        #A2
        #puts 1s for instances of roi, and 0s for errything else
        back.that.pee.up$a2_a1_barrier.pee<-ifelse(is.na(back.that.pee.up$a2roi),0,
                                                   ifelse((back.that.pee.up$a2roi=="a2_a1_barrier" & back.that.pee.up$a2.pee.vol>0),1,0))
        back.that.pee.up$a2_a1_barrier.pee<-cumsum(back.that.pee.up$a2_a1_barrier.pee)
        #
        back.that.pee.up$a2_a4_barrier.pee<-ifelse(is.na(back.that.pee.up$a2roi),0,
                                                   ifelse((back.that.pee.up$a2roi=="a2_a4_barrier" & back.that.pee.up$a2.pee.vol>0),1,0))
        back.that.pee.up$a2_a4_barrier.pee<-cumsum(back.that.pee.up$a2_a4_barrier.pee)
        #
        back.that.pee.up$a2_central_corner.pee<-ifelse(is.na(back.that.pee.up$a2roi),0,
                                                       ifelse((back.that.pee.up$a2roi=="a2_central_corner" & back.that.pee.up$a2.pee.vol>0),1,0))
        back.that.pee.up$a2_central_corner.pee<-cumsum(back.that.pee.up$a2_central_corner.pee)
        #
        back.that.pee.up$a2_water.pee<-ifelse(is.na(back.that.pee.up$a2roi),0,
                                              ifelse(((back.that.pee.up$a2roi=="a2_a1_water"|back.that.pee.up$a2roi=="a2_a4_water") & back.that.pee.up$a2.pee.vol>0),1,0))
        back.that.pee.up$a2_water.pee<-cumsum(back.that.pee.up$a2_water.pee)
        
        #A4
        #puts 1s for instances of roi, and 0s for errything else
        back.that.pee.up$a4_a2_barrier.pee<-ifelse(is.na(back.that.pee.up$a4roi),0,
                                                   ifelse((back.that.pee.up$a4roi=="a4_a2_barrier" & back.that.pee.up$a4.pee.vol>0),1,0))
        back.that.pee.up$a4_a2_barrier.pee<-cumsum(back.that.pee.up$a4_a2_barrier.pee)
        #
        back.that.pee.up$a4_a3_barrier.pee<-ifelse(is.na(back.that.pee.up$a4roi),0,
                                                   ifelse((back.that.pee.up$a4roi=="a4_a3_barrier"& back.that.pee.up$a4.pee.vol>0),1,0))
        back.that.pee.up$a4_a3_barrier.pee<-cumsum(back.that.pee.up$a4_a3_barrier.pee)
        #
        back.that.pee.up$a4_central_corner.pee<-ifelse(is.na(back.that.pee.up$a4roi),0,
                                                       ifelse((back.that.pee.up$a4roi=="a4_central_corner" & back.that.pee.up$a4.pee.vol>0),1,0))
        back.that.pee.up$a4_central_corner.pee<-cumsum(back.that.pee.up$a4_central_corner.pee)
        #
        back.that.pee.up$a4_water.pee<-ifelse(is.na(back.that.pee.up$a4roi),0,
                                              ifelse(((back.that.pee.up$a4roi=="a4_a3_water"|back.that.pee.up$a4roi=="a4_a2_water") & back.that.pee.up$a4.pee.vol>0),1,0))
        back.that.pee.up$a4_water.pee<-cumsum(back.that.pee.up$a4_water.pee)
        
        
        
      ###
      peeroller<-back.that.pee.up
      
      
      #Merge interpolated values from peeroller back into keeplongA13, keeplongA24
      ###################################

      
      BACKUP.A13<-keeplongA13;BACKUP.A24<-keeplongA24
      
      
      # keeplongA13[,"a1x0"]<-merge(keeplongA13[,c("dailysecond"),drop=FALSE],
      #                         peeroller[,c("dailysecond","a1x0")],by="dailysecond",
      #                         all.x=TRUE)[2]
      # keeplongA13[,"a1y0"]<-merge(keeplongA13[,c("dailysecond"),drop=FALSE],
      #                         peeroller[,c("dailysecond","a1y0")],by="dailysecond",
      #                         all.x=TRUE)[2]
      # keeplongA13[,"a3x0"]<-merge(keeplongA13[,c("dailysecond"),drop=FALSE],
      #                         peeroller[,c("dailysecond","a3x0")],by="dailysecond",
      #                         all.x=TRUE)[2]
      # keeplongA13[,"a3y0"]<-merge(keeplongA13[,c("dailysecond"),drop=FALSE],
      #                         peeroller[,c("dailysecond","a3y0")],by="dailysecond",
      #                         all.x=TRUE)[2]
      # ######3
      # keeplongA24[,"a2x0"]<-merge(keeplongA24[,c("dailysecond"),drop=FALSE],
      #                         peeroller[,c("dailysecond","a2x0")],by="dailysecond",
      #                         all.x=TRUE)[2]
      # keeplongA24[,"a2y0"]<-merge(keeplongA24[,c("dailysecond"),drop=FALSE],
      #                         peeroller[,c("dailysecond","a2y0")],by="dailysecond",
      #                         all.x=TRUE)[2]
      # keeplongA24[,"a4x0"]<-merge(keeplongA24[,c("dailysecond"),drop=FALSE],
      #                         peeroller[,c("dailysecond","a4x0")],by="dailysecond",
      #                         all.x=TRUE)[2]
      # keeplongA24[,"a4y0"]<-merge(keeplongA24[,c("dailysecond"),drop=FALSE],
      #                         peeroller[,c("dailysecond","a4y0")],by="dailysecond",
      #                         all.x=TRUE)[2]
      #########################################################################
      #Fix plotting locations regardless of whether plotting is to take place
      backup.A13<-keeplongA13
      backup.A24<-keeplongA24
      
      
      #keeplongA13<-backup.A13
      #keeplongA24<-backup.A24
      #
      # real.x.min<-min(c(keeplongA24$a4x0,keeplongA24$a2x0,keeplongA24$true.x),na.rm=TRUE)
      # xlimspacing<-real.x.min-10
      # plot(NA, xlim=c(xlimspacing,650),
      #      ylim=c(-500,500),
      #      xlab='',ylab='',
      #      axes=FALSE)
      # 
      # 
      # ############
      #a13 check
      # lines(A13.outline.xs,A13.outline.ys)
      # lines(keeplongA13$a1x0,keeplongA13$a1y0,col=add.alpha("blue",0.3))
      # lines(keeplongA13$a3x0,keeplongA13$a1y0,col=add.alpha("red",0.3))
      # 
      # justa1s<-keeplongA13[which(keeplongA13$quad=="a1"),]
      # justa3s<-keeplongA13[which(keeplongA13$quad=="a3"),]
      # 
      # points(justa1s$true.x,justa1s$true.y,col=add.alpha("darkblue",0.5))
      # points(justa3s$true.x,justa3s$true.y,col=add.alpha("darkred",0.75))
      # 
      # #a24 check
      # lines(A24.outline.xs,A24.outline.ys)
      # lines(keeplongA24$a2x0,keeplongA24$a2y0,col=add.alpha("green",0.3))
      # lines(keeplongA24$a4x0,keeplongA24$a4y0,col=add.alpha("orange",0.3))
      # 
      # justa2s<-keeplongA24[which(keeplongA24$quad=="a2"),]
      # justa4s<-keeplongA24[which(keeplongA24$quad=="a4"),]
      # 
      # points(justa2s$true.x,justa2s$true.y,col=add.alpha("darkgreen",0.5))
      # points(justa4s$true.x,justa4s$true.y,col=add.alpha("darkorange",0.75))

      pre.roi.lab.a13<-keeplongA13
      pre.roi.lab.a24<-keeplongA24
      
      keeplongA13.info<-peeroi(keeplongA13,trialspecific.ROI)
      keeplongA24.info<-peeroi(keeplongA24,trialspecific.ROI)
      
      a13.errors<-keeplongA13.info[[2]]
      a24.errors<-keeplongA24.info[[2]]
      
      errorframe<-rbind(a13.errors,a24.errors)
      
      if(roimismatchflag==1){
        roimismatchflag<-roimismatchflag+1
        roimismatchframe<-errorframe
      } else {
        roimismatchframe<-rbind(roimismatchframe,errorframe)
      }
      
      keeplongA13<-keeplongA13.info[[1]]
      keeplongA24<-keeplongA24.info[[1]]
      
      unadultered.coords.A13<-keeplongA13[,c("dailysecond","a1x0", "a1y0", "a3x0", "a3y0", "true.y", "true.x","quad")]
      unadultered.coords.A24<-keeplongA24[,c("dailysecond","a2x0", "a2y0", "a4x0", "a4y0", "true.y", "true.x","quad")]
        
        
      #A13
      #A13.outline.xs<-A13.reference[,c(10,8,14,16,6,12,10,8,6)]#old outline
      #A13.outline.ys<-A13.reference[,c(11,9,15,17,7,13,11,9,7)]#
      A13.outline.xs<-c(A13.left$X,A13.right$X)
      A13.outline.ys<-c(A13.left$Y,A13.right$Y)
      
      if(A13.reference$TopScreen=="inner") {
        keeplongA13$true.y<-keeplongA13$true.y*(-1)
        keeplongA13$a1y0<-keeplongA13$a1y0*(-1)
        keeplongA13$a3y0<-keeplongA13$a3y0*(-1)
        A13.outline.ys<-A13.outline.ys*(-1)
        
        add.back<-max(unlist(c(A13.outline.ys)),na.rm = TRUE)
        add.back<-(-1*add.back)
        
        keeplongA13$true.y<-keeplongA13$true.y+add.back
        keeplongA13$a1y0<-keeplongA13$a1y0+add.back
        keeplongA13$a3y0<-keeplongA13$a3y0+add.back
        A13.outline.ys<-A13.outline.ys+add.back
        
        A13.reference$quad.now.on.left<-A13.reference$ID.left
        A13.reference$quad.now.on.right<-A13.reference$ID.right
      }
      if(A13.reference$TopScreen=="outer") {
        keeplongA13$true.y<-keeplongA13$true.y-480
        keeplongA13$a1y0<-keeplongA13$a1y0-480
        keeplongA13$a3y0<-keeplongA13$a3y0-480
        A13.outline.ys<-A13.outline.ys-480
        
        keeplongA13$true.x<-keeplongA13$true.x*(-1)
        keeplongA13$a1x0<-keeplongA13$a1x0*(-1)
        keeplongA13$a3x0<-keeplongA13$a3x0*(-1)
        A13.outline.xs<-A13.outline.xs*(-1)
        
        add.back.x<-min(unlist(c(A13.outline.xs)),na.rm = TRUE)
        add.back.x<-(-1*add.back.x)
        
        add.back<-max(unlist(c(keeplongA13$true.y,keeplongA13$a1y0,keeplongA13$a3y0,A13.outline.ys)),na.rm = TRUE)
        add.back<-(-1*add.back)
        
        keeplongA13$true.y<-keeplongA13$true.y+add.back
        keeplongA13$a1y0<-keeplongA13$a1y0+add.back
        keeplongA13$a3y0<-keeplongA13$a3y0+add.back
        A13.outline.ys<-A13.outline.ys+add.back
        
        keeplongA13$true.x<-keeplongA13$true.x+add.back.x
        keeplongA13$a1x0<-keeplongA13$a1x0+add.back.x
        keeplongA13$a3x0<-keeplongA13$a3x0+add.back.x
        A13.outline.xs<-A13.outline.xs+add.back.x
        
        A13.reference$quad.now.on.left<-A13.reference$ID.right
        A13.reference$quad.now.on.right<-A13.reference$ID.left
        
        A13.reference$sex.left<-MetaTrial$sex[which(A13.reference$ID.right==MetaTrial$X8x8_quad)]
        A13.reference$sex.right<-MetaTrial$sex[which(A13.reference$ID.left==MetaTrial$X8x8_quad)]
      }
      ##############################################
      #A24
      #A24.outline.xs<-A24.reference[,c(10,8,14,16,6,12,10,8,6)]
      #A24.outline.ys<-A24.reference[,c(11,9,15,17,7,13,11,9,7)]
      A24.outline.xs<-c(A24.left$X,A24.right$X)
      A24.outline.ys<-c(A24.left$Y,A24.right$Y)
      
      
      if(A24.reference$TopScreen=="outer") {
        keeplongA24$true.y<-keeplongA24$true.y*(-1)
        keeplongA24$a2y0<-keeplongA24$a2y0*(-1)
        keeplongA24$a4y0<-keeplongA24$a4y0*(-1)
        A24.outline.ys<-A24.outline.ys*(-1)
        
        add.back<-min(unlist(c(A24.outline.ys)),na.rm = TRUE)
        add.back<-add.back*(-1)
        
        keeplongA24$true.y<-keeplongA24$true.y+add.back
        keeplongA24$a2y0<-keeplongA24$a2y0+add.back
        keeplongA24$a4y0<-keeplongA24$a4y0+add.back
        A24.outline.ys<-A24.outline.ys+add.back
        
        A24.reference$quad.now.on.left<-A24.reference$ID.left
        A24.reference$quad.now.on.right<-A24.reference$ID.right
      }
      if(A24.reference$TopScreen=="inner") {
        keeplongA24$true.y<-keeplongA24$true.y-480
        keeplongA24$a2y0<-keeplongA24$a2y0-480
        keeplongA24$a4y0<-keeplongA24$a4y0-480
        A24.outline.ys<-A24.outline.ys-480
        
        keeplongA24$true.x<-keeplongA24$true.x*(-1)
        keeplongA24$a2x0<-keeplongA24$a2x0*(-1)
        keeplongA24$a4x0<-keeplongA24$a4x0*(-1)
        A24.outline.xs<-A24.outline.xs*(-1)
        
        add.back.x<-min(unlist(c(A24.outline.xs)),na.rm = TRUE)
        add.back.x<-(-1*add.back.x)
        
        add.back<-min(unlist(c(A24.outline.ys)),na.rm = TRUE)
       # min(unlist(c(keeplongA24$true.y,keeplongA24$a2y0,keeplongA24$a4y0,A24.outline.ys)),na.rm = TRUE)
        add.back<-(-1*add.back)
        
        keeplongA24$true.y<-keeplongA24$true.y+add.back
        keeplongA24$a2y0<-keeplongA24$a2y0+add.back
        keeplongA24$a4y0<-keeplongA24$a4y0+add.back
        A24.outline.ys<-A24.outline.ys+add.back
        
        keeplongA24$true.x<-keeplongA24$true.x+add.back.x
        keeplongA24$a2x0<-keeplongA24$a2x0+add.back.x
        keeplongA24$a4x0<-keeplongA24$a4x0+add.back.x
        A24.outline.xs<-A24.outline.xs+add.back.x
        
      }

      #Centering
      if(centeradjust=="yes"){
        if(A13.reference$TopScreen=="inner"){
          a13x.c<-A13.reference$botdivide.x
        } else {
          a13x.c<-A13.reference$topdivide.x
        }
        
        if(A24.reference$TopScreen=="inner"){
          a24x.c<-A24.reference$topdivide.x
        } else {
          a24x.c<-A24.reference$botdivide.x
        }
        
        xadjust<-a24x.c-a13x.c
        keeplongA24$true.x<-keeplongA24$true.x-xadjust
        keeplongA24$a2x0<-keeplongA24$a2x0-xadjust
        keeplongA24$a4x0<-keeplongA24$a4x0-xadjust
        A24.outline.xs<-A24.outline.xs-xadjust
      }
      
      ###################################
      ###################################
      ###################################
      ###################################
      ###################################
      ###################################
      #############################################
      ######
      #
      ############
      ######
      ######################
      #############
      ########
      ##
      #
      keeplongA13<-merge(keeplongA13,IRdata,by.x="dailysecond",by.y="ds",all.x=TRUE)
      keeplongA24<-merge(keeplongA24,IRdata,by.x="dailysecond",by.y="ds",all.x=TRUE)
      peeroller<-merge(peeroller,IRdata,by.x="dailysecond",by.y="ds",all.x=TRUE)
      
      ########################
      #################
      ############
      #######
      thisday<-list(TrialID,p.day,keeplongA13,keeplongA24,peeroller,A13.reference,A24.reference,
                    A13.outline.xs,A13.outline.ys,A24.outline.xs,A24.outline.ys,unadultered.coords.A13,unadultered.coords.A24)
      names(thisday)<-c("TrialID","Day","keeplongA13","keeplongA24","peeroller","A13.reference","A24.reference",
                        "A13.outline.xs","A13.outline.ys","A24.outline.xs","A24.outline.ys",
                        "unadultered.coords.A13","unadultered.coords.A24")
      
      daylist[[dayflag]]<-thisday
      names(daylist)[dayflag]<-p.day
   
      # summary(A1data)
      # temp<-(as.data.frame(summary(A1data$a1roi)));colnames(temp)<-"n"
      # temp[order(temp$n),,drop=FALSE]
      # summary(keeplongA24)
 
      dayflag<-dayflag+1
      
      rm(Camera.24.day);rm(Camera.13.day)
    }#closes ifelse, if data for both cameras for this day NOT are missing
    
      
      
      
      
      if(exists("a1.track")){rm(a1.track)}
      if(exists("a2.track")){rm(a2.track)}
      if(exists("a3.track")){rm(a3.track)}
      if(exists("a4.track")){rm(a4.track)}
      
      if(exists("Camera.24.day")){rm(Camera.24.day)}
      if(exists("Camera.13.day")){rm(Camera.13.day)}
      
      if(exists("A1trj")){rm(A1trj)}
      if(exists("A2trj")){rm(A2trj)}
      if(exists("A3trj")){rm(A3trj)}
      if(exists("A4trj")){rm(A4trj)}


   } #day ---wrap
  alltrials[[trialflag]]<-daylist
  names(alltrials)[trialflag]<-TrialID
  trialflag<-trialflag+1
  
  save(daylist,file=paste(TrialID,"_2plot2analyze_relax.Rdata",sep=''))
  
  
}#trial ENDOFTRIALLOOP




library(summarytools)
view(dfSummary(roimismatchframe))

roimismatchframe->temper
temper<-temper[(grep("T001",temper$trial)),]
temper$trial<-factor(temper$trial)
view(dfSummary(temper))
temper$mismatches<-paste(temper$ralQuad,temper$roiID,sep=':')
hooke<-strsplit(as.character(temper$roiID),"_")
temper$mismatchtype<-paste(temper$ralQuad,unlist(lapply(hooke, `[[`, 1)),sep=':')



temper$mismatchtype<-factor(temper$mismatchtype)
summary(temper)
########################################

alltrials<-list()
completedprocessed<-list.files(pattern="2plot2analyze.Rdata")
for(p in 1:length(completedprocessed)){
  nomen<-AllTrialIDs[p]
  load(completedprocessed[p])
  alltrials[[p]]<-daylist
  rm(daylist)
}
names(alltrials)<-AllTrialIDs[c(1:17)]

####
dir<-"C:/Users/Rusty/Amazon Drive/MICE/Thermal/Next/8x8PeeData/PeeDetected"
setwd(dir)
#COOrdinates for polygons defining 'hard' boundaries for pee detection
coordinatefiles<-list.files("C:/Users/Rusty/Amazon Drive/MICE/Thermal/Arena_coordinates",full.names=TRUE,pattern="D1.")#creates list (2 items) with the megaframe, and its rowname file









if(makeplots=="yes"){
  library(svMisc)
  for(eachtrial in 1:length(alltrials)){
    
    alltrialdays<-alltrials[[eachtrial]]

    
    for(eachday in 1:length(alltrialdays)){
      
        ThisDayTrial<-alltrialdays[[eachday]]
      
        TrialID<-ThisDayTrial$TrialID
        p.day<-ThisDayTrial$Day
        peeroller<-ThisDayTrial$peeroller
        keeplongA13<-ThisDayTrial$keeplongA13
        keeplongA24<-ThisDayTrial$keeplongA24
        A13.reference<-ThisDayTrial$A13.reference
        A24.reference<-ThisDayTrial$A24.reference
        A13.outline.xs<-ThisDayTrial$A13.outline.xs
        A24.outline.xs<-ThisDayTrial$A24.outline.xs
        
        #A13.outline.ys<-ThisDayTrial$A13.outline.ys-max(ThisDayTrial$A13.outline.ys,na.rm = TRUE)
        A13.outline.ys<-ThisDayTrial$A13.outline.ys
       
        #A24.outline.ys<-  ThisDayTrial$A24.outline.ys-min(ThisDayTrial$A24.outline.ys,na.rm = TRUE)
        A24.outline.ys<-  ThisDayTrial$A24.outline.ys
        

        
        ############
        TheseCoords<-coordinatefiles[grepl(TrialID, coordinatefiles)]
        CoordsA13<-TheseCoords[grepl("A13", TheseCoords)]
        CoordsA24<-TheseCoords[grepl("A24", TheseCoords)]
        

        if(length(CoordsA13)>0){
          A13.left<-read.csv(CoordsA13[grep("L.c",CoordsA13)])
          A13.right<-read.csv(CoordsA13[grep("R.c",CoordsA13)])  
        } else {
          skipA13<-"yes"
        }
        
        if(length(CoordsA24)>0){
          A24.left<-read.csv(CoordsA24[grep("L.c",CoordsA24)])
          A24.right<-read.csv(CoordsA24[grep("R.c",CoordsA24)])
        } else {
          skipA24<-"yes"
        }
        
        ###################
        
        
        ###check/make directories for saving plot images
        whereplotgoes<-paste(TrialID,p.day,sep='.')
        newDir<-paste(vizdir,whereplotgoes,sep = "/")
        if (dir.exists(newDir)){
          setwd(newDir)
        } else {
          dir.create(newDir)
          setwd(newDir)
        }
        
        
        if(length(setdiff(keeplongA13$dailysecond,keeplongA24$dailysecond))==0){
          Just.Frames<-unique(keeplongA24$dailysecond)
        } else {
          print("Frames are not the same for A13 and A24")
        }
        
        
        ######################################################################################
        
        a1flag<-0;a1.pee.flag<-0;a1waterflag<-0
        a2flag<-0;a2.pee.flag<-0;a2waterflag<-0
        a3flag<-0;a3.pee.flag<-0;a3waterflag<-0
        a4flag<-0;a4.pee.flag<-0;a4waterflag<-0
        
        a1water.growing<-peeroller[1,];a1water.growing$dailysecond<-1
        a2water.growing<-peeroller[1,];a2water.growing$dailysecond<-1
        a3water.growing<-peeroller[1,];a3water.growing$dailysecond<-1
        a4water.growing<-peeroller[1,];a4water.growing$dailysecond<-1
        
        if(exists("a1.pee.build")){rm(a1.pee.build)};if(exists("a2.pee.build")){rm(a2.pee.build)};
        if(exists("a3.pee.build")){rm(a3.pee.build)};if(exists("a4.pee.build")){rm(a4.pee.build)}
        if(exists("time.series")){rm(time.series)}
        
        maxspeed<-max(peeroller[,c("a1.rollspeed","a2.rollspeed","a3.rollspeed","a4.rollspeed")],na.rm = TRUE)
        indiv.max.cum.disttrav<-data.frame(t(data.frame(colSums(peeroller[,c("a1.displacement","a2.displacement","a3.displacement","a4.displacement")],na.rm=TRUE))))
        maxdist<-max(indiv.max.cum.disttrav)
        
        
        maxcumpee<-max(colSums(peeroller[,c("a1.pee.vol","a2.pee.vol","a3.pee.vol","a4.pee.vol")],na.rm=TRUE)) 
        indiv.max.cum.pee<-data.frame(t(data.frame(colSums(peeroller[,c("a1.pee.vol","a2.pee.vol","a3.pee.vol","a4.pee.vol")],na.rm=TRUE))))
        indiv.max.instant.pee<-t(data.frame(c(max(peeroller[,c("a1.roll.pee")]),max(peeroller[,"a2.roll.pee"]),
                                              max(peeroller[,"a3.roll.pee"]),max(peeroller[,"a4.roll.pee"]))))
        colnames(indiv.max.instant.pee)<-c("a1.roll.pee","a2.roll.pee","a3.roll.pee","a4.roll.pee")
        indiv.max.instant.pee<-data.frame(indiv.max.instant.pee)
        
        xrange<-c(min(Just.Frames),max(Just.Frames))
        topout<-max(Just.Frames)
        teddyflag=1
        climbflag=1
        
        startvalue.13.1<-which((keeplongA13[!is.na(keeplongA13$a1x0),])[1,"dailysecond"]==Just.Frames)
        startvalue.13.3<-which((keeplongA13[!is.na(keeplongA13$a3x0),])[1,"dailysecond"]==Just.Frames)
        startvalue.13<-ifelse(startvalue.13.1<=startvalue.13.3,startvalue.13.1,startvalue.13.3)
        
        startvalue.24.2<-which((keeplongA24[!is.na(keeplongA24$a2x0),])[1,"dailysecond"]==Just.Frames)
        startvalue.24.4<-which((keeplongA24[!is.na(keeplongA24$a4x0),])[1,"dailysecond"]==Just.Frames)
        startvalue.24<-ifelse(startvalue.24.2<=startvalue.24.4,startvalue.24.2,startvalue.24.4)
        
        startvalue<-ifelse(startvalue.13<=startvalue.24,startvalue.13,startvalue.24)
        
        
        print(paste0((length(Just.Frames)-startvalue)," Frames to create for ",TrialID, ".",p.day))
        
        
        plotonlysummary<-"yes"
        
        if(plotonlysummary=="yes"){
          teddy=length(Just.Frames)
            #for(teddy in startvalue:(startvalue+10)){
            
            
           # txtProgressBar(min=startvalue,max=length(Just.Frames),initial=teddy,width = 80,style=3)
            
            
            
            
            current.roll<-peeroller[which(peeroller$dailysecond==Just.Frames[teddy]),]
            all.current.roll<-peeroller
            
            
            

            currentframe.A13<-keeplongA13[which(keeplongA13$dailysecond==Just.Frames[teddy]),]  
            #if all rows in current frame are equal (byproduct of cleaning data using location data)
            if(nrow(currentframe.A13)>1){
              #check if all columns in all rows of currentframe.A13 are the same, if so collapse
              if(all(apply(currentframe.A13, 2, function(x) length(unique(x)) == 1) == TRUE)) {
                currentframe.A13<-currentframe.A13[1,]
              }
            }
            currentframe.A24<-keeplongA24[which(keeplongA24$dailysecond==Just.Frames[teddy]),]  
            
            #if all rows in current frame are equal (byproduct of cleaning data using location data)
            if(nrow(currentframe.A24)>1){
              #check if all columns in all rows of currentframe.A13 are the same, if so collapse
              if(all(apply(currentframe.A24, 2, function(x) {length(unique(x))==1}) == TRUE)) {
                currentframe.A24<-currentframe.A24[1,]
              }
            }
            
            setwd(vizdir)
            
            png(paste("SUMMARY_",TrialID,".",p.day,".png",sep=''),bg="transparent",width=1200,height=800)#use if making transparent pngs
            
            #Set up plot locations/margins
            par(oma=c(1,1,1,1))
            par(mar = c(0, 0, 0, 0))
            
            layout(matrix(c(20,10,8,8,8,8,11,11,11,11,21,13,
                            20,10,9,9,9,9,12,12,12,12,21,13,
                            16,16,16,rep(1,6),17,17,17,
                            16,16,16,rep(1,6),17,17,17,
                            16,16,16,rep(1,6),17,17,17,
                            14,14,14,rep(1,6),15,15,15,
                            14,14,14,rep(1,6),15,15,15,
                            14,14,14,rep(1,6),15,15,15,
                            
                            18,4,2,2,2,2,5,5,5,5,19,7,
                            18,4,3,3,3,3,6,6,6,6,19,7), 10, 12, byrow = TRUE),
                   widths=rep(1,10), heights=rep(1,10))
            
            ######################
            # SETUP ARENA PLOT
            # SETUP ARENA PLOT
            real.x.min<-min(c(keeplongA24$a4x0,keeplongA24$a2x0,keeplongA24$true.x),na.rm=TRUE)
            xlimspacing<-real.x.min-30
            plot(NA, xlim=c(xlimspacing,670),
                 ylim=c(min(A13.outline.ys),max(A24.outline.ys)),
                 xlab='',ylab='',
                 axes=FALSE)
            
            
            ######A13
            #A13.left
            a13shape.left<-cbind(A13.outline.xs,A13.outline.ys)[c(1:nrow(A13.left)),];a13shape.left<-rbind(a13shape.left,a13shape.left[1,])
            a13shape.right<-cbind(A13.outline.xs,A13.outline.ys)[c((1+nrow(A13.left)):length(A13.outline.xs)),];a13shape.right<-rbind(a13shape.right,a13shape.right[1,])
            lines(a13shape.left)
            lines(a13shape.right)
            
            a13leftmousecoords<-getpositionalvalues(a13shape.left)
            a13rightmousecoords<-getpositionalvalues(a13shape.right)
            
            if(A13.reference$TopScreen=="inner"){
              text(min(a13shape.left[,1])-18,min(a13shape.left[,2]),A13.reference$quad.now.on.left,col="blue",cex=1.4)
              if(A13.reference$sex.left=="m"){
                text(min(a13shape.left[,1])-22,min(a13shape.left[,2])+60,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
              } else {
                text(min(a13shape.left[,1])-22,min(a13shape.left[,2])+60,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
              }
              
              
              text(max(a13shape.right[,1])+18,min(a13shape.right[,2]),A13.reference$quad.now.on.right,col="red",cex=1.4)
              if(A13.reference$sex.right=="m"){
                text(max(a13shape.right[,1])+22,min(a13shape.right[,2])+60,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
              } else {
                text(max(a13shape.right[,1])+22,min(a13shape.right[,2])+60,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
              }
            }
            if(A13.reference$TopScreen=="outer"){
              text(max(a13shape.left[,1])-18,max(a13shape.left[,2]),A13.reference$quad.now.on.left,col="blue",cex=1.4)
              if(A13.reference$sex.left=="m"){
                text(max(a13shape.left[,1])-22,max(a13shape.left[,2])+60,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
              } else {
                text(max(a13shape.left[,1])-22,max(a13shape.left[,2])+60,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
              }
              
              text(min(a13shape.right[,1])+18,max(a13shape.right[,2]),A13.reference$quad.now.on.right,col="red",cex=1.4)
              if(A13.reference$sex.right=="m"){
                text(min(a13shape.right[,1])+22,max(a13shape.right[,2])+60,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
              } else {
                text(min(a13shape.right[,1])+22,max(a13shape.right[,2])+60,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
              }
            }
            
            
            ####A24
            a24shape.left<-cbind(A24.outline.xs,A24.outline.ys)[c(1:nrow(A24.left)),];a24shape.left<-rbind(a24shape.left,a24shape.left[1,])
            a24shape.right<-cbind(A24.outline.xs,A24.outline.ys)[c((1+nrow(A24.left)):length(A24.outline.xs)),];a24shape.right<-rbind(a24shape.right,a24shape.right[1,])
            lines(a24shape.left)
            lines(a24shape.right)
            
            a24leftmousecoords<-getpositionalvalues(a24shape.left)
            a24rightmousecoords<-getpositionalvalues(a24shape.right)
            
            if(A24.reference$TopScreen=="outer"){
              text(min(a24shape.left[,1])-18,max(a24shape.left[,2]),A24.reference$quad.now.on.left,col="darkgreen",cex=1.4)
              if(A24.reference$sex.left=="m"){
                text(min(a24shape.left[,1])-18,max(a24shape.left[,2])-70,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
              } else {
                text(min(a24shape.left[,1])-18,max(a24shape.left[,2])-70,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
              }
              
              text(max(a24shape.right[,1])+18,max(a24shape.right[,2]),A24.reference$quad.now.on.right,col="darkgoldenrod3",cex=1.4)
              if(A24.reference$sex.right=="m"){
                text(max(a24shape.right[,1])+18,max(a24shape.right[,2])-70,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
              } else {
                text(max(a24shape.right[,1])+18,max(a24shape.right[,2])-70,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
              }
            }
            if(A24.reference$TopScreen=="inner"){
              text(max(a24shape.left[,1])-18,min(a24shape.left[,2]),A24.reference$quad.now.on.left,col="darkgreen",cex=1.4)
              if(A24.reference$sex.left=="m"){
                text(max(a24shape.left[,1])-18,min(a24shape.left[,2])-70,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
              } else {
                text(max(a24shape.left[,1])-18,min(a24shape.left[,2])-70,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
              }
              
              text(min(a24shape.right[,1])+18,min(a24shape.right[,2]),A24.reference$quad.now.on.right,col="darkgoldenrod3",cex=1.4)
              if(A24.reference$sex.right=="m"){
                text(min(a24shape.right[,1])+18,min(a24shape.right[,2])-70,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
              } else {
                text(min(a24shape.right[,1])+18,min(a24shape.right[,2])-70,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
              }
            }
            
            
            td <- seconds_to_period(Just.Frames[teddy])
            
            t.x<-(real.x.min-30)
            real.x.max<-max(c(keeplongA24$a4x0,keeplongA24$a2x0,keeplongA24$true.x),na.rm=TRUE)
            tt.x<-(real.x.max+10)
            
            text(t.x,-5,sprintf('%02d:%02d:%02d', td@hour, minute(td), second(td)),srt=90,cex = 2, col="black")
            
            text(tt.x,-5,Just.Frames[teddy],srt=90,cex = 1.8, col="darkgray")
            
            ####################################################################
            #A1
            ############
            for(vero in 1:1){
              
              check<-currentframe.A13[which(currentframe.A13$quad=="a1"),]
              a1.pee.build<-keeplongA13[which(keeplongA13$quad=="a1"),] 
              

                par(new=TRUE)
                plot(a1.pee.build$true.y~a1.pee.build$true.x,cex=log(a1.pee.build$Clustsize),
                     xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),
                     bty="o",
                     pch=16,
                     xlab='',ylab='',
                     axes=FALSE,
                     col=add.alpha("blue",0.09))  

                par(new=TRUE)

            }
            
            #A3
            ############
            for(vero in 1:1){
              
              check.3<-currentframe.A13[which(currentframe.A13$quad=="a3"),]
              a3.pee.build<-keeplongA13[which(keeplongA13$quad=="a3"),] 
  
                par(new=TRUE)
                plot(a3.pee.build$true.y~a3.pee.build$true.x,cex=log(a3.pee.build$Clustsize),
                     xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),bty="o",
                     pch=16,
                     xlab='',ylab='',
                     axes=FALSE,
                     col=add.alpha("darkred",0.09))  
            }
            
            #A2
            ############
            for(vero in 1:1){
              
              check.2<-currentframe.A24[which(currentframe.A24$quad=="a2"),]
              a2.pee.build<-keeplongA24[which(keeplongA24$quad=="a2"),] 
             
                par(new=TRUE)
                plot(a2.pee.build$true.y~a2.pee.build$true.x,cex=log(a2.pee.build$Clustsize),
                     xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),bty="o",
                     pch=16,
                     xlab='',ylab='',
                     axes=FALSE,
                     col=add.alpha("darkgreen",0.09))  
             }
            
            #A4
            ############
            for(vero in 1:1){
              
              check.4<-currentframe.A24[which(currentframe.A24$quad=="a4"),]
              a4.pee.build<-keeplongA24[which(keeplongA24$quad=="a4"),] 
                
                par(new=TRUE)
                plot(a4.pee.build$true.y~a4.pee.build$true.x,cex=log(a4.pee.build$Clustsize),
                     xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),bty="o",
                     pch=16,
                     xlab='',ylab='',
                     axes=FALSE,
                     col=add.alpha("darkgoldenrod3",0.09))  

            }
            
            ##############################################

              time.series<-peeroller
              

            
            ##################
            #1
            #A1 time-series
            #rolling
            par(mar = c(0, 1, 0, 1))
            
            plot(time.series$dailysecond,time.series$a1.roll.pee,
                 type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a1.roll.pee),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="dodgerblue")
            rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
            if(!is.na(current.roll$a1roi)){
              if(current.roll$a1roi=="a1_a3_water"|current.roll$a1roi=="a1_a2_water"){
                if(a1waterflag==0){
                  a1water.growing<-current.roll
                  a1waterflag<-a1waterflag+1
                } else {
                  a1water.growing<-rbind(a1water.growing,current.roll)
                }
              }
            }
            axis(side=3,at=a1water.growing$dailysecond,lwd=.8,tick=TRUE,labels=rep("",nrow(a1water.growing)),col="darkorchid1")
            par(new=TRUE)
            
            plot(time.series$dailysecond,time.series$a1.roll.pee,
                 type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a1.roll.pee),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="dodgerblue")
            par(new=TRUE)
            #cumulative
            plot(time.series$dailysecond,time.series$a1cumpee, 
                 type="l", xlim=xrange,ylim=c(0,indiv.max.cum.pee$a1.pee.vol),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,col="darkblue") 
            
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,indiv.max.cum.pee$a1.pee.vol),side=2,at=c(0,indiv.max.cum.pee$a1.pee.vol/2,indiv.max.cum.pee$a1.pee.vol),labels=c(0,0.5,1))
            title(ylab="Urine",line=-1.6,cex.lab=2)
            #title(xlab="Time",line=-1)
            
            ##################
            #VELOCITY/MOVEMENT
            par(mar = c(0, 1, 0, 1))
            plot(time.series$dailysecond,time.series$a1.rollspeed,
                 type="l", xlim=xrange,ylim=c(0,maxspeed),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="dodgerblue")
            par(new=TRUE)
            #cumulative
            plot(time.series$dailysecond,time.series$a1cumdist, 
                 type="l", xlim=xrange,ylim=c(0,indiv.max.cum.disttrav$a1.displacement),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,col="darkblue") 
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,indiv.max.cum.disttrav$a1.displacement),side=2,at=c(0,indiv.max.cum.disttrav$a1.displacement/2,indiv.max.cum.disttrav$a1.displacement),labels=c(0,0.5,1))
            title(ylab="Velocity",line=-1.6,cex.lab=2)
            title(xlab="Time",line=-1)
            
            par(mar = c(0, 2, 1, 1))
            barplot(current.roll$a1cumpee,ylim=c(0,maxcumpee),
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="blue")
            rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
            par(new=TRUE)
            barplot(current.roll$a1cumpee,ylim=c(0,maxcumpee),
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="blue")
            text(x=.75,y=maxcumpee/2,labels = "Total pee",cex=3,col="yellowgreen",srt=90)
            
            
            
            ###########################################################################
            #1
            #A3 time-series
            #rolling
            par(mar = c(0, 1, 0, 1))
            plot(time.series$dailysecond,time.series$a3.roll.pee,
                 type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a3.roll.pee),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="pink")
            rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
            if(!is.na(current.roll$a3roi)){
              if(current.roll$a3roi=="a3_a1_water"|current.roll$a3roi=="a3_a4_water"){
                if(a3waterflag==0){
                  a3water.growing<-current.roll
                  a3waterflag<-a3waterflag+1
                } else {
                  a3water.growing<-rbind(a3water.growing,current.roll)
                }
              }
            }           
            axis(side=3,at=a3water.growing$dailysecond,lwd=.8,tick=TRUE,labels=rep("",nrow(a3water.growing)),col="darkorchid1")
            
            
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a3.roll.pee,
                 type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a3.roll.pee),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="pink")
            par(new=TRUE)
            #cumulative
            plot(time.series$dailysecond,time.series$a3cumpee, 
                 type="l", xlim=xrange,ylim=c(0,indiv.max.cum.pee$a3.pee.vol),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,col="deeppink") 
            
            
            
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,indiv.max.cum.pee$a3.pee.vol),side=2,at=c(0,indiv.max.cum.pee$a3.pee.vol/2,indiv.max.cum.pee$a3.pee.vol),labels=c(0,0.5,1))
            title(ylab="Urine",line=-1.6,cex.lab=2)
            #title(xlab="Time",line=-1)
            
            ##################
            #VELOCITY/MOVEMENT
            par(mar = c(0, 1, 0, 1))
            plot(time.series$dailysecond,time.series$a3.rollspeed,
                 type="l", xlim=xrange,ylim=c(0,maxspeed),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="pink")
            par(new=TRUE)
            #cumulative
            plot(time.series$dailysecond,time.series$a3cumdist, 
                 type="l", xlim=xrange,ylim=c(0,indiv.max.cum.disttrav$a3.displacement),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,col="deeppink") 
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,indiv.max.cum.disttrav$a3.displacement),side=2,at=c(0,indiv.max.cum.disttrav$a3.displacement/2,indiv.max.cum.disttrav$a3.displacement),labels=c(0,0.5,1))
            title(ylab="Velocity",line=-1.6,cex.lab=2)
            title(xlab="Time",line=-1)
            
            par(mar = c(0, 2, 1, 1))
            barplot(current.roll$a3cumpee,ylim=c(0,maxcumpee),
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="darkred")
            rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
            par(new=TRUE)
            barplot(current.roll$a3cumpee,ylim=c(0,maxcumpee),
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="darkred")
            text(x=.75,y=maxcumpee/2,labels = "Total pee",cex=3,col="yellowgreen",srt=90)
            
            
            ###########################################################################
            #1
            #A2 time-series
            #rolling
            par(mar = c(0, 1, 0, 1))
            plot(time.series$dailysecond,time.series$a2.roll.pee,
                 type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a2.roll.pee),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="green")
            rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
            if(!is.na(current.roll$a2roi)){
              if(current.roll$a2roi=="a2_a1_water"|current.roll$a2roi=="a2_a4_water"){
                if(a2waterflag==0){
                  a2water.growing<-current.roll
                  a2waterflag<-a2waterflag+1
                } else {
                  a2water.growing<-rbind(a2water.growing,current.roll)
                }
              }
            }
            axis(side=3,at=a2water.growing$dailysecond,lwd=.8,tick=TRUE,labels=rep("",nrow(a2water.growing)),col="darkorchid1")
            
            
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a2.roll.pee,
                 type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a2.roll.pee),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="green")
            par(new=TRUE)
            #cumulative
            plot(time.series$dailysecond,time.series$a2cumpee, 
                 type="l", xlim=xrange,ylim=c(0,indiv.max.cum.pee$a2.pee.vol),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,col="darkgreen") 
            
            
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,indiv.max.cum.pee$a2.pee.vol),side=2,at=c(0,indiv.max.cum.pee$a2.pee.vol/2,indiv.max.cum.pee$a2.pee.vol),labels=c(0,0.5,1))
            title(ylab="Urine",line=-1.6,cex.lab=2)
            #title(xlab="Time",line=-1)
            
            ##################
            #VELOCITY/MOVEMENT
            par(mar = c(0, 1, 0, 1))
            plot(time.series$dailysecond,time.series$a2.rollspeed,
                 type="l", xlim=xrange,ylim=c(0,maxspeed),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="green")
            par(new=TRUE)
            #cumulative
            plot(time.series$dailysecond,time.series$a2cumdist, 
                 type="l", xlim=xrange,ylim=c(0,indiv.max.cum.disttrav$a2.displacement),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,col="darkgreen") 
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,indiv.max.cum.disttrav$a2.displacement),side=2,at=c(0,indiv.max.cum.disttrav$a2.displacement/2,indiv.max.cum.disttrav$a2.displacement),labels=c(0,0.5,1))
            title(ylab="Velocity",line=-1.6,cex.lab=2)
            title(xlab="Time",line=-1)
            
            par(mar = c(0, 2, 1, 1))
            barplot(current.roll$a2cumpee,ylim=c(0,maxcumpee),
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="darkgreen")
            rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
            par(new=TRUE)
            barplot(current.roll$a2cumpee,ylim=c(0,maxcumpee),
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="darkgreen")
            text(x=.75,y=maxcumpee/2,labels = "Total pee",cex=3,col="yellowgreen",srt=90)
            
            ###########################################################################
            #1
            #A4 time-series
            #rolling
            par(mar = c(0, 1, 0, 1))
            plot(time.series$dailysecond,time.series$a4.roll.pee,
                 type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a4.roll.pee),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="darkgoldenrod3")
            rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
            if(!is.na(current.roll$a4roi)){
              if(current.roll$a4roi=="a4_a3_water"|current.roll$a4roi=="a4_a2_water"){
                if(a4waterflag==0){
                  a4water.growing<-current.roll
                  a4waterflag<-a4waterflag+1
                } else {
                  a4water.growing<-rbind(a4water.growing,current.roll)
                }
              }
            }
            axis(side=3,at=a4water.growing$dailysecond,lwd=.8,tick=TRUE,labels=rep("",nrow(a4water.growing)),col="darkorchid1")
            
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a4.roll.pee,
                 type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a4.roll.pee),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="darkgoldenrod3")
            par(new=TRUE)
            #cumulative
            plot(time.series$dailysecond,time.series$a4cumpee, 
                 type="l", xlim=xrange,ylim=c(0,indiv.max.cum.pee$a4.pee.vol),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,col="darkorange") 
            
            
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,indiv.max.cum.pee$a4.pee.vol),side=2,at=c(0,indiv.max.cum.pee$a4.pee.vol/2,indiv.max.cum.pee$a4.pee.vol),labels=c(0,0.5,1))
            title(ylab="Urine",line=-1.6,cex.lab=2)
            #title(xlab="Time",line=-1)
            
            ##################
            #VELOCITY/MOVEMENT
            par(mar = c(0, 1, 0, 1)) #c(bottom, left, top, right)
            plot(time.series$dailysecond,time.series$a4.rollspeed,
                 type="l", xlim=xrange,ylim=c(0,maxspeed),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",lty=1,
                 lwd=1.5,col="darkgoldenrod3")
            par(new=TRUE)
            #cumulative
            plot(time.series$dailysecond,time.series$a4cumdist, 
                 type="l", xlim=xrange,ylim=c(0,indiv.max.cum.disttrav$a4.displacement),
                 pch=16,xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,col="darkorange") 
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,indiv.max.cum.disttrav$a4.displacement),side=2,at=c(0,indiv.max.cum.disttrav$a4.displacement/2,indiv.max.cum.disttrav$a4.displacement),labels=c(0,0.5,1))
            title(ylab="Velocity",line=-1.6,cex.lab=2)
            title(xlab="Time",line=-1)
            
            par(mar = c(0, 2, 1, 1))
            barplot(current.roll$a4cumpee,ylim=c(0,maxcumpee),
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="darkorange")
            rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
            par(new=TRUE)
            barplot(current.roll$a4cumpee,ylim=c(0,maxcumpee),
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="darkorange")
            text(x=.75,y=maxcumpee/2,labels = "Total pee",cex=3,col="yellowgreen",srt=90)
            
            
            
            if(wingplots=="stacked"){
              #Stacked Plots
              
              #a1
              h.a1<-as.matrix(c(current.roll$a1_water,
                                current.roll$a1_a3_barrier,
                                current.roll$a1_central_corner,
                                current.roll$a1_a2_barrier))
              H <- apply(h.a1, 2L, cumsum)
              H <- H - h.a1 / 2
              par(mar = c(2, 2, 1, 1))#c(bottom, left, top, right)
              x.a1<-barplot(h.a1,
                            col=c("darkorchid1","darkred","gray40","darkgreen"),
                            ylim=c(0,sum(c(current.roll$a1_a2_barrier,
                                           current.roll$a1_central_corner,
                                           current.roll$a1_a3_barrier,
                                           current.roll$a1_water))))
              text(rep(x.a1, each = nrow(H)), H, labels = c("water","a3","center","a2"),col="cornsilk")
              title(ylab="pee breakdown",line=-20,cex.lab=2)
              
              #a3
              h.a3<-as.matrix(c(current.roll$a3_water,
                                current.roll$a3_a1_barrier,
                                current.roll$a3_central_corner,
                                current.roll$a3_a4_barrier))
              H <- apply(h.a3, 2L, cumsum)
              H <- H - h.a3 / 2
              par(mar = c(2, 2, 1, 1))#c(bottom, left, top, right)
              x.a3<-barplot(h.a3,
                            col=c("darkorchid1","darkblue","gray40","darkorange"),
                            ylim=c(0,sum(c(current.roll$a3_a4_barrier,
                                           current.roll$a3_central_corner,
                                           current.roll$a3_a1_barrier,
                                           current.roll$a3_water))))
              text(rep(x.a3, each = nrow(H)), H, labels = c("water","a1","center","a4"),col="cornsilk")
              title(ylab="pee breakdown",line=-20,cex.lab=2)
              
              
              #a2
              h.a2<-as.matrix(c(current.roll$a2_a1_barrier,
                                current.roll$a2_central_corner,
                                current.roll$a2_a4_barrier,
                                current.roll$a2_water))
              H <- apply(h.a2, 2L, cumsum)
              H <- H - h.a2 / 2
              par(mar = c(0, 2, 2, 1))#c(bottom, left, top, right)
              x.a2<-barplot(h.a2,
                            col=c("darkblue","gray40","darkorange","darkorchid1"),
                            ylim=c(0,sum(c(current.roll$a2_a1_barrier,
                                           current.roll$a2_central_corner,
                                           current.roll$a2_a4_barrier,
                                           current.roll$a2_water))))
              text(rep(x.a2, each = nrow(H)), H, labels = c("a1","center","a4","water"),col="cornsilk")
              title(ylab="pee breakdown",line=-20,cex.lab=2)
              
              
              #a4
              h.a4<-as.matrix(c(current.roll$a4_a3_barrier,
                                current.roll$a4_central_corner,
                                current.roll$a4_a2_barrier,
                                current.roll$a4_water))
              H <- apply(h.a4, 2L, cumsum)
              H <- H - h.a4 / 2
              par(mar = c(0, 2, 2, 1))#c(bottom, left, top, right)
              
              x.a4<-barplot(h.a4,
                            col=c("darkred","gray40","darkgreen","darkorchid1"),
                            ylim=c(0,sum(c(current.roll$a4_a3_barrier,
                                           current.roll$a4_central_corner,
                                           current.roll$a4_a2_barrier,
                                           current.roll$a4_water))))
              text(rep(x.a4, each = nrow(H)), H, labels = c("a3","center","a2","water"),col="cornsilk")
              title(ylab="pee breakdown",line=-20,cex.lab=2)
            }
            if(wingplots=="climbinglines"){
              
              # if(climbflag==1){
              #   climbstart=teddy
              #   climbflag<-climbflag+1
              # }  
              # 
              # climb.roll<-peeroller[c(climbstart:teddy),]
              
              ######
              #A1
              a1.max.time.spend<-max(peeroller[,c("a1_a3_barrier","a1_a2_barrier","a1_central_corner","a1_water")])
              a1.max.pee.spend<-max(peeroller[,c("a1_a3_barrier.pee","a1_a2_barrier.pee","a1_central_corner.pee","a1_water.pee")])
              
              par(mar = c(1, 1, 1, 1))#c(bottom, left, top, right)
              plot(time.series$dailysecond,time.series$a1_a2_barrier,
                   xlim=xrange,
                   ylim=c(0,a1.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkgreen")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a1_a3_barrier,
                   xlim=xrange,
                   ylim=c(0,a1.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkred")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a1_central_corner,
                   xlim=xrange,
                   ylim=c(0,a1.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="gray40")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a1_water,
                   xlim=xrange,
                   ylim=c(0,a1.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkorchid1")
              par(new=TRUE)
              #################################
              plot(time.series$dailysecond,time.series$a1_a2_barrier.pee,
                   xlim=xrange,
                   ylim=c(0,a1.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkgreen")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a1_a3_barrier.pee,
                   xlim=xrange,
                   ylim=c(0,a1.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkred")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a1_central_corner.pee,
                   xlim=xrange,
                   ylim=c(0,a1.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="gray40")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a1_water.pee,
                   xlim=xrange,
                   ylim=c(0,a1.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkorchid1")
              
              
              Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
              Axis(x=c(0,a1.max.time.spend),side=2,at=c(0,a1.max.time.spend/2,a1.max.time.spend),labels=c(0,0.5,1))
              Axis(x=c(0,a1.max.pee.spend),side=4,at=c(0,a1.max.pee.spend/2,a1.max.pee.spend),labels=c(0,0.5,1))
              
              axis(4,at=a1.max.pee.spend/2,labels="x Pee in place x",pos=topout*.95,cex.axis=2,tick=FALSE)
              title(ylab="Time in place",line=-1.6,cex.lab=2)
              title(xlab="Time",line=-1)
              
              ######
              #A3
              a3.max.time.spend<-max(peeroller[,c("a3_a1_barrier","a3_a4_barrier","a3_central_corner","a3_water")])
              a3.max.pee.spend<-max(peeroller[,c("a3_a1_barrier.pee","a3_a4_barrier.pee","a3_central_corner.pee","a3_water.pee")])
              
              par(mar = c(1, 1, 1, 1))#c(bottom, left, top, right)
              plot(time.series$dailysecond,time.series$a3_a1_barrier,
                   xlim=xrange,
                   ylim=c(0,a3.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkblue")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a3_a4_barrier,
                   xlim=xrange,
                   ylim=c(0,a3.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkorange")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a3_central_corner,
                   xlim=xrange,
                   ylim=c(0,a3.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="gray40")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a3_water,
                   xlim=xrange,
                   ylim=c(0,a3.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkorchid1")
              par(new=TRUE)
              #################################
              plot(time.series$dailysecond,time.series$a3_a1_barrier.pee,
                   xlim=xrange,
                   ylim=c(0,a3.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkblue")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a3_a4_barrier.pee,
                   xlim=xrange,
                   ylim=c(0,a3.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkorange")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a3_central_corner.pee,
                   xlim=xrange,
                   ylim=c(0,a3.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="gray40")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a3_water.pee,
                   xlim=xrange,
                   ylim=c(0,a3.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkorchid1")
              
              
              Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
              Axis(x=c(0,a3.max.time.spend),side=2,at=c(0,a3.max.time.spend/2,a3.max.time.spend),labels=c(0,0.5,1))
              Axis(x=c(0,a3.max.pee.spend),side=4,at=c(0,a3.max.pee.spend/2,a3.max.pee.spend),labels=c(0,0.5,1))
              
              axis(4,at=a3.max.pee.spend/2,labels="x Pee in place x",pos=topout*.95,cex.axis=2,tick=FALSE)
              title(ylab="Time in place",line=-1.6,cex.lab=2)
              title(xlab="Time",line=-1)
              
              ######
              #A2
              
              a2.max.time.spend<-max(peeroller[,c("a2_a1_barrier","a2_a4_barrier","a2_central_corner","a2_water")])
              a2.max.pee.spend<-max(peeroller[,c("a2_a1_barrier.pee","a2_a4_barrier.pee","a2_central_corner.pee","a2_water.pee")])
              
              par(mar = c(1, 1, 1, 1))#c(bottom, left, top, right)
              plot(time.series$dailysecond,time.series$a2_a1_barrier,
                   xlim=xrange,
                   ylim=c(0,a2.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkblue")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a2_a4_barrier,
                   xlim=xrange,
                   ylim=c(0,a2.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkorange")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a2_central_corner,
                   xlim=xrange,
                   ylim=c(0,a2.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="gray40")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a2_water,
                   xlim=xrange,
                   ylim=c(0,a2.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkorchid1")
              par(new=TRUE)
              #################################
              plot(time.series$dailysecond,time.series$a2_a1_barrier.pee,
                   xlim=xrange,
                   ylim=c(0,a2.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkblue")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a2_a4_barrier.pee,
                   xlim=xrange,
                   ylim=c(0,a2.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkorange")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a2_central_corner.pee,
                   xlim=xrange,
                   ylim=c(0,a2.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="gray40")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a2_water.pee,
                   xlim=xrange,
                   ylim=c(0,a2.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkorchid1")
              
              
              Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
              Axis(x=c(0,a2.max.time.spend),side=2,at=c(0,a2.max.time.spend/2,a2.max.time.spend),labels=c(0,0.5,1))
              Axis(x=c(0,a2.max.pee.spend),side=4,at=c(0,a2.max.pee.spend/2,a2.max.pee.spend),labels=c(0,0.5,1))
              
              axis(4,at=a2.max.pee.spend/2,labels="x Pee in place x",pos=topout*.95,cex.axis=2,tick=FALSE)
              title(ylab="Time in place",line=-1.6,cex.lab=2)
              title(xlab="Time",line=-1)
              
              
              
              
              ######
              #A4
              a4.max.time.spend<-max(peeroller[,c("a4_a3_barrier","a4_a2_barrier","a4_central_corner","a4_water")])
              a4.max.pee.spend<-max(peeroller[,c("a4_a3_barrier.pee","a4_a2_barrier.pee","a4_central_corner.pee","a4_water.pee")])
              
              par(mar = c(1, 1, 1, 1))#c(bottom, left, top, right)
              plot(time.series$dailysecond,time.series$a4_a2_barrier,
                   xlim=xrange,
                   ylim=c(0,a4.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkgreen")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a4_a3_barrier,
                   xlim=xrange,
                   ylim=c(0,a4.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkred")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a4_central_corner,
                   xlim=xrange,
                   ylim=c(0,a4.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="gray40")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a4_water,
                   xlim=xrange,
                   ylim=c(0,a4.max.time.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   lwd=3,
                   type="l",col="darkorchid1")
              par(new=TRUE)
              #################################
              plot(time.series$dailysecond,time.series$a4_a2_barrier.pee,
                   xlim=xrange,
                   ylim=c(0,a4.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkgreen")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a4_a3_barrier.pee,
                   xlim=xrange,
                   ylim=c(0,a4.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkred")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a4_central_corner.pee,
                   xlim=xrange,
                   ylim=c(0,a4.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="gray40")
              par(new=TRUE)
              plot(time.series$dailysecond,time.series$a4_water.pee,
                   xlim=xrange,
                   ylim=c(0,a4.max.pee.spend),
                   xaxt='n',yaxt='n',
                   xlab="",ylab="",
                   pch=4,cex=.5,
                   #lwd=2,lty=3,type="l",
                   col="darkorchid1")
              
              
              Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
              Axis(x=c(0,a4.max.time.spend),side=2,at=c(0,a4.max.time.spend/2,a4.max.time.spend),labels=c(0,0.5,1))
              Axis(x=c(0,a4.max.pee.spend),side=4,at=c(0,a4.max.pee.spend/2,a4.max.pee.spend),labels=c(0,0.5,1))
              
              axis(4,at=a4.max.pee.spend/2,labels="x Pee in place x",pos=topout*.95,cex.axis=2,tick=FALSE)
              title(ylab="Time in place",line=-1.6,cex.lab=2)
              title(xlab="Time",line=-1)
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
              
            }
            
            ##############################
            par(mar = c(0, 2, 1, 1))#c(bottom, left, top, right)
            barplot(current.roll$a1cumdist,ylim=c(0,maxdist),
                    ylab="",yaxt='n',
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="royalblue4")
            text(x=.75,y=maxdist/2,labels = "Distance (m)",cex=2.5,col="magenta4",srt=90)
            Axis(x=c(0,maxdist),side=2,at=c(0,maxdist/2,maxdist),labels=(10*round(c(0,maxdist/2000,maxdist/1000))),cex.axis=2)
            
            par(mar = c(0, 2, 1, 1))#c(bottom, left, top, right)
            barplot(current.roll$a3cumdist,ylim=c(0,maxdist),
                    ylab="",yaxt='n',
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="firebrick3")
            text(x=.75,y=maxdist/2,labels = "Distance (m)",cex=2.5,col="magenta4",srt=90)
            Axis(x=c(0,maxdist),side=2,at=c(0,maxdist/2,maxdist),labels=(10*round(c(0,maxdist/2000,maxdist/1000))),cex.axis=2)
            
            par(mar = c(0, 2, 1, 1))#c(bottom, left, top, right)
            barplot(current.roll$a2cumdist,ylim=c(0,maxdist),
                    ylab="",yaxt='n',
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="springgreen4")
            text(x=.75,y=maxdist/2,labels = "Distance (m)",cex=2.5,col="magenta4",srt=90)
            Axis(x=c(0,maxdist),side=2,at=c(0,maxdist/2,maxdist),labels=(10*round(c(0,maxdist/2000,maxdist/1000))),cex.axis=2)
            
            par(mar = c(0, 2, 1, 1))#c(bottom, left, top, right)
            barplot(current.roll$a4cumdist,ylim=c(0,maxdist),
                    ylab="",yaxt='n',
                    #xlab="Total pee volume",
                    #main="Total pee",
                    col="chocolate1")
            text(x=.75,y=maxdist/2,labels = "Distance (m)",cex=2.5,col="magenta4",srt=90)
            Axis(x=c(0,maxdist),side=2,at=c(0,maxdist/2,maxdist),labels=(10*round(c(0,maxdist/2000,maxdist/1000))),cex.axis=2)
            
            
            
            
            
            dev.off()
            
          } else {
        
            for(teddy in startvalue:length(Just.Frames)){
        #for(teddy in startvalue:(startvalue+10)){
          
   
          txtProgressBar(min=startvalue,max=length(Just.Frames),initial=teddy,width = 80,style=3)


          
          
          current.roll<-peeroller[which(peeroller$dailysecond==Just.Frames[teddy]),]
          
          
          
          
          
          
          currentframe.A13<-keeplongA13[which(keeplongA13$dailysecond==Just.Frames[teddy]),]  
          #if all rows in current frame are equal (byproduct of cleaning data using location data)
          if(nrow(currentframe.A13)>1){
            #check if all columns in all rows of currentframe.A13 are the same, if so collapse
            if(all(apply(currentframe.A13, 2, function(x) length(unique(x)) == 1) == TRUE)) {
              currentframe.A13<-currentframe.A13[1,]
            }
          }
          currentframe.A24<-keeplongA24[which(keeplongA24$dailysecond==Just.Frames[teddy]),]  
          
          #if all rows in current frame are equal (byproduct of cleaning data using location data)
          if(nrow(currentframe.A24)>1){
            #check if all columns in all rows of currentframe.A13 are the same, if so collapse
            if(all(apply(currentframe.A24, 2, function(x) {length(unique(x))==1}) == TRUE)) {
              currentframe.A24<-currentframe.A24[1,]
            }
          }
          
          
          
          png(paste("q",Just.Frames[teddy],".png",sep=''),bg="transparent",width=1200,height=800)#use if making transparent pngs
          
          #Set up plot locations/margins
          par(oma=c(1,1,1,1))
          par(mar = c(0, 0, 0, 0))
          
          layout(matrix(c(20,10,8,8,8,8,11,11,11,11,21,13,
                          20,10,9,9,9,9,12,12,12,12,21,13,
                          16,16,16,rep(1,6),17,17,17,
                          16,16,16,rep(1,6),17,17,17,
                          16,16,16,rep(1,6),17,17,17,
                          14,14,14,rep(1,6),15,15,15,
                          14,14,14,rep(1,6),15,15,15,
                          14,14,14,rep(1,6),15,15,15,
                          
                          18,4,2,2,2,2,5,5,5,5,19,7,
                          18,4,3,3,3,3,6,6,6,6,19,7), 10, 12, byrow = TRUE),
                 widths=rep(1,10), heights=rep(1,10))
          
          ######################
          # SETUP ARENA PLOT
          # SETUP ARENA PLOT
          real.x.min<-min(c(keeplongA24$a4x0,keeplongA24$a2x0,keeplongA24$true.x),na.rm=TRUE)
          xlimspacing<-real.x.min-30
          plot(NA, xlim=c(xlimspacing,670),
               ylim=c(min(A13.outline.ys),max(A24.outline.ys)),
               xlab='',ylab='',
               axes=FALSE)
          

          ######A13
          #A13.left
          a13shape.left<-cbind(A13.outline.xs,A13.outline.ys)[c(1:nrow(A13.left)),];a13shape.left<-rbind(a13shape.left,a13shape.left[1,])
          a13shape.right<-cbind(A13.outline.xs,A13.outline.ys)[c((1+nrow(A13.left)):length(A13.outline.xs)),];a13shape.right<-rbind(a13shape.right,a13shape.right[1,])
          lines(a13shape.left)
          lines(a13shape.right)
          
          a13leftmousecoords<-getpositionalvalues(a13shape.left)
          a13rightmousecoords<-getpositionalvalues(a13shape.right)
          
          if(A13.reference$TopScreen=="inner"){
            text(min(a13shape.left[,1])-18,min(a13shape.left[,2]),A13.reference$quad.now.on.left,col="blue",cex=1.4)
            if(A13.reference$sex.left=="m"){
              text(min(a13shape.left[,1])-22,min(a13shape.left[,2])+60,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
            } else {
              text(min(a13shape.left[,1])-22,min(a13shape.left[,2])+60,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
            }
            

            text(max(a13shape.right[,1])+18,min(a13shape.right[,2]),A13.reference$quad.now.on.right,col="red",cex=1.4)
            if(A13.reference$sex.right=="m"){
              text(max(a13shape.right[,1])+22,min(a13shape.right[,2])+60,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
            } else {
              text(max(a13shape.right[,1])+22,min(a13shape.right[,2])+60,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
            }
          }
          if(A13.reference$TopScreen=="outer"){
            text(max(a13shape.left[,1])-18,max(a13shape.left[,2]),A13.reference$quad.now.on.left,col="blue",cex=1.4)
            if(A13.reference$sex.left=="m"){
              text(max(a13shape.left[,1])-22,max(a13shape.left[,2])+60,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
            } else {
              text(max(a13shape.left[,1])-22,max(a13shape.left[,2])+60,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
            }
            
            text(min(a13shape.right[,1])+18,max(a13shape.right[,2]),A13.reference$quad.now.on.right,col="red",cex=1.4)
            if(A13.reference$sex.right=="m"){
              text(min(a13shape.right[,1])+22,max(a13shape.right[,2])+60,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
            } else {
              text(min(a13shape.right[,1])+22,max(a13shape.right[,2])+60,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
            }
          }
          
          
          ####A24
          a24shape.left<-cbind(A24.outline.xs,A24.outline.ys)[c(1:nrow(A24.left)),];a24shape.left<-rbind(a24shape.left,a24shape.left[1,])
          a24shape.right<-cbind(A24.outline.xs,A24.outline.ys)[c((1+nrow(A24.left)):length(A24.outline.xs)),];a24shape.right<-rbind(a24shape.right,a24shape.right[1,])
          lines(a24shape.left)
          lines(a24shape.right)
          
          a24leftmousecoords<-getpositionalvalues(a24shape.left)
          a24rightmousecoords<-getpositionalvalues(a24shape.right)
          
          if(A24.reference$TopScreen=="outer"){
            text(min(a24shape.left[,1])-18,max(a24shape.left[,2]),A24.reference$quad.now.on.left,col="darkgreen",cex=1.4)
            if(A24.reference$sex.left=="m"){
              text(min(a24shape.left[,1])-18,max(a24shape.left[,2])-70,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
            } else {
              text(min(a24shape.left[,1])-18,max(a24shape.left[,2])-70,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
            }
            
            text(max(a24shape.right[,1])+18,max(a24shape.right[,2]),A24.reference$quad.now.on.right,col="darkgoldenrod3",cex=1.4)
            if(A24.reference$sex.right=="m"){
              text(max(a24shape.right[,1])+18,max(a24shape.right[,2])-70,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
            } else {
              text(max(a24shape.right[,1])+18,max(a24shape.right[,2])-70,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
            }
          }
          if(A24.reference$TopScreen=="inner"){
            text(max(a24shape.left[,1])-18,min(a24shape.left[,2]),A24.reference$quad.now.on.left,col="darkgreen",cex=1.4)
            if(A24.reference$sex.left=="m"){
              text(max(a24shape.left[,1])-18,min(a24shape.left[,2])-70,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
            } else {
              text(max(a24shape.left[,1])-18,min(a24shape.left[,2])-70,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
            }
            
            text(min(a24shape.right[,1])+18,min(a24shape.right[,2]),A24.reference$quad.now.on.right,col="darkgoldenrod3",cex=1.4)
            if(A24.reference$sex.right=="m"){
              text(min(a24shape.right[,1])+18,min(a24shape.right[,2])-70,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
            } else {
              text(min(a24shape.right[,1])+18,min(a24shape.right[,2])-70,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
            }
          }
          
          
          td <- seconds_to_period(Just.Frames[teddy])
          
          t.x<-(real.x.min-30)
          real.x.max<-max(c(keeplongA24$a4x0,keeplongA24$a2x0,keeplongA24$true.x),na.rm=TRUE)
          tt.x<-(real.x.max+10)
          
          text(t.x,-5,sprintf('%02d:%02d:%02d', td@hour, minute(td), second(td)),srt=90,cex = 2, col="black")
          
          text(tt.x,-5,Just.Frames[teddy],srt=90,cex = 1.8, col="darkgray")
          
          ####################################################################
          
          
          #A1
          ############
          #if there are track points to plot at this time interval
          if(!is.na(currentframe.A13$a1x0[1])){
            
            A13loc<-currentframe.A13[1,c("a1x0","a1y0")]
            
            if(a1flag==0){
              A13loc.growing<-A13loc
            } else {
              A13loc.growing<-rbind(A13loc.growing,A13loc)
            }
            sizeofA13loc<-nrow(A13loc.growing)
            
            if(a1flag>10){
              a1flag<-10
            } 
            
            a1.fromto.xs<-unlist(A13loc.growing[seq((sizeofA13loc-a1flag),sizeofA13loc,by=1),"a1x0"])  
            a1.fromto.ys<-unlist(A13loc.growing[seq((sizeofA13loc-a1flag),sizeofA13loc,by=1),"a1y0"])  
            
            lines(a1.fromto.xs,a1.fromto.ys,col=add.alpha("blue",.6),lwd=2)
            
            a1flag<-a1flag+1
            
          }
          
          #pee plots
          for(vero in 1:1){
            
            check<-currentframe.A13[which(currentframe.A13$quad=="a1"),]
            
            if(nrow(check)>0){
              
              if(a1.pee.flag==0){
                a1.pee.build<-check
              } else {
                a1.pee.build<-rbind(a1.pee.build,check)
              }
              par(new=TRUE)
              plot(a1.pee.build$true.y~a1.pee.build$true.x,cex=log(a1.pee.build$Clustsize),
                   xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),
                   bty="o",
                   pch=16,
                   xlab='',ylab='',
                   axes=FALSE,
                   col=add.alpha("blue",0.09))  
              
              
              
              par(new=TRUE)
              plot(check$true.y~check$true.x,cex=log(check$Clustsize),
                   xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),
                   bty="o",
                   pch=21,
                   xlab='',ylab='',
                   axes=FALSE,
                   bg="cyan",
                   col="cyan")
              par(new=TRUE)
              
              a1.pee.flag<-a1.pee.flag+1
            } else {
              if(exists("a1.pee.build")){
                par(new=TRUE)
                plot(a1.pee.build$true.y~a1.pee.build$true.x,cex=log(a1.pee.build$Clustsize),
                     xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),bty="o",
                     pch=16,
                     xlab='',ylab='',
                     axes=FALSE,
                     col=add.alpha("blue",0.09))  
              }
            }
            
          }
          
          #A3
          ############
          #if there are track points to plot at this time interval
          if(!is.na(currentframe.A13$a3x0[1])){
            
            A13loc.3<-currentframe.A13[1,c("a3x0","a3y0")]
            
            if(a3flag==0){
              A13loc.growing.3<-A13loc.3
            } else {
              A13loc.growing.3<-rbind(A13loc.growing.3,A13loc.3)
            }
            sizeofA13loc.3<-nrow(A13loc.growing.3)
            
            if(a3flag>10){
              a3flag<-10
            } 
            
            a3.fromto.xs<-unlist(A13loc.growing.3[seq((sizeofA13loc.3-a3flag),sizeofA13loc.3,by=1),"a3x0"])  
            a3.fromto.ys<-unlist(A13loc.growing.3[seq((sizeofA13loc.3-a3flag),sizeofA13loc.3,by=1),"a3y0"])  
            
            # print(Just.Frames[teddy])
            # print(cbind(a3.fromto.xs,a3.fromto.ys))
            # 
            #print(cbind(a3.fromto.xs,a3.fromto.ys))
            lines(a3.fromto.xs,a3.fromto.ys,col=add.alpha("red",.6),lwd=2)
            
            a3flag<-a3flag+1
            
          }
          
          
          for(vero in 1:1){
            
            check.3<-currentframe.A13[which(currentframe.A13$quad=="a3"),]
            
            if(nrow(check.3)>0){
              
              if(a3.pee.flag==0){
                a3.pee.build<-check.3
              } else {
                a3.pee.build<-rbind(a3.pee.build,check.3)
              }
              par(new=TRUE)
              plot(a3.pee.build$true.y~a3.pee.build$true.x,cex=log(a3.pee.build$Clustsize),
                   xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),bty="o",
                   pch=16,
                   xlab='',ylab='',
                   axes=FALSE,
                   col=add.alpha("darkred",0.09))  
              
              
              
              par(new=TRUE)
              plot(check.3$true.y~check.3$true.x,cex=log(check.3$Clustsize),
                   xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),
                   bty="o",
                   pch=21,
                   xlab='',ylab='',
                   axes=FALSE,
                   bg="deeppink",
                   col="deeppink")
              #par(new=TRUE)
              
              a3.pee.flag<-a3.pee.flag+1
            } else {
              if(exists("a3.pee.build")){
                par(new=TRUE)
                plot(a3.pee.build$true.y~a3.pee.build$true.x,cex=log(a3.pee.build$Clustsize),
                     xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),bty="o",
                     pch=16,
                     xlab='',ylab='',
                     axes=FALSE,
                     col=add.alpha("darkred",0.09))  
              }
            }
            
          }
          
          
          #A2
          ############
          #if there are track points to plot at this time interval
          if(!is.na(currentframe.A24$a2x0[1])){
            
            A24loc.2<-currentframe.A24[1,c("a2x0","a2y0")]
            
            if(a2flag==0){
              A24loc.growing.2<-A24loc.2
            } else {
              A24loc.growing.2<-rbind(A24loc.growing.2,A24loc.2)
            }
            sizeofA24loc.2<-nrow(A24loc.growing.2)
            
            if(a2flag>10){
              a2flag<-10
            } 
            
            a2.fromto.xs<-unlist(A24loc.growing.2[seq((sizeofA24loc.2-a2flag),sizeofA24loc.2,by=1),"a2x0"])  
            a2.fromto.ys<-unlist(A24loc.growing.2[seq((sizeofA24loc.2-a2flag),sizeofA24loc.2,by=1),"a2y0"]) 
            
            lines(a2.fromto.xs,a2.fromto.ys,col=add.alpha("green",.6),lwd=2)
            
            a2flag<-a2flag+1
            
          }
          
          
          for(vero in 1:1){
            
            check.2<-currentframe.A24[which(currentframe.A24$quad=="a2"),]
            
            if(nrow(check.2)>0){
              
              if(a2.pee.flag==0){
                a2.pee.build<-check.2
              } else {
                a2.pee.build<-rbind(a2.pee.build,check.2)
              }
              par(new=TRUE)
              plot(a2.pee.build$true.y~a2.pee.build$true.x,cex=log(a2.pee.build$Clustsize),
                   xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),bty="o",
                   pch=16,
                   xlab='',ylab='',
                   axes=FALSE,
                   col=add.alpha("darkgreen",0.09))  
              
              
              
              par(new=TRUE)
              plot(check.2$true.y~check.2$true.x,cex=log(check.2$Clustsize),
                   xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),
                   bty="o",
                   pch=21,
                   xlab='',ylab='',
                   axes=FALSE,
                   bg="green",
                   col="green")
              par(new=TRUE)
              
              a2.pee.flag<-a2.pee.flag+1
            } else {
              if(exists("a2.pee.build")){
                par(new=TRUE)
                plot(a2.pee.build$true.y~a2.pee.build$true.x,cex=log(a2.pee.build$Clustsize),
                     xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),bty="o",
                     pch=16,
                     xlab='',ylab='',
                     axes=FALSE,
                     col=add.alpha("darkgreen",0.09))  
              }
            }
            
          }
          
          #A4
          ############
          #if there are track points to plot at this time interval
          if(!is.na(currentframe.A24$a4x0[1])){
            
            A24loc.4<-currentframe.A24[1,c("a4x0","a4y0")]
            
            if(a4flag==0){
              A24loc.growing.4<-A24loc.4
            } else {
              A24loc.growing.4<-rbind(A24loc.growing.4,A24loc.4)
            }
            sizeofA24loc.4<-nrow(A24loc.growing.4)
            
            if(a4flag>10){
              a4flag<-10
            } 
            
            a4.fromto.xs<-unlist(A24loc.growing.4[seq((sizeofA24loc.4-a4flag),sizeofA24loc.4,by=1),"a4x0"])  
            a4.fromto.ys<-unlist(A24loc.growing.4[seq((sizeofA24loc.4-a4flag),sizeofA24loc.4,by=1),"a4y0"])  
            
            lines(a4.fromto.xs,a4.fromto.ys,col=add.alpha("darkorange",.6),lwd=2)

            a4flag<-a4flag+1
            
          }
          
          
          for(vero in 1:1){
            
            check.4<-currentframe.A24[which(currentframe.A24$quad=="a4"),]
            
            if(nrow(check.4)>0){
              
              if(a4.pee.flag==0){
                a4.pee.build<-check.4
              } else {
                a4.pee.build<-rbind(a4.pee.build,check.4)
              }
              par(new=TRUE)
              plot(a4.pee.build$true.y~a4.pee.build$true.x,cex=log(a4.pee.build$Clustsize),
                   xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),bty="o",
                   pch=16,
                   xlab='',ylab='',
                   axes=FALSE,
                   col=add.alpha("darkgoldenrod3",0.09))  
              
              
              
              par(new=TRUE)
              plot(check.4$true.y~check.4$true.x,cex=log(check.4$Clustsize),
                   xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),
                   bty="o",
                   pch=21,
                   xlab='',ylab='',
                   axes=FALSE,
                   bg="darkorange",
                   col="darkorange")
              #par(new=TRUE)
              
              a4.pee.flag<-a4.pee.flag+1
            } else {
              if(exists("a4.pee.build")){
                par(new=TRUE)
                plot(a4.pee.build$true.y~a4.pee.build$true.x,cex=log(a4.pee.build$Clustsize),
                     xlim=c(xlimspacing,650),ylim=c(min(A13.outline.ys),max(A24.outline.ys)),bty="o",
                     pch=16,
                     xlab='',ylab='',
                     axes=FALSE,
                     col=add.alpha("darkgoldenrod3",0.09))  
              }
            }
            
          }
          
          
          
          ##############################################
          if(teddyflag==1){
            time.series<-current.roll
            teddyflag<-teddyflag+1
          } else {
            time.series<-rbind(time.series,current.roll)
            
          }
          
          ##################
          #1
          #A1 time-series
          #rolling
          par(mar = c(0, 1, 0, 1))
          
          plot(time.series$dailysecond,time.series$a1.roll.pee,
               type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a1.roll.pee),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="dodgerblue")
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
          if(!is.na(current.roll$a1roi)){
            if(current.roll$a1roi=="a1_a3_water"|current.roll$a1roi=="a1_a2_water"){
              if(a1waterflag==0){
                a1water.growing<-current.roll
                a1waterflag<-a1waterflag+1
              } else {
                a1water.growing<-rbind(a1water.growing,current.roll)
              }
            }
          }
          axis(side=3,at=a1water.growing$dailysecond,lwd=.8,tick=TRUE,labels=rep("",nrow(a1water.growing)),col="darkorchid1")
          par(new=TRUE)
          
          plot(time.series$dailysecond,time.series$a1.roll.pee,
               type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a1.roll.pee),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="dodgerblue")
          par(new=TRUE)
          #cumulative
          plot(time.series$dailysecond,time.series$a1cumpee, 
               type="l", xlim=xrange,ylim=c(0,indiv.max.cum.pee$a1.pee.vol),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",
               lwd=3,col="darkblue") 
          
          
          Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
          Axis(x=c(0,indiv.max.cum.pee$a1.pee.vol),side=2,at=c(0,indiv.max.cum.pee$a1.pee.vol/2,indiv.max.cum.pee$a1.pee.vol),labels=c(0,0.5,1))
          title(ylab="Urine",line=-1.6,cex.lab=2)
          #title(xlab="Time",line=-1)
          
          ##################
          #VELOCITY/MOVEMENT
          par(mar = c(0, 1, 0, 1))
          plot(time.series$dailysecond,time.series$a1.rollspeed,
               type="l", xlim=xrange,ylim=c(0,maxspeed),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="dodgerblue")
          par(new=TRUE)
          #cumulative
          plot(time.series$dailysecond,time.series$a1cumdist, 
               type="l", xlim=xrange,ylim=c(0,indiv.max.cum.disttrav$a1.displacement),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",
               lwd=3,col="darkblue") 
          
          Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
          Axis(x=c(0,indiv.max.cum.disttrav$a1.displacement),side=2,at=c(0,indiv.max.cum.disttrav$a1.displacement/2,indiv.max.cum.disttrav$a1.displacement),labels=c(0,0.5,1))
          title(ylab="Velocity",line=-1.6,cex.lab=2)
          title(xlab="Time",line=-1)
          
          par(mar = c(0, 2, 1, 1))
          barplot(current.roll$a1cumpee,ylim=c(0,maxcumpee),
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="blue")
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
          par(new=TRUE)
          barplot(current.roll$a1cumpee,ylim=c(0,maxcumpee),
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="blue")
          text(x=.75,y=maxcumpee/2,labels = "Total pee",cex=3,col="yellowgreen",srt=90)
          
          
          
          ###########################################################################
          #1
          #A3 time-series
          #rolling
          par(mar = c(0, 1, 0, 1))
          plot(time.series$dailysecond,time.series$a3.roll.pee,
               type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a3.roll.pee),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="pink")
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
          if(!is.na(current.roll$a3roi)){
            if(current.roll$a3roi=="a3_a1_water"|current.roll$a3roi=="a3_a4_water"){
              if(a3waterflag==0){
                a3water.growing<-current.roll
                a3waterflag<-a3waterflag+1
              } else {
                a3water.growing<-rbind(a3water.growing,current.roll)
              }
            }
          }           
          axis(side=3,at=a3water.growing$dailysecond,lwd=.8,tick=TRUE,labels=rep("",nrow(a3water.growing)),col="darkorchid1")
          
          
          par(new=TRUE)
          plot(time.series$dailysecond,time.series$a3.roll.pee,
               type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a3.roll.pee),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="pink")
          par(new=TRUE)
          #cumulative
          plot(time.series$dailysecond,time.series$a3cumpee, 
               type="l", xlim=xrange,ylim=c(0,indiv.max.cum.pee$a3.pee.vol),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",
               lwd=3,col="deeppink") 
          
          
          
          
          Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
          Axis(x=c(0,indiv.max.cum.pee$a3.pee.vol),side=2,at=c(0,indiv.max.cum.pee$a3.pee.vol/2,indiv.max.cum.pee$a3.pee.vol),labels=c(0,0.5,1))
          title(ylab="Urine",line=-1.6,cex.lab=2)
          #title(xlab="Time",line=-1)
          
          ##################
          #VELOCITY/MOVEMENT
          par(mar = c(0, 1, 0, 1))
          plot(time.series$dailysecond,time.series$a3.rollspeed,
               type="l", xlim=xrange,ylim=c(0,maxspeed),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="pink")
          par(new=TRUE)
          #cumulative
          plot(time.series$dailysecond,time.series$a3cumdist, 
               type="l", xlim=xrange,ylim=c(0,indiv.max.cum.disttrav$a3.displacement),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",
               lwd=3,col="deeppink") 
          
          Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
          Axis(x=c(0,indiv.max.cum.disttrav$a3.displacement),side=2,at=c(0,indiv.max.cum.disttrav$a3.displacement/2,indiv.max.cum.disttrav$a3.displacement),labels=c(0,0.5,1))
          title(ylab="Velocity",line=-1.6,cex.lab=2)
          title(xlab="Time",line=-1)
          
          par(mar = c(0, 2, 1, 1))
          barplot(current.roll$a3cumpee,ylim=c(0,maxcumpee),
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="darkred")
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
          par(new=TRUE)
          barplot(current.roll$a3cumpee,ylim=c(0,maxcumpee),
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="darkred")
          text(x=.75,y=maxcumpee/2,labels = "Total pee",cex=3,col="yellowgreen",srt=90)
          
          
          ###########################################################################
          #1
          #A2 time-series
          #rolling
          par(mar = c(0, 1, 0, 1))
          plot(time.series$dailysecond,time.series$a2.roll.pee,
               type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a2.roll.pee),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="green")
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
          if(!is.na(current.roll$a2roi)){
            if(current.roll$a2roi=="a2_a1_water"|current.roll$a2roi=="a2_a4_water"){
              if(a2waterflag==0){
                a2water.growing<-current.roll
                a2waterflag<-a2waterflag+1
              } else {
                a2water.growing<-rbind(a2water.growing,current.roll)
              }
            }
          }
          axis(side=3,at=a2water.growing$dailysecond,lwd=.8,tick=TRUE,labels=rep("",nrow(a2water.growing)),col="darkorchid1")
          
          
          par(new=TRUE)
          plot(time.series$dailysecond,time.series$a2.roll.pee,
               type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a2.roll.pee),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="green")
          par(new=TRUE)
          #cumulative
          plot(time.series$dailysecond,time.series$a2cumpee, 
               type="l", xlim=xrange,ylim=c(0,indiv.max.cum.pee$a2.pee.vol),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",
               lwd=3,col="darkgreen") 
          
          
          
          Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
          Axis(x=c(0,indiv.max.cum.pee$a2.pee.vol),side=2,at=c(0,indiv.max.cum.pee$a2.pee.vol/2,indiv.max.cum.pee$a2.pee.vol),labels=c(0,0.5,1))
          title(ylab="Urine",line=-1.6,cex.lab=2)
          #title(xlab="Time",line=-1)
          
          ##################
          #VELOCITY/MOVEMENT
          par(mar = c(0, 1, 0, 1))
          plot(time.series$dailysecond,time.series$a2.rollspeed,
               type="l", xlim=xrange,ylim=c(0,maxspeed),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="green")
          par(new=TRUE)
          #cumulative
          plot(time.series$dailysecond,time.series$a2cumdist, 
               type="l", xlim=xrange,ylim=c(0,indiv.max.cum.disttrav$a2.displacement),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",
               lwd=3,col="darkgreen") 
          
          Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
          Axis(x=c(0,indiv.max.cum.disttrav$a2.displacement),side=2,at=c(0,indiv.max.cum.disttrav$a2.displacement/2,indiv.max.cum.disttrav$a2.displacement),labels=c(0,0.5,1))
          title(ylab="Velocity",line=-1.6,cex.lab=2)
          title(xlab="Time",line=-1)
          
          par(mar = c(0, 2, 1, 1))
          barplot(current.roll$a2cumpee,ylim=c(0,maxcumpee),
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="darkgreen")
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
          par(new=TRUE)
          barplot(current.roll$a2cumpee,ylim=c(0,maxcumpee),
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="darkgreen")
          text(x=.75,y=maxcumpee/2,labels = "Total pee",cex=3,col="yellowgreen",srt=90)
          
          ###########################################################################
          #1
          #A4 time-series
          #rolling
          par(mar = c(0, 1, 0, 1))
          plot(time.series$dailysecond,time.series$a4.roll.pee,
               type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a4.roll.pee),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="darkgoldenrod3")
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
          if(!is.na(current.roll$a4roi)){
            if(current.roll$a4roi=="a4_a3_water"|current.roll$a4roi=="a4_a2_water"){
              if(a4waterflag==0){
                a4water.growing<-current.roll
                a4waterflag<-a4waterflag+1
              } else {
                a4water.growing<-rbind(a4water.growing,current.roll)
              }
            }
          }
          axis(side=3,at=a4water.growing$dailysecond,lwd=.8,tick=TRUE,labels=rep("",nrow(a4water.growing)),col="darkorchid1")
          
          par(new=TRUE)
          plot(time.series$dailysecond,time.series$a4.roll.pee,
               type="l", xlim=xrange,ylim=c(0,indiv.max.instant.pee$a4.roll.pee),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="darkgoldenrod3")
          par(new=TRUE)
          #cumulative
          plot(time.series$dailysecond,time.series$a4cumpee, 
               type="l", xlim=xrange,ylim=c(0,indiv.max.cum.pee$a4.pee.vol),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",
               lwd=3,col="darkorange") 
          
          
          
          Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
          Axis(x=c(0,indiv.max.cum.pee$a4.pee.vol),side=2,at=c(0,indiv.max.cum.pee$a4.pee.vol/2,indiv.max.cum.pee$a4.pee.vol),labels=c(0,0.5,1))
          title(ylab="Urine",line=-1.6,cex.lab=2)
          #title(xlab="Time",line=-1)
          
          ##################
          #VELOCITY/MOVEMENT
          par(mar = c(0, 1, 0, 1)) #c(bottom, left, top, right)
          plot(time.series$dailysecond,time.series$a4.rollspeed,
               type="l", xlim=xrange,ylim=c(0,maxspeed),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",lty=1,
               lwd=1.5,col="darkgoldenrod3")
          par(new=TRUE)
          #cumulative
          plot(time.series$dailysecond,time.series$a4cumdist, 
               type="l", xlim=xrange,ylim=c(0,indiv.max.cum.disttrav$a4.displacement),
               pch=16,xaxt='n',yaxt='n',
               xlab="",ylab="",
               lwd=3,col="darkorange") 
          
          Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
          Axis(x=c(0,indiv.max.cum.disttrav$a4.displacement),side=2,at=c(0,indiv.max.cum.disttrav$a4.displacement/2,indiv.max.cum.disttrav$a4.displacement),labels=c(0,0.5,1))
          title(ylab="Velocity",line=-1.6,cex.lab=2)
          title(xlab="Time",line=-1)
          
          par(mar = c(0, 2, 1, 1))
          barplot(current.roll$a4cumpee,ylim=c(0,maxcumpee),
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="darkorange")
          rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "cornsilk1")
          par(new=TRUE)
          barplot(current.roll$a4cumpee,ylim=c(0,maxcumpee),
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="darkorange")
          text(x=.75,y=maxcumpee/2,labels = "Total pee",cex=3,col="yellowgreen",srt=90)
          
          
          
          if(wingplots=="stacked"){
            #Stacked Plots
            
            #a1
            h.a1<-as.matrix(c(current.roll$a1_water,
                              current.roll$a1_a3_barrier,
                              current.roll$a1_central_corner,
                              current.roll$a1_a2_barrier))
            H <- apply(h.a1, 2L, cumsum)
            H <- H - h.a1 / 2
            par(mar = c(2, 2, 1, 1))#c(bottom, left, top, right)
            x.a1<-barplot(h.a1,
                          col=c("darkorchid1","darkred","gray40","darkgreen"),
                          ylim=c(0,sum(c(current.roll$a1_a2_barrier,
                                         current.roll$a1_central_corner,
                                         current.roll$a1_a3_barrier,
                                         current.roll$a1_water))))
            text(rep(x.a1, each = nrow(H)), H, labels = c("water","a3","center","a2"),col="cornsilk")
            title(ylab="pee breakdown",line=-20,cex.lab=2)
            
            #a3
            h.a3<-as.matrix(c(current.roll$a3_water,
                              current.roll$a3_a1_barrier,
                              current.roll$a3_central_corner,
                              current.roll$a3_a4_barrier))
            H <- apply(h.a3, 2L, cumsum)
            H <- H - h.a3 / 2
            par(mar = c(2, 2, 1, 1))#c(bottom, left, top, right)
            x.a3<-barplot(h.a3,
                          col=c("darkorchid1","darkblue","gray40","darkorange"),
                          ylim=c(0,sum(c(current.roll$a3_a4_barrier,
                                         current.roll$a3_central_corner,
                                         current.roll$a3_a1_barrier,
                                         current.roll$a3_water))))
            text(rep(x.a3, each = nrow(H)), H, labels = c("water","a1","center","a4"),col="cornsilk")
            title(ylab="pee breakdown",line=-20,cex.lab=2)
            
            
            #a2
            h.a2<-as.matrix(c(current.roll$a2_a1_barrier,
                              current.roll$a2_central_corner,
                              current.roll$a2_a4_barrier,
                              current.roll$a2_water))
            H <- apply(h.a2, 2L, cumsum)
            H <- H - h.a2 / 2
            par(mar = c(0, 2, 2, 1))#c(bottom, left, top, right)
            x.a2<-barplot(h.a2,
                          col=c("darkblue","gray40","darkorange","darkorchid1"),
                          ylim=c(0,sum(c(current.roll$a2_a1_barrier,
                                         current.roll$a2_central_corner,
                                         current.roll$a2_a4_barrier,
                                         current.roll$a2_water))))
            text(rep(x.a2, each = nrow(H)), H, labels = c("a1","center","a4","water"),col="cornsilk")
            title(ylab="pee breakdown",line=-20,cex.lab=2)
            
            
            #a4
            h.a4<-as.matrix(c(current.roll$a4_a3_barrier,
                              current.roll$a4_central_corner,
                              current.roll$a4_a2_barrier,
                              current.roll$a4_water))
            H <- apply(h.a4, 2L, cumsum)
            H <- H - h.a4 / 2
            par(mar = c(0, 2, 2, 1))#c(bottom, left, top, right)
            
            x.a4<-barplot(h.a4,
                          col=c("darkred","gray40","darkgreen","darkorchid1"),
                          ylim=c(0,sum(c(current.roll$a4_a3_barrier,
                                         current.roll$a4_central_corner,
                                         current.roll$a4_a2_barrier,
                                         current.roll$a4_water))))
            text(rep(x.a4, each = nrow(H)), H, labels = c("a3","center","a2","water"),col="cornsilk")
            title(ylab="pee breakdown",line=-20,cex.lab=2)
          }
          if(wingplots=="climbinglines"){
            
            # if(climbflag==1){
            #   climbstart=teddy
            #   climbflag<-climbflag+1
            # }  
            # 
            # climb.roll<-peeroller[c(climbstart:teddy),]
            
            ######
            #A1
            a1.max.time.spend<-max(peeroller[,c("a1_a3_barrier","a1_a2_barrier","a1_central_corner","a1_water")])
            a1.max.pee.spend<-max(peeroller[,c("a1_a3_barrier.pee","a1_a2_barrier.pee","a1_central_corner.pee","a1_water.pee")])
            
            par(mar = c(1, 1, 1, 1))#c(bottom, left, top, right)
            plot(time.series$dailysecond,time.series$a1_a2_barrier,
                 xlim=xrange,
                 ylim=c(0,a1.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkgreen")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a1_a3_barrier,
                 xlim=xrange,
                 ylim=c(0,a1.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkred")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a1_central_corner,
                 xlim=xrange,
                 ylim=c(0,a1.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="gray40")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a1_water,
                 xlim=xrange,
                 ylim=c(0,a1.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkorchid1")
            par(new=TRUE)
            #################################
            plot(time.series$dailysecond,time.series$a1_a2_barrier.pee,
                 xlim=xrange,
                 ylim=c(0,a1.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkgreen")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a1_a3_barrier.pee,
                 xlim=xrange,
                 ylim=c(0,a1.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkred")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a1_central_corner.pee,
                 xlim=xrange,
                 ylim=c(0,a1.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="gray40")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a1_water.pee,
                 xlim=xrange,
                 ylim=c(0,a1.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkorchid1")
            
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,a1.max.time.spend),side=2,at=c(0,a1.max.time.spend/2,a1.max.time.spend),labels=c(0,0.5,1))
            Axis(x=c(0,a1.max.pee.spend),side=4,at=c(0,a1.max.pee.spend/2,a1.max.pee.spend),labels=c(0,0.5,1))
            
            axis(4,at=a1.max.pee.spend/2,labels="x Pee in place x",pos=topout*.95,cex.axis=2,tick=FALSE)
            title(ylab="Time in place",line=-1.6,cex.lab=2)
            title(xlab="Time",line=-1)
            
            ######
            #A3
            a3.max.time.spend<-max(peeroller[,c("a3_a1_barrier","a3_a4_barrier","a3_central_corner","a3_water")])
            a3.max.pee.spend<-max(peeroller[,c("a3_a1_barrier.pee","a3_a4_barrier.pee","a3_central_corner.pee","a3_water.pee")])
            
            par(mar = c(1, 1, 1, 1))#c(bottom, left, top, right)
            plot(time.series$dailysecond,time.series$a3_a1_barrier,
                 xlim=xrange,
                 ylim=c(0,a3.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkblue")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a3_a4_barrier,
                 xlim=xrange,
                 ylim=c(0,a3.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkorange")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a3_central_corner,
                 xlim=xrange,
                 ylim=c(0,a3.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="gray40")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a3_water,
                 xlim=xrange,
                 ylim=c(0,a3.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkorchid1")
            par(new=TRUE)
            #################################
            plot(time.series$dailysecond,time.series$a3_a1_barrier.pee,
                 xlim=xrange,
                 ylim=c(0,a3.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkblue")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a3_a4_barrier.pee,
                 xlim=xrange,
                 ylim=c(0,a3.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkorange")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a3_central_corner.pee,
                 xlim=xrange,
                 ylim=c(0,a3.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="gray40")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a3_water.pee,
                 xlim=xrange,
                 ylim=c(0,a3.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkorchid1")
            
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,a3.max.time.spend),side=2,at=c(0,a3.max.time.spend/2,a3.max.time.spend),labels=c(0,0.5,1))
            Axis(x=c(0,a3.max.pee.spend),side=4,at=c(0,a3.max.pee.spend/2,a3.max.pee.spend),labels=c(0,0.5,1))
            
            axis(4,at=a3.max.pee.spend/2,labels="x Pee in place x",pos=topout*.95,cex.axis=2,tick=FALSE)
            title(ylab="Time in place",line=-1.6,cex.lab=2)
            title(xlab="Time",line=-1)
            
            ######
            #A2
            
            a2.max.time.spend<-max(peeroller[,c("a2_a1_barrier","a2_a4_barrier","a2_central_corner","a2_water")])
            a2.max.pee.spend<-max(peeroller[,c("a2_a1_barrier.pee","a2_a4_barrier.pee","a2_central_corner.pee","a2_water.pee")])
            
            par(mar = c(1, 1, 1, 1))#c(bottom, left, top, right)
            plot(time.series$dailysecond,time.series$a2_a1_barrier,
                 xlim=xrange,
                 ylim=c(0,a2.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkblue")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a2_a4_barrier,
                 xlim=xrange,
                 ylim=c(0,a2.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkorange")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a2_central_corner,
                 xlim=xrange,
                 ylim=c(0,a2.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="gray40")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a2_water,
                 xlim=xrange,
                 ylim=c(0,a2.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkorchid1")
            par(new=TRUE)
            #################################
            plot(time.series$dailysecond,time.series$a2_a1_barrier.pee,
                 xlim=xrange,
                 ylim=c(0,a2.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkblue")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a2_a4_barrier.pee,
                 xlim=xrange,
                 ylim=c(0,a2.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkorange")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a2_central_corner.pee,
                 xlim=xrange,
                 ylim=c(0,a2.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="gray40")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a2_water.pee,
                 xlim=xrange,
                 ylim=c(0,a2.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkorchid1")
            
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,a2.max.time.spend),side=2,at=c(0,a2.max.time.spend/2,a2.max.time.spend),labels=c(0,0.5,1))
            Axis(x=c(0,a2.max.pee.spend),side=4,at=c(0,a2.max.pee.spend/2,a2.max.pee.spend),labels=c(0,0.5,1))
            
            axis(4,at=a2.max.pee.spend/2,labels="x Pee in place x",pos=topout*.95,cex.axis=2,tick=FALSE)
            title(ylab="Time in place",line=-1.6,cex.lab=2)
            title(xlab="Time",line=-1)
            
            
            
            
            ######
            #A4
            a4.max.time.spend<-max(peeroller[,c("a4_a3_barrier","a4_a2_barrier","a4_central_corner","a4_water")])
            a4.max.pee.spend<-max(peeroller[,c("a4_a3_barrier.pee","a4_a2_barrier.pee","a4_central_corner.pee","a4_water.pee")])
            
            par(mar = c(1, 1, 1, 1))#c(bottom, left, top, right)
            plot(time.series$dailysecond,time.series$a4_a2_barrier,
                 xlim=xrange,
                 ylim=c(0,a4.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkgreen")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a4_a3_barrier,
                 xlim=xrange,
                 ylim=c(0,a4.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkred")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a4_central_corner,
                 xlim=xrange,
                 ylim=c(0,a4.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="gray40")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a4_water,
                 xlim=xrange,
                 ylim=c(0,a4.max.time.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 lwd=3,
                 type="l",col="darkorchid1")
            par(new=TRUE)
            #################################
            plot(time.series$dailysecond,time.series$a4_a2_barrier.pee,
                 xlim=xrange,
                 ylim=c(0,a4.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkgreen")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a4_a3_barrier.pee,
                 xlim=xrange,
                 ylim=c(0,a4.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkred")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a4_central_corner.pee,
                 xlim=xrange,
                 ylim=c(0,a4.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="gray40")
            par(new=TRUE)
            plot(time.series$dailysecond,time.series$a4_water.pee,
                 xlim=xrange,
                 ylim=c(0,a4.max.pee.spend),
                 xaxt='n',yaxt='n',
                 xlab="",ylab="",
                 pch=4,cex=.5,
                 #lwd=2,lty=3,type="l",
                 col="darkorchid1")
            
            
            Axis(x=xrange,side=1,at=c(0,36001/4,36001/2,(36001*3)/4,36001),labels=c(0,0.25,0.5,0.75,1))
            Axis(x=c(0,a4.max.time.spend),side=2,at=c(0,a4.max.time.spend/2,a4.max.time.spend),labels=c(0,0.5,1))
            Axis(x=c(0,a4.max.pee.spend),side=4,at=c(0,a4.max.pee.spend/2,a4.max.pee.spend),labels=c(0,0.5,1))
            
            axis(4,at=a4.max.pee.spend/2,labels="x Pee in place x",pos=topout*.95,cex.axis=2,tick=FALSE)
            title(ylab="Time in place",line=-1.6,cex.lab=2)
            title(xlab="Time",line=-1)
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
          }
          
          ##############################
          par(mar = c(0, 2, 1, 1))#c(bottom, left, top, right)
          barplot(current.roll$a1cumdist,ylim=c(0,maxdist),
                  ylab="",yaxt='n',
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="royalblue4")
          text(x=.75,y=maxdist/2,labels = "Distance (m)",cex=2.5,col="magenta4",srt=90)
          Axis(x=c(0,maxdist),side=2,at=c(0,maxdist/2,maxdist),labels=(10*round(c(0,maxdist/2000,maxdist/1000))),cex.axis=2)
          
          par(mar = c(0, 2, 1, 1))#c(bottom, left, top, right)
          barplot(current.roll$a3cumdist,ylim=c(0,maxdist),
                  ylab="",yaxt='n',
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="firebrick3")
          text(x=.75,y=maxdist/2,labels = "Distance (m)",cex=2.5,col="magenta4",srt=90)
          Axis(x=c(0,maxdist),side=2,at=c(0,maxdist/2,maxdist),labels=(10*round(c(0,maxdist/2000,maxdist/1000))),cex.axis=2)
          
          par(mar = c(0, 2, 1, 1))#c(bottom, left, top, right)
          barplot(current.roll$a2cumdist,ylim=c(0,maxdist),
                  ylab="",yaxt='n',
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="springgreen4")
          text(x=.75,y=maxdist/2,labels = "Distance (m)",cex=2.5,col="magenta4",srt=90)
          Axis(x=c(0,maxdist),side=2,at=c(0,maxdist/2,maxdist),labels=(10*round(c(0,maxdist/2000,maxdist/1000))),cex.axis=2)
          
          par(mar = c(0, 2, 1, 1))#c(bottom, left, top, right)
          barplot(current.roll$a4cumdist,ylim=c(0,maxdist),
                  ylab="",yaxt='n',
                  #xlab="Total pee volume",
                  #main="Total pee",
                  col="chocolate1")
          text(x=.75,y=maxdist/2,labels = "Distance (m)",cex=2.5,col="magenta4",srt=90)
          Axis(x=c(0,maxdist),side=2,at=c(0,maxdist/2,maxdist),labels=(10*round(c(0,maxdist/2000,maxdist/1000))),cex.axis=2)
          
          
          
          
          
          dev.off()
          
        }
  
        }  
    }#each day
  }#each trial  
}#closes if(makeplots=="yes") bracket







##############################
#CODE TO PLOT FULL-FRAME FOR VIDEO OVERLAY
#########################
# EXPORT FULL SIZE, SINGLE CAMERA DATA
#####################################################################
#FULLFRAME

####################
#  keeplongA13.full<-keeplongA13
#  keeplongA24.full<-keeplongA24
#  
#  a13mult<-ifelse(max(keeplongA13[,c(3,6,11)],na.rm=TRUE)<0,1,-1)
#  a24mult<-ifelse(max(keeplongA24[,c(3,6,11)],na.rm=TRUE)<0,1,-1)
#  
#  keeplongA13.full[,c(3,6,11)]<-apply(keeplongA13.full[,c(3,6,11)],2,function(x){x*(a13mult)})
#  keeplongA24.full[,c(3,6,11)]<-apply(keeplongA24.full[,c(3,6,11)],2,function(x){x*(a24mult)})
#  
#  A13.outline.xs.full<-A13.outline.xs;
#  A13.outline.ys.full<-A13.outline.ys*a13mult;
#  A24.outline.xs.full<-A24.outline.xs;
#  A24.outline.ys.full<-A24.outline.ys*a24mult;
#  
#  whichcameratoplot<-"A13"
#  
#  a1flag<-0
#  a1.pee.flag<-0
#  a2flag<-0
#  a2.pee.flag<-0
#  a3flag<-0
#  a3.pee.flag<-0
#  a4flag<-0
#  a4.pee.flag<-0
#  
#  if(exists("a1.pee.build")){rm(a1.pee.build)};if(exists("a2.pee.build")){rm(a2.pee.build)};
#  if(exists("a3.pee.build")){rm(a3.pee.build)};if(exists("a4.pee.build")){rm(a4.pee.build)}
#  
#  for(teddy in 1:length(Just.Frames)){
# # for(teddy in 16030:16050){
#    
#    
#    currentframe.A13<-keeplongA13.full[which(keeplongA13.full$dailysecond==Just.Frames[teddy]),]  
#    #if all rows in current frame are equal (byproduct of cleaning data using location data)
#    if(nrow(currentframe.A13)>1){
#      #check if all columns in all rows of currentframe.A13 are the same, if so collapse
#      if(all(apply(currentframe.A13, 2, function(x) length(unique(x)) == 1) == TRUE)) {
#        currentframe.A13<-currentframe.A13[1,]
#      }
#    }
#    currentframe.A24<-keeplongA24.full[which(keeplongA24.full$dailysecond==Just.Frames[teddy]),]  
#    #if all rows in current frame are equal (byproduct of cleaning data using location data)
#    if(nrow(currentframe.A24)>1){
#      #check if all columns in all rows of currentframe.A13 are the same, if so collapse
#      if(all(apply(currentframe.A24, 2, function(x) length(unique(x)) == 1) == TRUE)) {
#        currentframe.A24<-currentframe.A24[1,]
#      }
#    }
#    
#    
#    if((sum(is.na(currentframe.A13[1,]))>11) & (sum(is.na(currentframe.A24[1,]))>11)){
#      #this ifelse basically lets us skip 'empty' time periods (e.g. at the beginning of Day1)
#    } else {
#      
#      png(paste(Just.Frames[teddy],".png",sep=''),bg="transparent",width=640,height=480)#use if making transparent pngs
#      
#      #Set up plot locations/margins
#      par(oma=c(0,0,0,0))
#      par(mar = c(0, 0, 0, 0))
#      
#      # layout(matrix(c(4,4,5,5,
#      #                 rep(1,32),
#      #                 2,2,3,3), 10, 4, byrow = TRUE),
#      #        widths=c(1,1,1,1), heights=rep(1,10))
#      # 
#      ###
#      # SETUP ARENA PLOT
#      plot(NA, xlim=c(1,640),ylim=c(-480,-1),
#           xlab='',ylab='',
#           axes=FALSE)
#      
#      
#      if(whichcameratoplot=="A13"){
#        lines(A13.outline.xs.full,A13.outline.ys.full)
#        
#        text(A13.outline.xs.full[,6]-18,A13.outline.ys.full[,6],A13.reference$quad.now.on.left,col="blue",cex=1.4)
#        if(A13.reference$sex.left=="m"){
#          text(A13.outline.xs.full[,6]-15,A13.outline.ys.full[,6]-20,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
#        } else {
#          text(A13.outline.xs.full[,6]-18,A13.outline.ys.full[,6]-20,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
#        }
#        
#        text(A13.outline.xs.full[,4]+18,A13.outline.ys.full[,4],A13.reference$quad.now.on.right,col="red",cex=1.4)
#        if(A13.reference$sex.right=="m"){
#          text(A13.outline.xs.full[,4]+15,A13.outline.ys.full[,4]-20,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
#        } else {
#          text(A13.outline.xs.full[,4]+18,A13.outline.ys.full[,4]-20,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
#        }
#      }
#      
#      if(whichcameratoplot=="A24"){
#        lines(A24.outline.xs.full,A24.outline.ys.full)
#        
#        text(A24.outline.xs.full[,1]-18,A24.outline.ys.full[,1],A24.reference$quad.now.on.left,col="darkgreen",cex=1.4)
#        if(A24.reference$sex.left=="m"){
#          text(A24.outline.xs.full[,1]-18,A24.outline.ys.full[,1]+20,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
#        } else {
#          text(A24.outline.xs.full[,1]-18,A24.outline.ys.full[,1]+20,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
#        }
#        
#        text(A24.outline.xs.full[,3]+18,A24.outline.ys.full[,3],A24.reference$quad.now.on.right,col="darkgoldenrod3",cex=1.4)
#        if(A24.reference$sex.right=="m"){
#          text(A24.outline.xs.full[,3]+18,A24.outline.ys.full[,3]+20,labels = "\\MA",vfont=c("sans serif","bold"),cex=3)
#        } else {
#          text(A24.outline.xs.full[,3]+18,A24.outline.ys.full[,3]+20,labels = "\\VE",vfont=c("sans serif","bold"),cex=3)
#        }
#      }
#      
#      td <- seconds_to_period(Just.Frames[teddy])
#      
#      text(-8,-30,sprintf('%02d:%02d:%02d', td@hour, minute(td), second(td)),srt=90,cex = 1.15, col="black")
#      
#      text(13,-30,Just.Frames[teddy],srt=90,cex = .95, col="darkgray")
#      
#      ####################################################################
#      
#      if(whichcameratoplot=="A13"){
#          #A1
#          ############
#          #if there are track points to plot at this time interval
#          if(!is.na(currentframe.A13$a1x0[1])){
#            
#            A13loc<-currentframe.A13[1,c("a1x0","a1y0")]
#            
#            if(a1flag==0){
#              A13loc.growing<-A13loc
#            } else {
#              A13loc.growing<-rbind(A13loc.growing,A13loc)
#            }
#            sizeofA13loc<-nrow(A13loc.growing)
#            
#            if(a1flag>10){
#              a1flag<-10
#            } 
#            
#            a1.fromto.xs<-A13loc.growing[seq((sizeofA13loc-a1flag),sizeofA13loc,by=1),"a1x0"]  
#            a1.fromto.ys<-A13loc.growing[seq((sizeofA13loc-a1flag),sizeofA13loc,by=1),"a1y0"]  
#            
#            lines(a1.fromto.xs,a1.fromto.ys,col=add.alpha("blue",.6),lwd=2)
#            
#            a1flag<-a1flag+1
#            
#          }
#          
#          
#          for(vero in 1:1){
#            
#            check<-currentframe.A13[which(currentframe.A13$quad=="a1"),]
#            
#            if(nrow(check)>0){
#              
#              if(a1.pee.flag==0){
#                a1.pee.build<-check
#              } else {
#                a1.pee.build<-rbind(a1.pee.build,check)
#              }
#              par(new=TRUE)
#              plot(a1.pee.build$true.y~a1.pee.build$true.x,cex=log(a1.pee.build$Clustsize),
#                   xlim=c(1,640),ylim=c(-480,-1),
#                   bty="o",
#                   pch=21,
#                   xlab='',ylab='',
#                   axes=FALSE,
#                   bg=add.alpha("cyan",0.00),
#                   col=add.alpha("cyan",0.1))  
#              
#              
#              
#              par(new=TRUE)
#              plot(check$true.y~check$true.x,cex=(1.4*log(check$Clustsize)),
#                   xlim=c(1,640),ylim=c(-480,-1),
#                   bty="o",
#                   pch=23,
#                   xlab='',ylab='',
#                   axes=FALSE,
#                   bg=add.alpha("cyan",0.00),
#                   col="blue")
#              par(new=TRUE)
#              
#              a1.pee.flag<-a1.pee.flag+1
#            } else {
#              if(exists("a1.pee.build")){
#                par(new=TRUE)
#                plot(a1.pee.build$true.y~a1.pee.build$true.x,cex=log(a1.pee.build$Clustsize),
#                     xlim=c(1,640),ylim=c(-480,-1),bty="o",
#                     pch=21,
#                     xlab='',ylab='',
#                     axes=FALSE,
#                     bg=add.alpha("cyan",0.00),
#                     col=add.alpha("cyan",0.1))  
#              }
#            }
#            
#          }
#          
#          #A3
#          ############
#          #if there are track points to plot at this time interval
#          if(!is.na(currentframe.A13$a3x0[1])){
#            
#            A13loc.3<-currentframe.A13[1,c("a3x0","a3y0")]
#            
#            if(a3flag==0){
#              A13loc.growing.3<-A13loc.3
#            } else {
#              A13loc.growing.3<-rbind(A13loc.growing.3,A13loc.3)
#            }
#            sizeofA13loc.3<-nrow(A13loc.growing.3)
#            
#            if(a3flag>10){
#              a3flag<-10
#            } 
#            
#            a3.fromto.xs<-A13loc.growing.3[seq((sizeofA13loc.3-a3flag),sizeofA13loc.3,by=1),"a3x0"]  
#            a3.fromto.ys<-A13loc.growing.3[seq((sizeofA13loc.3-a3flag),sizeofA13loc.3,by=1),"a3y0"]  
#            # 
#            # print(Just.Frames[teddy])
#            # print(cbind(a3.fromto.xs,a3.fromto.ys))
#            
#            #print(cbind(a3.fromto.xs,a3.fromto.ys))
#            lines(a3.fromto.xs,a3.fromto.ys,col=add.alpha("red",.6),lwd=2)
#            
#            a3flag<-a3flag+1
#            
#          }
#          
#          
#          for(vero in 1:1){
#            
#            check.3<-currentframe.A13[which(currentframe.A13$quad=="a3"),]
#            
#            if(nrow(check.3)>0){
#              
#              if(a3.pee.flag==0){
#                a3.pee.build<-check.3
#              } else {
#                a3.pee.build<-rbind(a3.pee.build,check.3)
#              }
#              par(new=TRUE)
#              plot(a3.pee.build$true.y~a3.pee.build$true.x,cex=log(a3.pee.build$Clustsize),
#                   xlim=c(1,640),ylim=c(-480,-1),bty="o",
#                   pch=21,
#                   xlab='',ylab='',
#                   axes=FALSE,
#                   bg=add.alpha("deeppink",0.00),
#                   col=add.alpha("deeppink",0.1))  
#              
#              
#              
#              par(new=TRUE)
#              plot(check.3$true.y~check.3$true.x,cex=(1.4*log(check.3$Clustsize)),
#                   xlim=c(1,640),ylim=c(-480,-1),
#                   bty="o",
#                   pch=23,
#                   xlab='',ylab='',
#                   axes=FALSE,
#                   bg=add.alpha("cyan",0.00),
#                   col="red")
#              par(new=TRUE)
#              
#              a3.pee.flag<-a3.pee.flag+1
#            } else {
#              if(exists("a3.pee.build")){
#                par(new=TRUE)
#                plot(a3.pee.build$true.y~a3.pee.build$true.x,cex=log(a3.pee.build$Clustsize),
#                     xlim=c(1,640),ylim=c(-480,-1),bty="o",
#                     pch=21,
#                     xlab='',ylab='',
#                     axes=FALSE,
#                     bg=add.alpha("deeppink",0.00),
#                     col=add.alpha("deeppink",0.1))  
#              }
#            }
#            
#          }
#      }    
#      
#      if(whichcameratoplot=="A24"){
#          #A2
#          ############
#          #if there are track points to plot at this time interval
#          if(!is.na(currentframe.A24$a2x0[1])){
#            
#            A24loc.2<-currentframe.A24[1,c("a2x0","a2y0")]
#            
#            if(a2flag==0){
#              A24loc.growing.2<-A24loc.2
#            } else {
#              A24loc.growing.2<-rbind(A24loc.growing.2,A24loc.2)
#            }
#            sizeofA24loc.2<-nrow(A24loc.growing.2)
#            
#            if(a2flag>10){
#              a2flag<-10
#            } 
#            
#            a2.fromto.xs<-A24loc.growing.2[seq((sizeofA24loc.2-a2flag),sizeofA24loc.2,by=1),"a2x0"]  
#            a2.fromto.ys<-A24loc.growing.2[seq((sizeofA24loc.2-a2flag),sizeofA24loc.2,by=1),"a2y0"]  
#            
#            lines(a2.fromto.xs,a2.fromto.ys,col=add.alpha("green",.6),lwd=2)
#            
#            a2flag<-a2flag+1
#            
#          }
#          
#          
#          for(vero in 1:1){
#            
#            check.2<-currentframe.A24[which(currentframe.A24$quad=="a2"),]
#            
#            if(nrow(check.2)>0){
#              
#              if(a2.pee.flag==0){
#                a2.pee.build<-check.2
#              } else {
#                a2.pee.build<-rbind(a2.pee.build,check.2)
#              }
#              par(new=TRUE)
#              plot(a2.pee.build$true.y~a2.pee.build$true.x,cex=log(a2.pee.build$Clustsize),
#                   xlim=c(1,640),ylim=c(-480,-1),bty="o",
#                   pch=16,
#                   xlab='',ylab='',
#                   axes=FALSE,
#                   col=add.alpha("darkgreen",0.09))  
#              
#              
#              
#              par(new=TRUE)
#              plot(check.2$true.y~check.2$true.x,cex=log(check.2$Clustsize),
#                   xlim=c(1,640),ylim=c(-480,-1),
#                   bty="o",
#                   pch=21,
#                   xlab='',ylab='',
#                   axes=FALSE,
#                   bg="green",
#                   col="green")
#              par(new=TRUE)
#              
#              a2.pee.flag<-a2.pee.flag+1
#            } else {
#              if(exists("a2.pee.build")){
#                par(new=TRUE)
#                plot(a2.pee.build$true.y~a2.pee.build$true.x,cex=log(a2.pee.build$Clustsize),
#                     xlim=c(1,640),ylim=c(-480,-1),bty="o",
#                     pch=16,
#                     xlab='',ylab='',
#                     axes=FALSE,
#                     col=add.alpha("darkgreen",0.09))  
#              }
#            }
#            
#          }
#          
#          #A4
#          ############
#          #if there are track points to plot at this time interval
#          if(!is.na(currentframe.A24$a4x0[1])){
#            
#            A24loc.4<-currentframe.A24[1,c("a4x0","a4y0")]
#            
#            if(a4flag==0){
#              A24loc.growing.4<-A24loc.4
#            } else {
#              A24loc.growing.4<-rbind(A24loc.growing.4,A24loc.4)
#            }
#            sizeofA24loc.4<-nrow(A24loc.growing.4)
#            
#            if(a4flag>10){
#              a4flag<-10
#            } 
#            
#            a4.fromto.xs<-A24loc.growing.4[seq((sizeofA24loc.4-a4flag),sizeofA24loc.4,by=1),"a4x0"]  
#            a4.fromto.ys<-A24loc.growing.4[seq((sizeofA24loc.4-a4flag),sizeofA24loc.4,by=1),"a4y0"]  
#            
#            lines(a4.fromto.xs,a4.fromto.ys,col=add.alpha("darkorange",.6),lwd=2)
#            
#            a4flag<-a4flag+1
#            
#          }
#      
#      
#          for(vero in 1:1){
#        
#        check.4<-currentframe.A24[which(currentframe.A24$quad=="a4"),]
#        
#        if(nrow(check.4)>0){
#          
#          if(a4.pee.flag==0){
#            a4.pee.build<-check.4
#          } else {
#            a4.pee.build<-rbind(a4.pee.build,check.4)
#          }
#          par(new=TRUE)
#          plot(a4.pee.build$true.y~a4.pee.build$true.x,cex=log(a4.pee.build$Clustsize),
#               xlim=c(1,640),ylim=c(-480,-1),bty="o",
#               pch=16,
#               xlab='',ylab='',
#               axes=FALSE,
#               col=add.alpha("darkgoldenrod3",0.09))  
#          
#          
#          
#          par(new=TRUE)
#          plot(check.4$true.y~check.4$true.x,cex=log(check.4$Clustsize),
#               xlim=c(1,640),ylim=c(-480,-1),
#               bty="o",
#               pch=21,
#               xlab='',ylab='',
#               axes=FALSE,
#               bg="darkorange",
#               col="darkorange")
#          par(new=TRUE)
#          
#          a4.pee.flag<-a4.pee.flag+1
#        } else {
#          if(exists("a4.pee.build")){
#            par(new=TRUE)
#            plot(a4.pee.build$true.y~a4.pee.build$true.x,cex=log(a4.pee.build$Clustsize),
#                 xlim=c(1,640),ylim=c(-480,-1),bty="o",
#                 pch=16,
#                 xlab='',ylab='',
#                 axes=FALSE,
#                 col=add.alpha("darkgoldenrod3",0.09))  
#          }
#        }
#        
#      }
#      
#      }
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      
#      dev.off()
#    }
#  }





