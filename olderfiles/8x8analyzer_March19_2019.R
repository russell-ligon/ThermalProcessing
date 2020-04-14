setwd("C:/Users/Rusty/Amazon Drive/MICE")
source("MouseFunctions.R", chdir = F)


#(p.a.x,p.a.y), (p.b.x,p.b.y), (p.c.x,p.c.y) each are a points with an x and y coordinate
#a = line point 1; b = line point 2; c = point to check against
#https://stackoverflow.com/questions/1560492/how-to-tell-whether-a-point-is-to-the-right-or-left-side-of-a-line
IsSide<-function(p.a.x,p.a.y, p.b.x,p.b.y, p.c.x, p.c.y){
  return ((p.b.x - p.a.x)*(p.c.y - p.a.y) - (p.b.y - p.a.y)*(p.c.x - p.a.x)) > 0;
}

dir<-"C:/Users/Rusty/Amazon Drive/MICE/Thermal/Next/8x8PeeData/Analyze"
setwd(dir)

HelperData<-read.csv("HelperData.csv")

folders<-list.dirs(dir,full.names=TRUE)# 
folders<-folders[-1]


#makeplots, not yet re-integrated with current pipeline
makeplots<-"yes"
makeplots<-"no"

trialflag<-1
#alltrials will be a list with n elements, corresponding to the # of folders analyzed (each corresponding to 1 trial)
# each element of alltrials will be a list, with m elements, corresponding to the # of days analyzed for that trial
alltrials<-list()
for(trial in 1:length(folders)){
  
  ClusterData<-list.files((folders[[trial]]),full.names=TRUE,pattern="(MASTER).*\\.Rdata$")# To date, I've named these megaframes "Summary_xxxxxx.csv"
  
  #Pull Trial-specific helper info
  TrialID<-strsplit(folders[[trial]],"/")[[1]][length(strsplit(folders[[trial]],"/")[[1]])]
  HelperTrial<-HelperData[which(HelperData$Trial==TrialID),]
  
  #Pull Trial-specific tracks
  trackfiles<-list.files(folders[[trial]],full.names = TRUE,pattern="(fixed).*\\.csv$")
  ROIfiles<-trackfiles[grepl("region",trackfiles)]
  trackfiles<-trackfiles[trackfiles!=ROIfiles]
  
  #Get directories with full-frame info
  FullFrameInfo<-list.dirs((folders[[trial]]),full.names=TRUE)# 
  FullFrameInfo<-FullFrameInfo[-1]
  A13.full.frame<-FullFrameInfo[grepl("13",FullFrameInfo)]
  A24.full.frame<-FullFrameInfo[grepl("24",FullFrameInfo)]
  full.frame.folder.list<-list(A13.full.frame,A24.full.frame)
  #Only run loop if folders containing BatchNames are in current folder
  if(length(full.frame.folder.list)>0){
    for(z in 1:(length(full.frame.folder.list))){
      batchfiles<-list.files(full.frame.folder.list[[z]],full.names=TRUE)
      for(z1 in 1:length(batchfiles)){
        frames<-read.csv(batchfiles[z1])
        frames<-frames[,2,drop=FALSE]
        if(z1==1){
          alldaframes<-frames
        } else {
          alldaframes<-rbind(alldaframes,frames)
        }
      }
      
      if(z==1){
        A13.frames<-alldaframes
      } else {
        A24.frames<-alldaframes
      }
    }
  }
  A13.frames$Name<-apply(A13.frames[,1,drop=FALSE],1,function(x) {strsplit(x,"ecord_2018-")[[1]][2]})#reference time
  A24.frames$Name<-apply(A24.frames[,1,drop=FALSE],1,function(x) {strsplit(x,"ecord_2018-")[[1]][2]})#reference time
  
  
  #Get names of Arena-specific, trial-specific files to read in
  A13s<-ClusterData[grepl("Arena1_3",ClusterData)]
  A24s<-ClusterData[grepl("Arena2_4",ClusterData)]
  Alist<-list(A13s,A24s)
 
  dayflag<-1
  daylist<-list()
    for(day in 1:length(A13s)){ #this length should equal the number of days
      #Set Day
      p.day<-ifelse(day==1,"D1",
                ifelse(day==2,"D2",
                    ifelse(day==3,"D3","D4")))
      
      
      #pull track,roi info
      alltracks<-trackfiles[grepl(p.day,trackfiles)]
      allrois<-ROIfiles[grepl(p.day,ROIfiles)]
      
      #Match up track info with full frame name info
      a1.track<-read.csv(alltracks[grepl("a1-",alltracks)]);a1.track$info<-as.character(a1.track$info);a1.track$info[3]<-as.character("a1")
        a1.start<-strsplit(a1.track$info[1],"ecord_2018-")[[1]][2]
        a1.track$position<-A13.frames$Name[c(which(A13.frames$Name==a1.start):(-1+nrow(a1.track)+which(A13.frames$Name==a1.start)))]
        if(length(allrois[grepl("a1-",allrois)])==1){
          a1.roi<-read.csv(allrois[grepl("a1-",allrois)]);colnames(a1.roi)[2]<-"a1.roi"
          a1.track$roi<-a1.roi$a1.roi
        } else {
          a1.track$roi<-NA
        }
        
      a2.track<-read.csv(alltracks[grepl("a2-",alltracks)]);a2.track$info<-as.character(a2.track$info);a2.track$info[3]<-as.character("a2")
        a2.start<-strsplit(a2.track$info[1],"ecord_2018-")[[1]][2]
        a2.track$position<-A24.frames$Name[c(which(A24.frames$Name==a2.start):(-1+nrow(a2.track)+which(A24.frames$Name==a2.start)))]
        if(length(allrois[grepl("a2-",allrois)])==1){
          a2.roi<-read.csv(allrois[grepl("a2-",allrois)]);colnames(a2.roi)[2]<-"a2.roi"
          a2.track$roi<-a2.roi$a2.roi
        } else {
          a2.track$roi<-NA
        }
        
      a3.track<-read.csv(alltracks[grepl("a3-",alltracks)]);a3.track$info<-as.character(a3.track$info);a3.track$info[3]<-as.character("a3")
       a3.start<-strsplit(a3.track$info[1],"ecord_2018-")[[1]][2]
       a3.track$position<-A13.frames$Name[c(which(A13.frames$Name==a3.start):(-1+nrow(a3.track)+which(A13.frames$Name==a3.start)))]
       if(length(allrois[grepl("a3-",allrois)])==1){
         a3.roi<-read.csv(allrois[grepl("a3-",allrois)]);colnames(a3.roi)[2]<-"a3.roi"
         a3.track$roi<-a3.roi$a3.roi
       } else {
         a3.track$roi<-NA
       }
       
      a4.track<-read.csv(alltracks[grepl("a4-",alltracks)]);a4.track$info<-as.character(a4.track$info);a4.track$info[3]<-as.character("a4")
        a4.start<-strsplit(a4.track$info[1],"ecord_2018-")[[1]][2]
        a4.track$position<-A24.frames$Name[c(which(A24.frames$Name==a4.start):(-1+nrow(a4.track)+which(A24.frames$Name==a4.start)))]
        if(length(allrois[grepl("a4-",allrois)])==1){
         a4.roi<-read.csv(allrois[grepl("a4-",allrois)]);colnames(a4.roi)[2]<-"a4.roi"
         a4.track$roi<-a4.roi$a4.roi
        } else {
          a4.track$roi<-NA
        }
      
      
      
      tracklist<-list(a1.track,a2.track,a3.track,a4.track)
      #pull track info to attach to pee info
      #creates dataframe called "combo_locations" with quad-specific location data
      for(TL in 1:length(tracklist)){
        focus<-tracklist[[TL]]
        trackerguide<-focus$info[c(1:3)]
        focus<-focus[,-4]
        focus<-focus[,c(2,3,1,4)]
        reftime<-strsplit(trackerguide[1],"ecord_2018-")[[1]][2]#reference time
        Htz<-trackerguide[2]
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
      
      
      #load cluster pee info
      Camera.13<-Alist[[1]]
        load(Camera.13[grepl(p.day,Camera.13)])
        Camera.13.day<-thisday;rm(thisday)
        Camera.13.day<-data.frame(Camera.13.day)
      Camera.24<-Alist[[2]]
        load(Camera.24[grepl(p.day,Camera.24)])
        Camera.24.day<-thisday;rm(thisday)
        Camera.24.day<-data.frame(Camera.24.day)
        
      
      #get frame info, camera13
      inf1<-(lapply(strsplit(rownames(Camera.13.day),"ecord_2018-"),function(x) {x[2]}))
      inf2<-(lapply(inf1,function(x){strsplit(x,"_")}))
      inf3<-(lapply(inf2,function(x) {x[[1]][2]}))
      hour<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][1]}))))
      minutes<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][2]}))))
      second<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][3]}))))
      dailyseconds<-second+(minutes*60)+(hour*60*60)
      frame.info<-paste(p.day,unlist(lapply(inf2,function(x) {x[[1]][2]})),sep=".")
      Camera.13.day$frameinfo<-frame.info
      Camera.13.day$dailysecond<-dailyseconds
      
      
      #get frame info, camera24
      inf1<-(lapply(strsplit(rownames(Camera.24.day),"ecord_2018-"),function(x) {x[2]}))
      inf2<-(lapply(inf1,function(x){strsplit(x,"_")}))
      inf3<-(lapply(inf2,function(x) {x[[1]][2]}))
      hour<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][1]}))))
      minutes<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][2]}))))
      second<-as.numeric(as.character(unlist(lapply(inf3,function(x) {(strsplit(x,"-"))[[1]][3]}))))
      dailyseconds<-second+(minutes*60)+(hour*60*60)
      frame.info<-paste(p.day,unlist(lapply(inf2,function(x) {x[[1]][2]})),sep=".")
      Camera.24.day$frameinfo<-frame.info
      Camera.24.day$dailysecond<-dailyseconds
      orig.singlecam.names<-colnames(Camera.24.day)
      
      ###################
      #unique colnames by camera
      colnames(Camera.13.day)<-paste(colnames(Camera.13.day),".A13",sep='');colnames(Camera.13.day)[103]<-"dailysecond"
      colnames(Camera.24.day)<-paste(colnames(Camera.24.day),".A24",sep='');colnames(Camera.24.day)[103]<-"dailysecond"
      
      #this merge temporally aligns the data from both cameras, then we can resplit and process separately
      wook<-merge(Camera.13.day,Camera.24.day,by="dailysecond",all=TRUE)
      keepframenum<-wook[,"dailysecond",drop=FALSE]
      
      #Full time (noon to 10pm)
      fulltimerunner<-as.data.frame(seq((12*60*60),(22*60*60),1));colnames(fulltimerunner)[1]<-"dailysecond"
      
      #merges wook to fulltimerunner
      wook2<-merge(fulltimerunner,wook,by="dailysecond",all.x=TRUE)
      
      #if row numbers don't match up, get rid of duplicates by taking larger clusters per time
      if(nrow(wook2)!=nrow(fulltimerunner)){
        A<-as.character(fulltimerunner$dailysecond)
        B<-as.character(wook2$dailysecond)
        C<-B[duplicated(B)]
        D<-wook2[wook2$dailysecond %in% C,]#dataframe of duplicated (timewise) rows
          for(q in 1:length(C)){
            q1<-D[which(D$dailysecond==C[q]),]
            a1<-q1$Clustsize.A13[1];b1<-q1$Clustsize.A13[2]
            
            if((is.na(a1) & is.na(b1)) | ((!is.na(a1) & is.na(b1)))){
              #if both are NA, or 1 is a num and 2 is na
              plugback<-q1[1,]
            } else {
              plugback<-q1[2,]
            }
            
            if(q==1){ #first iteration of loop, do this
              keeper<-plugback
            } else { #rbind on subsequent iterations
              keeper<-rbind(keeper,plugback)
            }
          }
        
        E<-wook2[!(wook2$dailysecond %in% C),] #full dataframe, minus duplicated times
        wook3<-rbind(E,keeper)
        wook4<-wook3[order(wook3$dailysecond),]
        wook2<-wook4
      }
      ####################################
      positionANDpee<-merge(wook2,combo_locations,by="dailysecond",all.x=TRUE)
      
      ####################################
      alignedA13<-positionANDpee[,c(1,206:208,212:214,103,2:102)]
      colnames(alignedA13)[c(9:109)]<-orig.singlecam.names[c(1:101)]
      alignedA24<-positionANDpee[,c(1,209:211,215:217,205,104:204)]
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
      bseq<-seq(9,108,4)
      for(b in 1:25){
        print(paste(b," of 25",sep=''))
        H<-bseq[b]+3
        if(b==1){
          longA13<- trim.alignedA13[,c(1:8,bseq[b]:H)]
          longA24<- trim.alignedA24[,c(1:8,bseq[b]:H)]
          
        } else {
          ta13<-trim.alignedA13[,c(1:8,bseq[b]:H)]
          colnames(ta13)<-colnames(longA13)
          longA13<-rbind(longA13,ta13)
          
          ta24<-trim.alignedA24[,c(1:8,bseq[b]:H)]
          colnames(ta24)<-colnames(longA24)
          longA24<-rbind(longA24,ta24)
        }
      }
      
      longA13<-longA13[order(longA13$dailysecond),]
      longA24<-longA24[order(longA24$dailysecond),]
      
      full.long.A13<-longA13 #keep as 'whole'
      full.long.A24<-longA24
      
      #droping by NAs in these named columns no longer works??
     # longA13.peeonly<-na.omit(longA13,cols=c("Clustsize","mTemp","true.y","true.x"))
     # longA24.peeonly<-na.omit(longA24,cols=c("Clustsize","mTemp","true.y","true.x"))
      
      longA13.peeonly<-longA13
      longA13.peeonly$nNas<-apply(longA13,1,function(x) sum(is.na(x)))
      longA13.peeonly<-longA13.peeonly[which(longA13.peeonly$nNas<5),]
      
      longA24.peeonly<-longA24
      longA24.peeonly$nNas<-apply(longA24,1,function(x) sum(is.na(x)))
      longA24.peeonly<-longA24.peeonly[which(longA24.peeonly$nNas<5),]
      
      
      longA13<-longA13.peeonly
      longA24<-longA24.peeonly
      
      # alignedA13[,seq(4,100,4)]<-(-1)*alignedA13[,seq(4,100,4)]#A13, make y-axis negative (0,-480)
      # alignedA24[,seq(4,100,4)]<-abs(alignedA24[,seq(4,100,4)]-480)#A24, flip y-axis, turn 1 to 479, and turn 480 to 0
      
      A13.reference<-HelperTrial[which(HelperTrial$Arena=="13"),]
      A24.reference<-HelperTrial[which(HelperTrial$Arena=="24"),]
      
      
      
      if(A13.reference$TopScreen=="inner"){
        longA13$quad<-ifelse(IsSide(A13.reference$topdivide.x,A13.reference$topdivide.y,
                             A13.reference$botdivide.x,A13.reference$botdivide.y,
                             longA13$true.x,longA13$true.y)<0,
                             "Right",
                             "Left")
        #keep below left top
        longA13$keep<-ifelse(longA13$quad=="Left",
                              ifelse(IsSide(A13.reference$toplefthoriz.x,A13.reference$toplefthoriz.y,
                                    A13.reference$topdivide.x,A13.reference$topdivide.y,
                                    longA13$true.x,longA13$true.y)>0,
                                    "in","out"),NA)
        #Keep below right top
        longA13$keep<-ifelse(longA13$quad=="Left",longA13$keep,
                              ifelse(IsSide(A13.reference$topdivide.x,A13.reference$topdivide.y,
                                            A13.reference$toprighthoriz.x,A13.reference$toprighthoriz.y,
                                            longA13$true.x,longA13$true.y)>0,"in","out"))
                              

        keeplongA13<-longA13[which(longA13$keep=="in"),]
        
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
        #keep below left top
        longA13$keep<-ifelse(longA13$quad=="Left",
                             ifelse(IsSide(A13.reference$botlefthoriz.x,A13.reference$botlefthoriz.y,
                                           A13.reference$botdivide.x,A13.reference$botdivide.y,
                                           longA13$true.x,longA13$true.y)>0,
                                    "out","in"),NA)
        #Keep below right top
        longA13$keep<-ifelse(longA13$quad=="Left",longA13$keep,
                             ifelse(IsSide(A13.reference$botdivide.x,A13.reference$botdivide.y,
                                           A13.reference$botrighthoriz.x,A13.reference$botrighthoriz.y,
                                           longA13$true.x,longA13$true.y)>0,"out","in"))
        
        
        keeplongA13<-longA13[which(longA13$keep=="in"),]
        
        keeplongA13$quad<-ifelse(keeplongA13$quad=="Left",
                                 as.character(A13.reference$ID.left),
                                 as.character(A13.reference$ID.right))
        
      }
      
      if(A24.reference$TopScreen=="inner"){
        longA24$quad<-ifelse(IsSide(A24.reference$topdivide.x,A24.reference$topdivide.y,
                                    A24.reference$botdivide.x,A24.reference$botdivide.y,
                                    longA24$true.x,longA24$true.y)<0,
                             "Right",
                             "Left")
        #keep below left top
        longA24$keep<-ifelse(longA24$quad=="Left",
                             ifelse(IsSide(A24.reference$toplefthoriz.x,A24.reference$toplefthoriz.y,
                                           A24.reference$topdivide.x,A24.reference$topdivide.y,
                                           longA24$true.x,longA24$true.y)>0,
                                    "in","out"),NA)
        #Keep below right top
        longA24$keep<-ifelse(longA24$quad=="Left",longA24$keep,
                             ifelse(IsSide(A24.reference$topdivide.x,A24.reference$topdivide.y,
                                           A24.reference$toprighthoriz.x,A24.reference$toprighthoriz.y,
                                           longA24$true.x,longA24$true.y)>0,"in","out"))
        
        
        keeplongA24<-longA24[which(longA24$keep=="in"),]
        
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
        #keep below left top
        longA24$keep<-ifelse(longA24$quad=="Left",
                             ifelse(IsSide(A24.reference$botlefthoriz.x,A24.reference$botlefthoriz.y,
                                           A24.reference$botdivide.x,A24.reference$botdivide.y,
                                           longA24$true.x,longA24$true.y)>0,
                                    "out","in"),NA)
        #Keep below right top
        longA24$keep<-ifelse(longA24$quad=="Left",longA24$keep,
                             ifelse(IsSide(A24.reference$botdivide.x,A24.reference$botdivide.y,
                                           A24.reference$botrighthoriz.x,A24.reference$botrighthoriz.y,
                                           longA24$true.x,longA24$true.y)>0,"out","in"))
        
        
        keeplongA24<-longA24[which(longA24$keep=="in"),]
        
        keeplongA24$quad<-ifelse(keeplongA24$quad=="Left",
                                 as.character(A24.reference$ID.left),
                                 as.character(A24.reference$ID.right))
        
      }
      
      #USE keeplongA13, and keeplongA24 here, but incorporate mouse location and roi earlier
      keeplongA13$distance<-ifelse(keeplongA13$quad=="A1",
                                   sqrt((keeplongA13$a1x0-keeplongA13$true.x)^2 + (keeplongA13$a1y0-keeplongA13$true.y)^2),
                                   sqrt((keeplongA13$a3x0-keeplongA13$true.x)^2 + (keeplongA13$a3y0-keeplongA13$true.y)^2))
    
      keeplongA24$distance<-ifelse(keeplongA24$quad=="A2",
                                   sqrt((keeplongA24$a2x0-keeplongA24$true.x)^2 + (keeplongA24$a2y0-keeplongA24$true.y)^2),
                                   sqrt((keeplongA24$a4x0-keeplongA24$true.x)^2 + (keeplongA24$a4y0-keeplongA24$true.y)^2))
      
      keeplongA13$quad<-as.factor(keeplongA13$quad)
      keeplongA24$quad<-as.factor(keeplongA24$quad)
      keeplongA13$a1roi<-as.factor(keeplongA13$a1roi)
      keeplongA13$a3roi<-as.factor(keeplongA13$a3roi)
      keeplongA24$a2roi<-as.factor(keeplongA24$a2roi)
      keeplongA24$a4roi<-as.factor(keeplongA24$a4roi)
      
      A1data<-subset(keeplongA13,quad=="A1")
      A3data<-subset(keeplongA13,quad=="A3")
      A2data<-subset(keeplongA24,quad=="A2")
      A4data<-subset(keeplongA24,quad=="A4")
      
      thisday<-list(A1data,A2data,A3data,A4data)
      names(thisday)<-c("A1data","A2data","A3data","A4data")
      
      daylist[[dayflag]]<-thisday
      names(daylist)[dayflag]<-p.day
   
      # summary(A1data)
      # temp<-(as.data.frame(summary(A1data$a1roi)));colnames(temp)<-"n"
      # temp[order(temp$n),,drop=FALSE]
      # summary(keeplongA24)
 
      dayflag<-dayflag+1
   } #day
  alltrials[[trialflag]]<-daylist
  names(alltrials)[trialflag]<-TrialID
  trialflag<-trialflag+1
}#trial
