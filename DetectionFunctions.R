



###############################################################################
#Get decimal level modal values
modefun2a<-function(x){
  ux = sort(unique(x))
  idx = match(x, ux)
  n = tabulate(idx, nbins=length(ux))
  df = data.frame(x=ux, n=n)
  ifelse(length(df[df$n ==  max(df$n),'x'])==1,df[df$n ==  max(df$n),'x'],df[df$n ==  max(df$n),'x'][1])
  #this *ifelse* fixes an earlier problem, where multiple values could be returned---effectively garbling all subsequent values
}


##
#clips matrix based on polygon coordinates, and subtracts side/polygon-specific median temp
clipmatrix<-function(inputlist,Lpolygoncoords,Rpolygoncoords,allposscombos){
  
  L.keepORout<-ifelse(point.in.polygon(allposscombos[,2],
                                       allposscombos[,1],
                                       Lpolygoncoords$X,
                                       Lpolygoncoords$Y)>0,1,0)
  
  R.keepORout<-ifelse(point.in.polygon(allposscombos[,2],
                                       allposscombos[,1],
                                       Rpolygoncoords$X,
                                       Rpolygoncoords$Y)>0,1,0)

  L.info<-cbind(allposscombos,L.keepORout);colnames(L.info)<-c("y","x","z")
  L.keepmatrix<- xtabs(z~y+x, data=L.info)
  L.keepmatrix2<-matrix(L.keepmatrix,nrow=480,ncol=640)
  
  R.info<-cbind(allposscombos,R.keepORout);colnames(R.info)<-c("y","x","z")
  R.keepmatrix<- xtabs(z~y+x, data=R.info)
  R.keepmatrix2<-matrix(R.keepmatrix,nrow=480,ncol=640)
  
  L.fullyclipped<-lapply(inputlist,function(x){
    big<-x
    big<-big*L.keepmatrix2
    big[big==0]<-NA
    med.L<-median(big,na.rm=TRUE) #subtract side/polygon-specific median temp
    big<-big-med.L
    return(big)})
  
  R.fullyclipped<-lapply(inputlist,function(x){
    big<-x
    big<-big*R.keepmatrix2
    big[big==0]<-NA
    med.R<-median(big,na.rm=TRUE) #subtract side/polygon-specific median temp
    big<-big-med.R
    return(big)})
  
  fullyclipped<-list(L.fullyclipped,R.fullyclipped)
  return(fullyclipped)
}



###############
#runs triple-step detector
# a) clusters *warm* pixels and labels the biggest cluster and/or any clusters over 195 as "mouse"
# b) identifies *warming* pixels (using frame[x]-frame[x-PrevFrame])
# c) then removes *warming* pixels that were also in the "mouse" cluster in the prev step
# d) then clusters the remaining pee pixels,and returns a 480x640 matrix with different numbers for each pee cluster

#triplestepdetector(L.list[[rory]],L.list[[prev]],L.list[[futuro]],150)
triplestepdetector<-function(thiswholeframe,backwholeframe,futurewholeframe,biggestexclude,PC1){
  #finds *warming* pixels in the last 20 seconds
  diffframe<-thiswholeframe-backwholeframe
  diffframe[(diffframe)<temprise]<-NA #drop all pixels that haven't warmed by at least 1 degree
  diffframe[!is.na(diffframe)]<-1 #sets all remaining, warming, pixels to 1
  
  #finds *warm* pixels, those that are more than 2 degrees higher than the frame median
  analyzewarm<-thiswholeframe
  analyzewarm[(analyzewarm)<2]<-NA
  analyzewarm[!is.na(analyzewarm)]<-1 #sets all remaining, warming, pixels to 1
  
  #finds *cool* pixels, that are more than "tempdrop" degrees BELOW than the frame median
  analyzecool<-futurewholeframe
  analyzecool[(analyzecool)>(tempdrop*-1)]<-NA
  analyzecool[!is.na(analyzecool)]<-1 #sets all remaining, cool, pixels to 1
  
  #if there are any warm pixels (2+), run with it, else fill with NAs
  if(sum(analyzewarm,na.rm=TRUE)>1){
    
    #warm step
    Coords.for.warm<-which(analyzewarm==1,arr.ind=TRUE)#Gets X/Y coordinates for pee pixels
    #Coords.for.warm<-cbind(Coords.for.warm,axualtemps)# Combines actual temps with X/Y coordinates for each pee pixel
    warm.dist.matrix<-dist(Coords.for.warm[,c(1,2)])#uses x/y coordinates to create distance matrix
    #Clusters pixels based on distance matrix
    warmcluster <-fastcluster::hclust(warm.dist.matrix, method="single", members=NULL) #Clusters pee based on distance matrix
    
    warmgroups<-cutree(warmcluster,h=specifdist)
    
    warmgroupinfo<-data.frame(table(warmgroups))
    #drop any warm clusters bigger than 190 pixels OR is the biggest cluster
    warmgroupinfo$mouse<-ifelse((warmgroupinfo$Freq>biggestexclude | warmgroupinfo$Freq==max(warmgroupinfo$Freq)),1,0)
    micecoords<-as.integer(warmgroupinfo$warmgroups[which(warmgroupinfo$mouse==1)])
    
    Coords.for.warm<-cbind(Coords.for.warm,warmgroups)#Determines cutoff for what is a cluster (based on 'specificDistance' based on folder title (MB1,MB2))
    colnames(Coords.for.warm)[3]<-"warmgroups"
    Coords.for.warm<-data.frame(Coords.for.warm)
    Coords.for.warm$mouse<-match(Coords.for.warm$warmgroups,micecoords)
    Coords.for.warm$mouse<-ifelse(is.na(Coords.for.warm$mouse),0,1)
    Coords.for.warm$mouse<-ifelse(Coords.for.warm$mouse==1,0,1)
    
    
    #create matrix with cells having zero
    mousezero<-Coords.for.warm[,c(1,2,4)];colnames(mousezero)<-c("y","x","z")
    fullup<-merge(PC1,mousezero,by.x=c("row","col"),by.y=c("y","x"),all.x=TRUE)
    colnames(fullup)[c(1:2)]<-c("y","x")
    fullup$z[is.na(fullup$z)]<-(-1)#have to make NAs something else for converting to correct matrix structure below
    mousematrix<- xtabs(z~y+x, data=fullup)
    mousematrix2<-matrix(mousematrix,nrow=480,ncol=640)
    mousematrix2[mousematrix2==(-1)]<-NA
  } else {
    mousematrix2<-matrix(NA,nrow=480,ncol=640)
  }
  
  
  #######################
  #mousematrix2 = pixels of warm huts or the mouse
  # diffframe = pixels that have warmed 1 degree since PastFrame
  
  peepixels<-mousematrix2*diffframe #multiply because any mouse/hut pixel in mousematrix2 is equal to 0 and therefore gets cancelled out by this multiplication
  peepixels[peepixels==0]<-NA
  
  
  #
  Coords.for.pee<-which(peepixels==1,arr.ind=TRUE)#Gets X/Y coordinates for pee pixels
  
  #IF EITHER the number of rows for Coords.for.pee is zero (no coords that were warm, and warmed that weren't mice/huts)
  # OR the sum of analyzecool (matrix with 1s for pixels that cooled) indicating NO future cooling (below tempdrop), 
  # then fill with NAs, ELSE analyze
  if(nrow(Coords.for.pee)==0 | sum(analyzecool,na.rm=TRUE)==0){
    peematrix2<-matrix(NA,nrow=480,ncol=640)
  } else {
    if(nrow(Coords.for.pee)>1){
      pee.dist.matrix<-dist(Coords.for.pee[,c(1,2)])#uses x/y coordinates to create distance matrix
      #Clusters pixels based on distance matrix
      peecluster <-fastcluster::hclust(pee.dist.matrix, method="single", members=NULL) #Clusters pee based on distance matrix
      peegroups<-cutree(peecluster,h=specifdist)
      peegroupinfo<-data.frame(table(peegroups))
      
      Coords.for.pee<-cbind(Coords.for.pee,peegroups);colnames(Coords.for.pee)<-c("y","x","z")
      
      
      ####################################
      ###############
      #cool step
      Coords.for.cool<-which(analyzecool==1,arr.ind=TRUE)#Gets X/Y coordinates for pee pixels
      cool.dist.matrix<-dist(Coords.for.cool[,c(1,2)])#uses x/y coordinates to create distance matrix
      #Clusters pixels based on distance matrix
      coolcluster <-fastcluster::hclust(cool.dist.matrix, method="single", members=NULL) #Clusters pee based on distance matrix
      
      coolgroups<-cutree(coolcluster,h=specifdist)#Determines cutoff for what is a cluster (based on 'specificDistance' based on folder title (MB1,MB2))
      
      coolgroupinfo<-data.frame(table(coolgroups))
      
      Coords.for.cool<-cbind(Coords.for.cool,coolgroups)
      Coords.for.cool<-data.frame(Coords.for.cool)
      
      
      #Combine Coords.for.pee and Coords.for.cool, which creates a new dataframe,
      # with y, x coordinates, z values (where 1 = warm spot that IS NOT a mouse)
      # and 'coolgroups' where cool-in-the-future pixels have a number, and non-cool-in-the-future pixels
      # have NA
      CompositeCoords<-merge(Coords.for.pee,Coords.for.cool,by.x=c("y","x"),by.y=c("row","col"),all=TRUE)
      
      #drop pixels/rows that don't have a 1 (non-mouse warming) AND a numeric value indicating cooling
      CompositeCoords<-na.omit(CompositeCoords)
      
      
      #if the NA dropped version of this matrix has 0 rows (i.e. there are no pixels that were warmed, then cooled),
      # just create an empty matrix, ELSE 
      if(nrow(CompositeCoords)==0){
        peematrix2<-matrix(NA,nrow=480,ncol=640)
      } else {
        
        warmthencool<-CompositeCoords[,c("y","x","z")]
        
        fullpee<-merge(PC1,warmthencool,by.x=c("row","col"),by.y=c("y","x"),all.x=TRUE)
        colnames(fullpee)[c(1:2)]<-c("y","x")
        fullpee$z[is.na(fullpee$z)]<-(-1)#have to make NAs something else for converting to correct matrix structure below
        peematrix<- xtabs(z~y+x, data=fullpee)
        peematrix2<-matrix(peematrix,nrow=480,ncol=640)
        peematrix2[peematrix2==(-1)]<-NA
      }
      
    }
    if(nrow(Coords.for.pee)==1){
      Coords.for.pee<-cbind(Coords.for.pee,1);colnames(Coords.for.pee)<-c("y","x","z")
      
      
      ###############
      #cool step
      Coords.for.cool<-which(analyzecool==1,arr.ind=TRUE)#Gets X/Y coordinates for pee pixels
      cool.dist.matrix<-dist(Coords.for.cool[,c(1,2)])#uses x/y coordinates to create distance matrix
      #Clusters pixels based on distance matrix
      coolcluster <-fastcluster::hclust(cool.dist.matrix, method="single", members=NULL) #Clusters pee based on distance matrix
      
      coolgroups<-cutree(coolcluster,h=specifdist)#Determines cutoff for what is a cluster (based on 'specificDistance' based on folder title (MB1,MB2))
      
      coolgroupinfo<-data.frame(table(coolgroups))
      
      Coords.for.cool<-cbind(Coords.for.cool,coolgroups)
      Coords.for.cool<-data.frame(Coords.for.cool)
      
      
      #Combine Coords.for.pee and Coords.for.cool, which creates a new dataframe,
      # with y, x coordinates, z values (where 1 = warm spot that IS NOT a mouse)
      # and 'coolgroups' where cool-in-the-future pixels have a number, and non-cool-in-the-future pixels
      # have NA
      CompositeCoords<-merge(Coords.for.pee,Coords.for.cool,by.x=c("y","x"),by.y=c("row","col"),all=TRUE)
      
      #drop pixels/rows that don't have a 1 (non-mouse warming) AND a numeric value indicating cooling
      CompositeCoords<-na.omit(CompositeCoords)
      
      
      #if the NA dropped version of this matrix has 0 rows (i.e. there are no pixels that were warmed, then cooled),
      # just create an empty matrix, ELSE 
      if(nrow(CompositeCoords)==0){
        peematrix2<-matrix(NA,nrow=480,ncol=640)
      } else {
        
        fullpee<-merge(PC1,CompositeCoords,by.x=c("row","col"),by.y=c("y","x"),all.x=TRUE)
        colnames(fullpee)[c(1:2)]<-c("y","x")
        fullpee$z[is.na(fullpee$z)]<-(-1)#have to make NAs something else for converting to correct matrix structure below
        peematrix<-xtabs(z~y+x, data=fullpee)
        peematrix2<-matrix(peematrix,nrow=480,ncol=640)
        peematrix2[peematrix2==(-1)]<-NA
      }
    }
  }
  
  
  
  
  
  
  #return(peematrix2)
  mm2<-mousematrix2
  pm2<-peematrix2
  
  #swaps (mice = 0, at first (for multiplying detections out), we convert these to 100)
  mousematrix2[(mousematrix2==0)]<-100
  mousematrix2[(mousematrix2==1)]<-0
  mousematrix2[is.na(mousematrix2)]<-0
  
  peematrix2[is.na(peematrix2)]<-0   
  
  newmatrix<-peematrix2+mousematrix2
  newmatrix[newmatrix==(0)]<-NA
  return(newmatrix)
}

############
#SingleWarmCluster(L.list[[rory]],L.list[[prev]],150)
# (looks for warm clusters, useful for identifying big warm clusters so that they can be marked as mice/huts)
SingleWarmCluster<-function(thiswholeframe,biggestexclude,PC2){
  
  #finds *warm* pixels, those that are more than 2 degrees higher than the frame median
  analyzewarm<-thiswholeframe
  analyzewarm[(analyzewarm)<2]<-NA
  analyzewarm[!is.na(analyzewarm)]<-1 #sets all remaining, warm, pixels to 1
  
  if(sum(analyzewarm,na.rm=TRUE)>1){
    
    #warm step
    Coords.for.Swarm<-which(analyzewarm==1,arr.ind=TRUE)#Gets X/Y coordinates for pee pixels
   
    
    #if there are more than 10,000 pixels that are 2 deg above median, then a person is in frame,
    # and we should ignore it, by assigning every pixel a value of 1, 
    # essentially saying that all pixels in frame are not to be used
    # ... this is a way to handle a person being in frame
    if(nrow(Coords.for.Swarm)>10000){
      mousematrix2<-matrix(1,nrow=480,ncol=640)
    } else {
           #Coords.for.Swarm<-cbind(Coords.for.Swarm,axualtemps)# Combines actual temps with X/Y coordinates for each pee pixel
    warm.dist.matrix<-dist(Coords.for.Swarm[,c(1,2)])#uses x/y coordinates to create distance matrix
    #Clusters pixels based on distance matrix
    warmcluster <-fastcluster::hclust(warm.dist.matrix, method="single", members=NULL) #Clusters pee based on distance matrix
    
    warmgroups<-cutree(warmcluster,h=specifdist)
    
    warmgroupinfo<-data.frame(table(warmgroups))
    #drop any warm clusters bigger than 190 pixels OR is the biggest cluster
    warmgroupinfo$mouse<-ifelse((warmgroupinfo$Freq>biggestexclude | warmgroupinfo$Freq==max(warmgroupinfo$Freq)),1,0)
    micecoords<-as.integer(warmgroupinfo$warmgroups[which(warmgroupinfo$mouse==1)])
    
    Coords.for.Swarm<-cbind(Coords.for.Swarm,warmgroups)#Determines cutoff for what is a cluster (based on 'specificDistance' based on folder title (MB1,MB2))
    colnames(Coords.for.Swarm)[3]<-"warmgroups"
    Coords.for.Swarm<-data.frame(Coords.for.Swarm)
    Coords.for.Swarm$mouse<-match(Coords.for.Swarm$warmgroups,micecoords)
    Coords.for.Swarm$mouse<-ifelse(is.na(Coords.for.Swarm$mouse),0,1)
    Coords.for.Swarm$mouse<-ifelse(Coords.for.Swarm$mouse==1,0,1)
    
    #create matrix with cells having zero
    mousezero<-Coords.for.Swarm[,c(1,2,4)];colnames(mousezero)<-c("y","x","z")
    fullup<-merge(PC2,mousezero,by.x=c("row","col"),by.y=c("y","x"),all.x=TRUE)
    colnames(fullup)[c(1:2)]<-c("y","x")
    fullup$z[is.na(fullup$z)]<-(-1)#have to make NAs something else for converting to correct matrix structure below
    mousematrix<- xtabs(z~y+x, data=fullup)
    mousematrix2<-matrix(mousematrix,nrow=480,ncol=640)
    mousematrix2[mousematrix2==(-1)]<-NA
    
    #######################
    #mousematrix2 = pixels of warm huts or the mouse
    
    mousematrix2[(mousematrix2==0)]<-100
    mousematrix2[(mousematrix2==1)]<-0
    mousematrix2[is.na(mousematrix2)]<-0
    }
    

  } else {
    mousematrix2<-matrix(NA,nrow=480,ncol=640)
  }
  return(mousematrix2)
}

############
#function to rearange vector into matrix, so we can apply over whole dataframe
thermalarrange<-function(arow,nr,nc){
  matrix(arow,nrow=nr,ncol=nc,byrow=TRUE)
}


############
# Recombine pee from L and R, but R is +20
LRCombine<-function(L1,R1){
  #L1<-Left.mega.pee.list[[pbj]]
  L1[is.na(L1)]<-0
  
  #R1<-Right.mega.pee.list[[pbj]]
  R1<-R1+20
  R1[is.na(R1)]<-0
  
  newmatrix<-L1+R1
  newmatrix[newmatrix==(0)]<-NA
  return(newmatrix)
}

############################
#########
#Function that incorporates most of the processing steps and functions for pee detection (run it all)
runitall<-function(submainframe,Lpolygoncoords,Rpolygoncoords){
  
  #Generates 307200 rows where the values in column 1 are possible y coordinates (1-480)
  #and the values in column 2 are possible x coordinates (1-640)
  possCombo<-expand.grid(seq(1,480,1),seq(1,640,1));colnames(possCombo)<-c("row","col")
  possCombo[,1]<-as.integer(possCombo[,1])
  possCombo[,2]<-as.integer(possCombo[,2])
  
  #turn into 3 matrices
  part.main.list<-list()
  for(fff in 1:3){
    part.main.list[[fff]]<-thermalarrange(submainframe[fff,],480,640)
  }
  rm(submainframe)
  gc()
  
  #clip matrix takes each frame in part.main.list, clips it to the correct dimensions, and subtracts the side-specific median temp
  LR.list<-clipmatrix(inputlist=part.main.list,Lpolygoncoords,Rpolygoncoords,allposscombos=possCombo)
  
  L.list<-LR.list[[1]]
  R.list<-LR.list[[2]]
  
  rm(LR.list)
 
  gc()


  ####################
  ##############################################################################################
  ##############################################################################################
  #FIND BIG WARM CLUSTERS BEFORE CLIPPING, assigns these pixels value of 300
  currentmedian<-median(part.main.list[[2]])
  First.Cluster.Frame<- SingleWarmCluster(thiswholeframe=(part.main.list[[2]]-currentmedian),biggestexclude=100,PC2=possCombo)
  
  First.Cluster.Frame[First.Cluster.Frame<100]<-NA
  First.Cluster.Frame[!is.na(First.Cluster.Frame)]<-300
  First.Cluster.Frame[is.na(First.Cluster.Frame)]<-0
  
  

  
  rm(part.main.list)
  ##############################################################################################
  #########################################################
  
  Left.mega.pee.list<- triplestepdetector(thiswholeframe=L.list[[2]],backwholeframe=L.list[[1]],
                                          futurewholeframe=L.list[[3]],
                                          biggestexclude=100,PC1=possCombo)
  Left.mega.pee.list[is.na(Left.mega.pee.list)]<-0
  Left.mega.pee.list<-Left.mega.pee.list+First.Cluster.Frame
  Left.mega.pee.list[Left.mega.pee.list<=0]<-NA
  Left.mega.pee.list[Left.mega.pee.list==400]<-100
  Left.mega.pee.list[Left.mega.pee.list==300]<-NA
  
  Right.mega.pee.list<- triplestepdetector(thiswholeframe=R.list[[2]],backwholeframe=R.list[[1]],
                                           futurewholeframe=R.list[[3]],
                                           biggestexclude=100,PC1=possCombo)
  Right.mega.pee.list[is.na(Right.mega.pee.list)]<-0
  Right.mega.pee.list<-Right.mega.pee.list+First.Cluster.Frame
  Right.mega.pee.list[Right.mega.pee.list<=0]<-NA
  Right.mega.pee.list[Right.mega.pee.list==400]<-100
  Right.mega.pee.list[Right.mega.pee.list==300]<-NA
  
  rm(L.list)
  rm(R.list)
  rm(First.Cluster.Frame)
  ##############################################################################################
  ##############################################################################################
  #####################################################
  #####################################################
  #####################################################
  #####################################################
  
  Mega.pee.list<- LRCombine(L1=Left.mega.pee.list,R1=Right.mega.pee.list)
  
  rm(Left.mega.pee.list)
  rm(Right.mega.pee.list)
  ##############################################################################################
  #######################################################################################
  #############################################################################
  ###############################################################
  #####################################################
  
  return(Mega.pee.list)
  gc()
}


#####
#listpeedetails 
# function that processes the lists (Left.mega.pee.list, Right.mega.pee.list) produced earlier
# returning a 1x101 vector per frame (list element), corresponding to nclusts (total for that frame),
# then clustsize and clust coordinates (3 values), repeating for all detected clusters
listpeedetails<-function(worklistelement){
  
  #a
  Coords.for.pee<-which(worklistelement>0 & worklistelement<50,arr.ind=TRUE)#Gets X/Y coordinates for pee pixels
  #If no pee clusters....
  if(nrow(Coords.for.pee)<1 ){
    peeclusts<-NA
  } else {#if some pee clusters...
    
    Coords.for.pee<-cbind(Coords.for.pee,worklistelement[which(worklistelement>0 & worklistelement<50)])
    colnames(Coords.for.pee)[3]<-c("peeID")
    
    gmon<-Coords.for.pee[,c(3,1,2),drop=FALSE]#combines cluster ID to coordinates/temp for each pixel
    
    if(nrow(gmon)>1)
      gmon<-gmon[order(gmon[,1]),]#orders s
    
    
    nclusts<-length(unique(gmon[,1]))#How many clusters of pee?
    clustIDs<-rev(order(c(table(gmon[,1])))) #orders cluster IDs by size (i.e. number of pixels)
    clustIDs<-unique(gmon[,1])[clustIDs]
    
    #Finds centroids of each pee-cluster and binds them into new df called 'peeclusts'
    # which is alternating x, y coordinates for each pee cluster
    for(vero in 1:nclusts){
      clust<-gmon[which(gmon[,1]==clustIDs[vero]),,drop=FALSE]
      centroid<-c(round(mean(clust[,2])),round(mean(clust[,3])))
      if(vero<2){
        peeclusts<-c(nclusts,nrow(clust),clust[1,1],centroid)
      } else {
        peeclusts<-c(peeclusts,nrow(clust),clust[1,1],centroid)
      }
    }
  }  
  
  #if peeclusts was successfully defined above, then it 'exists'
  if(exists("peeclusts")){
    
    #if it is of length 1, then it was defined above as having no clusters, or only a cluster with a single pixel of pee
    if(length(peeclusts)==1){
      peeclusts<-rep(NA,101)
    } else {
      #this is the 'typical' case, where peeclusts was defined above, and this section adds NAs so that it is 101 elements long
      #but the <100 just serves to take care of instances where there might be lots of pee detected, creating a negative value when
      #we subtract the length of peeclusts from 101, so we cap peeclusts at 101 elements
      if(length(peeclusts)<100){
        peeclusts<-c(peeclusts,rep(NA,101-length(peeclusts)))
      } else {
        peeclusts<-peeclusts[c(1:101)]
      }
    }
    
    
    
  } else {
    peeclusts<-rep(NA,101)
    #otherwise, use a vector of 101 NAs
  }
  peeclusts<-t(data.frame(peeclusts))
  colnames(peeclusts)<-c("Nclusts",rep(c("Clustsize","mTemp","x1","y1"),25))
  rownames(peeclusts) <- NULL
  return(peeclusts)
}


#COOL FUNCTIONS
#####################################################
##############################
newcooldetector<-function(backwholeframe,thiswholeframe,biggestexclude,PC1){
 
  
  
  
  #finds pixels that have changed
  # positive values indicate pixels that are NOW warmer than they were 180 (or PastFrame) seconds ago
  # negative values indicate pixels that have cooled (i.e. are cooler now than they were 180 sec ago)
  diffframe<-thiswholeframe-backwholeframe
  diffframe[(diffframe)>(tempdrop*-1)]<-NA #drop all pixels that haven't cooled by at least 1 degree
  diffframe[!is.na(diffframe)]<-1 #sets all remaining, cooled, pixels to 1
  
  #if there are any cool pixels (below tempdrop), run with it, else fill with NAs
  if(sum(diffframe,na.rm=TRUE)>1){
    Coords.for.cooled<-which(diffframe==1,arr.ind=TRUE)#Gets X/Y coordinates for cooled pixels
    cool.dist.matrix<-dist(Coords.for.cooled[,c(1,2)])#uses x/y coordinates to create distance matrix
    coolcluster <-fastcluster::hclust(cool.dist.matrix, method="single", members=NULL) #Clusters pee based on distance matrix
    coolgroups<-cutree(coolcluster,h=specifdist) #assigns groups based on 'specifdist' of clusters
    
   
    coolgroupinfo<-data.frame(table(coolgroups))
    
    #drop any cool clusters bigger than *biggestexclude* pixels OR is the biggest cluster
    coolgroupinfo$mouse<-ifelse((coolgroupinfo$Freq>biggestexclude | coolgroupinfo$Freq==max(coolgroupinfo$Freq)),1,0)
    
    micecoords<-as.integer(coolgroupinfo$coolgroups[which(coolgroupinfo$mouse==1)]) #what is the group ID of the mouse?
    
    
    Coords.for.cooled<-cbind(Coords.for.cooled,coolgroups) #assigns clusted ID as 'coolgroups' for all cooled pixels
    
    #creates new column (mouse), and if the group ID (coolgroup) of a given pixel matches the groupID for mice, then the value is 0 for the mouse column
    Coords.for.cooled<-data.frame(Coords.for.cooled)
    Coords.for.cooled$mouse<-match(Coords.for.cooled$coolgroups,micecoords)
    Coords.for.cooled$mouse<-ifelse(is.na(Coords.for.cooled$mouse),0,1)
    Coords.for.cooled$mouse<-ifelse(Coords.for.cooled$mouse==1,0,1)#mice pixels end up with '0' as value in 'mouse' column
    
    
    #create matrix with mouse cells having zero, all others equal NA
    mousezero<-Coords.for.cooled[,c(1,2,4)];colnames(mousezero)<-c("y","x","z")
    fullup<-merge(PC1,mousezero,by.x=c("row","col"),by.y=c("y","x"),all.x=TRUE)
    colnames(fullup)[c(1:2)]<-c("y","x")
    fullup$z[is.na(fullup$z)]<-(-1)#have to make NAs something else for converting to correct matrix structure below
    mousematrix<- xtabs(z~y+x, data=fullup)
    mousematrix2<-matrix(mousematrix,nrow=480,ncol=640)
    mousematrix2[mousematrix2==(-1)]<-NA
  
  } else {
    mousematrix2<-matrix(NA,nrow=480,ncol=640)
  }
  


  #######################
  #mousematrix2 = pixels of where the mouse was (i.e. and have gotten colder b/c the mouse is absent)
  # diffframe = pixels that have cooled 1 degree since PastFrame
  
  newcoldpixels<-mousematrix2*diffframe 
  #multiply together, which wipes out any newly cold pixels 
  #that are colder only because the mouse left because any mouse/hut pixel in mousematrix2 
  # is equal to 0 and therefore gets cancelled out by this multiplication
  
  newcoldpixels[newcoldpixels==0]<-NA #keeps only pixels with value 1, those that have cooled since PastFrame, and not b/c of the mouse moving
  
  
  #
  Coords.for.drip<-which(newcoldpixels==1,arr.ind=TRUE)#Gets X/Y coordinates for true, evaporatively cooled pixels

  return(newcoldpixels)
}


# Function to detect new cool clusters (evap, pee, dragging, trackking liquid, etc)
runitallcoolspotter<-function(submainframe,Lpolygoncoords,Rpolygoncoords){

  
  #Generates 307200 rows where the values in column 1 are possible y coordinates (1-480)
  #and the values in column 2 are possible x coordinates (1-640)
  possCombo<-expand.grid(seq(1,480,1),seq(1,640,1));colnames(possCombo)<-c("row","col")
  possCombo[,1]<-as.integer(possCombo[,1])
  possCombo[,2]<-as.integer(possCombo[,2])
  
  #turn into 2 matrices
  part.main.list<-list()
  for(fff in 1:2){
    part.main.list[[fff]]<-thermalarrange(submainframe[fff,],480,640)
  }
  rm(submainframe)
  gc()
  
  #clip matrix takes each frame in part.main.list, clips it to the correct dimensions, and subtracts the side-specific median temp
  LR.list<-clipmatrix(inputlist=part.main.list,Lpolygoncoords,Rpolygoncoords,allposscombos=possCombo)
  
  L.list<-LR.list[[1]]
  R.list<-LR.list[[2]]
  
  rm(LR.list)
  rm(part.main.list)
  gc()

  ##############################################################################################
  #########################################################
  
  Left.mega.pee.list<- newcooldetector(backwholeframe=L.list[[1]],
                                       thiswholeframe=L.list[[2]],
                                       biggestexclude=100,
                                        PC1=possCombo)

  
  Right.mega.pee.list<- newcooldetector(backwholeframe=R.list[[1]],
                                           thiswholeframe=R.list[[2]],
                                           biggestexclude=100,
                                           PC1=possCombo)

  
  rm(L.list)
  rm(R.list)
  ##############################################################################################
  ##############################################################################################
  #####################################################
  #####################################################
  #####################################################
  #####################################################
  
  Mega.pee.list<- LRCombine(L1=Left.mega.pee.list,R1=Right.mega.pee.list)
  
  rm(Left.mega.pee.list)
  rm(Right.mega.pee.list)
  ##############################################################################################
  #######################################################################################
  #############################################################################
  ###############################################################
  #####################################################
  
  return(Mega.pee.list)
  gc()
  
}