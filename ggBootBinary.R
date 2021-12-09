ggBootBinary<-function(data, boot=100, method="euclidean",mmd.size.cor=FALSE){
  
  #Required R packages
  require(MASS)
  require(ks)
  require(tidyverse)
  require(reshape2) 
  require(rgdal)
  require(rgeos)
  
  
  # 1.Basic info on the variables
  data<- data.frame(data)
  nvars<-ncol(data)-1
  data[,1]<-as.factor(data[,1])
  ngroups<-nlevels(data[,1])
  groupnames<-levels(data[,1])

  # 2. Frequency matrix
  ## 2.1 Matrix
  frequencies<-matrix(NA,ngroups*boot+ngroups,nvars)
  colnames(frequencies)<-colnames(data[,2:ncol(data)])

  ## 2.2 Assign row names (original+bootstrap)
  tmp.rownames<-rep(0,ngroups*boot+ngroups)
  tmp.rownames[1:ngroups]<-groupnames
  
  for(a in 1:ngroups){
    for(b in 1:boot){
      tmp.rownames[(a-1)*boot+ngroups+b]<-paste(groupnames[a],"Boot",b,sep="")
    }
  }
  row.names(frequencies)<-tmp.rownames
  
  ## 2.3 Calculate frequencies for the original data
  for(a in 1:ngroups){
    tmp.index<-which(data[,1]==groupnames[a])
    
    frequencies[a,]<-apply(data[tmp.index,2:ncol(data)],2,mean,na.rm=TRUE)
    
  }
  
  ## 2.4 Calculate the frequencies for each bootstrap
  ### First loop pools data by group
  for(a in 1:ngroups){
    group.index<-which(data[,1]==groupnames[a])
    tmp.data<-data[group.index,2:ncol(data)]
    
    ### Second loop removes NAS and counts the number of trait observations
    for(b in 1:nvars){
      tmp.obs.present<-which(is.na(tmp.data[,b])==FALSE)
      nobs<-length(tmp.obs.present)
      
      ### Third loop resamples data and calculates the mean
      for(c in 1:boot){
        boot.index<-sample.int(nobs,nobs,replace = TRUE)
        tmp.freq<-mean(tmp.data[tmp.obs.present[boot.index],b])
        frequencies[(a-1)*boot+ngroups+c,b]<-tmp.freq
        
      }
    }
  }
  
  # 3. Distance matrix
  ## 3.1 Euclidean Distance
  if(method== "euclidean"){
    distances<-dist(frequencies)
    distances<-as.matrix(distances)
    distances[which(distances==0)]<-0.00001
  }
  
  ## 3.2 Mean Measure of Divergence
  if(method=="MMD"){
    
    ### Sample sizes (original data)
    sample.sizes <-matrix(NA,ngroups*boot+ngroups,nvars)

    for(a in 1:ngroups){
      for(b in 1:nvars){
      group.index<-which(data[,1]==groupnames[a])
      tmp.data<-data[group.index,2:ncol(data)]
      
      obs.present<-which(is.na(tmp.data[,b])==FALSE)
      sample.sizes[a,b]<-length(obs.present)
      
      }
    }
    
    ### Sample sizes (bootstrap data)
    for(a in 1:ngroups){
      group.index<-which(data[,1]==groupnames[a])
      tmp.data<-data[group.index,2:ncol(data)]
      
      for(b in 1:nvars){
        tmp.obs.present<-which(is.na(tmp.data[,b])==FALSE)
        nobs<-length(tmp.obs.present)
        
        for(c in 1:boot){
          sample.sizes[(a-1)*boot+ngroups+c,b]<-nobs
        }
      }
    }
    
    ### Correct the observed frequencies according to Sjovold 1973. 
    frequencies[frequencies==1]<-1-(1/(4*sample.sizes[frequencies==1]))
    frequencies[frequencies==0]<-1/(4*sample.sizes[frequencies==0])
    
    ### MMD matrix
    distance.n <- ngroups*boot+ngroups
    distances<-matrix(0,distance.n,distance.n)
    colnames(distances)<-rownames(frequencies)
    rownames(distances)<-rownames(frequencies)
    
    ### Calculate MMD
    for (a in 1:distance.n){
      for (b in 1:distance.n){
        if (mmd.size.cor==TRUE){
          distances[a,b]<-sum((asin(1-2*frequencies[a,])-asin(1-2*frequencies[b,]))^2-((1/(sample.sizes[a,]))+(1/(sample.sizes[b,]))))/nvars
        } else { 
          distances[a,b]<-sum((asin(1-2*frequencies[a,])-asin(1-2*frequencies[b,]))^2)/nvars
        }
      }
    }
    
    #if values of MMD are negative (due to the size correction effect), replace them with 0
    distances[distances<0]<-0
  }

  # 4. Kruskal's Non-metric Multidimensional Scaling
  distances[which(distances==0)]<-0.00001
  MDSResultsOriginal<-isoMDS(as.dist(distances[1:ngroups,1:ngroups]))
  MDSResultsBoot<-isoMDS(as.dist(distances))
  
  # 5. Scatterplot
  ## 5.1 Vector with the categories that will divide the groups.
  tmp.rownames<-rep(0,ngroups*boot+ngroups)
  tmp.rownames[1:ngroups]<-groupnames
  for(a in 1:ngroups){
    for(b in 1:boot){
      tmp.rownames[(a-1)*boot+ngroups+b]<-paste(groupnames[a])
    }
  }
  
  ## 5.2 Get the MDS points.
  ggboot<-data.frame(MDSResultsBoot[["points"]])
  
  ## 5.3 Add the categories.
  ggboot[,3]<- tmp.rownames
  colnames(ggboot)<- c("x","y","name")
  
  ## 5.4 Separate the original distances from the bootstrap distances. 
  ### Original data.
  ggoriginal.index<- which(rownames(ggboot)==groupnames)
  ggoriginal<- ggboot[ggoriginal.index,]
  
  ### Bootstrapped distances.
  ggboot<-ggboot[-c(ggoriginal.index),]
  
  ## 5.5 Kernel densities and confidence intervals.
  ### Empty lists where we will store the Kernel estimations of each group.
  confidence.interval95<- list()
  confidence.interval85<- list()
  confidence.interval75<- list()
  confidence.interval50<- list()
  
  ### Calculate and store Kernel estimations for each confidence interval.
  for (a in 1:(ngroups)) {
    
    #### Subsets the data frames of each group.
    group.index<-which(ggboot[,3]==groupnames[a])
    df <- ggboot[group.index,c(1,2)] 
    kd <- kde(df, compute.cont=TRUE)
    
    #### Kernel estimation for a 95% confidence interval
    contour_95 <- with(kd, contourLines(x=eval.points[[1]],y=eval.points[[2]],z=estimate, 
                                        levels=cont["5%"])[[1]])
    contour_95 <- data.frame(contour_95) 
    contour_95$name<-rep(groupnames[a],nrow(contour_95))#Add a column with group name
    confidence.interval95[[a]]<-contour_95 #store it in the appropriate list
    
    #### Kernel estimation for a 85% confidence interval
    contour_85 <- with(kd, contourLines(x=eval.points[[1]],y=eval.points[[2]],z=estimate, 
                                        levels=cont["15%"])[[1]])
    contour_85 <- data.frame(contour_85)
    contour_85$name<-rep(groupnames[a],nrow(contour_85))
    confidence.interval85[[a]]<-contour_85
    
    #### Kernel estimation for a 75% confidence interval
    contour_75 <- with(kd, contourLines(x=eval.points[[1]],y=eval.points[[2]],z=estimate, 
                                        levels=cont["25%"])[[1]])
    contour_75 <- data.frame(contour_75)
    contour_75$name<-rep(groupnames[a],nrow(contour_75))
    confidence.interval75[[a]]<-contour_75
    
    #### Kernel estimation for a 50% confidence interval
    contour_50 <- with(kd, contourLines(x=eval.points[[1]],y=eval.points[[2]],z=estimate, 
                                        levels=cont["50%"])[[1]])
    contour_50 <- data.frame(contour_50)
    contour_50$name<-rep(groupnames[a],nrow(contour_50))
    confidence.interval50[[a]]<-contour_50
  }
  
  ### Merge the same confidence intervals of different groups into a data frame.
  conf.int95 = do.call(rbind, confidence.interval95)
  conf.int85 = do.call(rbind, confidence.interval85)
  conf.int75 = do.call(rbind, confidence.interval75)
  conf.int50 = do.call(rbind, confidence.interval50)
  
  ## 5.6 Estimate overlap between 95 confidence intervals
  overlap.matrix<- matrix(nrow = ngroups, 
                          ncol = ngroups,
                          dimnames = list(groupnames,groupnames))
  for (a in 1:ngroups-1) {
    for (b in 1:(ngroups)) {
      if(groupnames[a+1]==groupnames[b]){
        overlap.matrix[b,a+1]<- 1
      }else{
        group1 <- conf.int95[which(conf.int95$name==groupnames[a+1]),c(2,3)]
        group2 <- conf.int95[which(conf.int95$name==groupnames[b]),c(2,3)]
        
        group1_pol <- Polygons(list(Polygon(group1)), groupnames[a+1])
        group2_pol <- Polygons(list(Polygon(group2)), groupnames[b])
        shape <- SpatialPolygons(list(group1_pol, group2_pol))
        
        group1_area <- gArea(shape[groupnames[a+1]])
        intersections <- gIntersection(shape[groupnames[a+1]],
                                       shape[groupnames[b]])
        if(is.null(intersections)){
          overlap.matrix[b,a+1]<- 0
        }else{
          overlap.matrix[b,a+1]<-gArea(intersections)/group1_area
        }
      }
    }
  }

  p1<-
    ggplot()+
    geom_polygon(data=conf.int95,aes(x, y,fill=name),alpha=0.1)+
    geom_polygon(data=conf.int85,aes(x, y,fill=name),alpha=0.1)+
    geom_polygon(data=conf.int75,aes(x, y,fill=name),alpha=0.1)+
    geom_polygon(data=conf.int50,aes(x, y,fill=name),alpha=0.1)+
    geom_point(data=ggboot,aes(x=x,y=y, fill=name,color=name),alpha=0.15)+
    geom_text(data=ggoriginal,aes(x=x,y=y,label=name),vjust=-0.8,size=5)+
    geom_point(data=ggoriginal,aes(x=x,y=y,color=name,size=5))+
    labs(caption = paste("MDS Stress:", round(MDSResultsBoot[["stress"]],2),"%"),x="",y="")+
    theme_classic()+
    ylim(min(ggboot$y)+min(ggboot$y)/8,max(ggboot$y)+max(ggboot$y)/8)+
    xlim(min(ggboot$x)+min(ggboot$x)/8,max(ggboot$x)+max(ggboot$x)/8)+
    guides(size="none",color="none",alpha="none",fill="none")

  
  # 6. Heatmap - original distance matrix plot
  matrix.plot<- distances[1:ngroups,1:ngroups]
  matrix.plot[upper.tri(matrix.plot)]<-NA
  
  for (a in 1:ncol(matrix.plot)) {
    matrix.plot[a,a]<-NA
  }
  
  matrix.plot<- melt(matrix.plot,na.rm = TRUE)
  
  p2 <- ggplot(matrix.plot,aes(x=Var1,y=Var2,fill=value))+
    theme_classic()+
    geom_tile(color="white")+
    geom_text(aes(Var1,Var2,label= round(value,digits = 2)),color="black",size=4)+
    scale_fill_viridis_c(alpha = 0.8)+
    labs(fill=paste(method,"\n Distance"),x="",y="")+
    coord_fixed()+
    theme(axis.text.x = element_text(angle = 90))
  
  # 7. List of results
  results<-list(
    ## The distance matrix of the original data.
    OriginalDistances=as.matrix(distances)[1:ngroups,1:ngroups],
    #The distance matrix with all bootstraps.
    BootDistances=as.matrix(distances),
    #Frequency matrix showing how much each polygon overlaps
    polygon.overlap.matrix=overlap.matrix,
    #Scatterplot
    MDS.Plot=p1,
    #Heatmap
    Original.matrix=p2)
  
  return(results)
}
