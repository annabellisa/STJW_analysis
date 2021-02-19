# This function returns both richness and shannon's data frames in a list:

calc.div<-function(species.data, site.data){
  
  all.output<-list()
  
  rich.data<-list()
  shan.data<-list()
  simp.data<-list()
  invsimp.data<-list()
  
  for (i in 1:nrow(group_df)){
    
    name.thisrun<-as.character(group_df$group[i])
    vec.thisrun<-get(name.thisrun)
    
    data.thisrun<-species.data[,colnames(species.data) %in% vec.thisrun]
    head(data.thisrun,3); dim(data.thisrun)
    
    if(length(vec.thisrun)==1){
      data.thisrun<-data.frame(data.thisrun)
      colnames(data.thisrun)<-vec.thisrun
    } # close if
    
    rich.data[[i]]<-apply(data.thisrun,1,function(x)length(which(x>0)))
    
    # if there is only one species in a community (quadrat), shannon's diversity == 0, regardless of how many species are in functional group i. Thus, for functional groups with only one species, the value for all quadrats should be zero:
    if(length(vec.thisrun)==1) shan.data[[i]]<-rep(0, nrow(data.thisrun)) else shan.data[[i]]<-diversity(data.thisrun,index="shannon")
    
    if(length(vec.thisrun)==1) simp.data[[i]]<-rep(0, nrow(data.thisrun)) else simp.data[[i]]<-diversity(data.thisrun,index="simpson")
    
    if(length(vec.thisrun)==1) invsimp.data[[i]]<-rep(0, nrow(data.thisrun)) else invsimp.data[[i]]<-diversity(data.thisrun,index="invsimpson")
    
  } # close i for
  
  rich.res<-data.frame(do.call(cbind,rich.data))
  colnames(rich.res)<-group_df$group
  rich<-cbind(site.data,rich.res)
  
  shan.res<-data.frame(do.call(cbind,shan.data))
  colnames(shan.res)<-group_df$group
  shan<-cbind(site.data,shan.res)
  
  simp.res<-data.frame(do.call(cbind,simp.data))
  colnames(simp.res)<-group_df$group
  simp<-cbind(site.data,simp.res)
  
  invsimp.res<-data.frame(do.call(cbind,invsimp.data))
  colnames(invsimp.res)<-group_df$group
  invsimp<-cbind(site.data,invsimp.res)
  
  all.output$rich<-rich
  all.output$shan<-shan
  all.output$simp<-simp
  all.output$invsimp<-invsimp
  
  return(all.output)
  
} # close function