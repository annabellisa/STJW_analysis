# This function returns both richness and shannon's data frames in a list:

calc.div<-function(species.data, site.data){
  
  all.output<-list()
  
  rich.data<-list()
  shan.data<-list()
  
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
    
    # if there is only one species in a community (quadrat), shannon's diversity == 0, regardless of how many species are in functional group i. Thus, for functional groups with only one species, the value zero for all quadrats should be zero:
    if(length(vec.thisrun)==1) shan.data[[i]]<-rep(0, nrow(data.thisrun)) else shan.data[[i]]<-diversity(data.thisrun,index="shannon")
    
  } # close i for
  
  rich.res<-data.frame(do.call(cbind,rich.data))
  colnames(rich.res)<-group_df$group
  rich<-cbind(site.data,rich.res)
  
  shan.res<-data.frame(do.call(cbind,shan.data))
  colnames(shan.res)<-group_df$group
  shan<-cbind(site.data,shan.res)
  
  all.output$rich<-rich
  all.output$shan<-shan
  
  return(all.output)
  
} # close function