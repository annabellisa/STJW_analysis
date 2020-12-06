
# ------------------------------------ #
# ---------- STJW ANALYSIS  ---------- #
# ------------------------------------ #

### Analysis of chemical control experiment from STJW project
### Author: Annabel Smith & Raagini Muddaiah

# Load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

# Load libraries:
library(lme4); library(vegan); library(AICcmodavg); library(lmerTest)

# Load workspace
load("03_Workspaces/STJW_analysis.RData")

#  IMPORT & check data:    	# ----

data_dir<-"00_Data/Formatted_data/"
dir(data_dir)

# COUNT DATA:

ct17<-read.table(paste(data_dir,"count_2017.txt",sep=""),header=T)
ct18<-read.table(paste(data_dir,"count_2018.txt",sep=""),header=T)
ct19<-read.table(paste(data_dir,"count_2019.txt",sep=""),header=T)
rahead(ct17,3,7)
rahead(ct18,3,7)
rahead(ct19,3,7)

# PLANT INFO:
pinfo<-read.table(paste(data_dir,"plant_info.txt",sep="/"),header=T)
pinfo$Duration<-as.factor(pinfo$Duration)
pinfo$func_grp<-as.factor(pinfo$func_grp)
pinfo$Significance<-as.factor(pinfo$Significance)
pinfo$Status<-as.factor(pinfo$Status)
head(pinfo,3); dim(pinfo)
length(which(duplicated(pinfo$Sp))) # should be zero

# CHECK:

cnames17<-colnames(ct17)[which(colnames(ct17)=="Aca_ovi"):ncol(ct17)]
cnames18<-colnames(ct18)[which(colnames(ct18)=="Aca_ovi"):ncol(ct18)]
cnames19<-colnames(ct19)[which(colnames(ct19)=="Aca_ovi"):ncol(ct19)]

# these should all be TRUE:

# Are all names matching across years?
length(cnames17)==length(cnames18)
length(cnames17)==length(cnames19)
length(cnames18)==length(cnames19)

table(cnames17 %in% cnames18)
table(cnames17 %in% cnames19)
table(cnames18 %in% cnames19)

# Are all names in the data present in plant_info?
table(cnames17 %in% pinfo$Sp)
table(cnames18 %in% pinfo$Sp)
table(cnames19 %in% pinfo$Sp)

# close import data ----

#  FORMAT data:    	# ----

# replace date with year:
ct17$DATE<-2017
ct18$DATE<-2018
ct19$DATE<-2019

# QUADRAT_Direction is the treatment (A, B and C):
colnames(ct17)[which(colnames(ct17)=="QUADRAT_Direction")]<-"Treatment"
colnames(ct18)[which(colnames(ct18)=="QUADRAT_Direction")]<-"Treatment"
colnames(ct19)[which(colnames(ct19)=="QUADRAT_Direction")]<-"Treatment"

# create reserve variable from PLOT_ID:
ct17$reserve<-ct17$PLOT_ID
ct17$reserve[grep("J", ct17$reserve)]<-"J"
ct17$reserve[grep("M", ct17$reserve)]<-"M"

ct18$reserve<-ct18$PLOT_ID
ct18$reserve[grep("J", ct18$reserve)]<-"J"
ct18$reserve[grep("M", ct18$reserve)]<-"M"

ct19$reserve<-ct19$PLOT_ID
ct19$reserve[grep("J", ct19$reserve)]<-"J"
ct19$reserve[grep("M", ct19$reserve)]<-"M"

# re-arrange columns:
ct17<-ct17[,c(which(colnames(ct17) %in% c("DATE", "reserve")),which(!colnames(ct17) %in% c("DATE", "reserve")))]
ct18<-ct18[,c(which(colnames(ct18) %in% c("DATE", "reserve")),which(!colnames(ct18) %in% c("DATE", "reserve")))]
ct19<-ct19[,c(which(colnames(ct19) %in% c("DATE", "reserve")),which(!colnames(ct19) %in% c("DATE", "reserve")))]

# remove unwanted columns:
ct17$QUAD_ID<-NULL
ct18$QUAD_ID<-NULL
ct19$QUAD_ID<-NULL

# factorise variables which should be factors and make control the baseline level:
ct17$reserve<-as.factor(ct17$reserve)
ct17$PLOT_ID<-as.factor(ct17$PLOT_ID)
ct17$Treatment<-factor(ct17$Treatment,levels=c("C","A","B"))

ct18$reserve<-as.factor(ct18$reserve)
ct18$PLOT_ID<-as.factor(ct18$PLOT_ID)
ct18$Treatment<-factor(ct18$Treatment,levels=c("C","A","B"))

ct19$reserve<-as.factor(ct19$reserve)
ct19$PLOT_ID<-as.factor(ct19$PLOT_ID)
ct19$Treatment<-factor(ct19$Treatment,levels=c("C","A","B"))

# replace count categories with numbers:

# 16-50 	W	35
# 51-100	X	75
# >100	  Y	100

# Don't worry about NAs introduced by coercion, this'll be fixed when we sort out the NAs for Tri_ela and Ant_sca

ct17d<-as.matrix(ct17[,which(colnames(ct17)=="Aca_ovi"):ncol(ct17)])
ct17d[which(ct17d=="W")]<-35
ct17d[which(ct17d=="X")]<-75
ct17d[which(ct17d=="Y")]<-100
ct17d<-data.frame(apply(ct17d,2,as.numeric))
ct17<-data.frame(cbind(ct17[,1:which(colnames(ct17)=="Treatment")],ct17d))

ct18d<-as.matrix(ct18[,which(colnames(ct18)=="Aca_ovi"):ncol(ct18)])
ct18d[which(ct18d=="W")]<-35
ct18d[which(ct18d=="X")]<-75
ct18d[which(ct18d=="Y")]<-100
ct18d<-data.frame(apply(ct18d,2,as.numeric))
ct18<-data.frame(cbind(ct18[,1:which(colnames(ct18)=="Treatment")],ct18d))

ct19d<-as.matrix(ct19[,which(colnames(ct19)=="Aca_ovi"):ncol(ct19)])
ct19d[which(ct19d=="W")]<-35
ct19d[which(ct19d=="X")]<-75
ct19d[which(ct19d=="Y")]<-100
ct19d<-data.frame(apply(ct19d,2,as.numeric))
ct19<-data.frame(cbind(ct19[,1:which(colnames(ct19)=="Treatment")],ct19d))

# remove unidentified species (for now... but we need to discuss whether these should be included in any analyses):
unid<-pinfo$Sp[grep("Uni_", pinfo$Sp)]
pinfo<-pinfo[-which(pinfo$Sp %in% unid),]
pinfo<-tidy.df(pinfo)
head(pinfo, 3); dim(pinfo)

ct17<-ct17[,-which(colnames(ct17) %in% unid)]
ct18<-ct18[,-which(colnames(ct18) %in% unid)]
ct19<-ct19[,-which(colnames(ct19) %in% unid)]

rahead(ct17,3,7); dim(ct17)
rahead(ct18,3,7); dim(ct18)
rahead(ct19,3,7); dim(ct19)

# close format data ----

#  species DIVERSITY GROUPS:    	# ----

head(pinfo,3); dim(pinfo)
str(pinfo)

# Set-up vectors for grouping:

# All species:
all<-pinfo$Sp
native<-as.character(pinfo$Sp[which(pinfo$Status=="N")])
exotic<-as.character(pinfo$Sp[which(pinfo$Status=="E")])
(length(native)+length(exotic))==length(all) # should be TRUE

# Indicator all (A+B):
indic<-as.character(pinfo$Sp[which(pinfo$Indicator==1)])

# Significance:
sigA<-as.character(pinfo$Sp[which(pinfo$Significance=="A")])
sigB<-as.character(pinfo$Sp[which(pinfo$Significance=="B")])

# The length of indic should be equal to the length of significance A and B categories:
length(pinfo$Sp[pinfo$Significance=="A"])+length(pinfo$Sp[pinfo$Significance=="B"])==length(indic)

# For significance, X only has one species (Nas_tri). Combine with Y? This would make sigXY the "bad" weeds and sigZ the "not so bad" weeds:
table(pinfo$Significance)

sigC<-as.character(pinfo$Sp[which(pinfo$Significance=="C")])
sigXY<-as.character(pinfo$Sp[c(which(pinfo$Significance=="X"),which(pinfo$Significance=="Y"))])
sigZ<-as.character(pinfo$Sp[which(pinfo$Significance=="Z")])

# Herbs, including leguminous and non-leguminous herbs:
native_herb<-as.character(pinfo$Sp[which(pinfo$Herb==1 & pinfo$Status=="N")])
exotic_herb<-as.character(pinfo$Sp[which(pinfo$Herb==1 & pinfo$Status=="E")])

table(native_herb %in% exotic_herb) # should be FALSE
table(exotic_herb %in% native_herb) # should be FALSE

# Herb catgories:
exann_herb<-as.character(pinfo$Sp[which(pinfo$Herb==1 & pinfo$Status=="E" & pinfo$Duration=="Annual")])
exper_herb<-as.character(pinfo$Sp[which(pinfo$Herb==1 & pinfo$Status=="E" & pinfo$Duration=="Perennial")])
natann_herb<-as.character(pinfo$Sp[which(pinfo$Herb==1 & pinfo$Status=="N" & pinfo$Duration=="Annual")])
natper_herb<-as.character(pinfo$Sp[which(pinfo$Herb==1 & pinfo$Status=="N" & pinfo$Duration=="Perennial")])

# Non-leguminous Herbs:
native_nonlegherb<-as.character(pinfo$Sp[which(pinfo$Legume==0 & pinfo$Herb==1 & pinfo$Status=="N")])
exotic_nonlegherb<-as.character(pinfo$Sp[which(pinfo$Legume==0 & pinfo$Herb==1 & pinfo$Status=="E")])

# Leguminous Herbs:
native_legherb<-as.character(pinfo$Sp[which(pinfo$Legume==1 & pinfo$Herb==1 & pinfo$Status=="N")])
exotic_legherb<-as.character(pinfo$Sp[which(pinfo$Legume==1 & pinfo$Herb==1 & pinfo$Status=="E")])

head(pinfo,3); dim(pinfo)

# All grasses:
native_grass<-as.character(pinfo$Sp[which(pinfo$Grass==1 & pinfo$Status=="N")])
exotic_grass<-as.character(pinfo$Sp[which(pinfo$Grass==1 & pinfo$Status=="E")])

# Annual and perennial exotic Grasses 
# All native Grasses are perennial, so no need to split those
exotic_anngrass<-as.character(pinfo$Sp[which(pinfo$Duration=="Annual" & pinfo$Grass==1 & pinfo$Status=="E")])
exotic_perengrass<-as.character(pinfo$Sp[which(pinfo$Grass==1 & pinfo$Status=="E" & pinfo$Duration=="Perennial")])

# C3 and C4 Grasses:
# There are no exotic c4 grasses in this data set, so we won't look at that category
# All c4 grasses are native
c3_grass<-as.character(pinfo$Sp[which(pinfo$func_grp=="C3")])

native_c3<-as.character(pinfo$Sp[which(pinfo$func_grp=="C3" & pinfo$Status=="N")])
native_c4<-as.character(pinfo$Sp[which(pinfo$func_grp=="C4" & pinfo$Status=="N")])
exotic_c3<-as.character(pinfo$Sp[which(pinfo$func_grp=="C3" & pinfo$Status=="E")])

# Sedges and rushes:
sed_rus<-as.character(pinfo$Sp[which(pinfo$func_grp=="Sedge_Rush")])

# Categories to analyse
group_df<-data.frame(group=c("all","native","exotic","indic","sigA","sigB","sigC","sigXY","sigZ","native_herb","exotic_herb","exann_herb","exper_herb","natann_herb","natper_herb","native_nonlegherb","exotic_nonlegherb","native_legherb","exotic_legherb","native_grass","exotic_grass","exotic_anngrass","exotic_perengrass","c3_grass","native_c3","native_c4","exotic_c3","sed_rus"))

# close diversity groups ----

#  calculate species DIVERSITY:    	# ----

head(pinfo,3); dim(pinfo)
head(group_df)

rahead(ct17,3,7); dim(ct17)
rahead(ct18,3,7); dim(ct18)
rahead(ct19,3,7); dim(ct19)

# Add species richness and Shannon's index to count data:

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
    } 
    
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

# CALCULATE diversity for each year separately:

ct17_site.data<-ct17[,1:which(colnames(ct17)=="Treatment")]
div17<-calc.div(ct17, ct17_site.data)
rich17<-div17$rich
shan17<-div17$shan

ct18_site.data<-ct18[,1:which(colnames(ct18)=="Treatment")]
div18<-calc.div(ct18, ct18_site.data)
rich18<-div18$rich
shan18<-div18$shan

ct19_site.data<-ct19[,1:which(colnames(ct19)=="Treatment")]
div19<-calc.div(ct19, ct19_site.data)
rich19<-div19$rich
shan19<-div19$shan

# Then combine:

head(rich17[,1:10],3); dim(rich17)
head(rich18[,1:10],3); dim(rich18)
head(rich19[,1:10],3); dim(rich19)

rich<-rbind(rich17,rich18,rich19)
rahead(rich,4,7); dim(rich)

shan<-rbind(shan17, shan18, shan19)
rahead(shan,4,7); dim(shan)

# We need to decide the cut-off for analysis. We could say there must be more records (i.e. number of species recorded) than the number of quadrats (48)? Or the number of observations (within years?). Number of quadrats would cut out two responses (exotic_perengrass and sed_rus)
group_df$rich_records<-colSums(rich[,5:length(rich)])
gdf<-group_df
gdf<-tidy.df(gdf)
gdf

gdf$ylab<-c("All","Native","Exotic","Indicator","Significance A","Significance B","Common/Increaser","Significance X/Y","Significance Z","Native forb", "Exotic forb","Exotic annual forb","Exotic perennial forb","Native annual forb","Native perennial forb","Native non-leg. forb","Exotic non-leg. forb","Native leg. forb","Exotic leg. forb","Native grass","Exotic grass","Exotic annual grass","Exotic perennial grass","C3 grass","Native C3 grass","Native C4 grass","Exotic C3 grass","Sedge/Rush")

save.image("03_Workspaces/stjw_analysis.RData")

# close diversity calculation ----

#  VISUALISE raw data:    	# ----

head(gdf); dim(gdf)

rahead(rich,4,7); dim(rich)
rahead(shan,4,7); dim(shan)

# visualise raw data:

# year effects on species richness:

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(4,4,1,1), mgp=c(2.5,1,0))

for(i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  resp.dat<-rich[,which(colnames(rich)==resp.thisrun)]
  ylab.thisrun<-gdf$ylab[i]
  
  boxplot(resp.dat~rich$reserve+as.factor(rich$DATE), cex.axis=1, col=c("darkturquoise","darkolivegreen2"), ylab=ylab.thisrun, xlab="", las=1, xaxt="n")
  axis(side=1, at=c(1.5,3.5,5.5), labels=c(2017, 2018, 2019), cex.axis=0.8)
  arrows(c(2.5,4.5),0,c(2.5,4.5),40, length=0, col="grey70")
}
plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
legend(1,9,legend=c("Jerrabomberra","Mulangari"), col=c("darkturquoise","darkolivegreen2"), pch=15, bty="n", pt.cex = 3)

# treatement effects on species richness:

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(4,4,1,1), mgp=c(2.5,1,0))

for(i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  resp.dat<-rich[,which(colnames(rich)==resp.thisrun)]
  ylab.thisrun<-gdf$ylab[i]
  
  boxplot(resp.dat~rich$reserve+as.factor(rich$Treatment), cex.axis=1, col=c("darkturquoise","darkolivegreen2"), ylab=ylab.thisrun, xlab="", las=1, xaxt="n")
  axis(side=1, at=c(1.5,3.5,5.5), labels=c("C", "A", "B"), cex.axis=0.8)
  arrows(c(2.5,4.5),0,c(2.5,4.5),40, length=0, col="grey70")
}
plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
legend(1,9,legend=c("Jerrabomberra","Mulangari"), col=c("darkturquoise","darkolivegreen2"), pch=15, bty="n", pt.cex = 3)

# close visualise ----

#  ANALYSIS:    	# ----

head(gdf); dim(gdf)

rahead(rich,4,7); dim(rich)
rahead(shan,4,7); dim(shan)

# scale date only:
rich_sc<-rich
rich_sc$DATE<-rich_sc$DATE-2017
rahead(rich_sc,4,7); dim(rich_sc)

shan_sc<-shan
shan_sc$DATE<-shan_sc$DATE-2017
rahead(shan_sc,4,7); dim(shan_sc)

# new data for model estimates:
nd1<-data.frame(DATE=rep(c(0,1,2),rep(3,3)),Treatment=as.factor(c("C","A","B")))

# lists for storing coefficients, anova tables and model estimates:
coef.out<-list()
anova.out<-list()
preds.out<-list()

# run models and generate model estimates:
# update the data set (rich_sc or shan_sc) and re-run for each

data.set<-shan_sc

for (i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  data.thisrun<-data.set[,resp.thisrun]
  form.thisrun<-paste(resp.thisrun,"~Treatment*DATE+(1|reserve/PLOT_ID)", sep="")
  
  if(sum(data.thisrun,na.rm=T)==0){
  coef.out[[i]]<-NULL
  anova.out[[i]]<-NULL
  preds.out[[i]]<-NULL
  next
  }
  
  m1<-lmer(formula = form.thisrun, data=data.set)
  summary(m1)

  m1_coef<-coef.ext(m1)
  m1_anova<-anova.ext(m1)
  coef.out[[i]]<-m1_coef
  anova.out[[i]]<-m1_anova
  
  p1<-pred(model=m1, new.data = nd1,se.fit = T, type = "response")
  preds.out[[i]]<-p1

} # close models

## PLOT:

# Species Richness:

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(2,4,4,1), mgp=c(2.5,1,0))

for(i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  pred.thisrun<-preds.out[[i]]
  anova.thisrun<-anova.out[[i]]
  ylab.thisrun<-gdf$ylab[i]
  xofs<-0.2
  arrowlgth<-0.02
  
  plot(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$fit[pred.thisrun$Treatment=="C"], pch=15, ylim=c(min(pred.thisrun$lci), max(pred.thisrun$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
  axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
  arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$lci[pred.thisrun$Treatment=="C"],pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$uci[pred.thisrun$Treatment=="C"], code=3, angle=90, length=arrowlgth)
  points(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$fit[pred.thisrun$Treatment=="A"], pch=15, col="red")
  arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$lci[pred.thisrun$Treatment=="A"],pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$uci[pred.thisrun$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
  points(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$fit[pred.thisrun$Treatment=="B"], pch=15, col="blue")
  arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$lci[pred.thisrun$Treatment=="B"],pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$uci[pred.thisrun$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
  
  p.trt<-round(anova.thisrun$p[1],4)
  p.yr<-round(anova.thisrun$p[2],4)
  p.int<-round(anova.thisrun$p[3],4)
  
  if(p.trt>0.01) p.trt<-round(p.trt,2) else p.trt<-"<0.01"
  if(p.yr>0.01) p.yr<-round(p.yr,2) else p.yr<-"<0.01"
  if(p.int>0.01) p.int<-round(p.int,2) else p.int<-"<0.01"
  
  title(main=paste("P values: trt=",p.trt,"\nyr=",p.yr,"; int=",p.int), font.main=1, adj=0, cex=0.9, line=0.5)
  
}

plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
legend(1,9,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)

# Shannon's Diversity:

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(2,4,4,1), mgp=c(2.5,1,0))

for(i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  pred.thisrun<-preds.out[[i]]
  anova.thisrun<-anova.out[[i]]
  ylab.thisrun<-gdf$ylab[i]
  xofs<-0.2
  arrowlgth<-0.02
  
  if (is.null(pred.thisrun)){
    plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
    next
  }
  
  plot(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$fit[pred.thisrun$Treatment=="C"], pch=15, ylim=c(min(pred.thisrun$lci), max(pred.thisrun$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
  axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
  arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$lci[pred.thisrun$Treatment=="C"],pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$uci[pred.thisrun$Treatment=="C"], code=3, angle=90, length=arrowlgth)
  points(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$fit[pred.thisrun$Treatment=="A"], pch=15, col="red")
  arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$lci[pred.thisrun$Treatment=="A"],pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$uci[pred.thisrun$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
  points(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$fit[pred.thisrun$Treatment=="B"], pch=15, col="blue")
  arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$lci[pred.thisrun$Treatment=="B"],pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$uci[pred.thisrun$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
  
  p.trt<-round(anova.thisrun$p[1],4)
  p.yr<-round(anova.thisrun$p[2],4)
  p.int<-round(anova.thisrun$p[3],4)
  
  if(p.trt>0.01) p.trt<-round(p.trt,2) else p.trt<-"<0.01"
  if(p.yr>0.01) p.yr<-round(p.yr,2) else p.yr<-"<0.01"
  if(p.int>0.01) p.int<-round(p.int,2) else p.int<-"<0.01"
  
  title(main=paste("P values: trt=",p.trt,"\nyr=",p.yr,"; int=",p.int), font.main=1, adj=0, cex=0.9, line=0.5)
  
}

plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
legend(1,9,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)



# close analysis ----



