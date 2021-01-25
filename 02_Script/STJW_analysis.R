
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

#COVER DATA:
cv17<-read.table(paste(data_dir,"cover_2017.txt",sep=""),header=T)
cv18<-read.table(paste(data_dir,"cover_2018.txt",sep=""),header=T)
cv19<-read.table(paste(data_dir,"cover_2019.txt",sep=""),header=T)
rahead(cv17,3,7)
rahead(cv18,3,7)
rahead(cv19,3,7)


# PLANT INFO:

pinfo<-read.table(paste(data_dir,"plant_info.txt",sep="/"),header=T)
pinfo$Duration<-as.factor(pinfo$Duration)
pinfo$func_grp<-as.factor(pinfo$func_grp)
pinfo$Significance<-as.factor(pinfo$Significance)
pinfo$Status<-as.factor(pinfo$Status)
head(pinfo,3); dim(pinfo)
length(which(duplicated(pinfo$Sp))) # should be zero

# CHECK COUNT:

cnames17<-colnames(ct17)[which(colnames(ct17)=="Aca_ovi"):ncol(ct17)]
cnames18<-colnames(ct18)[which(colnames(ct18)=="Aca_ovi"):ncol(ct18)]
cnames19<-colnames(ct19)[which(colnames(ct19)=="Aca_ovi"):ncol(ct19)]

#CHECK COVER:

cvnames17<-colnames(cv17)[which(colnames(cv17)=="Aca_ovi"):ncol(cv17)]
cvnames18<-colnames(cv18)[which(colnames(cv17)=="Aca_ovi"):ncol(cv18)]
cvnames19<-colnames(cv19)[which(colnames(cv19)=="Aca_ovi"):ncol(cv19)]

# these should all be TRUE:

# Are all names matching across yearS? 
#COUNT:
length(cnames17)==length(cnames18)
length(cnames17)==length(cnames19)
length(cnames18)==length(cnames19)

table(cnames17 %in% cnames18)
table(cnames17 %in% cnames19)
table(cnames18 %in% cnames19)

#COVER:
length(cvnames17)==length(cvnames18)
length(cvnames17)==length(cvnames19)
length(cvnames18)==length(cvnames19)


table(cvnames17 %in% cvnames18)
table(cvnames17 %in% cvnames19)
table(cvnames18 %in% cvnames19)

#Are all names in the data present in plant_info?
#COUNT:
table(cnames17 %in% pinfo$Sp)
table(cnames18 %in% pinfo$Sp)
table(cnames19 %in% pinfo$Sp)

#COVER:
table(cvnames17 %in% pinfo$Sp)
table(cvnames18 %in% pinfo$Sp)
table(cvnames19 %in% pinfo$Sp)


# close import data ----

#  FORMAT data:    	# ----

#replace date with year:
#COUNT:
ct17$DATE<-2017
ct18$DATE<-2018
ct19$DATE<-2019

#COVER:
cv17$DATE<-2017
cv18$DATE<-2018
cv19$DATE<-2019

# QUADRAT_Direction is the treatment (A, B and C):
#COUNT:
colnames(ct17)[which(colnames(ct17)=="QUADRAT_Direction")]<-"Treatment"
colnames(ct18)[which(colnames(ct18)=="QUADRAT_Direction")]<-"Treatment"
colnames(ct19)[which(colnames(ct19)=="QUADRAT_Direction")]<-"Treatment"

#COVER:
colnames(cv17)[which(colnames(cv17)=="QUADRAT_Direction")]<-"Treatment"
colnames(cv18)[which(colnames(cv18)=="QUADRAT_Direction")]<-"Treatment"
colnames(cv19)[which(colnames(cv19)=="QUADRAT_Direction")]<-"Treatment"


# create reserve variable from PLOT_ID:
#COUNT:
ct17$reserve<-ct17$PLOT_ID
ct17$reserve[grep("J", ct17$reserve)]<-"J"
ct17$reserve[grep("M", ct17$reserve)]<-"M"

ct18$reserve<-ct18$PLOT_ID
ct18$reserve[grep("J", ct18$reserve)]<-"J"
ct18$reserve[grep("M", ct18$reserve)]<-"M"

ct19$reserve<-ct19$PLOT_ID
ct19$reserve[grep("J", ct19$reserve)]<-"J"
ct19$reserve[grep("M", ct19$reserve)]<-"M"

#COVER:
cv17$reserve<-cv17$PLOT_ID
cv17$reserve[grep("J", cv17$reserve)]<-"J"
cv17$reserve[grep("M", cv17$reserve)]<-"M"

cv18$reserve<-cv18$PLOT_ID
cv18$reserve[grep("J", cv18$reserve)]<-"J"
cv18$reserve[grep("M", cv18$reserve)]<-"M"

cv19$reserve<-cv19$PLOT_ID
cv19$reserve[grep("J", cv19$reserve)]<-"J"
cv19$reserve[grep("M", cv19$reserve)]<-"M"

# re-arrange columns:
#COUNT:
ct17<-ct17[,c(which(colnames(ct17) %in% c("DATE", "reserve")),which(!colnames(ct17) %in% c("DATE", "reserve")))]
ct18<-ct18[,c(which(colnames(ct18) %in% c("DATE", "reserve")),which(!colnames(ct18) %in% c("DATE", "reserve")))]
ct19<-ct19[,c(which(colnames(ct19) %in% c("DATE", "reserve")),which(!colnames(ct19) %in% c("DATE", "reserve")))]

#COVER:
cv17<-cv17[,c(which(colnames(cv17) %in% c("DATE", "reserve")),which(!colnames(cv17) %in% c("DATE", "reserve")))]
cv18<-cv18[,c(which(colnames(cv18) %in% c("DATE", "reserve")),which(!colnames(cv18) %in% c("DATE", "reserve")))]
cv19<-cv19[,c(which(colnames(cv19) %in% c("DATE", "reserve")),which(!colnames(cv19) %in% c("DATE", "reserve")))]


# remove unwanted columns:
#COUNT:
ct17$QUAD_ID<-NULL
ct18$QUAD_ID<-NULL
ct19$QUAD_ID<-NULL

#COVER:
cv17$QUAD_ID<-NULL
cv18$QUAD_ID<-NULL
cv19$QUAD_ID<-NULL

# factorise variables which should be factors and make control the baseline level:
#COUNT:
ct17$reserve<-as.factor(ct17$reserve)
ct17$PLOT_ID<-as.factor(ct17$PLOT_ID)
ct17$Treatment<-factor(ct17$Treatment,levels=c("C","A","B"))

ct18$reserve<-as.factor(ct18$reserve)
ct18$PLOT_ID<-as.factor(ct18$PLOT_ID)
ct18$Treatment<-factor(ct18$Treatment,levels=c("C","A","B"))

ct19$reserve<-as.factor(ct19$reserve)
ct19$PLOT_ID<-as.factor(ct19$PLOT_ID)
ct19$Treatment<-factor(ct19$Treatment,levels=c("C","A","B"))

#COVER:
cv17$reserve<-as.factor(cv17$reserve)
cv17$PLOT_ID<-as.factor(cv17$PLOT_ID)
cv17$Treatment<-factor(cv17$Treatment,levels=c("C","A","B"))


cv18$reserve<-as.factor(cv18$reserve)
cv18$PLOT_ID<-as.factor(cv18$PLOT_ID)
cv18$Treatment<-factor(cv18$Treatment,levels=c("C","A","B"))

cv19$reserve<-as.factor(cv19$reserve)
cv19$PLOT_ID<-as.factor(cv19$PLOT_ID)
cv19$Treatment<-factor(cv19$Treatment,levels=c("C","A","B"))


# replace COUNT categories with numbers:

# 16-50 	W	35
# 51-100	X	75
# >100	  Y	100


ct17d<-as.matrix(ct17[,which(colnames(ct17)=="Aca_ovi"):ncol(ct17)])
ct17d[which(ct17d=="W")]<-35
ct17d[which(ct17d=="X")]<-75
ct17d[which(ct17d=="Y")]<-100
ct17d<-data.frame(apply(ct17d,2,as.numeric)) #Warning messages: 1: In apply(ct17d, 2, as.numeric) : NAs introduced by coercion; 2: In apply(ct17d, 2, as.numeric) : NAs introduced by coercion

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

#replace COVER categories with numbers:

# >5%	A	5
# 5-10%	B	10
# 11-20%	C 	20
# 21-30%	D	30
# 31-40%	E	40
# 41-50%	F	50
# 51-60%	G	60
# 61-70%	H	70
# 71-80%	I	80
# 81-90%	J	90
# 91-100%	K	100

cv17d<-as.matrix(cv17[,which(colnames(cv17)=="Aca_ovi"):ncol(cv17)])
cv17d[which(cv17d=="A")]<-5
cv17d[which(cv17d=="B")]<-10
cv17d[which(cv17d=="C")]<-20
cv17d[which(cv17d=="D")]<-30
cv17d[which(cv17d=="E")]<-40
cv17d[which(cv17d=="F")]<-50
cv17d[which(cv17d=="G")]<-60
cv17d[which(cv17d=="H")]<-70
cv17d[which(cv17d=="I")]<-80
cv17d[which(cv17d=="J")]<-90
cv17d[which(cv17d=="K")]<-100
cv17d<-data.frame(apply(cv17d,2,as.numeric))
cv17d<-data.frame(cbind(cv17[,1:which(colnames(cv17)=="Treatment")],cv17d))

cv18d<-as.matrix(cv18[,which(colnames(cv18)=="Aca_ovi"):ncol(cv18)])
cv18d[which(cv18d=="A")]<-5
cv18d[which(cv18d=="B")]<-10
cv18d[which(cv18d=="C")]<-20
cv18d[which(cv18d=="D")]<-30
cv18d[which(cv18d=="E")]<-40
cv18d[which(cv18d=="F")]<-50
cv18d[which(cv18d=="G")]<-60
cv18d[which(cv18d=="H")]<-70
cv18d[which(cv18d=="I")]<-80
cv18d[which(cv18d=="J")]<-90
cv18d[which(cv18d=="K")]<-100
cv18d<-data.frame(apply(cv18d,2,as.numeric))
cv18d<-data.frame(cbind(cv18[,1:which(colnames(cv18)=="Treatment")],cv18d))

cv19d<-as.matrix(cv19[,which(colnames(cv19)=="Aca_ovi"):ncol(cv19)])
cv19d[which(cv19d=="A")]<-5
cv19d[which(cv19d=="B")]<-10
cv19d[which(cv19d=="C")]<-20
cv19d[which(cv19d=="D")]<-30
cv19d[which(cv19d=="E")]<-40
cv19d[which(cv19d=="F")]<-50
cv19d[which(cv19d=="G")]<-60
cv19d[which(cv19d=="H")]<-70
cv19d[which(cv19d=="I")]<-80
cv19d[which(cv19d=="J")]<-90
cv19d[which(cv19d=="K")]<-100
cv19d<-data.frame(apply(cv19d,2,as.numeric))
cv19d<-data.frame(cbind(cv19[,1:which(colnames(cv19)=="Treatment")],cv19d))

# remove unidentified species (for now... but we need to discuss whether these should be included in any analyses):
#COUNT:
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

#COVER:

cv17<-cv17[,-which(colnames(cv17) %in% unid)]
cv18<-cv18[,-which(colnames(cv18) %in% unid)]
cv19<-cv19[,-which(colnames(cv19) %in% unid)]

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

gdf$ylab<-c("All","Native","Exotic","Indicator","Significance A","Significance B","Common/Increaser","Significance X/Y","Significance Z","Native forb", "Exotic forb","Exotic annual forb","Exotic perennial forb","Native annual forb","Native perennial forb","Native non-leg. forb","Exotic non-leg. forb","Native leg. forb","Exotic leg. forb","Native grass","Exotic grass","Exotic annual grass","Exotic perennial grass","C3 grass","Native C3 grass","Native C4 grass","Exotic C3 grass","Sedge     /Rush")

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

# Notes on model selection to include in paper:

# We fit an interaction between treatment and date in all models to describe the before-after, control-impact nature of our design. 

# Our full model also included interactions between treatment and reserve , but we removed these if they were not significant, keeping only the treatment x date interaction and all main effects. 

# we did not expect the diversity differences to change over time, so we did not fit a year and reserve interaction

# Preliminary analysis showed no evidence of significant three way interactions between treatment, reserve and year. Thus, to avoid overly complex models we only fit first order interacions.

# Dcide the cut-off for analysis. There must be more records (i.e. number of species recorded) than the number of quadrats (48)? I.e. exclude two responses (exotic_perengrass and sed_rus):

gdf<-gdf[-which(gdf$rich_records<48),]
gdf<-tidy.df(gdf)
head(gdf); dim(gdf)
rahead(rich,4,7); dim(rich)
rahead(shan,4,7); dim(shan)

# three way interaction for RICHNESS:
data.set<-rich_sc

# scale date only:
rich_sc<-rich
rich_sc$DATE<-rich_sc$DATE-2017
rahead(rich_sc,4,7); dim(rich_sc)

# lists for storing model fits, coefficients, anova tables and model estimates:
fits.rich<-list()
coef.rich<-list()
anova.rich<-list()
preds.rich<-list()

# new data for model estimates (same for models with a date:treatment interaction only and models with a three way interaction; you can also use the same newdata frame for richness and shannon's):
nd1<-data.frame(DATE=rep(c(0,1,2),rep(3,3)),Treatment=as.factor(c("C","A","B")),reserve=c(rep("J",9),rep("M",9)))


for (i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  data.thisrun<-rich_sc[,resp.thisrun]
  form.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID)", sep="")
  
  # some functional groups have 0, its not going to run those functional groups here on
  if(sum(data.thisrun,na.rm=T)==0){
    fits.rich[[i]]<-NULL
    coef.rich[[i]]<-NULL
    anova.rich[[i]]<-NULL
    preds.rich[[i]]<-NULL
    next
  }
  
  m1<-lmer(formula = form.thisrun, data=data.set)
  summary(m1)
  anova(m1)
  
  m1_coef<-coef.ext(m1)
  m1_anova<-anova.ext(m1)
  
  # simplify model if the three way is not significant:
  if(m1_anova[which(m1_anova$term=="Treatment:DATE:reserve"),"p"]>0.05){
    
    # remove three way term from formula:
    form.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+(1|PLOT_ID)", sep="")
    
    # re-run model without three way:
    m1<-lmer(formula = form.thisrun, data=data.set)
    summary(m1)
    anova(m1)
    
    m1_coef<-coef.ext(m1)
    m1_anova<-anova.ext(m1)
    
  } # close if three way not signif
  
  fits.rich[[i]]<-m1
  coef.rich[[i]]<-m1_coef
  anova.rich[[i]]<-m1_anova
  
  # generate model predictions:
  
  m1_pr<-predictSE(mod=m1, newdata=nd1, se.fit = T)
  m1_pr<-data.frame(nd1,fit=m1_pr$fit,se=m1_pr$se.fit)
  m1_pr$lci<-m1_pr$fit-(m1_pr$se*1.96)
  m1_pr$uci<-m1_pr$fit+(m1_pr$se*1.96)
  head(m1_pr)
  
  preds.rich[[i]]<-m1_anova
  
} # close model

three.way.anovas.rich<-anova.rich #no significant three way interactions for richness except [19] (exotic_legherb), [26] (native_c4) and [5] (sigA)
  
# extract an individual model
summary(fits.rich[[24]])

# Three way interaction for DIVERSITY:
data.set<-shan_sc

#scale date only
shan_sc<-shan
shan_sc$DATE<-shan_sc$DATE-2017
rahead(shan_sc,4,7); dim(shan_sc)

for (i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  data.thisrun<-shan_sc[,resp.thisrun]
  form.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID)", sep="")
  
  # some functional groups have 0, its not going to run those functional groups here on
  if(sum(data.thisrun,na.rm=T)==0){
    coef.out[[i]]<-NULL
    anova.out[[i]]<-NULL
    preds.out[[i]]<-NULL
    next
  }
  
  m1<-lmer(formula = form.thisrun, data=data.set)
  summary(m1)
  anova(m1)
  
  m1_coef<-coef.ext(m1)
  m1_anova<-anova.ext(m1)
  coef.out[[i]]<-m1_coef
  anova.out[[i]]<-m1_anova
  
} # close model

three.way.anovas.div<-anova.out #no significant three way interaction for diversity except 19  (exotic_legherb)

# exotic leg herbs: (repeat for exotic_legherb diver, and for richness c4 and sigA)
# change data set and the response variable:
elh_mod<-lmer(exotic_legherb~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID), data=rich_sc)
summary(elh_mod)
elh_pr<-predictSE(elh_mod, nd1, se.fit = T)
elh_pr
elh_pr<-data.frame(nd1,fit=elh_pr$fit,se=elh_pr$se.fit)
elh_pr$lci<-elh_pr$fit-(elh_pr$se*1.96)
elh_pr$uci<-elh_pr$fit+(elh_pr$se*1.96)
head(elh_pr)

xofs<-0.2
arrowlgth<-0.02

dev.new(width=10,height=5,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(1,2), mar=c(2,4,4,1), mgp=c(2.5,1,0))

plot(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="C"]-xofs,elh_pr$fit[elh_pr$reserve=="J" & elh_pr$Treatment=="C"], pch=15, ylim=c(min(elh_pr$lci), max(elh_pr$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
arrows(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="C"]-xofs,elh_pr$lci[elh_pr$reserve=="J" & elh_pr$Treatment=="C"],elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="C"]-xofs,elh_pr$uci[elh_pr$reserve=="J" & elh_pr$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="A"],elh_pr$fit[elh_pr$reserve=="J" & elh_pr$Treatment=="A"], pch=15, col="red")
arrows(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="A"],elh_pr$lci[elh_pr$reserve=="J" & elh_pr$Treatment=="A"],elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="A"],elh_pr$uci[elh_pr$reserve=="J" & elh_pr$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="B"]+xofs,elh_pr$fit[elh_pr$reserve=="J" & elh_pr$Treatment=="B"], pch=15, col="blue")
arrows(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="B"]+xofs,elh_pr$lci[elh_pr$reserve=="J" & elh_pr$Treatment=="B"],elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="B"]+xofs,elh_pr$uci[elh_pr$reserve=="J" & elh_pr$Treatment=="B"], code=3, angle=90, length=arrowlgth)

legend(1.5,4,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)
mtext("Jerra",side=3, line=1)

# mulungarri:
plot(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="C"]-xofs,elh_pr$fit[elh_pr$reserve=="M" & elh_pr$Treatment=="C"], pch=15, ylim=c(min(elh_pr$lci), max(elh_pr$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
arrows(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="C"]-xofs,elh_pr$lci[elh_pr$reserve=="M" & elh_pr$Treatment=="C"],elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="C"]-xofs,elh_pr$uci[elh_pr$reserve=="M" & elh_pr$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="A"],elh_pr$fit[elh_pr$reserve=="M" & elh_pr$Treatment=="A"], pch=15, col="red")
arrows(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="A"],elh_pr$lci[elh_pr$reserve=="M" & elh_pr$Treatment=="A"],elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="A"],elh_pr$uci[elh_pr$reserve=="M" & elh_pr$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="B"]+xofs,elh_pr$fit[elh_pr$reserve=="M" & elh_pr$Treatment=="B"], pch=15, col="blue")
arrows(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="B"]+xofs,elh_pr$lci[elh_pr$reserve=="M" & elh_pr$Treatment=="B"],elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="B"]+xofs,elh_pr$uci[elh_pr$reserve=="M" & elh_pr$Treatment=="B"], code=3, angle=90, length=arrowlgth)
mtext("Mulangarri",side=3, line=1)



# new data for model estimates:
nd1<-data.frame(DATE=rep(c(0,1,2),rep(3,3)),Treatment=as.factor(c("C","A","B")),reserve=c(rep("J",9),rep("M",9)))

# lists for storing coefficients, anova tables and model estimates:
coef.out<-list()
anova.out<-list()
preds.out<-list()


  # SIMPLIFY MODEL: the interaction between treatment, date and reserve is not significant, remove:
  
if(which(rownames(anova(m1))=="Treatment:DATE:reserve")>0.05){
  form.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID)", sep="")
  m1<-lmer(formula = form.thisrun, data=data.set)
  summary(m1)
  anova(m1)
}

m1_coef<-coef.ext(m1)
m1_anova<-anova.ext(m1)
coef.out[[i]]<-m1_coef
anova.out[[i]]<-m1_anova

p1<-pred(model=m1, new.data = nd1,se.fit = T, type = "response")
preds.out[[i]]<-p1


   # close models

anova.out


## THREE WAY INTERACTION:
head(gdf); dim(gdf)
rahead(rich_sc,4,7); dim(rich_sc)
rahead(shan_sc,4,7); dim(shan_sc)

# Native leg. forb, exotic leg forb

resp.thisrun<-gdf$group[26]
form.threeway<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+Treatment:reserve+reserve:DATE+(1|PLOT_ID)", sep="")

nfmod2<-lmer(form.threeway,data=rich_sc)
summary(nfmod2)
anova(nfmod2)

form.simple<-paste(resp.thisrun,"~Treatment*DATE+reserve+(1|PLOT_ID)", sep="")
mod.simp<-lmer(form.simple,data=rich_sc)
summary(mod.simp)
anova(mod.simp); AIC(mod.simp)

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
  
  if (is.null(pred.thisrun)) next
  
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

save.image("03_Workspaces/STJW_analysis.RData")

# close analysis ----

#  COMPONENT 4 seed viability:    	# ----

sv<-read.table("00_Data/Formatted_data/seed_viability.txt", header=T)

rahead(rich,4,7); dim(rich)

# extract site data from richness data set
sdat<-rich[,2:4]
sdat$plot<-paste(sdat$PLOT_ID, sdat$Treatment, sep="")
sdat<-sdat[-which(duplicated(sdat$plot)),]
head(sdat); dim(sdat)

# combine with seed viability data:
sv$plot %in% sdat$plot
sv_site<-sdat[-which(!sdat$plot %in% sv$plot),]
sv_site$plot
head(sv_site)
head(sv); dim(sv)

sv<-merge(sv, sv_site, by="plot", all.x=T, all.y=F)
head(sv,3)

# fit models:
sv_mod1<-lmer(weight~Treatment+(1|reserve/PLOT_ID), data=sv)
summary(sv_mod1)
anova(sv_mod1)

sv_mod2<-lmer(germ7~Treatment+(1|reserve/PLOT_ID), data=sv)
summary(sv_mod2)
anova(sv_mod2)

sv_mod3<-lmer(germ21~Treatment+(1|reserve/PLOT_ID), data=sv)
summary(sv_mod3)
anova(sv_mod3)

sv_nd<-data.frame(Treatment=factor(levels(sv$Treatment),levels=c("C","A","B")))
sv_nd

# model estimates:
sv_pr1<-predictSE(sv_mod1, newdata = sv_nd, se.fit = T)
sv_pr1<-data.frame(sv_nd,fit=sv_pr1$fit,se=sv_pr1$se.fit)
sv_pr1$lci<-sv_pr1$fit-(sv_pr1$se*1.96)
sv_pr1$uci<-sv_pr1$fit+(sv_pr1$se*1.96)

sv_pr2<-predictSE(sv_mod2, newdata = sv_nd, se.fit = T)
sv_pr2<-data.frame(sv_nd,fit=sv_pr2$fit,se=sv_pr2$se.fit)
sv_pr2$lci<-sv_pr2$fit-(sv_pr2$se*1.96)
sv_pr2$uci<-sv_pr2$fit+(sv_pr2$se*1.96)

sv_pr3<-predictSE(sv_mod3, newdata = sv_nd, se.fit = T)
sv_pr3<-data.frame(sv_nd,fit=sv_pr3$fit,se=sv_pr3$se.fit)
sv_pr3$lci<-sv_pr3$fit-(sv_pr3$se*1.96)
sv_pr3$uci<-sv_pr3$fit+(sv_pr3$se*1.96)

# PLOT estimates:

dev.new(width=8,height=8,noRStudioGD = T,dpi=80, pointsize=16)
par(mfrow=c(2,2), mar=c(3,5,2,1), mgp=c(3.2,1,0))

plot(1:3, sv_pr1$fit, ylim=c(min(sv_pr1$lci),max(sv_pr1$uci)), las=1, type="p", xlim=c(0.75, 3.25), pch=20, xlab="", xaxt="n", ylab="Sample weight (g)")
arrows(1:3, sv_pr1$lci, 1:3, sv_pr1$uci, code=3, length=0.05, angle=90)
axis(side=1, at=1:3, labels=levels(sv_pr1$Treatment), xlab="")
title(xlab="Treatment", mgp=c(2,1,0))
text(0.75, max(sv_pr1$uci),paste("P = ",round(anova(sv_mod1)$"Pr(>F)",3),sep=""), adj=0)
mtext("A", side=3, line=0.5, cex=1, adj=0)

plot(1:3, sv_pr2$fit, ylim=c(min(sv_pr2$lci),max(sv_pr2$uci)), las=1, type="p", xlim=c(0.75, 3.25), pch=20, xlab="", xaxt="n", ylab="Germination at 7d / g")
arrows(1:3, sv_pr2$lci, 1:3, sv_pr2$uci, code=3, length=0.05, angle=90)
axis(side=1, at=1:3, labels=levels(sv_pr2$Treatment), xlab="")
title(xlab="Treatment", mgp=c(2,1,0))
text(0.75, max(sv_pr2$uci),paste("P = ",round(anova(sv_mod2)$"Pr(>F)",3),sep=""), adj=0)
mtext("B", side=3, line=0.5, cex=1, adj=0)

plot(1:3, sv_pr3$fit, ylim=c(min(sv_pr3$lci),max(sv_pr3$uci)), las=1, type="p", xlim=c(0.75, 3.25), pch=20, xlab="", xaxt="n", ylab="Germination at 21d / g")
arrows(1:3, sv_pr3$lci, 1:3, sv_pr3$uci, code=3, length=0.05, angle=90)
axis(side=1, at=1:3, labels=levels(sv_pr3$Treatment), xlab="")
title(xlab="Treatment", mgp=c(2,1,0))
text(0.75, max(sv_pr3$uci),paste("P = ",round(anova(sv_mod3)$"Pr(>F)",3),sep=""), adj=0)
mtext("C", side=3, line=0.5, cex=1, adj=0)

# PLOT raw data:

dev.new(width=8,height=8,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(2,2), mar=c(2,4,4,1), mgp=c(2.5,1,0))
plot(sv$Treatment, sv$weight, ylab="Weight")
plot(sv$Treatment, sv$germ7, ylab="Germ. 7")
plot(sv$Treatment, sv$germ21, ylab="Germ. 21")

# close component 4 ----

#  COMPONENT 2 spray drift:    	# ----

sdrift<-read.table("00_Data/Formatted_data/spray_drift.txt", header=T)
head(sdrift)

# fit model:
sd_mod1<-lm(percent_sprayed~treatment*method, data=sdrift)
summary(sd_mod1)
anova(sd_mod1)

# model estimates:

sd_nd<-data.frame(treatment=rep(unique(sdrift$treatment),rep(3,2)),method=unique(sdrift$method)[order(unique(sdrift$method))])

sd_pr<-predict(sd_mod1, newdata = sd_nd, se.fit = T)
sd_nd<-data.frame(sd_nd, fit=sd_pr$fit, se=sd_pr$se.fit)
sd_nd$lci<-sd_nd$fit-(sd_nd$se*1.96)
sd_nd$uci<-sd_nd$fit+(sd_nd$se*1.96)

# PLOT estimates:

dev.new(width=5,height=4,noRStudioGD = T,dpi=80, pointsize=14)
par(mfrow=c(1,1), mar=c(3.5,4,2.5,4.5), mgp=c(2.8,1,0))

plot(1:6, sd_nd$fit, ylim=c(min(sd_nd$lci),max(sd_nd$uci)), las=1, type="p", xlim=c(0.75, 6.25), pch=20, xlab="", xaxt="n", ylab="Percent sprayed",col=c("darkorange","darkturquoise","darkolivegreen2"))

arrows(1:6, sd_nd$lci, 1:6, sd_nd$uci, code=3, length=0.05, angle=90)
points(1:6, sd_nd$fit, pch=20, cex=2, col=c("darkorange","darkturquoise","darkolivegreen2"))
axis(side=1, at=c(2, 5), labels=c("A","B"))
title(xlab="Treatment", mgp=c(2.3,1,0))
arrows(c(3.5),0,c(3.5),1.5, length=0, col="grey70")

p.trt<-round(anova(sd_mod1)[1,5],3)
p.mth<-round(anova(sd_mod1)[2,5],3)
p.int<-round(anova(sd_mod1)[3,5],3)
p.mth<-"< 0.001"

title(main=paste("P values: Treatment =",p.trt,"\nMethod =",p.mth,"; Int. =",p.int), font.main=1, adj=0, cex.main=1, line=0.5)

par(xpd=NA)
legend(6.5,1,legend=c("Coarse","Fine","Spot"), col=c("darkorange","darkturquoise","darkolivegreen2"), pch=20, bty="n", pt.cex = 2)
par(xpd=T)

# PLOT raw data:

dev.new(width=5,height=4,noRStudioGD = T,dpi=80, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,6), mgp=c(2.8,1,0))
boxplot(sdrift$percent_sprayed~sdrift$method*as.factor(sdrift$treatment), col=c("darkorange","darkturquoise","darkolivegreen2"), ,las=2, xlab="", xaxt="n", ylab="Percent sprayed",at=c(0.7,1.7,2.7,4.3,5.3,6.3))
axis(side=1, at=c(2, 5), labels=c("A","B"))
title(xlab="Treatment", mgp=c(2.5,1,0))
arrows(c(3.5),0,c(3.5),1.5, length=0, col="grey70")
par(xpd=NA)
legend(7.5,1,legend=c("Coarse","Fine","Spot"), col=c("darkorange","darkturquoise","darkolivegreen2"), pch=15, bty="n", pt.cex = 3)
par(xpd=T)

# close component 2 ----





