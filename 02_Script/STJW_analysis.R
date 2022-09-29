
# ------------------------------------ #
# ---------- STJW ANALYSIS  ---------- #
# ------------------------------------ #

### Analysis of herbicide control experiment from STJW project
### Author: Annabel Smith & Raagini Muddaiah

# Load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""), function(x) source(x)))

# Load libraries:
library(lme4); library(vegan); library(AICcmodavg); library(lmerTest); library(glmmADMB); library(mgcv)

# Load workspace
load("03_Workspaces/stjw_analysis_R1.RData")

#  POST-ANALYSIS, summaries for paper:    	# ----

rahead(ct_dat,3,7); dim(ct_dat)
rahead(cv_dat,3,7); dim(cv_dat)

head(pinfo,3)

# Themeda = 1 species (The_tri)
# Austrostipa = 3 species (Aus_big, Aus_den, Aus_sca)
# Rytidosperma = 8 species ("Ryt_cae" "Ryt_car" "Ryt_lae" "Ryt_rac" "Ryt_sp." "Ryt_sp3" "Ryt_sp4" "Ryt_sp2")

aus_rty<-cv_dat[,c(1:4,c(grep("Ryt_", colnames(cv_dat)),grep("Aus_", colnames(cv_dat))))]
aus_rty$Ryt_Aus_ALL<-rowSums(aus_rty[,5:ncol(aus_rty)])
head(aus_rty,3); dim(aus_rty)

aus_rty$tag<-paste(aus_rty$DATE,aus_rty$reserve,aus_rty$PLOT_ID,aus_rty$Treatment, sep="_")
ar_all<-aus_rty[,c("Ryt_Aus_ALL", "tag")]
head(ar_all,3); dim(ar_all)

them<-cv_dat[,c(1:4,grep("The_", colnames(cv_dat)))]
them$tag<-paste(them$DATE,them$reserve,them$PLOT_ID,them$Treatment, sep="_")
head(them,3); dim(them)

table(ar_all$tag %in% them$tag)
table(them$tag %in% ar_all$tag)

mul_grass<-merge(them, ar_all, by="tag", all.x=T, all.y=F)
mul_grass<-mul_grass[which(mul_grass$reserve=="M"),]
mul_grass<-tidy.df(mul_grass)
mul_grass$RA_dom<-mul_grass$Ryt_Aus_ALL>mul_grass$The_tri
head(mul_grass,3); dim(mul_grass)

table(mul_grass$RA_dom, mul_grass$DATE)

# close data summaries ----

#  IMPORT, clean & transform data:    	# ----

data_dir<-"00_Data/Formatted_data/"
dir(data_dir)

# PLANT INFO:
pinfo<-read.table(paste(data_dir,"plant_info.txt",sep="/"),header=T)
pinfo$Duration<-as.factor(pinfo$Duration)
pinfo$func_grp<-as.factor(pinfo$func_grp)
pinfo$Significance<-as.factor(pinfo$Significance)
pinfo$Status<-as.factor(pinfo$Status)
head(pinfo,3); dim(pinfo)
length(which(duplicated(pinfo$Sp))) # should be zero

# remove unidentified species:
unid<-pinfo$Sp[grep("Uni_", pinfo$Sp)]
pinfo<-pinfo[-which(pinfo$Sp %in% unid),]
pinfo<-tidy.df(pinfo)
head(pinfo, 3); dim(pinfo)

# Re-classify Desmodium varians as a native leg forb
# Need to double check this is OK with RM. 
# PlantNet classifies it as a "Prostrate or climbing herb"
# And it is the ONLY species in the Legume category
table(pinfo$func_grp)
pinfo[which(pinfo$Sp=="Des_var"),]$func_grp<-"Leguminous_herb"
pinfo<-tidy.df(pinfo)
table(pinfo$func_grp)
head(pinfo,2); dim(pinfo)

# ALL data:
a17<-read.table(paste(data_dir,"all_2017.txt",sep="/"),header=F)
a18<-read.table(paste(data_dir,"all_2018.txt",sep="/"),header=F)
a19<-read.table(paste(data_dir,"all_2019.txt",sep="/"),header=F)

# Data column order for each quadrat:
# 1 = cover
# 2 = count (main, unshaded, clump)
# 3 = count (second column, shaded for Tricoryne and Lomandra only, tuft)

# For count we will use the main, unshaded clump count, so remove the third column from all data files:
# Warning: run only once:

a17<-a17[,-seq(4,ncol(a17), by=3)]
a18<-a18[,-seq(4,ncol(a18), by=3)]
a19<-a19[,-seq(4,ncol(a19), by=3)]

rahead(a17,6,6); dim(a17)
rahead(a18,6,6); dim(a18)
rahead(a19,6,6); dim(a19)

# Separate count and cover:

ct17<-a17[,c(1,seq(3,ncol(a17),by=2))]
cv17<-a17[,c(1,seq(2,ncol(a17),by=2))]

ct18<-a18[,c(1,seq(3,ncol(a18),by=2))]
cv18<-a18[,c(1,seq(2,ncol(a18),by=2))]

ct19<-a19[,c(1,seq(3,ncol(a19),by=2))]
cv19<-a19[,c(1,seq(2,ncol(a19),by=2))]

rahead(ct17,6,6); dim(ct17)
rahead(ct18,6,6); dim(ct18)
rahead(ct19,6,6); dim(ct19)

rahead(cv17,6,6); dim(cv17)
rahead(cv18,6,6); dim(cv18)
rahead(cv19,6,6); dim(cv19)

# replace COUNT categories with numbers:

# 16-50 	W	35
# 51-100	X	75
# >100	  Y	100

# ** WARNING: numeric subsets
ct17d<-as.matrix(ct17[4:nrow(ct17),2:ncol(ct17)])
ct17d[which(ct17d=="W")]<-35
ct17d[which(ct17d=="X")]<-75
ct17d[which(ct17d=="Y")]<-100
ct17site<-ct17[1:3,]
ct17site[,1:10]
ct17<-data.frame(cbind(ct17[4:nrow(ct17),1]),ct17d)
names(ct17)<-names(ct17site)
ct17<-data.frame(rbind(ct17site, ct17))

ct18d<-as.matrix(ct18[4:nrow(ct18),2:ncol(ct18)])
ct18d[which(ct18d=="W")]<-35
ct18d[which(ct18d=="X")]<-75
ct18d[which(ct18d=="Y")]<-100
ct18site<-ct18[1:3,]
ct18site[,1:10]
ct18<-data.frame(cbind(ct18[4:nrow(ct18),1]),ct18d)
names(ct18)<-names(ct18site)
ct18<-data.frame(rbind(ct18site, ct18))

ct19d<-as.matrix(ct19[4:nrow(ct19),2:ncol(ct19)])
ct19d[which(ct19d=="W")]<-35
ct19d[which(ct19d=="X")]<-75
ct19d[which(ct19d=="Y")]<-100
ct19site<-ct19[1:3,]
ct19site[,1:10]
ct19<-data.frame(cbind(ct19[4:nrow(ct19),1]),ct19d)
names(ct19)<-names(ct19site)
ct19<-data.frame(rbind(ct19site, ct19))

# replace COVER categories with numbers:

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

# ** WARNING: numeric subsets
cv17d<-as.matrix(cv17[4:nrow(cv17),2:ncol(cv17)])
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
cv17site<-cv17[1:3,]
cv17site[,1:10]
cv17<-data.frame(cbind(cv17[4:nrow(cv17),1]),cv17d)
names(cv17)<-names(cv17site)
cv17<-data.frame(rbind(cv17site, cv17))

cv18d<-as.matrix(cv18[4:nrow(cv18),2:ncol(cv18)])
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
cv18site<-cv18[1:3,]
cv18site[,1:10]
cv18<-data.frame(cbind(cv18[4:nrow(cv18),1]),cv18d)
names(cv18)<-names(cv18site)
cv18<-data.frame(rbind(cv18site, cv18))

cv19d<-as.matrix(cv19[4:nrow(cv19),2:ncol(cv19)])
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
cv19site<-cv19[1:3,]
cv19site[,1:10]
cv19<-data.frame(cbind(cv19[4:nrow(cv19),1]),cv19d)
names(cv19)<-names(cv19site)
cv19<-data.frame(rbind(cv19site, cv19))

# preserve species names and transpose:

## COUNT
ct17sp<-ct17$V1[4:nrow(ct17)]
ct17t<-t(ct17)
colnames(ct17t)<-ct17t[1,]
ct17t<-ct17t[2:nrow(ct17t),]
ct17<-data.frame(ct17t)
ct17<-tidy.df(ct17)

ct18sp<-ct18$V1[4:nrow(ct18)]
ct18t<-t(ct18)
colnames(ct18t)<-ct18t[1,]
ct18t<-ct18t[2:nrow(ct18t),]
ct18<-data.frame(ct18t)
ct18<-tidy.df(ct18)

ct19sp<-ct19$V1[4:nrow(ct19)]
ct19t<-t(ct19)
colnames(ct19t)<-ct19t[1,]
ct19t<-ct19t[2:nrow(ct19t),]
ct19<-data.frame(ct19t)
ct19<-tidy.df(ct19)

## COVER
cv17sp<-cv17$V1[4:nrow(cv17)]
cv17t<-t(cv17)
colnames(cv17t)<-cv17t[1,]
cv17t<-cv17t[2:nrow(cv17t),]
cv17<-data.frame(cv17t)
cv17<-tidy.df(cv17)

cv18sp<-cv18$V1[4:nrow(cv18)]
cv18t<-t(cv18)
colnames(cv18t)<-cv18t[1,]
cv18t<-cv18t[2:nrow(cv18t),]
cv18<-data.frame(cv18t)
cv18<-tidy.df(cv18)

cv19sp<-cv19$V1[4:nrow(cv19)]
cv19t<-t(cv19)
colnames(cv19t)<-cv19t[1,]
cv19t<-cv19t[2:nrow(cv19t),]
cv19<-data.frame(cv19t)
cv19<-tidy.df(cv19)

# Make numeric: 
ct17[,4:ncol(ct17)]<-apply(ct17[,4:ncol(ct17)],2,as.numeric)
ct18[,4:ncol(ct18)]<-apply(ct18[,4:ncol(ct18)],2,as.numeric)
ct19[,4:ncol(ct19)]<-apply(ct19[,4:ncol(ct19)],2,as.numeric)

cv17[,4:ncol(cv17)]<-apply(cv17[,4:ncol(cv17)],2,as.numeric)
cv18[,4:ncol(cv18)]<-apply(cv18[,4:ncol(cv18)],2,as.numeric)
cv19[,4:ncol(cv19)]<-apply(cv19[,4:ncol(cv19)],2,as.numeric)

# Rytidosperma sp (ALL) is only in ct19, not the other years, so we will remove this column:
ct19sp[which(!ct19sp %in% ct17sp)]
cv19sp[which(!cv19sp %in% cv17sp)]

ct19<-ct19[,-which(colnames(ct19)=="Rytidosperma.sp..ALL.")]
ct19sp<-ct19sp[-which(ct19sp=="Rytidosperma sp (ALL)")]

cv19<-cv19[,-which(colnames(cv19)=="Rytidosperma.sp..ALL.")]
cv19sp<-cv19sp[-which(cv19sp=="Rytidosperma sp (ALL)")]

# Check that all species columns are the same:
table(ct17sp==ct18sp)
table(ct17sp==ct19sp)
table(ct18sp==ct19sp)

table(cv17sp==cv18sp)
table(cv17sp==cv19sp)
table(cv18sp==cv19sp)

# Find out which species have no data in ANY of the years:

## COUNT

sum_ct17<-colSums(ct17[4:ncol(ct17)],na.rm = T)
sum_ct18<-colSums(ct18[4:ncol(ct18)],na.rm = T)
sum_ct19<-colSums(ct19[4:ncol(ct19)],na.rm = T)

sum_ctdat<-data.frame(name17=names(sum_ct17),sum_ct17=sum_ct17,name18=names(sum_ct18),sum_ct18=sum_ct18,name19=names(sum_ct19),sum_ct19=sum_ct19)
sum_ctdat<-tidy.df(sum_ctdat)
table(sum_ctdat$name17==sum_ctdat$name18)
table(sum_ctdat$name17==sum_ctdat$name19)
sum_ctdat$name18<-NULL
sum_ctdat$name19<-NULL
sum_ctdat$all_yrs<-rowSums(sum_ctdat[,c("sum_ct17","sum_ct18","sum_ct19")])
sum_ctdat$name_orig<-ct17sp
check.rows(sum_ctdat[,c("name17","name_orig")])
check.rows(sum_ctdat)
head(sum_ctdat)

## COVER

sum_cv17<-colSums(cv17[4:ncol(cv17)],na.rm = T)
sum_cv18<-colSums(cv18[4:ncol(cv18)],na.rm = T)
sum_cv19<-colSums(cv19[4:ncol(cv19)],na.rm = T)

sum_cvdat<-data.frame(name17=names(sum_cv17),sum_cv17=sum_cv17,name18=names(sum_cv18),sum_cv18=sum_cv18,name19=names(sum_cv19),sum_cv19=sum_cv19)
sum_cvdat<-tidy.df(sum_cvdat)
table(sum_cvdat$name17==sum_cvdat$name18)
table(sum_cvdat$name17==sum_cvdat$name19)
sum_cvdat$name18<-NULL
sum_cvdat$name19<-NULL
sum_cvdat$all_yrs<-rowSums(sum_cvdat[,c("sum_cv17","sum_cv18","sum_cv19")])
sum_cvdat$name_orig<-cv17sp
check.rows(sum_cvdat[,c("name17","name_orig")])
check.rows(sum_cvdat)
head(sum_cvdat)

# These are the species with no data in any year:
nodat_ct<-sum_ctdat[which(sum_ctdat$all_yrs==0),]
nodat_ct_sp<-nodat_ct$name17
head(nodat_ct); dim(nodat_ct)

nodat_cv<-sum_cvdat[which(sum_cvdat$all_yrs==0),]
nodat_cv_sp<-nodat_cv$name17
head(nodat_cv); dim(nodat_cv)

# Make sure that the species with no data are the same for count and cover:
table(nodat_ct_sp==nodat_cv_sp)

# Remove species with no data:
# Reduces from 322 columns to 126 (123 species):
ct17<-ct17[,-which(colnames(ct17) %in% nodat_ct_sp)]
ct18<-ct18[,-which(colnames(ct18) %in% nodat_ct_sp)]
ct19<-ct19[,-which(colnames(ct19) %in% nodat_ct_sp)]

cv17<-cv17[,-which(colnames(cv17) %in% nodat_cv_sp)]
cv18<-cv18[,-which(colnames(cv18) %in% nodat_cv_sp)]
cv19<-cv19[,-which(colnames(cv19) %in% nodat_cv_sp)]

# Remove unidentified species:
# Reduces to 120 species:
ct17<-ct17[,-grep("Unidentified",colnames(ct17))]
ct18<-ct18[,-grep("Unidentified",colnames(ct18))]
ct19<-ct19[,-grep("Unidentified",colnames(ct19))]

cv17<-cv17[,-grep("Unidentified",colnames(cv17))]
cv18<-cv18[,-grep("Unidentified",colnames(cv18))]
cv19<-cv19[,-grep("Unidentified",colnames(cv19))]

# Make sure count and cover colnames are the same
table(colnames(ct17)==colnames(cv17))
table(colnames(ct18)==colnames(cv18))
table(colnames(ct19)==colnames(cv19))

# Create data set of plant names in new, filtered data
# Since the colnames are lining up precisely, this can be done on only one dataset and applied to all:
spdat<-sum_ctdat[-which(sum_ctdat$name17 %in% nodat_ct_sp),]
# Should include only those with data:
which(spdat$all_yrs==0)

# remove unidentified species from name data:
spdat<-spdat[-grep("Unidentified", spdat$name_orig),]
spdat<-tidy.df(spdat)
head(spdat); dim(spdat)

# Duplicated in the data:
# Rytidosperma sp.
# Wahlenbergia sp.
spdat[which(duplicated(spdat$name_orig)),]
spdat[which(spdat$name_orig %in% c("Rytidosperma sp.","Wahlenbergia sp.")),]

# "Rytidosperma.sp." and "Wahlenbergia sp." are duplicated in the data. R re-names duplictae colnames with ".1", so the second records are coming up as "Rytidosperma.sp..1" and "Wahlenbergia.sp..1"

# FIX Rytidosperma sp.
# In 2017 and 2019 the first Rytidosperma.sp. col was used, while in 2018 the second was used. 
# For 2018 data, the second Rytidosperma sp. column is the correct one (called "Rytidosperma.sp..1" in R notation). Take this column and overwirte the first Rytidosperma sp. column (called "Rytidosperma.sp."). Then delete the second one. 
# For 2017 and 2019 simply delete the second Rytidosperma sp. ("Rytidosperma.sp..1")

# Rytidosperma 2018
# Overwrite the first col using data from the second
ct18[,which(colnames(ct18)=="Rytidosperma.sp.")]<-ct18[,which(colnames(ct18)=="Rytidosperma.sp..1")]
cv18[,which(colnames(cv18)=="Rytidosperma.sp.")]<-cv18[,which(colnames(cv18)=="Rytidosperma.sp..1")]

# Check it worked:
ct18[,which(colnames(ct18)=="Rytidosperma.sp.")]
ct18[,which(colnames(ct18)=="Rytidosperma.sp..1")]

cv18[,which(colnames(cv18)=="Rytidosperma.sp.")]
cv18[,which(colnames(cv18)=="Rytidosperma.sp..1")]

# Then delete the second column:
ct18[,which(colnames(ct18)=="Rytidosperma.sp..1")]<-NULL
cv18[,which(colnames(cv18)=="Rytidosperma.sp..1")]<-NULL

# Rytidosperma 2017 and 2019
# Delete the second Rytidosperma col ("Rytidosperma.sp..1")
ct17[,which(colnames(ct17)=="Rytidosperma.sp.")]
ct17[,which(colnames(ct17)=="Rytidosperma.sp..1")]
ct19[,which(colnames(ct19)=="Rytidosperma.sp.")]
ct19[,which(colnames(ct19)=="Rytidosperma.sp..1")]

cv17[,which(colnames(cv17)=="Rytidosperma.sp.")]
cv17[,which(colnames(cv17)=="Rytidosperma.sp..1")]
cv19[,which(colnames(cv19)=="Rytidosperma.sp.")]
cv19[,which(colnames(cv19)=="Rytidosperma.sp..1")]

ct17[,which(colnames(ct17)=="Rytidosperma.sp..1")]<-NULL
ct19[,which(colnames(ct19)=="Rytidosperma.sp..1")]<-NULL
cv17[,which(colnames(cv17)=="Rytidosperma.sp..1")]<-NULL
cv19[,which(colnames(cv19)=="Rytidosperma.sp..1")]<-NULL

# 119 species:
rahead(ct18,6,6); dim(ct18)
rahead(ct17,6,6); dim(ct17)
rahead(ct19,6,6); dim(ct19)

rahead(cv18,6,6); dim(cv18)
rahead(cv17,6,6); dim(cv17)
rahead(cv19,6,6); dim(cv19)

# Also remove from spdat:
spdat[grep("Rytid",spdat$name_orig),]
spdat<-spdat[-which(spdat$name17=="Rytidosperma.sp..1"),]
spdat<-tidy.df(spdat)
head(spdat); dim(spdat)

# FIX Whalenbergia
# There are 2 x "Whalenbergia sp." in the data. The first has data in all years and the second has data in 2017 and 2019
# In this case, keep the first one in all data files and change the second one to "Wahlenbergia. sp. 1"
# Since R has already done this in its automated column re-naming, I've added the second Wahlenbergia ("Wahlenbergia. sp. 1") to the master pinfo, making it 119 species
# Also need to update it in spdat:
head(spdat); dim(spdat)
spdat[grep("Wahl",spdat$name_orig),]
spdat$name_orig[which(spdat$name17=="Wahlenbergia.sp..1")]<-"Wahlenbergia sp. 1"
spdat<-tidy.df(spdat)
head(spdat); dim(spdat)

# Find which data set has duplicates:

rahead(ct17,6,6); dim(ct17)
rahead(ct18,6,6); dim(ct18)
rahead(ct19,6,6); dim(ct19)

which(duplicated(colnames(ct17[,4:ncol(ct17)])))
which(duplicated(colnames(ct18[,4:ncol(ct18)])))
which(duplicated(colnames(ct19[,4:ncol(ct19)])))

# Make data set of species names that don't appear in pinfo:
head(pinfo,2); dim(pinfo)
problem_names<-data.frame(dat_name=spdat$name_orig[which(!spdat$name_orig %in% pinfo$Species)])

# problem_names have been updated directly in pinfo so that the names match exactly the names in the original data. 
# For some species, it was simply a missing full stop after the sp.: Rytidosperma sp., Dichelachne sp., Microtis sp. 
# Others had different notation for the names in pinfo:
# "Hordeum (Critesion) sp."
# "Anthosachne scaber (Elymus)" changed to "Anthosachne scaber" # Lomandra coriacea filiformis changed to Lomandra filiformis coriacea
# "Linaria sp. (arvensis?)" changed to "Linaria spp (arvensis?)"
# "Linum trigynum" changed to "Linum trigynum ? (yellow)". Also updated pinfo using data from plantnet. Perennial changed to annual. func_grp changed to Herb
# Trailing spaced removed after Rumex brownii in original data, not pinfo
# Rytidosperma sp2 changed to Rytidosperma sp. 2
# "Sonchus oleruceus" changed to "Sonchus oleraceus"
# "Vittadinia green cuneata" changed to "Vittadinia cuneata (green)"
# "Zornia dyctocarpa" changed to Zornia dyctiocarpa

# this should now be zero:
nrow(problem_names)

# Use names in spdat and pinfo to replace the column names with species codes:
head(pinfo,2); dim(pinfo)

table(spdat$name_orig %in% pinfo$Species)
pin2<-pinfo[,c("Species","Sp")]

spdat<-merge(spdat, pin2, by.x="name_orig", by.y="Species", all.x=T, all.y=F)

# spdat is sorted by the name_orig col
# the data are in their orignal order and are unsorted
is.unsorted(spdat$name_orig)
is.unsorted(colnames(ct17)[4:ncol(ct17)])

cnames17<-data.frame(cname=colnames(ct17)[4:ncol(ct17)], index=1:length(colnames(ct17)[4:ncol(ct17)]))
head(cnames17); dim(cnames17)

# Order spdat so it's in the same order as the data:
spdat<-merge(spdat, cnames17, by.x="name17", by.y="cname", all.x=T, all.y=F)
spdat<-spdat[order(spdat$index),]
spdat<-tidy.df(spdat)
head(spdat); dim(spdat)

# Check that the new order lines up with all data files:
table(spdat$name17==colnames(ct17)[4:ncol(ct17)])
table(spdat$name17==colnames(ct18)[4:ncol(ct18)])
table(spdat$name17==colnames(ct19)[4:ncol(ct19)])

table(spdat$name17==colnames(cv17)[4:ncol(cv17)])
table(spdat$name17==colnames(cv18)[4:ncol(cv18)])
table(spdat$name17==colnames(cv19)[4:ncol(cv19)])

# Replace the colnames with species codes
colnames(ct17)[4:ncol(ct17)]<-spdat$Sp
colnames(ct18)[4:ncol(ct18)]<-spdat$Sp
colnames(ct19)[4:ncol(ct19)]<-spdat$Sp

colnames(cv17)[4:ncol(cv17)]<-spdat$Sp
colnames(cv18)[4:ncol(cv18)]<-spdat$Sp
colnames(cv19)[4:ncol(cv19)]<-spdat$Sp

# Make sure they're all numeric:
table(apply(ct17[2:nrow(ct17),4:ncol(ct17)],2,function(x)is.numeric(x)))
table(apply(ct18[2:nrow(ct17),4:ncol(ct18)],2,function(x)is.numeric(x)))
table(apply(ct19[2:nrow(ct17),4:ncol(ct19)],2,function(x)is.numeric(x)))

# Make sure there's no NAs:
# The only NAs should be six plots in 2017 for Tri_ela
table(unlist(apply(ct17[2:nrow(ct17),4:ncol(ct17)],2,function(x) which(is.na(x)))))
table(unlist(apply(ct18[2:nrow(ct18),4:ncol(ct18)],2,function(x) which(is.na(x)))))
table(unlist(apply(ct19[2:nrow(ct19),4:ncol(ct19)],2,function(x) which(is.na(x)))))

# these should all be TRUE:

# Are all names matching across yearS? 
#COUNT:
table(colnames(ct17)[4:ncol(ct17)]==colnames(ct18)[4:ncol(ct18)])
table(colnames(ct17)[4:ncol(ct17)]==colnames(ct19)[4:ncol(ct19)])
table(colnames(ct18)[4:ncol(ct18)]==colnames(ct19)[4:ncol(ct19)])

table(colnames(ct17)[4:ncol(ct17)] %in% colnames(ct18)[4:ncol(ct18)])
table(colnames(ct17)[4:ncol(ct17)] %in% colnames(ct19)[4:ncol(ct19)])
table(colnames(ct18)[4:ncol(ct18)] %in% colnames(ct19)[4:ncol(ct19)])

#COVER:
table(colnames(cv17)[4:ncol(cv17)]==colnames(cv18)[4:ncol(cv18)])
table(colnames(cv17)[4:ncol(cv17)]==colnames(cv19)[4:ncol(cv19)])
table(colnames(cv18)[4:ncol(cv18)]==colnames(cv19)[4:ncol(cv19)])

table(colnames(cv17)[4:ncol(cv17)] %in% colnames(cv18)[4:ncol(cv18)])
table(colnames(cv17)[4:ncol(cv17)] %in% colnames(cv19)[4:ncol(cv19)])
table(colnames(cv18)[4:ncol(cv18)] %in% colnames(cv19)[4:ncol(cv19)])

rahead(ct17,6,6); dim(ct17)
rahead(ct18,6,6); dim(ct18)
rahead(ct19,6,6); dim(ct19)

# save.image("03_Workspaces/stjw_analysis_R1.RData")

# close import data ----

#  FORMAT data:    	# ----

# Add year:
# COUNT:
ct17$DATE<-2017
ct18$DATE<-2018
ct19$DATE<-2019

# COVER:
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
# COUNT:
ct17$reserve<-ct17$PLOT_ID
ct17$reserve[grep("J", ct17$reserve)]<-"J"
ct17$reserve[grep("M", ct17$reserve)]<-"M"

ct18$reserve<-ct18$PLOT_ID
ct18$reserve[grep("J", ct18$reserve)]<-"J"
ct18$reserve[grep("M", ct18$reserve)]<-"M"

ct19$reserve<-ct19$PLOT_ID
ct19$reserve[grep("J", ct19$reserve)]<-"J"
ct19$reserve[grep("M", ct19$reserve)]<-"M"

# COVER:
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
# COUNT:
ct17<-ct17[,c(which(colnames(ct17) %in% c("DATE", "reserve")),which(!colnames(ct17) %in% c("DATE", "reserve")))]
ct18<-ct18[,c(which(colnames(ct18) %in% c("DATE", "reserve")),which(!colnames(ct18) %in% c("DATE", "reserve")))]
ct19<-ct19[,c(which(colnames(ct19) %in% c("DATE", "reserve")),which(!colnames(ct19) %in% c("DATE", "reserve")))]

# COVER:
cv17<-cv17[,c(which(colnames(cv17) %in% c("DATE", "reserve")),which(!colnames(cv17) %in% c("DATE", "reserve")))]
cv18<-cv18[,c(which(colnames(cv18) %in% c("DATE", "reserve")),which(!colnames(cv18) %in% c("DATE", "reserve")))]
cv19<-cv19[,c(which(colnames(cv19) %in% c("DATE", "reserve")),which(!colnames(cv19) %in% c("DATE", "reserve")))]

# remove unwanted columns:
# COUNT:
ct17$QUAD_ID<-NULL
ct18$QUAD_ID<-NULL
ct19$QUAD_ID<-NULL

# COVER:
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

# IMPUTE Tricoryne

# For Tricoryne, there were six plots in 2017, where only a tuft count (shaded column) and no clump was recorded
# We will use data from the other clump/tuft observations to impute clump counts at these sites
# M1A, M1B, M1C, M4A, J4C, J6A

# Tri_ela data:
te<-read.table(paste(data_dir, "Tri_ela_dat.txt", sep="/"), header=T)
te$reserve<-te$PLOT_ID 
te$reserve<- substr(te$reserve,1,1)
te$reserve<-as.factor(te$reserve)
te<-te[-which(te$type=="cover"),]
te<-tidy.df(te)

# Divide clump and tuft
te_clump<-te[te$type=="clump",]
te_tuft<-te[te$type=="tuft",]
head(te_clump,3); dim(te_clump)
head(te_tuft,3); dim(te_tuft)

# Combine so that clump and tuft have their own columns:
te<-te_clump
colnames(te)[which(colnames(te)=="Tri_ela")]<-"clump"
te$type<-NULL

# Make sure quadrats are in the same order:
table(te$QUAD_ID==te_tuft$QUAD_ID)

# Then combine and re-name columns:
te<-cbind(te, te_tuft$Tri_ela)
colnames(te)[which(colnames(te)=="te_tuft$Tri_ela")]<-"tuft"

# Add yr
te$yr<-te$year-min(te$year)

# Correlation:
# plot(te$clump, te$tuft)
cor.test(te$clump, te$tuft)

# Remove outlier:
te<-te[-which(te$clump>25),]
head(te,3); dim(te)

# Fit models:
temod1<-lm(clump~tuft*yr+reserve,data=te) 
temod2<-lm(clump~tuft+yr+reserve,data=te) 
temod3<-lm(clump~tuft*yr,data=te) 
temod4<-lm(clump~tuft+yr,data=te) 
summary(temod1)
summary(temod2)
AIC(temod1); AIC(temod2); AIC(temod3); AIC(temod4)

# Predict from model:
te.nd<-data.frame(tuft=seq(min(te$tuft),max(te$tuft),length.out = 50), yr=rep(seq(min(te$yr),max(te$yr)),rep(50,3)),reserve=rep(c("J","M"),rep(150,2)))
te.pr<-predict(temod1, newdata = te.nd)
te.nd$predicted_clump<-te.pr

# Plot predictions:
# plot(te.nd$tuft[te.nd$yr==0 & te.nd$reserve=="J"], te.nd$predicted_clump[te.nd$yr==0  & te.nd$reserve=="J"], type="l", col="red")
# lines(te.nd$tuft[te.nd$yr==0 & te.nd$reserve=="M"], te.nd$predicted_clump[te.nd$yr==0  & te.nd$reserve=="M"], lty=2, col="red")
# lines(te.nd$tuft[te.nd$yr==1 & te.nd$reserve=="J"], te.nd$predicted_clump[te.nd$yr==1  & te.nd$reserve=="J"], lty=1, col="orange")
# lines(te.nd$tuft[te.nd$yr==1 & te.nd$reserve=="M"], te.nd$predicted_clump[te.nd$yr==1  & te.nd$reserve=="M"], lty=2, col="orange")
# lines(te.nd$tuft[te.nd$yr==2 & te.nd$reserve=="J"], te.nd$predicted_clump[te.nd$yr==2  & te.nd$reserve=="J"], lty=1, col="yellow")
# lines(te.nd$tuft[te.nd$yr==2 & te.nd$reserve=="M"], te.nd$predicted_clump[te.nd$yr==2  & te.nd$reserve=="M"], lty=2, col="yellow")

# Use yr zero only:
head(te.nd)
te17<-te.nd[te.nd$yr==0,]
te17<-tidy.df(te17)

# PREDICTED DATA, round variables:
te17$tuft<-round(te17$tuft,0)
te17$predicted_clump<-round(te17$predicted_clump,0)
head(te17); dim(te17)

# Original Tri_ela data, including tuft counts:
te_dat17<-te[which(te$year==2017),]
te_dat17<-tidy.df(te_dat17)
imp_dat<-te_dat17[which(is.na(te_dat17$clump)),]

# Merge original data with predicted data:
imp_dat<-merge(imp_dat, te17, by.x=c("reserve","tuft"), by.y=c("reserve","tuft"), all.x=T, all.y=F)

# check it:
head(imp_dat); dim(imp_dat)
te17[which(te17$tuft==113),]
te17[which(te17$tuft==51),]
te17[which(te17$tuft==98),]
te17[which(te17$tuft==113),]

# There's no value for 74, so use 75
imp_dat$predicted_clump[which(imp_dat$tuft==74)]<-te17$predicted_clump[which(te17$tuft==75 & te17$reserve=="M")]
imp_dat<-imp_dat[,c("reserve","PLOT_ID","QUADRAT_Direction", "predicted_clump")]
colnames(imp_dat)[which(colnames(imp_dat)=="QUADRAT_Direction")]<-"Treatment"

# Impute missing data in 2017:
te.missing<-ct17[,c(1:4,which(colnames(ct17)=="Tri_ela"))]
te.missing[which(is.na(te.missing$Tri_ela)),]

te.missing<-merge(te.missing, imp_dat, by=c("reserve","PLOT_ID","Treatment"), all.x=T, all.y=F)
te.missing$Tri_ela[which(is.na(te.missing$Tri_ela))]<-te.missing$predicted_clump[which(is.na(te.missing$Tri_ela))]
te.missing$predicted_clump<-NULL
te.missing$DATE<-NULL
head(te.missing); dim(te.missing)

# Merge with original data
# This is going to change the row order of the main data frame AND the column order:
ct17$Tri_ela<-NULL
ct17<-merge(ct17, te.missing, by=c("reserve","PLOT_ID","Treatment"))

# Re-order ALL data frames:
ct17<-ct17[,c("DATE","reserve","PLOT_ID","Treatment",colnames(ct17)[5:ncol(ct17)][order(colnames(ct17)[5:ncol(ct17)])])]
ct18<-ct18[,c("DATE","reserve","PLOT_ID","Treatment",colnames(ct18)[5:ncol(ct18)][order(colnames(ct18)[5:ncol(ct18)])])]
ct19<-ct19[,c("DATE","reserve","PLOT_ID","Treatment",colnames(ct19)[5:ncol(ct19)][order(colnames(ct19)[5:ncol(ct19)])])]

cv17<-cv17[,c("DATE","reserve","PLOT_ID","Treatment",colnames(cv17)[5:ncol(cv17)][order(colnames(cv17)[5:ncol(cv17)])])]
cv18<-cv18[,c("DATE","reserve","PLOT_ID","Treatment",colnames(cv18)[5:ncol(cv18)][order(colnames(cv18)[5:ncol(cv18)])])]
cv19<-cv19[,c("DATE","reserve","PLOT_ID","Treatment",colnames(cv19)[5:ncol(cv19)][order(colnames(cv19)[5:ncol(cv19)])])]

# Re-order rows:
ct17<-ct17[order(ct17$reserve,ct17$PLOT_ID,ct17$Treatment),]
ct18<-ct18[order(ct18$reserve,ct18$PLOT_ID,ct18$Treatment),]
ct19<-ct19[order(ct19$reserve,ct19$PLOT_ID,ct19$Treatment),]

cv17<-cv17[order(cv17$reserve,cv17$PLOT_ID,cv17$Treatment),]
cv18<-cv18[order(cv18$reserve,cv18$PLOT_ID,cv18$Treatment),]
cv19<-cv19[order(cv19$reserve,cv19$PLOT_ID,cv19$Treatment),]

# TIDY:
ct17<-tidy.df(ct17)
ct18<-tidy.df(ct18)
ct19<-tidy.df(ct19)

cv17<-tidy.df(cv17)
cv18<-tidy.df(cv18)
cv19<-tidy.df(cv19)

rahead(ct17,4,10); dim(ct17)
rahead(ct18,4,10); dim(ct18)
rahead(ct19,4,10); dim(ct19)

rahead(cv17,4,10); dim(cv17)
rahead(cv18,4,10); dim(cv18)
rahead(cv19,4,10); dim(cv19)

# save.image("03_Workspaces/stjw_analysis_diversity.RData")

# Combine count and cover data:

# COUNT:
rahead(ct17,4,10); dim(ct17)
rahead(ct18,4,10); dim(ct18)
rahead(ct19,4,10); dim(ct19)

# Check colname order:
table(colnames(ct17)==colnames(ct18))
table(colnames(ct17)==colnames(ct19))
table(colnames(ct18)==colnames(ct19))
colnames(ct17)[which(!colnames(ct17)==colnames(ct19))]
colnames(ct19)[which(!colnames(ct19)==colnames(ct17))]

ct_dat<-rbind(ct17, ct18, ct19)
rahead(ct_dat,3,7); dim(ct_dat)

# COVER
rahead(cv17,3,7); dim(ct17)
rahead(cv18,3,7); dim(ct18)
rahead(cv19,3,7); dim(ct19)

# Check colname order:
table(colnames(cv17)==colnames(cv18))
table(colnames(cv17)==colnames(cv19))
table(colnames(cv18)==colnames(cv19))

cv_dat<-rbind(cv17, cv18, cv19)
rahead(cv_dat,3,7); dim(cv_dat)
rahead(ct_dat,3,7); dim(ct_dat)

# Check that all species have been recorded
head(spdat); dim(spdat)
spcounts<-spdat[,c("Sp","sum_ct17","sum_ct18","sum_ct19")]
spcounts<-tidy.df(spcounts)
which(rowSums(spcounts[,2:4])==0)
head(pinfo,3); dim(pinfo)
table(pinfo$Status)

# Calculate the number of species for each group:
head(gdf,2)
no_sp<-gdf[,c("group","ylabn")]
no_sp$no_sp<-unlist(lapply(gdf$group,FUN=function(x) length(get(x))))
# write.table(no_sp,"no_sp.txt",quote=F, row.names = F, sep="\t")

# save.image("03_Workspaces/stjw_analysis_R1.RData")

# close format data ----

#  species DIVERSITY GROUPS:    	# ----

head(pinfo,3); dim(pinfo)
str(pinfo)
table(pinfo$func_grp)
exotic_legherb
native_legherb

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
# Indicator is actually a direct combination of sigA and sigB and should not be included in the analysis. I haven't updated it in the code, I've just removed it from the MS. 
group_df<-data.frame(group=c("all","native","exotic","indic","sigA","sigB","sigC","sigXY","sigZ","native_herb","exotic_herb","exann_herb","exper_herb","natann_herb","natper_herb","native_nonlegherb","exotic_nonlegherb","native_legherb","exotic_legherb","native_grass","exotic_grass","exotic_anngrass","exotic_perengrass","c3_grass","native_c3","native_c4","exotic_c3","sed_rus"))

# save.image("03_Workspaces/stjw_analysis_R1.RData")

# close diversity groups ----

#  calculate species DIVERSITY:    	# ----

head(pinfo,3); dim(pinfo)
head(group_df)

rahead(ct17,3,7); dim(ct17)
rahead(ct18,3,7); dim(ct18)
rahead(ct19,3,7); dim(ct19)

# Add species richness and Shannon's index to count data:

# CALCULATE diversity for each year separately:

ct17_site.data<-ct17[,1:which(colnames(ct17)=="Treatment")]
div17<-calc.div(ct17, ct17_site.data)
rich17<-div17$rich
shan17<-div17$shan
simp17<-div17$simp
invsimp17<-div17$invsimp

ct18_site.data<-ct18[,1:which(colnames(ct18)=="Treatment")]
div18<-calc.div(ct18, ct18_site.data)
rich18<-div18$rich
shan18<-div18$shan
simp18<-div18$simp
invsimp18<-div18$invsimp

ct19_site.data<-ct19[,1:which(colnames(ct19)=="Treatment")]
div19<-calc.div(ct19, ct19_site.data)
rich19<-div19$rich
shan19<-div19$shan
simp19<-div19$simp
invsimp19<-div19$invsimp

# Then combine:

head(rich17[,1:10],3); dim(rich17)
head(rich18[,1:10],3); dim(rich18)
head(rich19[,1:10],3); dim(rich19)

rich<-rbind(rich17,rich18,rich19)
rahead(rich,4,7); dim(rich)

shan<-rbind(shan17, shan18, shan19)
rahead(shan,4,7); dim(shan)

simp<-rbind(simp17, simp18, simp19)
rahead(simp,4,7); dim(simp)

# When richness is zero, shan==0, simp==1 and invsimp==Inf
# When richness is one, shan==0, simp==0, and invsimp==1

# Representing a single species as zero inflates the zeros in the data, so shan and simp are not ideal; simp also represents zero as one which is not ideal. 

# But, if we replace the 'Inf' with zero in invsimp, then zero will be zero and one will be one. 

invsimp<-rbind(invsimp17, invsimp18, invsimp19)
invsimp<-cbind(invsimp[1:4],apply(invsimp[,5:ncol(invsimp)],2,function(x) ifelse(x=="Inf",0,x)))
rahead(invsimp,4,7); dim(invsimp)

rahead(rich,3,6); dim(rich)
rahead(shan,3,6); dim(shan)
rahead(simp,3,6); dim(simp)
rahead(invsimp,3,6); dim(invsimp)

# plot r.ships between diversity metrics
dev.new(width=6, height=6, dpi=80, pointsize=16,noRStudioGD = T)
par(mfrow=c(2,2),mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(rich$exotic,shan$exotic, pch=20, xlab="richness", ylab="shannon", las=1)
plot(rich$exotic,simp$exotic, pch=20, xlab="richness", ylab="simpsons", las=1)
plot(rich$exotic,invsimp$exotic, pch=20, xlab="richness", ylab="invsimpsons", las=1)

# Create group data frame and data to assist with deciding cut-offs for analysis. 
group_df$rich_records<-colSums(rich[,5:length(rich)])
gdf<-group_df 
gdf<-tidy.df(gdf)
gdf

gdf$ylab<-c("All","Native","Exotic","Indicator","Significance A","Significance B","Common/Increaser","Significance X/Y","Significance Z","Native forb", "Exotic forb","Exotic annual forb","Exotic perennial forb","Native annual forb","Native perennial forb","Native non-leg. forb","Exotic non-leg. forb","Native leg. forb","Exotic leg. forb","Native grass","Exotic grass","Exotic annual grass","Exotic perennial grass","C3 grass","Native C3 grass","Native C4 grass","Exotic C3 grass","Sedge/Rush")

# save.image("03_Workspaces/stjw_analysis_R1.RData")

# **** SCALE DATE:

rich_sc<-rich
rich_sc$DATE<-rich_sc$DATE-2017

shan_sc<-shan
shan_sc$DATE<-shan_sc$DATE-2017

simp_sc<-simp
simp_sc$DATE<-simp_sc$DATE-2017

invsimp_sc<-invsimp
invsimp_sc$DATE<-invsimp_sc$DATE-2017

rahead(rich_sc,4,7); dim(rich_sc)
rahead(shan_sc,4,7); dim(shan_sc)
rahead(simp_sc,4,7); dim(simp_sc)
rahead(invsimp_sc,4,7); dim(invsimp_sc)

# Examine data distributions:

# RICHNESS:
rahead(rich_sc,6,6)

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(2,4,4,1), mgp=c(2.5,1,0))

data.update<-rich_sc[,5:ncol(rich_sc)]

for (i in 1:ncol(data.update)){
  sp.thisrun<-colnames(data.update)[i]
  data.thisrun<-data.update[,sp.thisrun]
  hist(data.thisrun, main=sp.thisrun)
} # close

# Shannon's (zero inflated)

rahead(shan_sc,6,6)

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(2,4,4,1), mgp=c(2.5,1,0))

data.update<-shan_sc[,5:ncol(shan_sc)]

for (i in 1:ncol(data.update)){
  sp.thisrun<-colnames(data.update)[i]
  data.thisrun<-data.update[,sp.thisrun]
  hist(data.thisrun, main=sp.thisrun)
} # close

# Simpsons

rahead(simp_sc,6,6)

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(2,4,4,1), mgp=c(2.5,1,0))

data.update<-simp_sc[,5:ncol(simp_sc)]

for (i in 1:ncol(data.update)){
  sp.thisrun<-colnames(data.update)[i]
  data.thisrun<-data.update[,sp.thisrun]
  hist(data.thisrun, main=sp.thisrun)
} # close

# Inverse Simpsons

rahead(invsimp_sc,6,6)

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(2,4,4,1), mgp=c(2.5,1,0))

data.update<-invsimp_sc[,5:ncol(invsimp_sc)]

for (i in 1:ncol(data.update)){
  sp.thisrun<-colnames(data.update)[i]
  data.thisrun<-data.update[,sp.thisrun]
  hist(data.thisrun, main=sp.thisrun)
} # close

# save.image("03_Workspaces/stjw_analysis_R1.RData")

# close diversity calculation ----

#  Spray drift ANALYSIS:    	# ----

sdrift<-read.table("00_Data/Formatted_data/spray_drift.txt", header=T)
sdrift$method<-factor(sdrift$method, levels=c("Spot", "Fine","Coarse"))
sdrift$treatment<-as.factor(sdrift$treatment)
head(sdrift,3); dim(sdrift)

# fit beta regression model in mgcv:
# The Intercept is treatmentA (no STJW) and Spot spray
sd_mod2a<-gam(percent_sprayed~treatment*method, family=betar,data=sdrift)
summary(sd_mod2a)
anova(sd_mod2a)

### ----- CONTRASTS

summary(sd_mod2a)$p.table

# For spray drift, there are two treatments (A, B, relating to no and mod/high STJW density) and three methods (coarse, fine, spot) making 6 treatments total and 15 contrasts. 
choose(6,2)

### Create a "unique model matrix":

spray.z.mat<-lm(percent_sprayed~treatment*method,data=sdrift,x=T)$x
head(spray.z.mat)
spray.uzm<-unique(spray.z.mat)
rownames(spray.uzm)<-1:nrow(spray.uzm)
spray.uzm; dim(spray.uzm)

# Re-order the uzm so that the baseline (i.e. the rows with zeros) come out first - this will make naming and interpretation easier
spray.uzm<-spray.uzm[c(2,1,3,5,4,6),]
rownames(spray.uzm)<-1:nrow(spray.uzm)
spray.uzm; dim(spray.uzm)

### Create a "difference matrix"

# each row must be a vector with a length equal to the number of rows in the coefficient table. I.e. 6 numbers for this spray drift model. 15 contrasts=15 rows. Each row will specify ONE contrast. 

summary(sd_mod2a)$p.table

spray.z.diff<-rbind(
  # first year treatment differences:
  c(-1,1,0,0,0,0),
  c(-1,0,1,0,0,0),
  c(-1,0,0,1,0,0),
  c(-1,0,0,0,1,0),
  c(-1,0,0,0,0,1),
  c(0,-1,1,0,0,0),
  c(0,-1,0,1,0,0),
  c(0,-1,0,0,1,0),
  c(0,-1,0,0,0,1),
  c(0,0,-1,1,0,0),
  c(0,0,-1,0,1,0),
  c(0,0,-1,0,0,1),
  c(0,0,0,-1,1,0),
  c(0,0,0,-1,0,1),
  c(0,0,0,0,-1,1)
)

# Now we have our unique model matrix
spray.uzm

# And our difference matrix:
spray.z.diff

# Create names for the contrasts:
# This has to relate to how the uzm and the z.diff are set up, so make sure they are re-ordered in the same way?

spray.lt<-c(paste(levels(sdrift$method),levels(sdrift$treatment)[1],sep=""),paste(levels(sdrift$method),levels(sdrift$treatment)[2],sep=""))
spray.lt<-paste(combn(spray.lt,2)[1,],combn(spray.lt,2)[2,],sep=" vs ")

# Difference estimates, SE and CI:
spray.nd1<-data.frame(diff.est(sd_mod2a,spray.uzm, spray.z.diff),Contrast=spray.lt)
spray.nd1$diff<-ifelse(sign(spray.nd1$lci)==sign(spray.nd1$uci),1,0)
summary(sd_mod2a)$p.table

# beta reg model estimates:

# order new data as: spot, fine, coarse, so that it's in the same order as the plant results:
sd_nd<-data.frame(treatment=rep(unique(sdrift$treatment),rep(3,2)),method=unique(sdrift$method)[c(2,1,3)])

sd_pr<-predict(sd_mod2a, newdata = sd_nd, se.fit = T, type="response")
sd_nd<-data.frame(sd_nd, fit=sd_pr$fit, se=sd_pr$se.fit)
sd_nd$lci<-sd_nd$fit-(sd_nd$se*1.96)
sd_nd$uci<-sd_nd$fit+(sd_nd$se*1.96)
head(sd_nd)

# PLOT estimates:

dev.new(width=6, height=4, dpi=100, pointsize=16,noRStudioGD = T)
par(mfrow=c(1,1),mar=c(4,4,1,1), oma=c(0,0,0,6), mgp=c(2.5,1,0))

plot(1:6, sd_nd$fit, ylim=c(min(sd_nd$lci),max(sd_nd$uci)), las=1, type="p", xlim=c(0.75, 6.25),pch=15, xlab="", xaxt="n", ylab="Proportion sprayed",col=c("red","blue","cornflowerblue"))

arrows(1:6, sd_nd$lci, 1:6, sd_nd$uci, code=3, length=0.05, angle=90,col=c("red","blue","cornflowerblue"))
points(1:6, sd_nd$fit, pch=15, cex=1, col=c("red","blue","cornflowerblue"))
axis(side=1, at=c(2, 5), labels=c("None","Mod-high"))
title(xlab=bquote(italic(H.~perforatum)~density), mgp=c(2.3,1,0))
arrows(c(3.5),0,c(3.5),1.5, length=0, col="grey70")

p.trt<-round(anova(sd_mod1)[1,5],3)
p.mth<-round(anova(sd_mod1)[2,5],3)
p.int<-round(anova(sd_mod1)[3,5],3)
p.mth<-"< 0.001"

par(xpd=NA)
legend(6.5,0.45,legend=c("Spot spray","Fine boom","Coarse boom"), col=c("red","blue","cornflowerblue"), pch=15, bty="n", pt.cex = 2.7)
par(xpd=F)

# Use the contrasts to manually add grouping numbers
spray.nd1[,c("Contrast","diff")]
par(xpd=NA)
text(1:6,1.07,labels = rep(c("a","b","b"),2), pos=2,offset=-0.2)
par(xpd=F)

# PLOT raw data:

dev.new(width=5,height=4,noRStudioGD = T,dpi=80, pointsize=14)
par(mfrow=c(1,1), mar=c(4,4,1,6), mgp=c(2.8,1,0))
boxplot(sdrift$percent_sprayed~sdrift$method*as.factor(sdrift$treatment), col=c("darkorange","darkturquoise","darkolivegreen2"),las=2, xlab="", xaxt="n", ylab="Percent sprayed",at=c(0.7,1.7,2.7,4.3,5.3,6.3))
axis(side=1, at=c(2, 5), labels=c("A","B"))
title(xlab="Treatment", mgp=c(2.5,1,0))
arrows(c(3.5),0,c(3.5),1.5, length=0, col="grey70")
par(xpd=NA)
legend(7.5,1,legend=c("Coarse boom","Fine boom","Spot spray"), col=c("darkorange","darkturquoise","darkolivegreen2"), pch=15, bty="n", pt.cex = 3)
par(xpd=T)

# close spray drift ----

#  Richness & Diversity ANALYSIS:    	# ----

# MODELLING SUMMARY:

# We used the inverse Simpson's Diversity Index as a measure of biodiversity because it avoided inflating the zeros in the data set; richness values of 1 were equal to one, unlike other metrics where richness values of one are equal to zero. 

# Functional groups with a high proportion of zeros also had very little data overall (Proportion_ZERO.pdf), creating a risk of missing important results because of data limitations. For these groups, we therefore fit a binomial GLMM to the data, to model the probability of detecting any plant in that group. 

# We then modelled each functional group using a negative binomial distribution for richness (count data) and a gamma distribution for diversity (continuous data, positive values only), to accommodate for heterogeneity in the data (i.e. overdispersion). We excluded exotic_perengrass and sed_rus from these analyses because they had fewer total records than the number of quadrats in the study (n=5 and n=27, respectively). Although zeros were not included in the diversity analysis, a binomial model was fitted for any group with a high proportion (> 20%) of zeros.  

# https://www.umass.edu/landeco/teaching/ecodata/schedule/distributions.pdf
# The Gamma is an excellent choice when only positive real numbers are possible (note the normal is not lower bounded like Gamma).
# The Gamma is similarly used in phenomenological fashion with continuous, positive data having too much variance and a right skew; in other words, an overdispersed normal distribution with a right skew.
# The Gamma distribution is the continuous counterpart of the negative binomial, which recall is used to describe overdispersed Poisson with heterogeneous data

# Overdispersed count data:
# The negative binomial distribution ... is more often used phenomenologically to describe a patchy or clustered distribution with no intrinsic upper limit that has more variance thanthe Poisson; i.e., an overdispersed Poisson distribution in which the variance is greater than themean

rahead(invsimp_sc,3,6); dim(invsimp_sc)
apply(invsimp_sc[,5:ncol(invsimp_sc)],2,range)
head(gdf); dim(gdf)

gdf$no_rich0<-apply(rich[,5:ncol(rich)], 2, function(x) length(which(x==0)))
gdf$no_invsimp0<-apply(invsimp[,5:ncol(invsimp)], 2, function(x) length(which(x==0)))

# The zeros all match (this should be zero):
nrow(gdf[!(gdf$no_rich0==0)==(gdf$no_invsimp0==0),])

# Plot cut-off for models with binomial only (lots of zeros, little data), two part models (enough data but zero inflated), and positive values only (lots of data, few zeros)

dev.new(width=6, height=6, dpi=80, pointsize=16,noRStudioGD = T)
par(mar=c(4,4,1,1))
plot(gdf$rich_records,gdf$no_rich0/144, pch=20, xlab="", ylab="Proportion zeros", las=1, xlim=c(0,3100))
arrows(0,0.2,max(gdf$rich_records),0.2, col="red",code=0)
# text(gdf$rich_records,gdf$no_rich0/144,labels=gdf$group, adj=0, pos=3, offset=0.5,cex=0.5, srt=50)
title(xlab="Total number of observations", mgp=c(2.5,1,0))
cor.test(gdf$rich_records,gdf$no_rich0)

# BINOMIAL CUT-OFF # 1: Groups with 20% or more zeros (i.e. 29 out of 144) will ONLY get a binomial model

gdf[,c("group","no_rich0","no_shan0", "fit_bin")]
gdf$fit_bin<-ifelse(gdf$no_rich0>=29, "yes", "no")
head(gdf,3); dim(gdf)

# POSITIVE CUT-OFF # 2: don't fit models for +ve data to groups where the total count is less than the number of quadrats. 
rahead(rich,3,6)

# This excludes exotic_perengrass and sed_rus which have almost no data:
gdf$fit_pos<-ifelse(gdf$rich_records<=length(unique(paste(rich$PLOT_ID,rich$Treatment,sep=""))), "no", "yes")

# new data for model estimates (same for models with a date:treatment interaction only and models with a three way interaction; you can also use the same newdata frame for richness and shannon's):
nd1<-data.frame(DATE=rep(c(0,1,2),rep(3,3)),Treatment=factor(c("C","A","B"),levels=c("C","A","B")),reserve=factor(c(rep("J",9),rep("M",9)),levels=c("J","M")))

# lists for storing model fits, coefficients, anova tables and model estimates:
fits.binom<-list()
coef.binom<-list()
preds.binom<-list()
anova.binom<-list()

fits.rich<-list()
coef.rich<-list()
preds.rich<-list()
anova.rich<-list()

fits.invsimp<-list()
coef.invsimp<-list()
preds.invsimp<-list()
anova.invsimp<-list()

# Add summary data to gdf:
head(gdf,3); dim(gdf)

gdf$bin_3wayP<-NA
gdf$bin_2wayP<-NA

gdf$rich_3wayP<-NA
gdf$rich_2wayP<-NA

gdf$invsimp_3wayP<-NA
gdf$invsimp_2wayP<-NA

head(gdf,3); dim(gdf)

rahead(rich_sc,4,7); dim(rich_sc)
rahead(invsimp_sc,4,7); dim(invsimp_sc)

# save workspace:
# save.image("03_Workspaces/stjw_analysis_R1.RData")

## *** Nov 2021 update:
# Re-running with correct levels for nd1 c("C","A","B"). The binomial model, using predictSE, automatically put the levels at the same as the fitted model. 
# The problem was in the glmmadmb models where we were using the custom pred() function with hand-calculated back transformations; this used the original levels of nd1 which was causing problems with plotting. 

# save.image("03_Workspaces/stjw_analysis_R1.RData")

## RUN FULL ANALYSIS LOOP:

## TAKES approx 15 min; i==27 and 28 have errors (all C3 grass is n)

for (i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  bin.thisrun<-gdf$fit_bin[i]
  pos.thisrun<-gdf$fit_pos[i]
  
  # Positive models are only run for species where the count is greater than the number of quadrats:
  if (pos.thisrun=="no"){
   
    fits.rich[[i]]<-NULL
    coef.rich[[i]]<-NULL
    preds.rich[[i]]<-NULL
    anova.rich[[i]]<-NULL
    
    fits.invsimp[[i]]<-NULL
    coef.invsimp[[i]]<-NULL
    preds.invsimp[[i]]<-NULL
    anova.invsimp[[i]]<-NULL
    
  } # close no binomial
  
  if (bin.thisrun=="yes"){
    
    data.set<-rich_sc
    
    data.set[,resp.thisrun]<-ifelse(data.set[,resp.thisrun]==0,0,1)
    data.thisrun<-data.set[,resp.thisrun]
    rahead(data.set,6,6); dim(data.set)
    
    form.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID)", sep="")
    
    # some functional groups have 0, its not going to run those functional groups here on
    if(sum(data.thisrun,na.rm=T)==0){
      fits.shan[[i]]<-NULL
      coef.shan[[i]]<-NULL
      anova.shan[[i]]<-NULL
      preds.shan[[i]]<-NULL
      next
    }
    
    m1<-glmer(formula = form.thisrun, family="binomial", data=data.set)
    summary(m1)
    anova(m1)
    
    # run model without three-way:
    form.twoway<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+(1|PLOT_ID)", sep="")
    mod_twoway<-glmer(formula = form.twoway, family="binomial", data=data.set)
    
    m1_coef<-coef.ext(m1)
    m1_anova<-anova(mod_twoway, m1)
    
    p_anova<-m1_anova[2,which(colnames(m1_anova)=="Pr(>Chisq)")]
    gdf$bin_3wayP[i]<-p_anova
    head(gdf,3); dim(gdf)
    
    anova.binom[[i]]<-m1_anova
    
    # simplify model if the three way is not significant:
    if(p_anova>0.05){
      
      # remove three way term from formula:
      form.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+(1|PLOT_ID)", sep="")
      
      # re-run model without three way:
      m1<-glmer(formula = form.thisrun, family="binomial", data=data.set)
      summary(m1)
      
      # run model without any interaction (to figure out if the two way interaction is significant):
      form.noint<-paste(resp.thisrun,"~Treatment+DATE+reserve+(1|PLOT_ID)", sep="")
      mod_noint<-glmer(formula = form.noint, family="binomial", data=data.set)
      
      int_term_anova<-anova(mod_noint, m1)
      int_term_p<-int_term_anova[2,which(colnames(m1_anova)=="Pr(>Chisq)")]
      anova.binom[[i]]<-int_term_anova
      
      m1_coef<-coef.ext(m1)
      
      gdf$bin_2wayP[i]<-int_term_p
      head(gdf,3); dim(gdf)
      
    } # close if three way not signif
    
    fits.binom[[i]]<-m1
    coef.binom[[i]]<-m1_coef
    
    # generate model predictions:
    
    m1_pr<-predictSE(mod=m1,newdata=nd1,type="response", se.fit = T)
    m1_pr<-data.frame(nd1,fit=m1_pr$fit,se=m1_pr$se.fit)
    m1_pr$lci<-m1_pr$fit-(m1_pr$se*1.96)
    m1_pr$uci<-m1_pr$fit+(m1_pr$se*1.96)
    head(m1_pr)
    
    preds.binom[[i]]<-m1_pr
    
  } # close binomial model
  
  if (pos.thisrun=="yes"){
    
    rich.thisrun<-rich_sc[,resp.thisrun]
    invsimp.thisrun<-invsimp_sc[,resp.thisrun]
    
    ## FIT three-way models:
    
    threeway.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID)", sep="")
    
    m3way_rich<-glmmadmb(as.formula(threeway.thisrun), family="nbinom", data=rich_sc)
    summary(m3way_rich)
    
    m3way_invsimp<-glmmadmb(as.formula(threeway.thisrun), family="gamma", data=invsimp_sc[which(invsimp.thisrun>0),])
    summary(m3way_invsimp)
    
    # STORE fits and coefficients from three-way
    
    fits.rich[[i]]<-m3way_rich
    fits.invsimp[[i]]<-m3way_invsimp
    
    coef.rich[[i]]<-coef.ext(m3way_rich)
    coef.invsimp[[i]]<-coef.ext(m3way_invsimp)
  
    ## FIT two-way models:
    
    twoway.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+(1|PLOT_ID)", sep="")
    
    m2way_rich<-glmmadmb(as.formula(twoway.thisrun), family="nbinom", data=rich_sc)
    summary(m2way_rich)
    
    m2way_invsimp<-glmmadmb(as.formula(twoway.thisrun), family="gamma", data=invsimp_sc[which(invsimp.thisrun>0),])
    summary(m2way_invsimp)
    
    # Get p-value for three way:
    
    aov3way_rich<-anova(m2way_rich,m3way_rich)
    aov3way_invsimp<-anova(m2way_invsimp,m3way_invsimp)
    
    head(gdf,3)
    
    anova.rich[[i]]<-aov3way_rich
    anova.invsimp[[i]]<-aov3way_invsimp
    
    p3way_rich<-round(aov3way_rich[2,"Pr(>Chi)"],4)
    p3way_invsimp<-round(aov3way_invsimp[2,"Pr(>Chi)"],4)
    
    gdf$rich_3wayP[i]<-p3way_rich
    gdf$invsimp_3wayP[i]<-p3way_invsimp
    
    # Store current model:
    m1_rich<-m3way_rich
    m1_invsimp<-m3way_invsimp
    
    # simplify RICHNESS if the three way is not significant:
    if(p3way_rich>0.05){
      
      # Overwrite fits and coefs:
      fits.rich[[i]]<-m2way_rich
      coef.rich[[i]]<-coef.ext(m2way_rich)

      # run model without any interaction (to get p value for two-way interaction):
      
      form.noint<-paste(resp.thisrun,"~Treatment+DATE+reserve+(1|PLOT_ID)", sep="")
      
      rich_noint<-glmmadmb(as.formula(form.noint), family="nbinom", data=data.set)
      
      aov2way_rich<-anova(rich_noint,m2way_rich)
      p2way_rich<-round(aov2way_rich[2,"Pr(>Chi)"],4)
      gdf$rich_2wayP[i]<-p2way_rich      
    
      # Overwrite anova table:
      anova.rich[[i]]<-aov2way_rich
      
      # Update current model:
      m1_rich<-m2way_rich
      
      } # close if TWO-WAY richness
    
    # simplify INVSIMP if the three way is not significant:
    if(p3way_invsimp>0.05){
      
      # Overwrite fits and coefs:
      fits.invsimp[[i]]<-m2way_invsimp
      coef.invsimp[[i]]<-coef.ext(m2way_invsimp)
      
      # run model without any interaction (to get p value for two-way interaction):
      
      form.noint<-paste(resp.thisrun,"~Treatment+DATE+reserve+(1|PLOT_ID)", sep="")
      
      invsimp_noint<-glmmadmb(as.formula(form.noint), family="gamma", data=invsimp_sc[which(invsimp.thisrun>0),])
      
      aov2way_invsimp<-anova(invsimp_noint,m2way_invsimp)
      p2way_invsimp<-round(aov2way_invsimp[2,"Pr(>Chi)"],4)
      gdf$invsimp_2wayP[i]<-p2way_invsimp   
      
      # Overwrite anova table:
      anova.invsimp[[i]]<-aov2way_invsimp
      
      # Update current model:
      m1_invsimp<-m2way_invsimp
      
    } # close if TWO-WAY invsimp
    
    # generate model predictions:
    
    rich_pr<-pred(model = m1_rich,new.data=nd1,se.fit=T)
    invsimp_pr<-pred(model = m1_invsimp,new.data=nd1,se.fit=T)
    
    head(rich_pr,3)
    head(invsimp_pr,3)
    
    preds.rich[[i]]<-rich_pr
    preds.invsimp[[i]]<-invsimp_pr
    
  } # close positive models
  
} # close ANALYSIS loop

# save workspace:
# save.image("03_Workspaces/stjw_analysis_R1.RData")

anova.binom
fits.binom
coef.binom
preds.binom

anova.rich
fits.rich
coef.rich
preds.rich

anova.invsimp
fits.invsimp
coef.invsimp
preds.invsimp

head(gdf,4); dim(gdf)

# Binomial models:

# 18 native leg herb three way
# 19 exotic leg herb two way

gdf[which(gdf$bin_3wayP<0.05),]
gdf[which(gdf$bin_2wayP<0.05),]

head(preds.binom[[19]])
head(preds.binom[[18]])

# Richness

gdf[which(gdf$rich_3wayP<0.05),]
gdf[which(gdf$rich_2wayP<0.05),]

# 12 exotic annual herb THREE way
# 3 exotic TWO way
# 11 exotic herb TWO way

# NO richness trends, the only one was exotic, which came up in the two way results

gdf[which(gdf$rich_3wayP>0.05 & gdf$rich_3wayP<0.1),]
gdf[which(gdf$rich_2wayP>0.05 & gdf$rich_2wayP<0.1),]

# Inverse Simpsons

gdf[which(gdf$invsimp_3wayP<0.05),]
gdf[which(gdf$invsimp_2wayP<0.05),]
summary(fits.invsimp[[6]])
anova.invsimp[[6]]

# 6 sig B THREE way
# 12 exotic annual herb THREE WAY
# 19 exotic leg herb TWO way

# FIVE groups have TWO way TRENDS for inverse simpsons: native herb, exotic herb, native annual herb, c3 grass and native c3 (10, 11, 14, 24, 25):
gdf[which(gdf$invsimp_3wayP>0.05 & gdf$invsimp_3wayP<0.1),]
gdf[which(gdf$invsimp_2wayP>0.05 & gdf$invsimp_2wayP<0.1),]

# Summarise reserve effects:

# P values and coefficient for binomial models
# The coefficient for M is shown for reserve effects; if it's positive there is a higher probability of detecting that group at Mulanggari
head(gdf,3); dim(gdf)

gdf$bin_resP<-NA
gdf$bin_resM_coef<-NA
gdf$bin_resM_se<-NA

gdf$bin_resP<-ifelse(gdf$fit_bin=="yes",round(unlist(lapply(coef.binom,FUN=function(x) x[which(x$term=="reserveM"),"P"])),4),gdf$bin_resP)
gdf$bin_resM_coef<-ifelse(gdf$fit_bin=="yes",round(unlist(lapply(coef.binom,FUN=function(x) x[which(x$term=="reserveM"),"est"])),4),gdf$bin_resP)
gdf$bin_resM_se<-ifelse(gdf$fit_bin=="yes",round(unlist(lapply(coef.binom,FUN=function(x) x[which(x$term=="reserveM"),"se"])),4),gdf$bin_resP)

gdf[which(gdf$fit_bin=="yes"),]

# P values and coefficient for richness models (two-way models only)
# The P value and coefficient are not shown for models with three way interactions, because the reserve effect is captured in the three way term

gdf$rich2w_resP<-NA
gdf$rich2w_resM_coef<-NA
gdf$rich2w_resM_se<-NA

head(gdf,5)
gdf$rich2w_resP<-ifelse(gdf$rich_3wayP>0.05,round(unlist(lapply(coef.rich,FUN=function(x) x[which(x$term=="reserveM"),"P"])),4),gdf$rich2w_resP)
gdf$rich2w_resM_coef<-ifelse(gdf$rich_3wayP>0.05,round(unlist(lapply(coef.rich,FUN=function(x) x[which(x$term=="reserveM"),"est"])),4),gdf$rich2w_resP)
gdf$rich2w_resM_se<-ifelse(gdf$rich_3wayP>0.05,round(unlist(lapply(coef.rich,FUN=function(x) x[which(x$term=="reserveM"),"se"])),4),gdf$rich2w_resP)

# P values and coefficient for diversity models (two-way models only)
# The P value and coefficient are not shown for models with three way interactions, because the reserve effect is captured in the three way term

gdf$invsimp2w_resP<-NA
gdf$invsimp2w_resM_coef<-NA
gdf$invsimp2w_resM_se<-NA

gdf$invsimp2w_resP<-ifelse(gdf$invsimp_3wayP>0.05,round(unlist(lapply(coef.invsimp,FUN=function(x) x[which(x$term=="reserveM"),"P"])),4),gdf$invsimp2w_resP)
gdf$invsimp2w_resM_coef<-ifelse(gdf$invsimp_3wayP>0.05,round(unlist(lapply(coef.invsimp,FUN=function(x) x[which(x$term=="reserveM"),"est"])),4),gdf$invsimp2w_resP)
gdf$invsimp2w_resM_se<-ifelse(gdf$invsimp_3wayP>0.05,round(unlist(lapply(coef.invsimp,FUN=function(x) x[which(x$term=="reserveM"),"se"])),4),gdf$invsimp2w_resP)

# Model results are exactly the same after update, the only difference should be in the predictions:

# write.table(gdf, "div_sum.txt", sep="\t", quote=F, row.names = F)

# save.image("03_Workspaces/stjw_analysis_R1.RData")

# close richness & diversity ----

# If we need to re-do the contrasts to deal with the boundary issue, we could try either fitting the binomial models in glmmADMB or try doing the contrasts with multicomp
# https://cran.r-project.org/web/packages/multcomp/multcomp.pdf

#  CONTRASTS (Richness & Diversity):    	# ----

# BINOMIAL MODELS (two signif models)
# Native legume forb (index 18), three-way
# Exotic legume forb (index 19), two-way

head(gdf,2)
gdf[18:19,]
coef.binom[[18]]
coef.binom[[19]]
fits.binom[[18]]
fits.binom[[19]]

# For all models, there are three treatments (C, A, B), one date date variable with three years (0:3, linear numeric), and two reserves (J, M).

# For BOTH the two-way AND the three-way interaction, this will equal a total of 45 contrasts - i.e. 15 contrasts per year, bewteen all treatements.

# This means we can use the same difference matrix for all of the calculations. 

# We need two different uzm matrices for the three-way and two-way models. 

### Create a "unique model matrix" THREE-WAY:

int3way.z.mat<-lm(native_legherb~Treatment + DATE + reserve + Treatment:DATE +Treatment:DATE:reserve,data=rich_sc,x=T)$x
int3way.uzm<-unique(int3way.z.mat)
rownames(int3way.uzm)<-1:nrow(int3way.uzm)
int3way.uzm; dim(int3way.uzm)

# The natural order for this difference matrix is:
# C, A, B, yr 0, Jerra 
# C, A, B, yr 0, Mulang 
# C, A, B, yr 1, Jerra 
# C, A, B, yr 1, Mulang 
# C, A, B, yr 2, Jerra 
# C, A, B, yr 2, Mulang 

# So there is no need to re-order
head(int3way.uzm); dim(int3way.uzm)

### Create a "difference matrix"

# each row must be a vector with a length equal to the number of rows in the uzm matrix, i.e. 18; 45 contrasts=45 rows. Each row will specify ONE contrast. 

coef.binom[[18]]

int.z.diff<-rbind(
  # year 0 treatment differences:
  c(rep(0,0),-1,1,0,0,0,0,rep(0,12)),
  c(rep(0,0),-1,0,1,0,0,0,rep(0,12)),
  c(rep(0,0),-1,0,0,1,0,0,rep(0,12)),
  c(rep(0,0),-1,0,0,0,1,0,rep(0,12)),
  c(rep(0,0),-1,0,0,0,0,1,rep(0,12)),
  c(rep(0,0),0,-1,1,0,0,0,rep(0,12)),
  c(rep(0,0),0,-1,0,1,0,0,rep(0,12)),
  c(rep(0,0),0,-1,0,0,1,0,rep(0,12)),
  c(rep(0,0),0,-1,0,0,0,1,rep(0,12)),
  c(rep(0,0),0,0,-1,1,0,0,rep(0,12)),
  c(rep(0,0),0,0,-1,0,1,0,rep(0,12)),
  c(rep(0,0),0,0,-1,0,0,1,rep(0,12)),
  c(rep(0,0),0,0,0,-1,1,0,rep(0,12)),
  c(rep(0,0),0,0,0,-1,0,1,rep(0,12)),
  c(rep(0,0),0,0,0,0,-1,1,rep(0,12)),

  # year 1 treatment differences:
  c(rep(0,6),-1,1,0,0,0,0,rep(0,6)),
  c(rep(0,6),-1,0,1,0,0,0,rep(0,6)),
  c(rep(0,6),-1,0,0,1,0,0,rep(0,6)),
  c(rep(0,6),-1,0,0,0,1,0,rep(0,6)),
  c(rep(0,6),-1,0,0,0,0,1,rep(0,6)),
  c(rep(0,6),0,-1,1,0,0,0,rep(0,6)),
  c(rep(0,6),0,-1,0,1,0,0,rep(0,6)),
  c(rep(0,6),0,-1,0,0,1,0,rep(0,6)),
  c(rep(0,6),0,-1,0,0,0,1,rep(0,6)),
  c(rep(0,6),0,0,-1,1,0,0,rep(0,6)),
  c(rep(0,6),0,0,-1,0,1,0,rep(0,6)),
  c(rep(0,6),0,0,-1,0,0,1,rep(0,6)),
  c(rep(0,6),0,0,0,-1,1,0,rep(0,6)),
  c(rep(0,6),0,0,0,-1,0,1,rep(0,6)),
  c(rep(0,6),0,0,0,0,-1,1,rep(0,6)),
  
  # year 2 treatment differences:
  c(rep(0,12),-1,1,0,0,0,0,rep(0,0)),
  c(rep(0,12),-1,0,1,0,0,0,rep(0,0)),
  c(rep(0,12),-1,0,0,1,0,0,rep(0,0)),
  c(rep(0,12),-1,0,0,0,1,0,rep(0,0)),
  c(rep(0,12),-1,0,0,0,0,1,rep(0,0)),
  c(rep(0,12),0,-1,1,0,0,0,rep(0,0)),
  c(rep(0,12),0,-1,0,1,0,0,rep(0,0)),
  c(rep(0,12),0,-1,0,0,1,0,rep(0,0)),
  c(rep(0,12),0,-1,0,0,0,1,rep(0,0)),
  c(rep(0,12),0,0,-1,1,0,0,rep(0,0)),
  c(rep(0,12),0,0,-1,0,1,0,rep(0,0)),
  c(rep(0,12),0,0,-1,0,0,1,rep(0,0)),
  c(rep(0,12),0,0,0,-1,1,0,rep(0,0)),
  c(rep(0,12),0,0,0,-1,0,1,rep(0,0)),
  c(rep(0,12),0,0,0,0,-1,1,rep(0,0))
)

# Now we have our unique model matrix
int3way.uzm

# And our difference matrix:
int.z.diff

# Create names for the contrasts:

# This has to relate to how the uzm and the z.diff are set up, so make sure they are re-ordered in the same way
# C, A, B, yr 0, Jerra 
# C, A, B, yr 0, Mulang 
# C, A, B, yr 1, Jerra 
# C, A, B, yr 1, Mulang 
# C, A, B, yr 2, Jerra 
# C, A, B, yr 2, Mulang 

coef.binom[[18]]

# For each year, there are 15 contrasts
# Yr 0: CJer vs AJer
# Yr 0: CJer vs BJer
# Yr 0: CJer vs CMul
# Yr 0: CJer vs AMul
# Yr 0: CJer vs BMul

# Yr 0: AJer vs BJer
# Yr 0: AJer vs CMul
# Yr 0: AJer vs AMul
# Yr 0: AJer vs BMul

# Yr 0: BJer vs CMul
# Yr 0: BJer vs AMul
# Yr 0: BJer vs BMul

# Yr 0: CMul vs AMul
# Yr 0: CMul vs BMul
# Yr 0: AMul vs BMul

int.names<-data.frame(Year=c(rep(0,15),rep(1,15),rep(2,15)),Contrast=c("CJer vs AJer","CJer vs BJer","CJer vs CMul","CJer vs AMul","CJer vs BMul","AJer vs BJer","AJer vs CMul","AJer vs AMul","AJer vs BMul","BJer vs CMul","BJer vs AMul","BJer vs BMul","CMul vs AMul","CMul vs BMul","AMul vs BMul"),res=c("Jer:Jer","Jer:Jer","Jer:Mul","Jer:Mul","Jer:Mul","Jer:Jer","Jer:Mul","Jer:Mul","Jer:Mul","Jer:Mul","Jer:Mul","Jer:Mul","Mul:Mul","Mul:Mul","Mul:Mul"),trt=c("C:A","C:B","C:C","C:A","C:B","A:B","A:C","A:A","A:B","B:C","B:A","B:B","C:A","C:B","A:B"))

# Difference estimates, SE and CI:
# Native legume forb (index 18), three-way BINOMIAL

natlegforb.diff<-data.frame(diff.est(fits.binom[[18]],int3way.uzm, int.z.diff),int.names)
natlegforb.diff$diff<-ifelse(sign(natlegforb.diff$lci)==sign(natlegforb.diff$uci),1,0)
natlegforb.diff

# Compare differences estimates with inbuilt contrasts function in lsmeans:
natlegforb<-fits.binom[[18]]
summary(natlegforb)$coefficients
library(lsmeans)
xx<-emmeans(natlegforb, pairwise~Treatment:reserve, by="DATE", adjust="tukey")
xx

# Cannot figure out how to use lsmeans to get the per year contrasts. It's not subsetting by year because DATE is numeric; it estimates at the mean of the numeric variable (i.e. year 1). However, the main reason I wanted to use this was to explore whether my hand calculations were correct, so I can do this by comparing the year 1 contrasts from both methods:

# The following shows they are exactly the same - my method is correct:

# Year 1 contrasts from lsmeans:
xx$contrasts

# Year 1 contrasts from my method:
binom.nd1[16:30,c(1:4,length(binom.nd1))]

# For the two-way interaction, this gives 7 coef lines and 21 contrasts
choose(7,2)

### Create a "unique model matrix" TWO-WAY:
coef.binom[[19]]
summary(fits.binom[[19]])

int2way.z.mat<-lm(exotic_legherb~Treatment + DATE + reserve + Treatment:DATE ,data=rich_sc,x=T)$x
int2way.uzm<-unique(int2way.z.mat)
rownames(int2way.uzm)<-1:nrow(int2way.uzm)
int2way.uzm; dim(int2way.uzm)

# The natural order for this difference matrix is:
# C, A, B, yr 0, Jerra 
# C, A, B, yr 0, Mulang 
# C, A, B, yr 1, Jerra 
# C, A, B, yr 1, Mulang 
# C, A, B, yr 2, Jerra 
# C, A, B, yr 2, Mulang 

# So there is no need to re-order
head(int2way.uzm); dim(int2way.uzm)

# Now we have our unique model matrix
int2way.uzm

# Our difference matrix (from above):
int.z.diff

# And our names (from above):
int.names

# Difference estimates, SE and CI:
# Exotic legume forb (index 19), two-way (binomial)

exlegforb.diff<-data.frame(diff.est(fits.binom[[19]],int2way.uzm, int.z.diff),int.names)
exlegforb.diff$diff<-ifelse(sign(exlegforb.diff$lci)==sign(exlegforb.diff$uci),1,0)
exlegforb.diff

# RICHNESS CONTRASTS
coef.rich[[3]] # All exotic two-way
coef.rich[[11]] # Exotic forb two-way
coef.rich[[12]] # Exotic annual forb three-way

# All exotic RICHNESS two-way
allex.diff<-data.frame(diff.est(fits.rich[[3]],int2way.uzm, int.z.diff),int.names)
allex.diff$diff<-ifelse(sign(allex.diff$lci)==sign(allex.diff$uci),1,0)
allex.diff

# Exotic forb RICHNESS two-way
exforb.diff<-data.frame(diff.est(fits.rich[[11]],int2way.uzm, int.z.diff),int.names)
exforb.diff$diff<-ifelse(sign(exforb.diff$lci)==sign(exforb.diff$uci),1,0)
exforb.diff

# Exotic annual forb RICHNESS two-way
exannforb.diff<-data.frame(diff.est(fits.rich[[12]],int3way.uzm, int.z.diff),int.names)
exannforb.diff$diff<-ifelse(sign(exannforb.diff$lci)==sign(exannforb.diff$uci),1,0)
exannforb.diff

# DIVERSITY CONTRASTS
coef.invsimp[[6]] # Sig B three-way
coef.invsimp[[12]] # Exotic annual forb three-way
coef.invsimp[[19]] # Exotic leg. forb two-way

# Significance B inv.simp three-way
sigB.diff<-data.frame(diff.est(fits.invsimp[[6]],int3way.uzm, int.z.diff),int.names)
sigB.diff$diff<-ifelse(sign(sigB.diff$lci)==sign(sigB.diff$uci),1,0)
sigB.diff

# Exotic annual forb inv.simp three-way
exanfbDIV.diff<-data.frame(diff.est(fits.invsimp[[12]],int3way.uzm, int.z.diff),int.names)
exanfbDIV.diff$diff<-ifelse(sign(exanfbDIV.diff$lci)==sign(exanfbDIV.diff$uci),1,0)
exanfbDIV.diff

# Exotic legume forb inv.simp two-way
exlegfbDIV.diff<-data.frame(diff.est(fits.invsimp[[19]],int2way.uzm, int.z.diff),int.names)
exlegfbDIV.diff$diff<-ifelse(sign(exlegfbDIV.diff$lci)==sign(exlegfbDIV.diff$uci),1,0)
exlegfbDIV.diff

# Reduce to within reserve contrasts:
# binomial
nlf.diff<-natlegforb.diff[-which(natlegforb.diff$res=="Jer:Mul"),]
exlf.diff<-exlegforb.diff[-which(exlegforb.diff$res=="Jer:Mul"),]
# richness
allex.diff<-allex.diff[-which(allex.diff$res=="Jer:Mul"),]
exfb.diff<-exforb.diff[-which(exforb.diff$res=="Jer:Mul"),]
exaf.diff<-exannforb.diff[-which(exannforb.diff$res=="Jer:Mul"),]
# inv.simp
sigB.diff2<-sigB.diff[-which(sigB.diff$res=="Jer:Mul"),]
exanfbDIV.diff2<-exanfbDIV.diff[-which(exanfbDIV.diff$res=="Jer:Mul"),]
exlegfbDIV.diff2<-exlegfbDIV.diff[-which(exlegfbDIV.diff$res=="Jer:Mul"),]


nlf.diff<-tidy.df(nlf.diff)
exlf.diff<-tidy.df(exlf.diff)

allex.diff<-tidy.df(allex.diff)
exfb.diff<-tidy.df(exfb.diff)
exaf.diff<-tidy.df(exaf.diff)

sigB.diff2<-tidy.df(sigB.diff2)
exanfbDIV.diff2<-tidy.df(exanfbDIV.diff2)
exlegfbDIV.diff2<-tidy.df(exlegfbDIV.diff2)

# close contrasts ----

#  PLOT ESTIMATES (Richness & Diversity):    	# ----

# Plot significant binomial models:

# update y labels for manuscript:

head(gdf,2)
gdf$ylabn<-gdf$ylab
gdf$ylabn[gdf$group=="native_legherb"]<-"Native legume forb"
gdf$ylabn[gdf$group=="exotic_legherb"]<-"Exotic legume forb"

# Contrasts
# The only significant contrast for the native leg forb model is for C:A at Mul in yr 2
nlf.diff
coef.binom[[18]]

# HOWEVER, boundary issues have caused a problems computing confidence intervals because many of the plots had complete cover of native leg herbs; therefore it would be worth indicating this on the plot

# nat leg herb
nlh_dat<-rich_sc[,c(1:4,which(colnames(rich_sc)=="native_legherb"))]
nlh_dat$nlh<-nlh_dat$native_legherb
nlh_dat$nlh<- ifelse(nlh_dat$nlh>0,1,0)
head(nlh_dat)
table(nlh_dat$nlh>0, nlh_dat$reserve,nlh_dat$DATE,nlh_dat)
nlh_sum<-aggregate(nlh~Treatment:DATE:reserve, data=nlh_dat,FUN=function(x) length(which(x>0)))
nlh_sum$prop<-nlh_sum$nlh/8

panel.size<-2

dev.new(width=panel.size*3.5,height=panel.size*3,noRStudioGD = T,dpi=80, pointsize=(panel.size*3.5)*2)
par(mfrow=c(2,2), mar=c(4,4,3,2), oma=c(0,0,0,6), mgp=c(2.5,1,0))

for(i in 18:19){
  
  resp.thisrun<-gdf$group[i]
  pred.thisrun<-preds.binom[[i]]
  anova.thisrun<-anova.binom[[i]]
  coef.thisrun<-coef.binom[[i]]
  ylab.thisrun<-gdf$ylabn[i]
  meta.thisrun<-gdf[i,]
  xofs<-0.2
  arrowlgth<-0.02
  if(i==18) codes.thisrun<-letters[1:2] else codes.thisrun<-letters[3:4]
  
  head(pred.thisrun)
  
  mul.preds<-pred.thisrun[which(pred.thisrun$reserve=="M"),]
  jerra.preds<-pred.thisrun[which(pred.thisrun$reserve=="J"),]
  
  # CIs for control in 2017 are all 1, and cannot draw the arrow. It's actually running, just sending out a warning that it doesn't plot the very small ones. 
  
  # Mulanggari
  
  plot(mul.preds$DATE[mul.preds$Treatment=="C"]-xofs,mul.preds$fit[mul.preds$Treatment=="C"], pch=15, ylim=c(min(mul.preds$lci), max(mul.preds$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=ylab.thisrun, las=1)
  axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
  
  arrows(mul.preds$DATE[mul.preds$Treatment=="C"]-xofs,mul.preds$lci[mul.preds$Treatment=="C"],mul.preds$DATE[mul.preds$Treatment=="C"]-xofs,mul.preds$uci[mul.preds$Treatment=="C"], code=3, angle=90, length=arrowlgth)
  
  points(mul.preds$DATE[mul.preds$Treatment=="A"],mul.preds$fit[mul.preds$Treatment=="A"], pch=15, col="red")
  arrows(mul.preds$DATE[mul.preds$Treatment=="A"],mul.preds$lci[mul.preds$Treatment=="A"],mul.preds$DATE[mul.preds$Treatment=="A"],mul.preds$uci[mul.preds$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
  points(mul.preds$DATE[mul.preds$Treatment=="B"]+xofs,mul.preds$fit[mul.preds$Treatment=="B"], pch=15, col="blue")
  arrows(mul.preds$DATE[mul.preds$Treatment=="B"]+xofs,mul.preds$lci[mul.preds$Treatment=="B"],mul.preds$DATE[mul.preds$Treatment=="B"]+xofs,mul.preds$uci[mul.preds$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
  mtext(paste("(",codes.thisrun[1],") ", "Mulanggari",sep=""), side=3, line=1.2, adj=0, cex=0.9)
  
  # Add star for significant differences
  nlf.diff
  exlf.diff

  # A is different to control in yr 1 and yr 2
  if (i==18) points(mul.preds$DATE[mul.preds$Treatment=="A"][c(3)],mul.preds$uci[mul.preds$Treatment=="A"][c(3)]+0.1, pch="*", cex=1.5,col="black")
  if (i==18) points((mul.preds$DATE[mul.preds$Treatment=="C"]-xofs)[c(2)],(mul.preds$lci[mul.preds$Treatment=="C"])[c(2)]-0.15, pch="†", cex=1,col="black")
  
  if (i==19) points(mul.preds$DATE[mul.preds$Treatment=="A"][c(2,3)],mul.preds$uci[mul.preds$Treatment=="A"][c(2,3)]+0.1, pch="*", cex=1.5,col="black")
  
  # JERRA
  
  plot(jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$fit[jerra.preds$Treatment=="C"], pch=15, ylim=c(min(jerra.preds$lci), max(jerra.preds$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=ylab.thisrun, las=1)
  axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
  arrows(jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$lci[jerra.preds$Treatment=="C"],jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$uci[jerra.preds$Treatment=="C"], code=3, angle=90, length=arrowlgth)
  points(jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$fit[jerra.preds$Treatment=="A"], pch=15, col="red")
  arrows(jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$lci[jerra.preds$Treatment=="A"],jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$uci[jerra.preds$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
  points(jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$fit[jerra.preds$Treatment=="B"], pch=15, col="blue")
  arrows(jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$lci[jerra.preds$Treatment=="B"],jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$uci[jerra.preds$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
  
  mtext(paste("(",codes.thisrun[2],") ", "Jerrabomberra", sep=""), side=3, line=1.2, adj=0, cex=0.9)
  
  # Add star for only significant difference
  if (i==18) points((jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs)[c(2,3)],(jerra.preds$uci[jerra.preds$Treatment=="C"])[c(2,3)]-0.15, pch="†", cex=1,col="black")
  
  exlf.diff
  if (i==19) points(jerra.preds$DATE[jerra.preds$Treatment=="A"][c(2,3)],jerra.preds$uci[jerra.preds$Treatment=="A"][c(2,2)]+0.1, pch="*", cex=1.5,col="black")
  
  if (gdf$bin_3wayP[i]<0.05){
    
    p.trt_yr_int<-round(anova.thisrun[2,which(colnames(anova.thisrun)=="Pr(>Chisq)")],3)
    title(main=bquote(Three-way~int.~italic(P)~"="~.(p.trt_yr_int)), font.main=1, adj=0, cex.main=1, line=0.5)
    
  } # close if 3-way signif
  
  if(is.na(gdf$bin_2wayP[i])) next
  
  if (gdf$bin_2wayP[i]<0.05){
    
    p.trt_yr_int<-round(anova.thisrun[2,which(colnames(anova.thisrun)=="Pr(>Chisq)")],3)
    title(main=bquote(Two-way~int.~italic(P)~"="~.(p.trt_yr_int)), font.main=1.4, adj=0, cex.main=1, line=0.5)
    
  } # close if 2-way signif
  
} # close plot binomial

par(xpd=NA)
legend(2.6,0.6,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)
par(xpd=F)

# Plot significant richness models:

# update y labels for manuscript:

head(gdf,2)
gdf$ylabn<-gdf$ylab
gdf$ylabn[gdf$group=="exotic"]<-"All exotic plants"

allex.diff
exfb.diff
exaf.diff

dev.new(width=panel.size*3.5,height=panel.size*4,noRStudioGD = T,dpi=80, pointsize=(panel.size*4)*2)
par(mfrow=c(3,2), mar=c(4,4,3,2), oma=c(0,0,0,6), mgp=c(2.5,1,0))

for(i in c(3,11,12)){
  
  resp.thisrun<-gdf$group[i]
  pred.thisrun<-preds.rich[[i]]
  anova.thisrun<-anova.rich[[i]]
  coef.thisrun<-coef.rich[[i]]
  ylab.thisrun<-gdf$ylabn[i]
  meta.thisrun<-gdf[i,]
  xofs<-0.2
  arrowlgth<-0.02
  
  if(i==3) codes.thisrun<-letters[1:2] 
  if(i==11) codes.thisrun<-letters[3:4] 
  if(i==12) codes.thisrun<-letters[5:6] 
  
  head(pred.thisrun)
  
  mul.preds<-pred.thisrun[which(pred.thisrun$reserve=="M"),]
  jerra.preds<-pred.thisrun[which(pred.thisrun$reserve=="J"),]
  
  # Mulanggari
  
  plot(mul.preds$DATE[mul.preds$Treatment=="C"]-xofs,mul.preds$fit.resp[mul.preds$Treatment=="C"], pch=15, ylim=c(min(mul.preds$lci.resp), max(mul.preds$uci.resp)+1), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=ylab.thisrun, las=1)
  axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
  
  arrows(mul.preds$DATE[mul.preds$Treatment=="C"]-xofs,mul.preds$lci.resp[mul.preds$Treatment=="C"],mul.preds$DATE[mul.preds$Treatment=="C"]-xofs,mul.preds$uci.resp[mul.preds$Treatment=="C"], code=3, angle=90, length=arrowlgth)
  
  points(mul.preds$DATE[mul.preds$Treatment=="A"],mul.preds$fit.resp[mul.preds$Treatment=="A"], pch=15, col="red")
  arrows(mul.preds$DATE[mul.preds$Treatment=="A"],mul.preds$lci.resp[mul.preds$Treatment=="A"],mul.preds$DATE[mul.preds$Treatment=="A"],mul.preds$uci.resp[mul.preds$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
  points(mul.preds$DATE[mul.preds$Treatment=="B"]+xofs,mul.preds$fit.resp[mul.preds$Treatment=="B"], pch=15, col="blue")
  arrows(mul.preds$DATE[mul.preds$Treatment=="B"]+xofs,mul.preds$lci.resp[mul.preds$Treatment=="B"],mul.preds$DATE[mul.preds$Treatment=="B"]+xofs,mul.preds$uci.resp[mul.preds$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
  mtext(paste("(",codes.thisrun[1],") ", "Mulanggari",sep=""), side=3, line=1.2, adj=0, cex=0.7)
  
  # Add star for significant differences
  allex.diff[allex.diff$res=="Mul:Mul",]  # [[3]] # All exotic two-way
  exfb.diff[allex.diff$res=="Mul:Mul",] # [[11]] # Exotic forb two-way
  exaf.diff[allex.diff$res=="Mul:Mul",] # [[12]] # Exotic annual forb three-way
  
  # A is different to control in yr 1 and yr 2
  if (i==3) points(mul.preds$DATE[mul.preds$Treatment=="A"][c(2,3)],mul.preds$uci.resp[mul.preds$Treatment=="A"][c(2,3)]+1, pch="*", cex=1.5,col="black")
  
  if (i==11) points(mul.preds$DATE[mul.preds$Treatment=="A"][c(3)],mul.preds$uci.resp[mul.preds$Treatment=="A"][c(3)]+1, pch="*", cex=1.5,col="black")
  
  # JERRA
  
  if(i==12) jerra.maxlim<-max(jerra.preds$uci.resp)+1 else jerra.maxlim<-max(jerra.preds$uci.resp) 
  plot(jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$fit.resp[jerra.preds$Treatment=="C"], pch=15, ylim=c(min(jerra.preds$lci.resp), jerra.maxlim), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=ylab.thisrun, las=1)
  axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
  arrows(jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$lci.resp[jerra.preds$Treatment=="C"],jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$uci.resp[jerra.preds$Treatment=="C"], code=3, angle=90, length=arrowlgth)
  points(jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$fit.resp[jerra.preds$Treatment=="A"], pch=15, col="red")
  arrows(jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$lci.resp[jerra.preds$Treatment=="A"],jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$uci.resp[jerra.preds$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
  points(jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$fit.resp[jerra.preds$Treatment=="B"], pch=15, col="blue")
  arrows(jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$lci.resp[jerra.preds$Treatment=="B"],jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$uci.resp[jerra.preds$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
  
  mtext(paste("(",codes.thisrun[2],") ", "Jerrabomberra", sep=""), side=3, line=1.2, adj=0, cex=0.7)
  
  # Add star for significant differences
  allex.diff[allex.diff$res=="Jer:Jer",] # [[3]] # All exotic two-way
  exfb.diff[allex.diff$res=="Jer:Jer",] # [[11]] # Exotic forb two-way
  exaf.diff[allex.diff$res=="Jer:Jer",] # [[12]] # Exotic annual forb three-way
  
  # A is different to control in yr 1 and yr 2
  if (i==3) points(jerra.preds$DATE[jerra.preds$Treatment=="A"][c(2,3)],jerra.preds$uci.resp[jerra.preds$Treatment=="A"][c(2,3)]+1, pch="*", cex=1.5,col="black")
  
  if (i==11) points(jerra.preds$DATE[jerra.preds$Treatment=="A"][c(3)],jerra.preds$uci.resp[jerra.preds$Treatment=="A"][c(3)]+1, pch="*", cex=1.5,col="black")
  
  if (i==12) points(jerra.preds$DATE[jerra.preds$Treatment=="A"][c(2,3)],jerra.preds$uci.resp[jerra.preds$Treatment=="A"][c(2,3)]+1, pch="*", cex=1.5,col="black")
  
  if (gdf$rich_3wayP[i]<0.05){
    
    p.trt_yr_int<-round(anova.thisrun[2,which(colnames(anova.thisrun)=="Pr(>Chi)")],3)
    title(main=bquote(Three-way~int.~italic(P)~"="~.(p.trt_yr_int)), font.main=1, adj=0, cex.main=1, line=0.5)
    
  } # close if 3-way signif
  
  if(is.na(gdf$rich_2wayP[i])) next
  
  if (gdf$rich_2wayP[i]<0.05){
    
    p.trt_yr_int<-round(anova.thisrun[2,which(colnames(anova.thisrun)=="Pr(>Chi)")],3)
    if (p.trt_yr_int<0.001) title(main=bquote(Two-way~int.~italic(P)~"<"~0.001), font.main=1, adj=0, cex.main=1, line=0.5) else title(main=bquote(Two-way~int.~italic(P)~"="~.(p.trt_yr_int)), font.main=1.4, adj=0, cex.main=1, line=0.5)
    
  } # close if 2-way signif
  
} # close plot richness

par(xpd=NA)
legend(2.6,3,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)
par(xpd=F)

# Plot significant invsimp models:

# update y labels for manuscript:

head(gdf,2)
gdf$ylabn<-gdf$ylab
gdf$ylabn[gdf$group=="sigB"]<-"Indicator B"
gdf$ylabn[gdf$group=="exotic_legherb"]<-"Exotic legume forb"

# For sigB and exotic leg forb, none of the within year, within reserve contrasts were significant. The interactions were driven by differences between reserves
sigB.diff2 # 6
exanfbDIV.diff2 # 12
exlegfbDIV.diff2 # 19

sigB.diff # e.g. CMul vs CJer
exlegfbDIV.diff # e.g. CMul vs AJer 

dev.new(width=panel.size*3.5,height=panel.size*4,noRStudioGD = T,dpi=80, pointsize=(panel.size*4)*2)
par(mfrow=c(3,2), mar=c(4,4,3,2), oma=c(0,0,0,6), mgp=c(2.5,1,0))

for(i in c(6,12,19)){
  
  resp.thisrun<-gdf$group[i]
  pred.thisrun<-preds.invsimp[[i]]
  anova.thisrun<-anova.invsimp[[i]]
  coef.thisrun<-coef.invsimp[[i]]
  ylab.thisrun<-gdf$ylabn[i]
  meta.thisrun<-gdf[i,]
  xofs<-0.2
  arrowlgth<-0.02
  
  if(i==6) codes.thisrun<-letters[1:2] 
  if(i==12) codes.thisrun<-letters[3:4] 
  if(i==19) codes.thisrun<-letters[5:6] 
  
  head(pred.thisrun)
  
  mul.preds<-pred.thisrun[which(pred.thisrun$reserve=="M"),]
  jerra.preds<-pred.thisrun[which(pred.thisrun$reserve=="J"),]
  
  # Mulanggari
  
  plot(mul.preds$DATE[mul.preds$Treatment=="C"]-xofs,mul.preds$fit.resp[mul.preds$Treatment=="C"], pch=15, ylim=c(min(mul.preds$lci.resp), max(mul.preds$uci.resp)), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=ylab.thisrun, las=1)
  axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
  
  arrows(mul.preds$DATE[mul.preds$Treatment=="C"]-xofs,mul.preds$lci.resp[mul.preds$Treatment=="C"],mul.preds$DATE[mul.preds$Treatment=="C"]-xofs,mul.preds$uci.resp[mul.preds$Treatment=="C"], code=3, angle=90, length=arrowlgth)
  
  points(mul.preds$DATE[mul.preds$Treatment=="A"],mul.preds$fit.resp[mul.preds$Treatment=="A"], pch=15, col="red")
  arrows(mul.preds$DATE[mul.preds$Treatment=="A"],mul.preds$lci.resp[mul.preds$Treatment=="A"],mul.preds$DATE[mul.preds$Treatment=="A"],mul.preds$uci.resp[mul.preds$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
  points(mul.preds$DATE[mul.preds$Treatment=="B"]+xofs,mul.preds$fit.resp[mul.preds$Treatment=="B"], pch=15, col="blue")
  arrows(mul.preds$DATE[mul.preds$Treatment=="B"]+xofs,mul.preds$lci.resp[mul.preds$Treatment=="B"],mul.preds$DATE[mul.preds$Treatment=="B"]+xofs,mul.preds$uci.resp[mul.preds$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
  mtext(paste("(",codes.thisrun[1],") ", "Mulanggari",sep=""), side=3, line=1.2, adj=0, cex=0.7)
  
  # JERRA
  
  if(i==12) jerra.maxlim<-max(jerra.preds$uci.resp)+1 else jerra.maxlim<-max(jerra.preds$uci.resp) 
  plot(jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$fit.resp[jerra.preds$Treatment=="C"], pch=15, ylim=c(min(jerra.preds$lci.resp), jerra.maxlim), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=ylab.thisrun, las=1)
  axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
  arrows(jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$lci.resp[jerra.preds$Treatment=="C"],jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$uci.resp[jerra.preds$Treatment=="C"], code=3, angle=90, length=arrowlgth)
  points(jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$fit.resp[jerra.preds$Treatment=="A"], pch=15, col="red")
  arrows(jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$lci.resp[jerra.preds$Treatment=="A"],jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$uci.resp[jerra.preds$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
  points(jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$fit.resp[jerra.preds$Treatment=="B"], pch=15, col="blue")
  arrows(jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$lci.resp[jerra.preds$Treatment=="B"],jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$uci.resp[jerra.preds$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
  
  mtext(paste("(",codes.thisrun[2],") ", "Jerrabomberra", sep=""), side=3, line=1.2, adj=0, cex=0.7)
  
  exanfbDIV.diff2 
  
  if (i==12) points(jerra.preds$DATE[jerra.preds$Treatment=="A"][c(3)],jerra.preds$uci.resp[jerra.preds$Treatment=="A"][c(3)]+1, pch="*", cex=1.5,col="black")
  if (i==12) points((jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs)[c(3)],jerra.preds$uci.resp[jerra.preds$Treatment=="B"][c(3)]+1, pch="*", cex=1.5,col="black")
  
  if (gdf$invsimp_3wayP[i]<0.05){
    
    p.trt_yr_int<-round(anova.thisrun[2,which(colnames(anova.thisrun)=="Pr(>Chi)")],3)
    title(main=bquote(Three-way~int.~italic(P)~"="~.(p.trt_yr_int)), font.main=1, adj=0, cex.main=1, line=0.5)
    
  } # close if 3-way signif
  
  if(is.na(gdf$invsimp_2wayP[i])) next
  
  if (gdf$invsimp_2wayP[i]<0.05){
    
    p.trt_yr_int<-round(anova.thisrun[2,which(colnames(anova.thisrun)=="Pr(>Chi)")],3)
    if (p.trt_yr_int<0.001) title(main=bquote(Two-way~int.~italic(P)~"<"~0.001), font.main=1, adj=0, cex.main=1, line=0.5) else title(main=bquote(Two-way~int.~italic(P)~"="~.(p.trt_yr_int)), font.main=1.4, adj=0, cex.main=1, line=0.5)
    
  } # close if 2-way signif
  
} # close plot invsimp

par(xpd=NA)
legend(2.6,1.5,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)
par(xpd=F)

# close plot Richness & Diversity ----

# For individual species:
# load("03_Workspaces/stjw_analysis.RData")

#  INDIVIDUAL SPECIES data set-up:    	# ----

# Need to check with RM: why was Desmodium varians classified as a legume rather than a native leg forb?
# Is it OK for us to classify it as native leg forb (I have done this above in the 'format data' section)?
# PlantNet classifies it as a "Prostrate or climbing herb"
# And it is the ONLY species in the Legume category
table(pinfo$func_grp)

# CHECK DATA

# Make sure all species data columns match:

rahead(ct_dat,3,7); dim(ct_dat)
rahead(cv_dat,3,7); dim(cv_dat)

# The species data columns in cover and count are not in the same order (but we can work around this):
table(colnames(ct_dat)[5:ncol(ct_dat)]==colnames(cv_dat)[5:ncol(cv_dat)])

# The same species are present in both data sets:
table(colnames(ct_dat)[5:ncol(ct_dat)] %in% colnames(cv_dat)[5:ncol(cv_dat)])
table(colnames(cv_dat)[5:ncol(cv_dat)] %in% colnames(ct_dat)[5:ncol(ct_dat)])

# Make sure that the cover and count data have data in the same plots:

sp.totest<-colnames(ct_dat)[5:ncol(ct_dat)]

sp.testout<-data.frame(Sp=sp.totest, lineup=NA, mm17="-", mm18="-", mm19="-")

for (i in 1:length(sp.totest)){
  
  sp.thisrun<-sp.totest[i]
  ct.thisrun<-ct_dat[,sp.thisrun]
  cv.thisrun<-cv_dat[,sp.thisrun]
  
  test.thisrun<-table((ct.thisrun>0)==(cv.thisrun>0))
  res.thisrun<-unlist(attr(test.thisrun,"dimnames"))
  
  # If false is in the result, it means that some of them are lining up. So if this if statement is true, they are NOT lining up. If the if statement is false, they ARE lining up:
  if("FALSE" %in% res.thisrun) sp.testout$lineup[i]<-"no" else sp.testout$lineup[i]<-"yes"
  
  if("FALSE" %in% res.thisrun){
    
    test.dat<-ct_dat[,c(1:4,which(colnames(ct_dat)==sp.thisrun))]
    colnames(test.dat)[ncol(test.dat)]<-"count_data"
    test.dat<-cbind(test.dat, cv_dat[,which(colnames(cv_dat)==sp.thisrun)])
    colnames(test.dat)[ncol(test.dat)]<-"cover_data"
    head(test.dat)
    
    test.dat<-test.dat[which(!(test.dat$count_data>0)==(test.dat$cover_data>0)),]
    dates.thisrun<-unique(test.dat$DATE)
    if("2017" %in% dates.thisrun) sp.testout$mm17[i]<-"yes"
    if("2018" %in% dates.thisrun) sp.testout$mm18[i]<-"yes"
    if("2019" %in% dates.thisrun) sp.testout$mm19[i]<-"yes"
    
  } # close if
  
} # close for

# NO more problems in re-formatted data. 
sp.testout

# Individual species to model:
head(pinfo,3); dim(pinfo)

# In total, 14 native forb species were used to define the study plots:
ind.sp<-data.frame(ind_species=c("Eryngium ovinum", "Chrysocephalum apiculatum", "Arthropodium fimbriatum", "Wurmbea dioica", "Desmodium varians", "Plantago varia", "Tricoryne elatior", "Triptilodiscus pygmaeus","Lomandra filiformis coriacea", "Lomandra bracteata", "Lomandra filiformis","Lomandra multiflora","Glycine clandestina","Glycine tabacina"))

# Typically how many of each of the 14 species were present in the study plots at the start of the experiment?
rahead(ct_dat,3,7); dim(ct_dat)
unique(ct_dat$DATE)
indsp
head(ind.sp,2)

ct_exp<-ct_dat[ct_dat$DATE==2017,c(1:4,which(colnames(ct_dat) %in% indsp))]
rahead(ct_exp,3,7); dim(ct_exp)
colSums(ct_exp[,5:length(ct_exp)])
range(apply(ct_exp[,5:length(ct_exp)], 1, function(x) length(which(x>0))))
summary(rowSums(ct_exp[,5:length(ct_exp)]))
# colSums(ct_dat[,5:length(ct_dat)])

# What if we only include the five species from RM's methods doc (that's also excluding Gly_cla which had only 7 records at the start of the study (or the 8 species that we  modelled):
# Or even adding Gly_cla back in?
indsp5<-c("Chr_api","Ery_ovi", "Tri_ela", "Des_var", "Gly_tab")
indsp6<-c("Chr_api","Ery_ovi", "Tri_ela", "Des_var", "Gly_tab","Gly_cla")
indsp8<-c(indsp5,"Lom_cor", "Pla_var" , "Tri_pyg")

ct_exp5<-ct_dat[ct_dat$DATE==2017,c(1:4,which(colnames(ct_dat) %in% indsp5))]
rahead(ct_exp5,3,7); dim(ct_exp5)
range(apply(ct_exp5[,5:length(ct_exp5)], 1, function(x) length(which(x>0))))
summary(apply(ct_exp5[,5:length(ct_exp5)], 1, function(x) length(which(x>0))))

# Adding Gly_cla back in doesn't change the minimum number:
ct_exp6<-ct_dat[ct_dat$DATE==2017,c(1:4,which(colnames(ct_dat) %in% indsp6))]
rahead(ct_exp6,3,7); dim(ct_exp6)
range(apply(ct_exp6[,5:length(ct_exp6)], 1, function(x) length(which(x>0))))
summary(range(apply(ct_exp6[,5:length(ct_exp6)], 1, function(x) length(which(x>0)))))

ct_exp8<-ct_dat[ct_dat$DATE==2017,c(1:4,which(colnames(ct_dat) %in% indsp8))]
rahead(ct_exp8,3,7); dim(ct_exp8)
range(apply(ct_exp8[,5:length(ct_exp8)], 1, function(x) length(which(x>0))))

# These should all be TRUE:
table(ind.sp$ind_species %in% pinfo$Species)

ind.sp<-merge(ind.sp, pinfo, by.x="ind_species", by.y="Species")
ind.sp
indsp<-ind.sp$Sp

# Add proportion of zeros column
# exclude species that have more than 95% zeros
ind.sp$prop0<-NA
ind.sp$prop0[which(ind.sp$Sp=="Chr_api")]<-table(ind_po$Chr_api)[1]/sum(table(ind_po$Chr_api))
ind.sp$prop0[which(ind.sp$Sp=="Gly_cla")]<-table(ind_po$Gly_cla)[1]/sum(table(ind_po$Gly_cla))
ind.sp$prop0[which(ind.sp$Sp=="Art_fim")]<-table(ind_po$Art_fim)[1]/sum(table(ind_po$Art_fim))
ind.sp$prop0[which(ind.sp$Sp=="Des_var")]<-table(ind_po$Des_var)[1]/sum(table(ind_po$Des_var))
ind.sp$prop0[which(ind.sp$Sp=="Ery_ovi")]<-table(ind_po$Ery_ovi)[1]/sum(table(ind_po$Ery_ovi))
ind.sp$prop0[which(ind.sp$Sp=="Gly_tab")]<-table(ind_po$Gly_tab)[1]/sum(table(ind_po$Gly_tab))
ind.sp$prop0[which(ind.sp$Sp=="Lom_bra")]<-table(ind_po$Lom_bra)[1]/sum(table(ind_po$Lom_bra))
ind.sp$prop0[which(ind.sp$Sp=="Lom_fil")]<-table(ind_po$Lom_fil)[1]/sum(table(ind_po$Lom_fil))
ind.sp$prop0[which(ind.sp$Sp=="Lom_cor")]<-table(ind_po$Lom_cor)[1]/sum(table(ind_po$Lom_cor))
ind.sp$prop0[which(ind.sp$Sp=="Lom_mul")]<-table(ind_po$Lom_mul)[1]/sum(table(ind_po$Lom_mul))
ind.sp$prop0[which(ind.sp$Sp=="Pla_var")]<-table(ind_po$Pla_var)[1]/sum(table(ind_po$Pla_var))
ind.sp$prop0[which(ind.sp$Sp=="Tri_ela")]<-table(ind_po$Tri_ela)[1]/sum(table(ind_po$Tri_ela))
ind.sp$prop0[which(ind.sp$Sp=="Tri_pyg")]<-table(ind_po$Tri_pyg)[1]/sum(table(ind_po$Tri_pyg))
ind.sp$prop0[which(ind.sp$Sp=="Wur_dio")]<-table(ind_po$Wur_dio)[1]/sum(table(ind_po$Wur_dio))

# Add count column
ind.sp$count<-NA
ind.sp$count[which(ind.sp$Sp=="Chr_api")]<-sum(ct_ind$Chr_api)
ind.sp$count[which(ind.sp$Sp=="Wur_dio")]<-sum(ct_ind$Wur_dio)
ind.sp$count[which(ind.sp$Sp=="Art_fim")]<-sum(ct_ind$Art_fim)
ind.sp$count[which(ind.sp$Sp=="Des_var")]<-sum(ct_ind$Des_var)
ind.sp$count[which(ind.sp$Sp=="Ery_ovi")]<-sum(ct_ind$Ery_ovi)
ind.sp$count[which(ind.sp$Sp=="Gly_cla")]<-sum(ct_ind$Gly_cla) 
ind.sp$count[which(ind.sp$Sp=="Gly_tab")]<-sum(ct_ind$Gly_tab)
ind.sp$count[which(ind.sp$Sp=="Lom_bra")]<-sum(ct_ind$Lom_bra)
ind.sp$count[which(ind.sp$Sp=="Lom_fil")]<-sum(ct_ind$Lom_fil)
ind.sp$count[which(ind.sp$Sp=="Lom_cor")]<-sum(ct_ind$Lom_cor)
ind.sp$count[which(ind.sp$Sp=="Lom_mul")]<-sum(ct_ind$Lom_mul)
ind.sp$count[which(ind.sp$Sp=="Pla_var")]<-sum(ct_ind$Pla_var)
ind.sp$count[which(ind.sp$Sp=="Tri_ela")]<-sum(ct_ind$Tri_ela)
ind.sp$count[which(ind.sp$Sp=="Tri_pyg")]<-sum(ct_ind$Tri_pyg)

# Establish minimum data thresholds for analysis

# For binomial models, exclude sp. that have more than 90 % zeros Art_fim; Gly_cla; Lom_bra; Lom_fil; Lom_mul; Wur_dio

# The data for Art_fim, although abundant, are spatially clustered and not related to treatment (i.e. all treatments in M1 in each year); thus, we only have a single quadrat and couldn't model that

head(ind.sp)

ind.sp$fit_binom<-"yes"
ind.sp$fit_binom[which(ind.sp$prop0> 0.9)]<-"no"
ind.sp<-tidy.df(ind.sp)
head(ind.sp); dim(ind.sp)

# close individual species data ----

# DATA TO MODEL (ind_po) ----

ct_ind<-ct_dat[,c(1:4, which(colnames(ct_dat) %in% indsp))]
cv_ind<-cv_dat[,c(1:4, which(colnames(cv_dat) %in% indsp))]

# Make sure columns match:
table(colnames(ct_ind)==colnames(cv_ind))

# Model individual species using count data, as presence/absence:
rahead(ct_ind,6,6); dim(ct_ind)
rahead(cv_ind,6,6)

# Format date for modelling:
ind_po<-ct_ind
ind_po$DATE<-ind_po$DATE-min(ind_po$DATE)
rahead(ind_po,6,6); dim(ind_po)

# close data to model ----

# This is the workspace with the individual species:
# save.image("03_Workspaces/stjw_analysis.RData")

#  INDIVIDUAL SPECIES BINOMIAL MODELS:    	# ----

# Fit binomial models for eight species:
ind.sp$Sp[which(ind.sp$fit_binom=="yes")]

# Chr_api 
ind_po$Chr_api<-ifelse(ind_po$Chr_api>0,1,0)
three_way_Chr_api<-glmer(Chr_api~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(three_way_Chr_api)

two_way_Chr_api<-glmer(Chr_api~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID), family="binomial", data=ind_po)
noint_Chr_api<-glmer(Chr_api~DATE+Treatment+reserve+(1|PLOT_ID), family="binomial", data=ind_po)

summary(two_way_Chr_api) #final model
Chr_api3way<-anova(three_way_Chr_api,two_way_Chr_api) #p value for 3way int, not significant
Chr_api2way<-anova(two_way_Chr_api, noint_Chr_api) #p value for 2way, no effect of spraying

# Des_var - three and two way
ind_po$Des_var<-ifelse(ind_po$Des_var>0,1,0)
ind_po$DATE<-ind_po$DATE-min(ind_po$DATE)
three_way_Des_var<-glmer(Des_var~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(three_way_Des_var)

two_way_Des_var<-glmer(Des_var~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID), family="binomial", data=ind_po)
noint_Des_var<-glmer(Des_var~DATE+Treatment+reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(two_way_Des_var) #final model
anova(three_way_Des_var,two_way_Des_var)
anova(two_way_Des_var, noint_Des_var) 

Des_var3way<-anova(three_way_Des_var,two_way_Des_var)
Des_var2way<-anova(two_way_Des_var, noint_Des_var) 

# Ery_ovi 
ind_po$Ery_ovi<-ifelse(ind_po$Ery_ovi>0,1,0)
ind_po$DATE<-ind_po$DATE-min(ind_po$DATE)
three_way_Ery_ovi<-glmer(Ery_ovi~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(three_way_Ery_ovi)

two_way_Ery_ovi<-glmer(Ery_ovi~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID), family="binomial", data=ind_po)
noint_Ery_ovi<-glmer(Ery_ovi~DATE+Treatment+reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(two_way_Ery_ovi) #final model
anova(three_way_Ery_ovi,two_way_Ery_ovi) 
anova(two_way_Ery_ovi, noint_Ery_ovi) 

Ery_ovi3way<-anova(three_way_Ery_ovi,two_way_Ery_ovi) 
Ery_ovi2way<-anova(two_way_Ery_ovi, noint_Ery_ovi) 

# Gly_tab
ind_po$Gly_tab<-ifelse(ind_po$Gly_tab>0,1,0)
ind_po$DATE<-ind_po$DATE-min(ind_po$DATE)
three_way_Gly_tab<-glmer(Gly_tab~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(three_way_Gly_tab)

two_way_Gly_tab<-glmer(Gly_tab~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID), family="binomial", data=ind_po)
noint_Gly_tab<-glmer(Gly_tab~DATE+Treatment+reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(two_way_Gly_tab) #final model
anova(three_way_Gly_tab,two_way_Gly_tab)
anova(two_way_Gly_tab, noint_Gly_tab) 

Gly_tab3way<-anova(three_way_Gly_tab,two_way_Gly_tab)
Gly_tab2way<-anova(two_way_Gly_tab, noint_Gly_tab) 

# Lom_cor 
ind_po$Lom_cor<-ifelse(ind_po$Lom_cor>0,1,0)
ind_po$DATE<-ind_po$DATE-min(ind_po$DATE)
three_way_Lom_cor<-glmer(Lom_cor~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(three_way_Lom_cor) 

two_way_Lom_cor<-glmer(Lom_cor~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID), family="binomial", data=ind_po)
noint_Lom_cor<-glmer(Lom_cor~DATE+Treatment+reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(two_way_Lom_cor) #final model
anova(three_way_Lom_cor,two_way_Lom_cor) 
anova(two_way_Lom_cor, noint_Lom_cor) 

Lom_cor3way<-anova(three_way_Lom_cor,two_way_Lom_cor)
Lom_cor2way<-anova(two_way_Lom_cor, noint_Lom_cor) 

# Pla_var
ind_po$Pla_var<-ifelse(ind_po$Pla_var>0,1,0)
ind_po$DATE<-ind_po$DATE-min(ind_po$DATE)
three_way_Pla_var<-glmer(Pla_var~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(three_way_Pla_var)

two_way_Pla_var<-glmer(Pla_var~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID), family="binomial", data=ind_po)
noint_Pla_var<-glmer(Pla_var~DATE+Treatment+reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(two_way_Pla_var) #final model
anova(three_way_Pla_var,two_way_Pla_var) 
anova(two_way_Pla_var, noint_Pla_var) 

Pla_var3way<-anova(three_way_Pla_var,two_way_Pla_var) 
Pla_var2way<-anova(two_way_Pla_var, noint_Pla_var) 

# Tri_ela
ind_po$Tri_ela<-ifelse(ind_po$Tri_ela>0,1,0)
ind_po$DATE<-ind_po$DATE-min(ind_po$DATE)
three_way_Tri_ela<-glmer(Tri_ela~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(three_way_Tri_ela)

two_way_Tri_ela<-glmer(Tri_ela~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID), family="binomial", data=ind_po)
noint_Tri_ela<-glmer(Tri_ela~DATE+Treatment+reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(two_way_Tri_ela) #final model
anova(three_way_Tri_ela,two_way_Tri_ela) 
anova(two_way_Tri_ela, noint_Tri_ela) 

Tri_ela3way<-anova(three_way_Tri_ela,two_way_Tri_ela) 
Tri_ela2way<-anova(two_way_Tri_ela, noint_Tri_ela)

# Tri_pyg
ind_po$Tri_pyg<-ifelse(ind_po$Tri_pyg>0,1,0)
ind_po$DATE<-ind_po$DATE-min(ind_po$DATE)
three_way_Tri_pyg<-glmer(Tri_pyg~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(three_way_Tri_pyg)

two_way_Tri_pyg<-glmer(Tri_pyg~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID), family="binomial", data=ind_po)
noint_Tri_pyg<-glmer(Tri_pyg~DATE+Treatment+reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(two_way_Tri_pyg) #final model
anova(three_way_Tri_pyg,two_way_Tri_pyg) 
anova(two_way_Tri_pyg, noint_Tri_pyg) 

Tri_pyg3way<-anova(three_way_Tri_pyg,two_way_Tri_pyg) 
Tri_pyg2way<-anova(two_way_Tri_pyg, noint_Tri_pyg) 

# add p value columns
ind.sp$binom_3wayP<-NA
ind.sp$binom_2wayP<-NA

ind.sp$binom_3wayP[which(ind.sp$Sp=="Chr_api")]<-Chr_api3way[2,which(colnames(Chr_api3way)=="Pr(>Chisq)")]
ind.sp$binom_2wayP[which(ind.sp$Sp=="Chr_api")]<-Chr_api2way[2,which(colnames(Chr_api2way)=="Pr(>Chisq)")]

ind.sp$binom_3wayP[which(ind.sp$Sp=="Des_var")]<-Des_var3way[2,which(colnames(Des_var3way)=="Pr(>Chisq)")]
ind.sp$binom_2wayP[which(ind.sp$Sp=="Des_var")]<-Des_var2way[2,8]

ind.sp$binom_3wayP[which(ind.sp$Sp=="Ery_ovi")]<-Ery_ovi3way[2,8]
ind.sp$binom_2wayP[which(ind.sp$Sp=="Ery_ovi")]<-Ery_ovi2way[2,8]

ind.sp$binom_2wayP[which(ind.sp$Sp=="Gly_tab")]<-Gly_tab3way[2,which(colnames(Gly_tab3way)=="Pr(>Chisq)")]
ind.sp$binom_3wayP[which(ind.sp$Sp=="Gly_tab")]<-Gly_tab2way[2,which(colnames(Gly_tab2way)=="Pr(>Chisq)")]

ind.sp$binom_3wayP[which(ind.sp$Sp=="Lom_cor")]<-Lom_cor3way[2,8]
ind.sp$binom_2wayP[which(ind.sp$Sp=="Lom_cor")]<-Lom_cor2way[2,8]

ind.sp$binom_3wayP[which(ind.sp$Sp=="Pla_var")]<-Pla_var3way[2,8]
ind.sp$binom_2wayP[which(ind.sp$Sp=="Pla_var")]<-Pla_var2way[2,8]

ind.sp$binom_3wayP[which(ind.sp$Sp=="Tri_ela")]<-Tri_ela3way[2,8]
ind.sp$binom_2wayP[which(ind.sp$Sp=="Tri_ela")]<-Tri_ela2way[2,8]

ind.sp$binom_3wayP[which(ind.sp$Sp=="Tri_pyg")]<-Tri_pyg3way[2,8]
ind.sp$binom_2wayP[which(ind.sp$Sp=="Tri_pyg")]<-Tri_pyg2way[2,8]

# save.image("03_Workspaces/stjw_analysis.RData")

# predictions  
# use nd1 from above so the levels are correct:
nd1<-data.frame(DATE=rep(c(0,1,2),rep(3,3)),Treatment=factor(c("C","A","B"),levels=c("C","A","B")),reserve=factor(c(rep("J",9),rep("M",9)),levels=c("J","M")))
str(nd1)

Chr_api_pr<-pred(model=two_way_Chr_api,new.data=nd1,se.fit=T,type="response")
Des_var_pr<-pred(model=two_way_Des_var, new.data=nd1, se.fit=T, type="response")
Ery_ovi_pr<-pred(model=two_way_Ery_ovi, new.data=nd1, se.fit=T, type="response") #check again
Gly_tab_pr<-pred(model=two_way_Gly_tab, new.data=nd1, se.fit=T, type="response")
Lom_cor_pr<-pred(model=two_way_Lom_cor, new.data=nd1, se.fit=T,type="response")
Pla_var_pr<-pred(model=two_way_Pla_var, new.data=nd1, se.fit=T,type="response")
Tri_ela_pr<-pred(model=two_way_Tri_ela, new.data=nd1, se.fit=T,type="response")
Tri_pyg_pr<-pred(model=two_way_Tri_pyg, new.data=nd1, se.fit=T,type="response")

# save.image("03_Workspaces/stjw_analysis.RData")

# close binomial models ----

#  INDIVIDUAL SPECIES ABUNDANCE MODELS:    	# ----

# Start with species that passed the binomial threshold:
sp.tomodel<-ind.sp$Sp[ind.sp$fit_binom=="yes"]

# For the count data with zeroes removed, use truncated negative binomial model - in glmmADMB package

# Abundance data to model:
ct_ind2<-ct_ind[,c(1:4,which(colnames(ct_ind)%in% sp.tomodel))]
unique(ct_ind2$DATE)
ct_ind2$DATE<-ct_ind2$DATE-min(unique(ct_ind2$DATE))
rahead(ct_ind2,6,6);dim(ct_ind2)

# histograms of count data

dev.new(width=8, height=8, dpi=80, pointsize=16,noRStudioGD = T)
par(mfrow=c(3,3),mar=c(4,4,2,1), oma=c(0,0,0,0), mgp=c(2.5,1,0))

for (i in 1:length(sp.tomodel)){
  sp.thisrun<-sp.tomodel[i]
  data.thisrun<-ct_ind2[,sp.thisrun]
  hist(data.thisrun[which(data.thisrun>0)], xlab="", ylab="", main=sp.thisrun, font.main=1)
} # close for

# How many data points do we have for the truncated nbin models?
abund_data<-data.frame(Sp=names(apply(ct_ind2[,5:ncol(ct_ind2)], 2, function(x) length(which(x>0)))),abund=apply(ct_ind2[,5:ncol(ct_ind2)], 2, function(x) length(which(x>0))))
abund_data<-tidy.df(abund_data)

ind.sp<-merge(ind.sp, abund_data, by="Sp", all.x=T, all.y=F)
ind.sp$fit_abund<-ind.sp$fit_binom

# Fit abundance models for species with greater abundance than the number of quadrats (48). Anything less is not enough data because we have such a complicated model structure. 

ind.sp$fit_abund[which(ind.sp$abund<48)]<-"no"

# Remove Gly_tab from the analysis as even the simple models would converge

ind.sp$fit_abund[which(ind.sp$Sp=="Gly_tab")]<-"no"

# This leaves five species for which we can (attemp to) fit abundance models:
# Some may still not have enough data for complex three way interactions
abund.tomodel<-ind.sp$Sp[ind.sp$fit_abund=="yes"]

# Add indicator to whether three way converged or not. This can be our indication of whether to look at the three way:

ind.sp$converge_3way<-NA

# save.image("03_Workspaces/stjw_analysis.RData")

# use anova to determine the significance of the three way (anova(3way, 2way)) and the significance of the two way (anova(2way, null_model))

# Chr_api
three_way_ca_nb<-glmmadmb(Chr_api~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID),family="truncnbinom1", data=ct_ind2[ct_ind2$Chr_api>0,])

ind.sp$converge_3way[ind.sp$Sp=="Chr_api"]<-"yes"
summary(three_way_ca_nb) #significant 

two_way_ca_nb<-glmmadmb(Chr_api~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID), family="truncnbinom1", data=ct_ind2[ct_ind2$Chr_api>0,])
summary(two_way_ca_nb) #not significant

noint_ca_nb<-glmmadmb(Chr_api~DATE+Treatment+reserve+(1|PLOT_ID), family="truncnbinom1", data=ct_ind2[ct_ind2$Chr_api>0,])
anova(three_way_ca_nb,two_way_ca_nb) #significant 3 way 
anova(two_way_ca_nb,noint_ca_nb) #significant 2 way
AIC(three_way_ca_nb); AIC(two_way_ca_nb)

Chrapi_3waynb<-anova(three_way_ca_nb,two_way_ca_nb)
Chrapi_2waynb<-anova(two_way_ca_nb,noint_ca_nb)

# Des_var 
three_way_dv_nb<-glmmadmb(Des_var~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID),family="truncnbinom1", data=ct_ind2[ct_ind2$Des_var>0,])
summary(three_way_dv_nb) #not significant

ind.sp$converge_3way[ind.sp$Sp=="Des_var"]<-"no"

two_way_dv_nb<-glmmadmb(Des_var~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID), family="truncnbinom1", data=ct_ind2[ct_ind2$Des_var>0,])
summary(two_way_dv_nb) #not significant
noint_dv_nb<-glmmadmb(Des_var~DATE+Treatment+reserve+(1|PLOT_ID), family="truncnbinom1", data=ct_ind2[ct_ind2$Des_var>0,])

# 3 way model did not converge, do not compare models; take the two-way to be the final model:
# anova(three_way_dv_nb,two_way_dv_nb) #not signigicant 
# but note there are fitting problems with this 2 way, possibly not enough abundance data (70)
anova(two_way_dv_nb, noint_dv_nb) # significant 2 way

Desvar_3waynb<-NA
Desvar_2waynb<-anova(two_way_dv_nb, noint_dv_nb)

# Ery_ovi 
three_way_eo_nb<-glmmadmb(Ery_ovi~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID),family="truncnbinom1",data=ct_ind2[ct_ind2$Ery_ovi>0,])
summary(three_way_eo_nb) # not significant

ind.sp$converge_3way[ind.sp$Sp=="Ery_ovi"]<-"yes"

two_way_eo_nb<-glmmadmb(Ery_ovi~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID),family="truncnbinom1",data=ct_ind2[ct_ind2$Ery_ovi>0,])
summary(two_way_eo_nb) #not significant
noint_eo_nb<-glmmadmb(Ery_ovi~DATE+Treatment+reserve+(1|PLOT_ID), family="truncnbinom1", data=ct_ind2[ct_ind2$Ery_ovi>0,])

anova(three_way_eo_nb,two_way_eo_nb) # significant 3way
anova(two_way_eo_nb, noint_eo_nb) # not significant
AIC(three_way_eo_nb); AIC(two_way_eo_nb)

Eryovi_3waynb<-anova(three_way_eo_nb,two_way_eo_nb)
Eryovi_2waynb<-anova(two_way_eo_nb, noint_eo_nb)

# Lom_cor 

three_way_lc_nb<-glmmadmb(Lom_cor~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID),family="truncnbinom1", data=ct_ind2[ct_ind2$Lom_cor>0,])
summary(three_way_lc_nb) # not significant

ind.sp$converge_3way[ind.sp$Sp=="Lom_cor"]<-"no"

two_way_lc_nb<-glmmadmb(Lom_cor~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID),family="truncnbinom1",data=ct_ind2[ct_ind2$Lom_cor>0,])
summary(two_way_lc_nb) # not significant
noint_lc_nb<-glmmadmb(Lom_cor~DATE+Treatment+reserve+(1|PLOT_ID), family="truncnbinom1", data=ct_ind2[ct_ind2$Lom_cor>0,]) 

# 3 way model did not converge, do not compare models; take the two-way to be the final model:
# anova(three_way_lc_nb,two_way_lc_nb) #very highly significant 3way- check

anova(two_way_lc_nb, noint_lc_nb) # not significant

Lomcor_3waynb<-NA
Lomcor_2waynb<-anova(two_way_lc_nb, noint_lc_nb)

# Tri_ela 
three_way_te_nb<-glmmadmb(Tri_ela~DATE+Treatment+reserve+DATE:Treatment+DATE:Treatment:reserve+(1|PLOT_ID), family="truncnbinom1", data=ct_ind2[ct_ind2$Tri_ela>0,])
summary(three_way_te_nb) # not significant

ind.sp$converge_3way[ind.sp$Sp=="Tri_ela"]<-"yes"

two_way_te_nb<-glmmadmb(Tri_ela~DATE+Treatment+reserve+DATE:Treatment+(1|PLOT_ID),family="truncnbinom1",data=ct_ind2[ct_ind2$Tri_ela>0,])
summary(two_way_te_nb) #not significant
noint_te_nb<-glmmadmb(Tri_ela~DATE+Treatment+reserve+(1|PLOT_ID), family="truncnbinom1", data=ct_ind2[ct_ind2$Tri_ela>0,]) 

anova(three_way_te_nb,two_way_te_nb) # not significant
anova(two_way_te_nb, noint_te_nb) # not significant
AIC(three_way_te_nb); AIC(two_way_te_nb)

Triela_3waynb<-anova(three_way_te_nb,two_way_te_nb)
Triela_2waynb<-anova(two_way_te_nb, noint_te_nb)

# add columns to ind.sp table
ind.sp$abund.3p<-NA
ind.sp$abund.2p<-NA

ind.sp$abund.3p[which(ind.sp$Sp=="Chr_api")]<-round(Chrapi_3waynb[2,5],4)
ind.sp$abund.2p[which(ind.sp$Sp=="Chr_api")]<-round(Chrapi_2waynb[2,5],4)

ind.sp$abund.3p[which(ind.sp2$Sp=="Des_var")]<-Desvar_3waynb
ind.sp$abund.2p[which(ind.sp$Sp=="Des_var")]<-Desvar_2waynb[2,5]

ind.sp$abund.3p[which(ind.sp$Sp=="Ery_ovi")]<-Eryovi_3waynb[2,5]
ind.sp$abund.2p[which(ind.sp$Sp=="Ery_ovi")]<-Eryovi_2waynb[2,5]

ind.sp$abund.3p[which(ind.sp$Sp=="Lom_cor")]<-Lomcor_3waynb
ind.sp$abund.2p[which(ind.sp$Sp=="Lom_cor")]<-Lomcor_2waynb[2,5]

ind.sp$abund.3p[which(ind.sp$Sp=="Tri_ela")]<-Triela_3waynb[2,5]
ind.sp$abund.2p[which(ind.sp$Sp=="Tri_ela")]<-Triela_2waynb[2,5]

ind.sp[,c(1,which(colnames(ind.sp)=="abund"):ncol(ind.sp))]

# predictions  

ind.sp[,c(1,which(colnames(ind.sp)=="abund"):ncol(ind.sp))]
head(ind.sp)

# when using the pred function on glmmadmb (i.e. not glmer, for which the pred function uses predictSE), we need to make sure the levels in the new data frame match the levels in the data that were modelled. 

nd1<-data.frame(DATE=rep(c(0,1,2),rep(3,3)),Treatment=factor(c("C","A","B"), levels=c("C","A","B")),reserve=factor(c(rep("J",9),rep("M",9)),levels=c("J","M")))

#Chr_api
model=three_way_ca_nb
new.data=nd1
se.fit=T
type="response"
summary(two_way_ca_nb)$coefficients
Chr_api_nbpr<-pred(model=three_way_ca_nb, new.data = nd1,se.fit = T,type="response") 

Chr_api_nbpr_2way<-pred(model=two_way_ca_nb, new.data = nd1,se.fit = T,type="response") 

# Ery_ovi
model=three_way_af_nb
new.data=nd1
se.fit=T
type="response"
summary(three_way_eo_nb)$coefficients
Ery_ovi_nbpr<-pred(model=three_way_eo_nb,new.data = nd1, se.fit=T, type="response")

#Des_var (predictions not reliable; too many zeros > 50%)
model=two_way_dv_nb
new.data=nd1
se.fit=T
type="response"
Des_var_nbpr<-pred(model=two_way_dv_nb,new.data=nd1,type="response",se.fit=T)
head(Des_var_nbpr)

# save.image("03_Workspaces/stjw_analysis.RData")

# close abundance models ----

#  INDIVIDUAL SPECIES CONTRASTS:    	# ----

# Final models to plot:
summary(two_way_dv_nb) # Des_var two-way (binomial)
summary(three_way_ca_nb) # Chr_api three-way (abund)
summary(three_way_eo_nb) # Ery_ovi three-way (abund)

# Des_var two-way (binomial)
dv.diff<-data.frame(diff.est(two_way_dv_nb,int2way.uzm, int.z.diff),int.names)
dv.diff$diff<-ifelse(sign(dv.diff$lci)==sign(dv.diff$uci),1,0)
dv.diff<-dv.diff[-which(dv.diff$res=="Jer:Mul"),]
dv.diff

# Chr_api three-way (abund)
ca.diff<-data.frame(diff.est(three_way_ca_nb,int3way.uzm, int.z.diff),int.names)
ca.diff$diff<-ifelse(sign(ca.diff$lci)==sign(ca.diff$uci),1,0)
ca.diff<-ca.diff[-which(ca.diff$res=="Jer:Mul"),]
ca.diff

# Ery_ovi three-way (abund)
eo.diff<-data.frame(diff.est(three_way_eo_nb,int3way.uzm, int.z.diff),int.names)
eo.diff$diff<-ifelse(sign(eo.diff$lci)==sign(eo.diff$uci),1,0)
eo.diff<-eo.diff[-which(eo.diff$res=="Jer:Mul"),]
eo.diff

# close CONTRASTS ----

#  PLOT INDIVIDUAL SPECIES:    	# ----

# PLOT all on the same page (binomial and abundance)

# BINOMIAL

# For binomial, there were not significant 3 way interactions and only one significant 2 way (Des_var)
ind.sp[which(ind.sp$binom_3wayP<0.05),]
ind.sp[which(ind.sp$binom_2wayP<0.05),]

# Des_var data
Des_var_pr_M<-Des_var_pr[which(Des_var_pr$reserve=="M"),]
Des_var_pr_J<-Des_var_pr[which(Des_var_pr$reserve=="J"),]
DesvarP<-round(Des_var2way$"Pr(>Chisq)"[2],3)
summary(two_way_Des_var) #final model
head(Des_var_pr,3)

# ABUNDANCE

ind.sp[which(ind.sp$abund.3p<0.05),]
ind.sp[which(ind.sp$abund.2p<0.05),]

# Since Chr_api, Ery_ovi have significant 3 way interations, plot effects for both Mulangarri and Jerra. 
# Des_var has significant 2-way interaction but the model has a poor fit - see abundance estimates

summary(three_way_ca_nb) #final model
summary(three_way_eo_nb) #final model
summary(two_way_dv_nb) #final model

head(Chr_api_nbpr,3)
head(Ery_ovi_nbpr,3)
head(Des_var_nbpr,3)

# Subset data by location:

ca_nbpr_M<-Chr_api_nbpr[which(Chr_api_nbpr$reserve=="M"),]
ca_nbpr_J<-Chr_api_nbpr[which(Chr_api_nbpr$reserve=="J"),]
eo_nbpr_M<-Ery_ovi_nbpr[which(Ery_ovi_nbpr$reserve=="M"),]
eo_nbpr_J<-Ery_ovi_nbpr[which(Ery_ovi_nbpr$reserve=="J"),]
dv_nbpr_M<-Des_var_nbpr[which(Des_var_nbpr$reserve=="M"),]
dv_nbpr_J<-Des_var_nbpr[which(Des_var_nbpr$reserve=="J"),]

# save.image("03_Workspaces/stjw_analysis.RData")

# CONTRASTS
dv.diff
ca.diff
eo.diff

# PLOT

xofs<-0.2
arrowlgth<-0.02

panel.size<-2

dev.new(width=panel.size*3.5,height=panel.size*3,noRStudioGD = T,dpi=80, pointsize=(panel.size*3.5)*2)
par(mfrow=c(2,2), mar=c(4,4,3,2), oma=c(0,0,0,6), mgp=c(2.5,1,0))

## Des_var (binomial)

# Mulanggari
plot(Des_var_pr_M$DATE[Des_var_pr_M$Treatment=="C" ]-xofs,Des_var_pr_M$fit[Des_var_pr_M$Treatment=="C"], pch=15, ylim=c(min(Des_var_pr_M$lci), max(Des_var_pr_M$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=bquote(italic("D. varians ")~.("occurence")), las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))

arrows(Des_var_pr_M$DATE[Des_var_pr_M$Treatment=="C"]-xofs,Des_var_pr_M$lci[Des_var_pr_M$Treatment=="C"],Des_var_pr_M$DATE[Des_var_pr_M$Treatment=="C"]-xofs,Des_var_pr_M$uci[Des_var_pr_M$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(Des_var_pr_M$DATE[Des_var_pr_M$Treatment=="A"],Des_var_pr_M$fit[Des_var_pr_M$Treatment=="A"], pch=15, col="red")
arrows(Des_var_pr_M$DATE[Des_var_pr_M$Treatment=="A"],Des_var_pr_M$lci[Des_var_pr_M$Treatment=="A"],Des_var_pr_M$DATE[Des_var_pr_M$Treatment=="A"],Des_var_pr_M$uci[Des_var_pr_M$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
points(Des_var_pr_M$DATE[Des_var_pr_M$Treatment=="B"]+xofs,Des_var_pr_M$fit[Des_var_pr_M$Treatment=="B"], pch=15, col="blue")
arrows(Des_var_pr_M$DATE[Des_var_pr_M$Treatment=="B"]+xofs,Des_var_pr_M$lci[Des_var_pr_M$Treatment=="B"],Des_var_pr_M$DATE[Des_var_pr_M$Treatment=="B"]+xofs,Des_var_pr_M$uci[Des_var_pr_M$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")

mtext(paste("(",letters[1],") ", "Mulanggari",sep=""), side=3, line=1.2, adj=0, cex=0.9)

# Add star for significant differences
dv.diff[dv.diff$res=="Mul:Mul",]

points(Des_var_pr_M$DATE[Des_var_pr_M$Treatment=="A"][c(2,3)],Des_var_pr_M$uci[Des_var_pr_M$Treatment=="A"][c(2,3)]+0.1, pch="*", cex=1.5,col="black")

# JERRA
plot(Des_var_pr_J$DATE[Des_var_pr_J$Treatment=="C" ]-xofs,Des_var_pr_J$fit[Des_var_pr_J$Treatment=="C"], pch=15, ylim=c(min(Des_var_pr_J$lci), max(Des_var_pr_J$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=bquote(italic("D. varians ")~.("occurence")), las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))

arrows(Des_var_pr_J$DATE[Des_var_pr_J$Treatment=="C"]-xofs,Des_var_pr_J$lci[Des_var_pr_J$Treatment=="C"],Des_var_pr_J$DATE[Des_var_pr_J$Treatment=="C"]-xofs,Des_var_pr_J$uci[Des_var_pr_J$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(Des_var_pr_J$DATE[Des_var_pr_J$Treatment=="A"],Des_var_pr_J$fit[Des_var_pr_J$Treatment=="A"], pch=15, col="red")
arrows(Des_var_pr_J$DATE[Des_var_pr_J$Treatment=="A"],Des_var_pr_J$lci[Des_var_pr_J$Treatment=="A"],Des_var_pr_J$DATE[Des_var_pr_J$Treatment=="A"],Des_var_pr_J$uci[Des_var_pr_J$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
points(Des_var_pr_J$DATE[Des_var_pr_J$Treatment=="B"]+xofs,Des_var_pr_J$fit[Des_var_pr_J$Treatment=="B"], pch=15, col="blue")
arrows(Des_var_pr_J$DATE[Des_var_pr_J$Treatment=="B"]+xofs,Des_var_pr_J$lci[Des_var_pr_J$Treatment=="B"],Des_var_pr_J$DATE[Des_var_pr_J$Treatment=="B"]+xofs,Des_var_pr_J$uci[Des_var_pr_J$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
mtext(paste("(",letters[2],") ", "Jerrabomberra", sep=""), side=3, line=1.2, adj=0, cex=0.9)

if (DesvarP<0.001) title(main=bquote(Two-way~int.~italic(P)~"<"~0.001), font.main=1, adj=0, cex.main=1, line=0.5) else title(main=bquote(Two-way~int.~italic(P)~"="~.(DesvarP)), font.main=1.4, adj=0, cex.main=1, line=0.5)

# Add star for significant differences
dv.diff[dv.diff$res=="Jer:Jer",]

points(Des_var_pr_J$DATE[Des_var_pr_J$Treatment=="A"][c(2,3)],Des_var_pr_J$uci[Des_var_pr_J$Treatment=="A"][c(2,3)]+0.1, pch="*", cex=1.5,col="black")

## Chr_api Mulanggari

plot(ca_nbpr_M$DATE[ca_nbpr_M$Treatment=="C" ]-xofs,ca_nbpr_M$fit.resp[ca_nbpr_M$Treatment=="C"], pch=15, ylim=c(min(ca_nbpr_M$lci.resp), max(ca_nbpr_M$uci.resp)), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=bquote(italic("C. apiculatum")~.("abundance")), las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))

arrows(ca_nbpr_M$DATE[ca_nbpr_M$Treatment=="C"]-xofs,ca_nbpr_M$lci.resp[ca_nbpr_M$Treatment=="C"],ca_nbpr_M$DATE[ca_nbpr_M$Treatment=="C"]-xofs,ca_nbpr_M$uci.resp[ca_nbpr_M$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(ca_nbpr_M$DATE[ca_nbpr_M$Treatment=="A"],ca_nbpr_M$fit.resp[ca_nbpr_M$Treatment=="A"], pch=15, col="red")
arrows(ca_nbpr_M$DATE[ca_nbpr_M$Treatment=="A"],ca_nbpr_M$lci.resp[ca_nbpr_M$Treatment=="A"],ca_nbpr_M$DATE[ca_nbpr_M$Treatment=="A"],ca_nbpr_M$uci.resp[ca_nbpr_M$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(ca_nbpr_M$DATE[ca_nbpr_M$Treatment=="B"]+xofs,ca_nbpr_M$fit.resp[ca_nbpr_M$Treatment=="B"], pch=15, col="blue")
arrows(ca_nbpr_M$DATE[ca_nbpr_M$Treatment=="B"]+xofs,ca_nbpr_M$lci.resp[ca_nbpr_M$Treatment=="B"],ca_nbpr_M$DATE[ca_nbpr_M$Treatment=="B"]+xofs,ca_nbpr_M$uci.resp[ca_nbpr_M$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")

mtext(paste("(",letters[3],") ", "Mulanggari",sep=""), side=3, line=1.2, adj=0, cex=0.9)

# Add star for significant differences (none for this one)
ca.diff[ca.diff$res=="Mul:Mul",]

## Chr_api Jerrabomberra
plot(ca_nbpr_J$DATE[ca_nbpr_J$Treatment=="C" ]-xofs,ca_nbpr_J$fit.resp[ca_nbpr_J$Treatment=="C"], pch=15, ylim=c(min(ca_nbpr_J$lci.resp), max(ca_nbpr_J$uci.resp)), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=bquote(italic("C. apiculatum")~.("abundance")), las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))

arrows(ca_nbpr_J$DATE[ca_nbpr_J$Treatment=="C"]-xofs,ca_nbpr_J$lci.resp[ca_nbpr_J$Treatment=="C"],ca_nbpr_J$DATE[ca_nbpr_J$Treatment=="C"]-xofs,ca_nbpr_J$uci.resp[ca_nbpr_J$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(ca_nbpr_J$DATE[ca_nbpr_J$Treatment=="A"],ca_nbpr_J$fit.resp[ca_nbpr_J$Treatment=="A"], pch=15, col="red")
arrows(ca_nbpr_J$DATE[ca_nbpr_J$Treatment=="A"],ca_nbpr_J$lci.resp[ca_nbpr_J$Treatment=="A"],ca_nbpr_J$DATE[ca_nbpr_J$Treatment=="A"],ca_nbpr_J$uci.resp[ca_nbpr_J$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(ca_nbpr_J$DATE[ca_nbpr_J$Treatment=="B"]+xofs,ca_nbpr_J$fit.resp[ca_nbpr_J$Treatment=="B"], pch=15, col="blue")
arrows(ca_nbpr_J$DATE[ca_nbpr_J$Treatment=="B"]+xofs,ca_nbpr_J$lci.resp[ca_nbpr_J$Treatment=="B"],ca_nbpr_J$DATE[ca_nbpr_J$Treatment=="B"]+xofs,ca_nbpr_J$uci.resp[ca_nbpr_J$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")

mtext(paste("(",letters[4],") ", "Jerrabomberra", sep=""), side=3, line=1.2, adj=0, cex=0.9)

# Add star for significant differences
ca.diff[ca.diff$res=="Jer:Jer",]

points(ca_nbpr_J$DATE[ca_nbpr_J$Treatment=="A"][c(3)],ca_nbpr_J$uci.resp[ca_nbpr_J$Treatment=="A"][c(3)]+3, pch="*", cex=1.5,col="black")

p.Chrapi_3waynb<-round(Chrapi_3waynb[2,which(colnames(Chrapi_3waynb)=="Pr(>Chi)")],3)

if (p.Chrapi_3waynb<0.001) title(main=bquote(Three-way~int.~italic(P)~"<"~0.001), font.main=1, adj=0, cex.main=1, line=0.5) else title(main=bquote(Three-way~int.~italic(P)~"="~.(p.Chrapi_3waynb)), font.main=1.4, adj=0, cex.main=1, line=0.5)

par(xpd=NA)
legend(2.6,16,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)
par(xpd=F)




## Ery_ovi for SI

xofs<-0.2
arrowlgth<-0.02

panel.size<-2

dev.new(width=9,height=4,noRStudioGD = T,dpi=100, pointsize=14)
par(mfrow=c(1,2), mar=c(4,4,3,2), oma=c(0,0,0,6), mgp=c(2.5,1,0))

## Ery_ovi Mulanggari 

plot(eo_nbpr_M$DATE[eo_nbpr_M$Treatment=="C" ]-xofs,eo_nbpr_M$fit.resp[eo_nbpr_M$Treatment=="C"], pch=15, ylim=c(min(eo_nbpr_M$lci.resp), max(eo_nbpr_M$uci.resp)+10), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=bquote(italic("E. ovinum")~.("abundance")), las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))

arrows(eo_nbpr_M$DATE[eo_nbpr_M$Treatment=="C"]-xofs,eo_nbpr_M$lci.resp[eo_nbpr_M$Treatment=="C"],eo_nbpr_M$DATE[eo_nbpr_M$Treatment=="C"]-xofs,eo_nbpr_M$uci.resp[eo_nbpr_M$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(eo_nbpr_M$DATE[eo_nbpr_M$Treatment=="A"],eo_nbpr_M$fit.resp[eo_nbpr_M$Treatment=="A"], pch=15, col="red")
arrows(eo_nbpr_M$DATE[eo_nbpr_M$Treatment=="A"],eo_nbpr_M$lci.resp[eo_nbpr_M$Treatment=="A"],eo_nbpr_M$DATE[eo_nbpr_M$Treatment=="A"],eo_nbpr_M$uci.resp[eo_nbpr_M$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
points(eo_nbpr_M$DATE[eo_nbpr_M$Treatment=="B"]+xofs,eo_nbpr_M$fit.resp[eo_nbpr_M$Treatment=="B"], pch=15, col="blue")
arrows(eo_nbpr_M$DATE[eo_nbpr_M$Treatment=="B"]+xofs,eo_nbpr_M$lci.resp[eo_nbpr_M$Treatment=="B"],eo_nbpr_M$DATE[eo_nbpr_M$Treatment=="B"]+xofs,eo_nbpr_M$uci.resp[eo_nbpr_M$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")

mtext(paste("(",letters[1],") ", "Mulanggari",sep=""), side=3, line=1.2, adj=0, cex=1)

# Add star for significant differences
eo.diff[eo.diff$res=="Mul:Mul",]

points((eo_nbpr_M$DATE[eo_nbpr_M$Treatment=="B"]+xofs)[c(1,2,3)],eo_nbpr_M$uci.resp[eo_nbpr_M$Treatment=="B"][c(1,2,3)]+5, pch="*", cex=1.5,col="black")

## Ery_ovi Jerrabomberra
plot(eo_nbpr_J$DATE[eo_nbpr_J$Treatment=="C" ]-xofs,eo_nbpr_J$fit.resp[eo_nbpr_J$Treatment=="C"], pch=15, ylim=c(min(eo_nbpr_J$lci.resp), max(eo_nbpr_J$uci.resp)), xlim=c(-0.3,2.3), xaxt="n", xlab="Year", ylab=bquote(italic("E. ovinum")~.("abundance")), las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))

arrows(eo_nbpr_J$DATE[eo_nbpr_J$Treatment=="C"]-xofs,eo_nbpr_J$lci.resp[eo_nbpr_J$Treatment=="C"],eo_nbpr_J$DATE[eo_nbpr_J$Treatment=="C"]-xofs,eo_nbpr_J$uci.resp[eo_nbpr_J$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(eo_nbpr_J$DATE[eo_nbpr_J$Treatment=="A"],eo_nbpr_J$fit.resp[eo_nbpr_J$Treatment=="A"], pch=15, col="red")
arrows(eo_nbpr_J$DATE[eo_nbpr_J$Treatment=="A"],eo_nbpr_J$lci.resp[eo_nbpr_J$Treatment=="A"],eo_nbpr_J$DATE[eo_nbpr_J$Treatment=="A"],eo_nbpr_J$uci.resp[eo_nbpr_J$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
points(eo_nbpr_J$DATE[eo_nbpr_J$Treatment=="B"]+xofs,eo_nbpr_J$fit.resp[eo_nbpr_J$Treatment=="B"], pch=15, col="blue")
arrows(eo_nbpr_J$DATE[eo_nbpr_J$Treatment=="B"]+xofs,eo_nbpr_J$lci.resp[eo_nbpr_J$Treatment=="B"],eo_nbpr_J$DATE[eo_nbpr_J$Treatment=="B"]+xofs,eo_nbpr_J$uci.resp[eo_nbpr_J$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")

mtext(paste("(",letters[2],") ", "Jerrabomberra", sep=""), side=3, line=1.2, adj=0, cex=1)

# Add star for significant differences
eo.diff[eo.diff$res=="Jer:Jer",]

points((eo_nbpr_J$DATE[eo_nbpr_J$Treatment=="B"]+xofs)[c(1)],eo_nbpr_J$uci.resp[eo_nbpr_J$Treatment=="B"][c(1)]+5, pch="*", cex=1.5,col="black")

p.Eryovi_3waynb<-round(Eryovi_3waynb[2,which(colnames(Eryovi_3waynb)=="Pr(>Chi)")],3)

if (p.Eryovi_3waynb<0.001) title(main=bquote(Three-way~int.~italic(P)~"<"~0.001), font.main=1, adj=0, cex.main=1, line=0.5) else title(main=bquote(Three-way~int.~italic(P)~"="~.(p.Eryovi_3waynb)), font.main=1.4, adj=0, cex.main=1, line=0.5)

par(xpd=NA)
legend(2.6,16,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)
par(xpd=F)

# close invidividual species plots ----




