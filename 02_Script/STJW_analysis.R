
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
load("03_Workspaces/stjw_analysis.RData")

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

# save.image("03_Workspaces/stjw_analysis.RData")

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

rahead(ct17,4,10); dim(ct17)
rahead(ct18,4,10); dim(ct18)
rahead(ct19,4,10); dim(ct19)

rahead(cv17,4,10); dim(cv17)
rahead(cv18,4,10); dim(cv18)
rahead(cv19,4,10); dim(cv19)

# save.image("03_Workspaces/stjw_analysis.RData")

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

# save.image("03_Workspaces/stjw_analysis.RData")

# close format data ----

#  species DIVERSITY GROUPS:    	# ----

head(pinfo,3); dim(pinfo)
str(pinfo)
table(pinfo$func_grp)

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

# save.image("03_Workspaces/stjw_analysis.RData")

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

#  ANALYSIS (COMPONENT 1):    	# ----

# Notes on model selection to include in paper:

# We fit an interaction between treatment and date in all models to describe the before-after, control-impact nature of our design. 

# Our full model included an interaction between date, treatment and reserve, but we removed this term if it was not significant (< 0.05)

# we did not expect the diversity differences to change over time, so we did not fit a year and reserve interaction

# Dcide the cut-off for analysis. There must be more records (i.e. number of species recorded) than the number of quadrats (48)? I.e. exclude two responses (exotic_perengrass and sed_rus):

gdf<-gdf[-which(gdf$rich_records<48),]
gdf<-tidy.df(gdf)
head(gdf); dim(gdf)
rahead(rich,4,7); dim(rich)
rahead(shan,4,7); dim(shan)

# save workspace:
# save.image("03_Workspaces/stjw_analysis.RData")

# **** RICHNESS:

# scale date only:
rich_sc<-rich
rich_sc$DATE<-rich_sc$DATE-2017
rahead(rich_sc,4,7); dim(rich_sc)

# lists for storing model fits, coefficients, anova tables and model estimates:
fits.rich<-list()
coef.rich<-list()
preds.rich<-list()

# This is a different kind of anova table than for diversity because we are just using it to test the significance of the two-way date:treatment interaction in the final models. We can get the other p values from the coefficient table:
anova.rich<-list()

# new data for model estimates (same for models with a date:treatment interaction only and models with a three way interaction; you can also use the same newdata frame for richness and shannon's):
nd1<-data.frame(DATE=rep(c(0,1,2),rep(3,3)),Treatment=as.factor(c("C","A","B")),reserve=c(rep("J",9),rep("M",9)))

# THE ACTUAL ANALYSIS CODE:
# RUN MODELLING LOOP:

# **** RICHNESS:

# follow 'best bet' instructions for installation:
# https://glmmadmb.r-forge.r-project.org/
library("glmmADMB")

for (i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  data.set<-rich_sc
  data.thisrun<-rich_sc[,resp.thisrun]
  form.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID)", sep="")
  
  # some functional groups have 0, its not going to run those functional groups here on
  if(sum(data.thisrun,na.rm=T)==0){
    fits.rich[[i]]<-NULL
    coef.rich[[i]]<-NULL
    preds.rich[[i]]<-NULL
    anova.rich[[i]]<-NULL
    next
  }
  
  m1<-glmmadmb(as.formula(form.thisrun), family="poisson", data=data.set)
  
  # run model without three-way:
  form.twoway<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+(1|PLOT_ID)", sep="")
  mod_twoway<-glmmadmb(as.formula(form.twoway), family="poisson", data=data.set)
  
  summary(m1)
  
  m1_coef<-coef.ext(m1)
  m1_anova<-anova(mod_twoway, m1)
  
  p_anova<-m1_anova[2,5]
  anova.rich[[i]]<-m1_anova
  
  # simplify model if the three way is not significant:
  if(p_anova>0.05){
    
    # remove three way term from formula:
    form.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+(1|PLOT_ID)", sep="")
    
    # re-run model without three way:
    m1<-glmmadmb(as.formula(form.thisrun), family="poisson", data=data.set)
    summary(m1)
    
    m1_coef<-coef.ext(m1)
    
    # run model without any interaction (to figure out if the two way interaction is significant):
    form.noint<-paste(resp.thisrun,"~Treatment+DATE+reserve+(1|PLOT_ID)", sep="")
    mod_noint<-glmmadmb(as.formula(form.noint), family="poisson", data=data.set)
    
    int_term_anova<-anova(mod_noint, m1)
    anova.rich[[i]]<-int_term_anova
    
  } # close if three way not signif
  
  fits.rich[[i]]<-m1
  coef.rich[[i]]<-m1_coef
  
  # generate model predictions:
  
  m1_pr<-predict(object=m1,type="link",newdata=nd1, se.fit = T)
  m1_pr<-data.frame(nd1,fit=m1_pr$fit,se=m1_pr$se.fit)
  m1_pr$lci<-m1_pr$fit-(m1_pr$se*1.96)
  m1_pr$uci<-m1_pr$fit+(m1_pr$se*1.96)
  m1_pr$fit.resp<-round(exp(m1_pr$fit),4)
  m1_pr$lci.resp<-round(exp(m1_pr$lci),4)
  m1_pr$uci.resp<-round(exp(m1_pr$uci),4)
  
  head(m1_pr)
  
  preds.rich[[i]]<-m1_pr
  
} # close model

# Only one three-way significant for richness (12, exotic annual forb):

summary(fits.rich[[1]])
coef.rich
anova.rich[[12]]
preds.rich

# Summarise results (RICHNESS):

gdf.rich<-gdf
head(gdf.rich); dim(gdf.rich)

# The anova table for richness (Poisson models) was a Likelihood Ratio Test comparing the full 3-way interacion model and the nested model without the 3-way. Where P was > 0.05 the 3-way was removed. When the 3-way was removed, we repeated the LRT test on the interaction model, compared with a nested model without any interactions. Thus, the anova table gives a P value for the 3 way interaction if it was kept, and a P value for the 2 way interaction for all other models:
rich.anova.ps<-unlist(lapply(anova.rich,function(x)x[2,5]))

# Was the three-way significant?

# Three-way interactions have a coefficient table with 10 rows, while those without have only seven rows. Thus we can use the length of the coef table to designate whether the three way was significant:

gdf.rich$sig.3way<-ifelse(unlist(lapply(coef.rich,nrow))==10,"yes","no")
gdf.rich$P.3way<-NA
gdf.rich$P.3way[which(gdf.rich$sig.3way=="yes")]<-rich.anova.ps[which(gdf.rich$sig.3way=="yes")]

# For models without a three way, was the two-way significant?

gdf.rich$P.2way<-rich.anova.ps
gdf.rich$P.2way[which(gdf.rich$sig.3way=="yes")]<-NA
gdf.rich$sig.2way<-ifelse(gdf.rich$P.2way<0.05,"yes","no")

# For models without a three way, was the reserve main effect significant?

# This can come directly from the coefficient table since there are only two levels in this factor:

gdf.rich$P.res<-unlist(lapply(coef.rich,function(x)x[x$term=="reserveM",4]))
gdf.rich$sig.res<-ifelse(gdf.rich$P.res<0.05,"yes","no")
gdf.rich$P.res[which(gdf.rich$sig.3way=="yes")]<-NA
gdf.rich$sig.res[which(gdf.rich$sig.3way=="yes")]<-NA
head(gdf.rich); dim(gdf.rich)

# save.image("03_Workspaces/stjw_analysis.RData")

# **** DIVERSITY:

#scale date only
shan_sc<-shan
shan_sc$DATE<-shan_sc$DATE-2017
rahead(shan_sc,4,7); dim(shan_sc)

# lists for storing model fits, coefficients, anova tables and model estimates:
fits.shan<-list()
coef.shan<-list()
anova.shan<-list()
preds.shan<-list()

# new data for model estimates (same for models with a date:treatment interaction only and models with a three way interaction; you can also use the same newdata frame for richness and shannon's):
nd1<-data.frame(DATE=rep(c(0,1,2),rep(3,3)),Treatment=as.factor(c("C","A","B")),reserve=c(rep("J",9),rep("M",9)))

# RUN MODELLING LOOP:
# **** DIVERSITY:

for (i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  data.set<-shan_sc
  data.thisrun<-shan_sc[,resp.thisrun]
  form.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID)", sep="")
  
  # some functional groups have 0, its not going to run those functional groups here on
  if(sum(data.thisrun,na.rm=T)==0){
    fits.shan[[i]]<-NULL
    coef.shan[[i]]<-NULL
    anova.shan[[i]]<-NULL
    preds.shan[[i]]<-NULL
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
  
  fits.shan[[i]]<-m1
  coef.shan[[i]]<-m1_coef
  anova.shan[[i]]<-m1_anova
  
  # generate model predictions:
  
  m1_pr<-predictSE(mod=m1,newdata=nd1, se.fit = T)
  m1_pr<-data.frame(nd1,fit=m1_pr$fit,se=m1_pr$se.fit)
  m1_pr$lci<-m1_pr$fit-(m1_pr$se*1.96)
  m1_pr$uci<-m1_pr$fit+(m1_pr$se*1.96)
  head(m1_pr)
  
  preds.shan[[i]]<-m1_pr
  
} # close model

summary(fits.shan[[1]])
coef.shan
anova.shan[[19]] # 19 and 6 significant three-way
preds.shan

rahead(shan_sc,3,5)

# save.image("03_Workspaces/stjw_analysis.RData")

# Summarise results (DIVERISTY):

gdf.shan<-gdf
head(gdf.shan); dim(gdf.shan)

# The anova table for diversity (normal models) gives P values for all terms directly, rather than the LRT format of the Poisson models. 

# Was the three-way significant?

# Three-way interactions have a anova table with five rows, while those without have only four rows. Use the length of the anova table to designate whether the three way was significant:

gdf.shan$sig.3way<-ifelse(unlist(lapply(anova.shan,nrow))==5,"yes","no")
gdf.shan$P.3way<-NA
gdf.shan$P.3way[which(gdf.shan$sig.3way=="yes")]<-unlist(lapply(anova.shan,function(x) x[which(x$term=="Treatment:DATE:reserve"),"p"]))
head(gdf.shan); dim(gdf.shan)

# For models without a three way, was the two-way significant?

gdf.shan$P.2way<-unlist(lapply(anova.shan,function(x) x[which(x$term=="Treatment:DATE"),"p"]))
gdf.shan$P.2way[which(gdf.shan$sig.3way=="yes")]<-NA
gdf.shan$sig.2way<-ifelse(gdf.shan$P.2way<0.05,"yes","no")

# For models without a three way, was the reserve main effect significant?

# This can come directly from the coefficient table since there are only two levels in this factor:

gdf.shan$P.res<-unlist(lapply(coef.shan,function(x)x[x$term=="reserveM",4]))
gdf.shan$sig.res<-ifelse(gdf.shan$P.res<0.05,"yes","no")
gdf.shan$P.res[which(gdf.shan$sig.3way=="yes")]<-NA
gdf.shan$sig.res[which(gdf.shan$sig.3way=="yes")]<-NA
head(gdf.shan); dim(gdf.shan)

# save.image("03_Workspaces/stjw_analysis.RData")

# what's going on with native_legherb?
nlh<-lmer(native_legherb~Treatment+DATE+reserve+Treatment:DATE+(1|PLOT_ID), data=shan_sc)
summary(nlh)
# hist(shan_sc$native_legherb) # there's very little data in this group and there isvery high proportion of zeros
gdf.shan
rahead(shan_sc,6,6)
range(shan_sc$native_legherb)
range(rich_sc$native_legherb)
range(shan_sc[2:nrow(shan_sc),5:ncol(shan_sc)])
range(rich_sc[2:nrow(rich_sc),5:ncol(rich_sc)])

# what's going on with exann_herb?
exanh_mod<-lmer(exann_herb~Treatment+DATE+reserve+Treatment:DATE+(1|PLOT_ID), data=shan_sc)
summary(exanh_mod)
anova(exanh_mod)

# close analysis ----

#  PLOT ESTIMATES (COMPONENT 1):    	# ----

# Species Richness:

# Twenty-six fitted models
# One significant 3-way (Exotic annual forb) 
# Eight significant 2-way
# Seventeen significant reserve effects

# For all models, plot the treatment by year effect, even if the two way wasn't significant. This is so we have an overall visual of the results. 
# This will give 26 treatment:year plots + one 3-way = 27 panels
# For all models without a three-way, plot MULANGARRI only. We will plot the reserve effects separately, on another panel. For the three way, plot both reserves

head(gdf.rich); dim(gdf.rich)
fits.rich
coef.rich
anova.rich
preds.rich

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(2,4,4,1), mgp=c(2.5,1,0))

for(i in 1:nrow(gdf.rich)){
  
  resp.thisrun<-gdf.rich$group[i]
  pred.thisrun<-preds.rich[[i]]
  pred.thisrun<-pred.thisrun[which(pred.thisrun$reserve=="M"),]
  anova.thisrun<-anova.rich[[i]]
  coef.thisrun<-coef.rich[[i]]
  ylab.thisrun<-gdf.rich$ylab[i]
  meta.thisrun<-gdf.rich[i,]
  xofs<-0.2
  arrowlgth<-0.02
  
  head(pred.thisrun)
  
  if (is.null(pred.thisrun)) next
  
  if(meta.thisrun$sig.3way=="yes"){
    
    jerra.preds<-preds.rich[[i]]
    jerra.preds<-jerra.preds[which(jerra.preds$reserve=="J"),]
    
    # MULANGARRI
    
    plot(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$fit.resp[pred.thisrun$Treatment=="C"], pch=15, ylim=c(min(pred.thisrun$lci.resp), max(pred.thisrun$uci.resp)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
    axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
    arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$lci.resp[pred.thisrun$Treatment=="C"],pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$uci.resp[pred.thisrun$Treatment=="C"], code=3, angle=90, length=arrowlgth)
    points(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$fit.resp[pred.thisrun$Treatment=="A"], pch=15, col="red")
    arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$lci.resp[pred.thisrun$Treatment=="A"],pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$uci.resp[pred.thisrun$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
    points(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$fit.resp[pred.thisrun$Treatment=="B"], pch=15, col="blue")
    arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$lci.resp[pred.thisrun$Treatment=="B"],pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$uci.resp[pred.thisrun$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
    
    p.trt_yr_int<-round(anova.thisrun$"Pr(>Chi)"[2],4)
    title(main=paste("Three-way int, ","P=",p.trt_yr_int, sep=""), font.main=1, adj=0, cex=0.9, line=0.5)
    
    mtext("Mulangarri", side=3, line=1.5, adj=0)
    
    # JERRA
    
    plot(jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$fit.resp[jerra.preds$Treatment=="C"], pch=15, ylim=c(min(jerra.preds$lci.resp), max(jerra.preds$uci.resp)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
    axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
    arrows(jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$lci.resp[jerra.preds$Treatment=="C"],jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$uci.resp[jerra.preds$Treatment=="C"], code=3, angle=90, length=arrowlgth)
    points(jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$fit.resp[jerra.preds$Treatment=="A"], pch=15, col="red")
    arrows(jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$lci.resp[jerra.preds$Treatment=="A"],jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$uci.resp[jerra.preds$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
    points(jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$fit.resp[jerra.preds$Treatment=="B"], pch=15, col="blue")
    arrows(jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$lci.resp[jerra.preds$Treatment=="B"],jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$uci.resp[jerra.preds$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
    
    p.trt_yr_int<-round(anova.thisrun$"Pr(>Chi)"[2],4)
    title(main=paste("Three-way int, ","P=",p.trt_yr_int, sep=""), font.main=1, adj=0, cex=0.9, line=0.5)
    
    mtext("Jerrabomberra", side=3, line=1.5, adj=0)
    
  } # close three way plot
  
  if(meta.thisrun$sig.3way=="no"){
    
  plot(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$fit.resp[pred.thisrun$Treatment=="C"], pch=15, ylim=c(min(pred.thisrun$lci.resp), max(pred.thisrun$uci.resp)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
  axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
  arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$lci.resp[pred.thisrun$Treatment=="C"],pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$uci.resp[pred.thisrun$Treatment=="C"], code=3, angle=90, length=arrowlgth)
  points(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$fit.resp[pred.thisrun$Treatment=="A"], pch=15, col="red")
  arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$lci.resp[pred.thisrun$Treatment=="A"],pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$uci.resp[pred.thisrun$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
  points(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$fit.resp[pred.thisrun$Treatment=="B"], pch=15, col="blue")
  arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$lci.resp[pred.thisrun$Treatment=="B"],pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$uci.resp[pred.thisrun$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
  
  p.res<-round(coef.thisrun$P[which(coef.thisrun$term=="reserveM")],4)
  p.yr<-round(coef.thisrun$P[which(coef.thisrun$term=="DATE")],4)
  p.trt_yr_int<-round(anova.thisrun$"Pr(>Chi)"[2],4)
  
  if(p.res>0.01) p.res<-round(p.res,2) else p.res<-"<0.01"
  if(p.yr>0.01) p.yr<-round(p.yr,2) else p.yr<-"<0.01"
  if(p.trt_yr_int>0.01) p.trt_yr_int<-round(p.trt_yr_int,2) else p.trt_yr_int<-"<0.01"
  
  title(main=paste("P values: res=",p.res,"\nyr=",p.yr,"; trt x yr=",p.trt_yr_int, sep=""), font.main=1, adj=0, cex=0.9, line=0.5)
  
  } # close plot 2-way
  
} # close richness plot

plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
legend(1,9,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)

## PLOT RESERVE EFFECTS (richness)

res.rich<-gdf.rich[which(gdf.rich$sig.res=="yes"),]
res.rich<-tidy.df(res.rich)
head(res.rich); dim(res.rich)

dev.new(width=11.69,height=(8.27/5)*3,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(3,6), mar=c(2,4,2,1), mgp=c(2.5,1,0))

for (i in 1:nrow(res.rich)){
  
  resp.thisrun<-res.rich$group[i]
  index.ingdf<-which(gdf.rich$group==resp.thisrun)
  pred.thisrun<-preds.rich[[index.ingdf]]
  pred.thisrun<-pred.thisrun[pred.thisrun$DATE==0 & pred.thisrun$Treatment=="C",]
  coef.thisrun<-coef.rich[[index.ingdf]]
  ylab.thisrun<-res.rich$ylab[i]

  plot(1:2,pred.thisrun$fit.resp, xlim=c(0.5, 2.5), ylim=c(min(pred.thisrun$lci.resp), max(pred.thisrun$uci.resp)), pch=20, xaxt="n", ylab=ylab.thisrun, las=1)
  arrows(1:2, pred.thisrun$lci.resp, 1:2, pred.thisrun$uci.resp, code=3, angle=90, length=0.02)
  axis(side=1,at = 1:2,labels=c("Jerra","Mulangarri"))
  p.res<-round(coef.thisrun$P[which(coef.thisrun$term=="reserveM")],4)
  if(p.res>0.001) p.res<-round(p.res,3) else p.res<-"<0.001"
  
  title(main=paste("P value: res=",p.res,sep=""), font.main=1, adj=0, cex=0.9, line=0.5)
  
} # close plot reserve

plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
text(1,8,labels="Estimates of reserve\neffects displayed for\ncontrol in 2017", adj=0)

# Diversity:

# Twenty-six fitted models
# Two significant 3-way (Exotic leg. forb and SigB) 
# One significant 2-way (Exotic forb)
# Thirteen significant reserve effects

# For all models, plot the treatment by year effect, even if the two way wasn't significant. This is so we have an overall visual of the results. 
# This will give 24 treatment:year plots + two 3-way = 28 panels
# For all models without a three-way, plot MULANGARRI only. We will plot the reserve effects separately, on another panel. For the three way, plot both reserves

head(gdf.shan); dim(gdf.shan)
fits.shan
coef.shan
anova.shan
preds.shan

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(2,4,4,1), mgp=c(2.5,1,0))

for(i in 1:nrow(gdf.shan)){
  
  resp.thisrun<-gdf.shan$group[i]
  pred.thisrun<-preds.shan[[i]]
  pred.thisrun<-pred.thisrun[which(pred.thisrun$reserve=="M"),]
  anova.thisrun<-anova.shan[[i]]
  coef.thisrun<-coef.shan[[i]]
  ylab.thisrun<-gdf.shan$ylab[i]
  meta.thisrun<-gdf.shan[i,]
  xofs<-0.2
  arrowlgth<-0.02
  
  if (is.null(pred.thisrun)) next
  
  if(meta.thisrun$sig.3way=="yes"){
    
    jerra.preds<-preds.shan[[i]]
    jerra.preds<-jerra.preds[which(jerra.preds$reserve=="J"),]
    
    # MULANGARRI
    
    plot(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$fit[pred.thisrun$Treatment=="C"], pch=15, ylim=c(min(pred.thisrun$lci), max(pred.thisrun$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
    axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
    arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$lci[pred.thisrun$Treatment=="C"],pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$uci[pred.thisrun$Treatment=="C"], code=3, angle=90, length=arrowlgth)
    points(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$fit[pred.thisrun$Treatment=="A"], pch=15, col="red")
    arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$lci[pred.thisrun$Treatment=="A"],pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$uci[pred.thisrun$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
    points(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$fit[pred.thisrun$Treatment=="B"], pch=15, col="blue")
    arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$lci[pred.thisrun$Treatment=="B"],pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$uci[pred.thisrun$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
    
    p.trt_yr_int<-round(anova.thisrun$p[which(anova.thisrun$term=="Treatment:DATE:reserve")],4)
    title(main=paste("Three-way int, ","P=",p.trt_yr_int, sep=""), font.main=1, adj=0, cex=0.9, line=0.5)
    
    mtext("Mulangarri", side=3, line=1.5, adj=0)
    
    # JERRA
    
    plot(jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$fit[jerra.preds$Treatment=="C"], pch=15, ylim=c(min(jerra.preds$lci), max(jerra.preds$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
    axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
    arrows(jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$lci[jerra.preds$Treatment=="C"],jerra.preds$DATE[jerra.preds$Treatment=="C"]-xofs,jerra.preds$uci[jerra.preds$Treatment=="C"], code=3, angle=90, length=arrowlgth)
    points(jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$fit[jerra.preds$Treatment=="A"], pch=15, col="red")
    arrows(jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$lci[jerra.preds$Treatment=="A"],jerra.preds$DATE[jerra.preds$Treatment=="A"],jerra.preds$uci[jerra.preds$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
    points(jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$fit[jerra.preds$Treatment=="B"], pch=15, col="blue")
    arrows(jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$lci[jerra.preds$Treatment=="B"],jerra.preds$DATE[jerra.preds$Treatment=="B"]+xofs,jerra.preds$uci[jerra.preds$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
    
    p.trt_yr_int<-round(anova.thisrun$p[which(anova.thisrun$term=="Treatment:DATE:reserve")],4)
    title(main=paste("Three-way int, ","P=",p.trt_yr_int, sep=""), font.main=1, adj=0, cex=0.9, line=0.5)
    
    mtext("Jerrabomberra", side=3, line=1.5, adj=0)
    
  } # close three way plot
  
  if(meta.thisrun$sig.3way=="no"){
    
    plot(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$fit[pred.thisrun$Treatment=="C"], pch=15, ylim=c(min(pred.thisrun$lci), max(pred.thisrun$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
    axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
    arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$lci[pred.thisrun$Treatment=="C"],pred.thisrun$DATE[pred.thisrun$Treatment=="C"]-xofs,pred.thisrun$uci[pred.thisrun$Treatment=="C"], code=3, angle=90, length=arrowlgth)
    points(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$fit[pred.thisrun$Treatment=="A"], pch=15, col="red")
    arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$lci[pred.thisrun$Treatment=="A"],pred.thisrun$DATE[pred.thisrun$Treatment=="A"],pred.thisrun$uci[pred.thisrun$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")
    points(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$fit[pred.thisrun$Treatment=="B"], pch=15, col="blue")
    arrows(pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$lci[pred.thisrun$Treatment=="B"],pred.thisrun$DATE[pred.thisrun$Treatment=="B"]+xofs,pred.thisrun$uci[pred.thisrun$Treatment=="B"], code=3, angle=90, length=arrowlgth, col="blue")
    
    p.res<-round(coef.thisrun$P[which(coef.thisrun$term=="reserveM")],4)
    p.yr<-round(coef.thisrun$P[which(coef.thisrun$term=="DATE")],4)
    p.trt_yr_int<-round(anova.thisrun$p[which(anova.thisrun$term=="Treatment:DATE")],4)
    
    if(p.res>0.01) p.res<-round(p.res,2) else p.res<-"<0.01"
    if(p.yr>0.01) p.yr<-round(p.yr,2) else p.yr<-"<0.01"
    if(p.trt_yr_int>0.001) p.trt_yr_int<-round(p.trt_yr_int,3) else p.trt_yr_int<-"<0.001"
    
    title(main=paste("P values: res=",p.res,"\nyr=",p.yr,"; trt x yr=",p.trt_yr_int, sep=""), font.main=1, adj=0, cex=0.9, line=0.5)
    
  } # close plot 2-way
  
} # close diversity plot

plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
legend(1,9,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)

## PLOT RESERVE EFFECTS (diversity)

res.shan<-gdf.shan[which(gdf.shan$sig.res=="yes"),]
res.shan<-tidy.df(res.shan)
head(res.shan); dim(res.shan)

dev.new(width=(11.69/6)*5,height=(8.27/5)*3,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(3,5), mar=c(2,4,2,1), mgp=c(2.5,1,0))

for (i in 1:nrow(res.shan)){
  
  resp.thisrun<-res.shan$group[i]
  index.ingdf<-which(gdf.shan$group==resp.thisrun)
  pred.thisrun<-preds.shan[[index.ingdf]]
  pred.thisrun<-pred.thisrun[pred.thisrun$DATE==0 & pred.thisrun$Treatment=="C",]
  coef.thisrun<-coef.shan[[index.ingdf]]
  ylab.thisrun<-res.shan$ylab[i]
  
  plot(1:2,pred.thisrun$fit, xlim=c(0.5, 2.5), ylim=c(min(pred.thisrun$lci), max(pred.thisrun$uci)), pch=20, xaxt="n", ylab=ylab.thisrun, las=1)
  arrows(1:2, pred.thisrun$lci, 1:2, pred.thisrun$uci, code=3, angle=90, length=0.02)
  axis(side=1,at = 1:2,labels=c("Jerra","Mulangarri"))
  p.res<-round(coef.thisrun$P[which(coef.thisrun$term=="reserveM")],4)
  if(p.res>0.001) p.res<-round(p.res,3) else p.res<-"<0.001"
  
  title(main=paste("P value: res=",p.res,sep=""), font.main=1, adj=0, cex=0.9, line=0.5)
  
} # close plot reserve

plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
text(1,8,labels="Estimates of reserve\neffects displayed for\ncontrol in 2017", adj=0)

# close plot component 1 ----

#  INDIVIDUAL SPECIES:    	# ----

# Need to check with RM: why was Desmodium varians classified as a legume rather than a native leg forb?
# Is it OK for us to classify it as native leg forb (I have done this above in the 'format data' section)?
# PlantNet classifies it as a "Prostrate or climbing herb"
# And it is the ONLY species in the Legume category
table(pinfo$func_grp)

# save.image("03_Workspaces/stjw_analysis.RData")

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

ind.sp<-data.frame(ind_species=c("Eryngium ovinum", "Chrysocephalum apiculatum", "Arthropodium fimbriatum", "Wurmbea dioica", "Desmodium varians", "Plantago varia", "Tricoryne elatior", "Triptilodiscus pygmaeus","Lomandra filiformis coriacea", "Lomandra bracteata", "Lomandra filiformis","Lomandra multiflora","Glycine clandestina","Glycine tabacina"))

# These should all be TRUE:
table(ind.sp$ind_species %in% pinfo$Species)

ind.sp<-merge(ind.sp, pinfo, by.x="ind_species", by.y="Species")
ind.sp
indsp<-ind.sp$Sp

ct_ind<-ct_dat[,c(1:4, which(colnames(ct_dat) %in% indsp))]
cv_ind<-cv_dat[,c(1:4, which(colnames(cv_dat) %in% indsp))]

# Make sure columns match:
table(colnames(ct_ind)==colnames(cv_ind))

# Model individual species using these data sets:
rahead(ct_ind,6,6)
rahead(cv_ind,6,6)

ind_po<-ct_ind
rahead(ind_po,6,6); dim(ind_po)

ind_po$Chr_api<-ifelse(ind_po$Chr_api>0,1,0)
ind_po$DATE<-ind_po$DATE-min(ind_po$DATE)

test_mod<-glmer(Chr_api~DATE*Treatment*reserve+(1|PLOT_ID), family="binomial", data=ind_po)
summary(test_mod)


# close indiv species ----

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





