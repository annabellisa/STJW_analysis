
# ------------------------------------ #
# ---------- STJW ANALYSIS  ---------- #
# ------------------------------------ #

### Analysis of chemical control experiment from STJW project
### Author: Annabel Smith & Raagini Muddaiah

# Load functions:
invisible(lapply(paste("01_Functions/",dir("01_Functions"),sep=""),function(x) source(x)))

# Load libraries:
library(lme4)

# Load workspace
load("03_workspaces/xx.RData")

#  IMPORT data:    	# ----

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
head(pinfo,3)
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


rahead(ct17,3,7)
rahead(ct18,3,7)
rahead(ct19,3,7)





