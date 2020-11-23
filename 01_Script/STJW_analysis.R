# updating script

Data_dir<-"Data/Final_datasets/Formatted_data"
Data_dir
dir(Data_dir)


pinfo<-read.table(paste(Data_dir,"plant_info.txt",sep="/"),header=T)
head(pinfo)
which(duplicated(pinfo$Sp))
levels(pinfo$Duration)
str(pinfo)
unique(pinfo$Duration)


c17<-read.table(paste(Data_dir,"2017_count.txt",sep="/"),header=T)
head(c17[,1:10],3)
c17$QUAD_ID
length(c17$QUAD_ID)
which(duplicated(c17$QUAD_ID))
c17sp<-colnames(c17)[6:length(c17)]
dim(c17)
c17sp[order(c17sp)]

#Nov 14th - updated count and cover sheets
Data_dir<-"Data/Final_datasets/copies to edit"
Data_dir
dir(Data_dir)

ct17<-read.table(paste(Data_dir,"count_2017.txt",sep="/"),header=T)
head(ct17[,1:10])
which(duplicated(colnames(ct17))) #check for duplicated  sp names
nchar(colnames(ct17)) #no of characters - 7

pinfo<-read.table(paste(Data_dir,"plant_info.txt",sep="/"),header=T)
head(pinfo,3)

which(duplicated(pinfo$Sp)) 
levels(pinfo$Duration) #
str(pinfo)
pinfo$Duration<-as.factor(pinfo$Duration) #turning ch variable into a factor
colnames(ct17)[5:ncol(ct17)] #col names of all plants
#to check theure all in pinfo
colnames(ct17[,5:ncol(ct17)]) 
colnames(ct17[,5:ncol(ct17)])%in% pinfo$Sp
which(!colnames(ct17[,5:ncol(ct17)])%in% pinfo$Sp) 
colnames(ct17)[which(!colnames(ct17[,5:ncol(ct17)])%in% pinfo$Sp)] #wihch sp was not (!) in there
head(ct17[,1:10])

ct17d<-ct17[,5:ncol(ct17)]
head(ct17d[,1:10])
grep(LETTERS[23:25],ct17d[,1])
LETTERS[23:25]
length(which(ct17d[,1]%in%LETTERS [23:25])) #which has letters in them
apply(ct17d,2,function(x)length(which(x%in%LETTERS [23:25]))) #which sp. has letters like W and need to be replaced






#old stuff:

pn<-read.table("Data/Plant_info.txt",header=T)
s17<-read.table("Data/Survey_data_observed_sp_2017.txt",header=F)
s18<-read.table("Data/Survey_data_observed_sp_2018.txt",header=F)
s19<-read.table("Data/Survey_data_observed_sp_2019.txt",header=F)
d1<-read.table("Data/DAT.txt",header=T)

head(pn); dim(pn)
head(s17[,1:10],8); dim(s17)
head(s18[,1:10],8); dim(s18)
head(s19[,1:10],8); dim(s19)
head(c17[,1:10],8); dim(c17)
head(c18[,1:10],8); dim(c18)
head(c19[,1:10],8); dim(c19)


head(count17[,1:6])
dir("Data")

pdat<-read.table("Data/Plant_info.txt", header=T)
head(pdat)

cnames<-colnames(count17[,which(colnames(count17)=="Air_sp"):ncol(count17)])
head(cnames)
length(cnames)

head(pdat, 3)
head(count17[,1:6])

# you want this line to give you nothing:
which(!cnames %in% pdat$Sp)

which(pdat$Sp=="Air_sp")
dim(count17)
311+6

head(count17[,1:10],3)
count17[,length(count17)]

colSums(count17[,which(colnames(count17)=="Air_sp"):ncol(count17)])




