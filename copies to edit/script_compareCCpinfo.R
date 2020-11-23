
#Nov 14th - updated count and cover sheets
Data_dir<-"Data/Final_datasets/copies to edit"
Data_dir
dir(Data_dir)

ct17<-read.table(paste(Data_dir,"count_2017.txt",sep="/"),header=T)
head(ct17[,1:10])
which(duplicated(colnames(ct17))) #check for duplicated  sp names
nchar(colnames(ct17)) #no of characters - 7; do the same for all count and cover sheets

#compare pinfo with count&cover sheets to make sure all sp are present  

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

colnames(ct19)[which(!colnames(ct19[,5:ncol(ct19)])%in% pinfo$Sp)] #wihch sp was not (!) in there
head(ct17[,1:10])

#23rd nov
which(colnames(ct17)=="Dic_cri")
which(!colnames(ct17[,5:ncol(ct17)])%in% pinfo$Sp)
which(pinfo$Sp=="Dic_cri")


ct17d<-ct17[,5:ncol(ct17)]
head(ct17d[,1:10])
grep(LETTERS[23:25],ct17d[,1])
LETTERS[23:25]
length(which(ct17d[,1]%in%LETTERS [23:25])) #which has letters in them
apply(ct17d,2,function(x)length(which(x%in%LETTERS [23:25]))) #which sp. has letters like W and need to be replaced





#Nov 19th - checking count 2018
Data_dir<-"copies to edit/count_data"
Data_dir
dir(Data_dir)
ct18<-read.table(paste(Data_dir,"count_2018.txt",sep="/"),header=T)
ct18<-read.table(paste(Data_dir,"count_2018.txt",sep="/"),header=T)
head(ct18[,1:10])
which(duplicated(colnames(ct18)))
nchar(colnames(ct18))

#compare pinfo with count & cover sheets - 2018 count
pinfo<-read.table(paste(Data_dir,"plant_info.txt",sep="/"),header=T)
head(pinfo,5)
which(duplicated(pinfo$Sp))
str(pinfo)
head(pinfo,3)
head(ct18[,1:10])
colnames(ct18)[10:ncol(ct18)]
colnames(ct18[,2:ncol(ct18)])

colnames(ct18[,5:ncol(ct18)])%in% pinfo$Sp
which(colnames(ct18[,5:ncol(ct18)])%in% pinfo$Sp)
which(!colnames(ct18[,5:ncol(ct18)])%in% pinfo$Sp)
which(!colnames(ct18[,5:ncol(ct18)])%in% pinfo$Sp)
colnames(ct18[,which(!colnames(ct18[,5:ncol(ct18)])%in% pinfo$Sp)])
colnames(ct18)[which(!colnames(ct18[,5:ncol(ct18)])%in% pinfo$Sp)]
head(ct18[,1:10])
ct18d<-ct18[,1:ncol(ct18)]
head(ct18d)
head(ct18d[,1:10])
ct18d<-ct18[.5:ncol(ct18d)]
ct18d
head(ct18d[,5:10])





#checking 2019 count
Data_dir<-"copies to edit/count_data"
Data_dir
dir(Data_dir)
ct19<-read.table(paste(Data_dir, "count_2019.txt", sep "/"), header=T)
ct19<-read.table(paste(Data_dir,"count_2019.txt",sep="/"),header=T)
head(ct19[,1:10])
which(duplicated(colnames(ct19)))
nchar(colnames(ct19))

#compare pinfo with count2019
Data_dir<-"copies to edit/plant_info"
dir(Data_dir)
pinfo<-read.table(paste(Data_dir,"plant_info.txt",sep="/"),header=T)
head(pinfo,5)
which(duplicated(pinfo$Sp))
str(pinfo)
pinfo$Duration<-as.factor(levels(pinfo$Duration)) #turning ch variable into a factor
pinfo$Duration<-as.factor(pinfo$Duration)
str(pinfo)
levels(pinfo$Duration)
head(pinfo,3)
head(ct19[,1:10])
colnames(ct19)[5:ncol(ct19)]
colnames(ct19[,5:ncol(ct19)])
colnames(ct19[,5:ncol(ct19)]) %in% pinfo$Sp
which(colnames(ct19[,5:ncol(ct19)]) %in% pinfo$Sp)
which(!colnames(ct19[,5:ncol(ct19)]) %in% pinfo$Sp)
colnames(ct19[,which(!colnames(ct19[,5:ncol(ct19)])%in% pinfo$Sp)])
colnames(ct19)[which(!colnames(ct19[,5:ncol(ct19)])%in% pinfo$Sp)]
head(ct19)[.1:10]
#tried changing Dic_cri; running again to check 
colnames(ct19)[which(!colnames(ct19[,5:ncol(ct19)])%in% pinfo$Sp)]
colnames(ct17)[which(!colnames(ct17[,5:ncol(ct17)])%in% pinfo$Sp)]


head(ct19[,1:10])
ct19d<-ct19[,1:ncol(ct19)]
head(ct19d)
head(ct19d[,1:10])





#checking cover2017
Data_dir<-"copies to edit/cover_data"
Data_dir
dir(Data_dir)
cv17<-read.table(paste(Data_dir,"cover_2017.txt",sep="/"),header=T)
head(cv17[,1:10])
which(duplicated(colnames(cv17)))
nchar(colnames(cv17))

#compare pinfo with cover2017
Data_dir<-"copies to edit/plant_info"
dir(Data_dir)
pinfo<-read.table(paste(Data_dir,"plant_info.txt",sep="/"),header=T)
head(pinfo,5)
which(duplicated(pinfo$Sp))
str(pinfo)
pinfo$Duration<-as.factor(levels(pinfo$Duration))

#turning ch variable into a factor
pinfo$Duration<-as.factor(pinfo$Duration)
str(pinfo)
levels(pinfo$Duration)
head(pinfo,3)
head(cv17[,1:10])
colnames(cv17)[6:ncol(cv17)]
colnames(cv17[,5:ncol(cv17)])
colnames(cv17[,5:ncol(cv17)]) %in% pinfo$Sp
colnames(cv17[,7:ncol(cv17)]) %in% pinfo$Sp
which(colnames(cv17[,5:ncol(cv17)]) %in% pinfo$Sp)
which(!colnames(ct19[,5:ncol(ct19)]) %in% pinfo$Sp)
colnames(cv17[,which(!colnames(cv17[,5:ncol(cv17)])%in% pinfo$Sp)])
colnames(cv17)[which(!colnames(cv17[,5:ncol(cv17)])%in% pinfo$Sp)]
head(ct19)[.1:10]
#col names correspond b/w pinfo and cover2017

head(cv17[,1:10])
cv17d<-cv17[,1:ncol(cv17)]
head(cv17)
head(cv17d[,1:10])





#checking cover 2018
Data_dir<-"copies to edit/cover_data"
Data_dir
dir(Data_dir)
cv18<-read.table(paste(Data_dir,"cover_2018.txt",sep="/"),header=T)
head(cv18[,1:10])
which(duplicated(colnames(cv18)))
nchar(colnames(cv18))

#compare pinfo with cover2018
Data_dir<-"copies to edit/plant_info"
dir(Data_dir)
pinfo<-read.table(paste(Data_dir, "plant_info.txt", sep="/"), header=T)
head(pinfo,5)
which(duplicated(pinfo$Sp))
str(pinfo)
#changing ch variable into factor
pinfo$Duration<-as.factor(pinfo$Duration)
str(pinfo)
levels(pinfo$Duration)
head(pinfo,3)
head(cv18[,6:12])
colnames(cv18[,6:ncol(cv18)])
colnames(cv18[,5:ncol(cv18)])
colnames(cv18[,6:ncol(cv18)]) %in% pinfo$Sp
colnames(cv17[,5:ncol(cv17)]) %in% pinfo$Sp
colnames(cv18[,6:ncol(cv18)]) %in% pinfo$Sp
colnames(cv18[,4:ncol(cv18)]) %in% pinfo$Sp
colnames(cv18[,6:ncol(cv18)]) %in% pinfo$Sp
#col names correspond b/w pinfo and cover2018





#checking Cover2019
Data_dir<-"copies to edit/cover_data"
dir(Data_dir)
cv19<-read.table(paste(Data_dir,"cover_2019.txt",sep="/"), header=T)
head(cv19[,6:ncol(cv19)])
which(duplicated(colnames(cv19)))
nchar(colnames(cv19))

#compare pinfo with cover2019
Data_dir<-"copies to edit/plant_info"
dir(Data_dir)
pinfo<-read.table(paste(Data_dir,"plant_info.txt",sep="/"), header=T)
head(pinfo,5)
which(duplicated(pinfo$Sp))
str(pinfo)
pinfo$Duration<-as.factor(pinfo$Duration)
str(pinfo)
levels(pinfo$Duration)
head(pinfo,3)
head(cv19,4)
colnames(cv19[,6:ncol(cv19)])
colnames(cv19[,6:ncol(cv19)]) %in% pinfo$Sp
colnames(cv19[,1:ncol(cv19)]) %in% pinfo$Sp #sp. names start only after column 6
#col names correspond b/w pinfo and cover 2019




#23d nov

ct17d<-ct17[,5:ncol(ct17)]
head(ct17d[,1:10])
grep(LETTERS[23:25],ct17d[,1])
LETTERS[23:25]
length(which(ct17d[,1]%in%LETTERS [23:25])) #which has letters in them
apply(ct17d,2,function(x)length(which(x%in%LETTERS [23:25]))) #which sp. has letters like W and need to be replaced
ct17d$Air_sp.
replace(ct17d$Air_sp.,list = c("W","X","Y"),values=c(35, 75, 100))
replace(ct17d$Air_sp., lapply())
repdf<-data.frame(oldval=c("W","X","Y"), newval=c(35, 75, 100))
repdf        
replace(ct17d$Air_sp., which(repdf$oldval==