# updating script

dir()
dir("Data")


pn<-read.table("Data/Plant_info.txt",header=T)
s17<-read.table("Data/Survey_data_observed_sp_2017.txt",header=F)
s18<-read.table("Data/Survey_data_observed_sp_2018.txt",header=F)
s19<-read.table("Data/Survey_data_observed_sp_2019.txt",header=F)
c17<-read.table("Data/Count_and_cover_sheet_2017.txt",header=F)
c18<-read.table("Data/Count_and_cover_sheet_2018.txt",header=F)
c19<-read.table("Data/Count_and_cover_sheet_2019.txt",header=F)
count17<-read.table("Data/Not_fixed_count_sheet_2017.txt",header=T)
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
