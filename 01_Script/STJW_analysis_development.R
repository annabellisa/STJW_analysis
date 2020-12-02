# updating script

dir()
dir("Data")


pn<-read.table("Data/Plant_names.txt",header=T)
pnex<-read.table("Data/plant_info_example.txt",header=T)

head(pn); dim(pn)
head(pnex,2); dim(pnex)

#old script

# Lacey is Emily's dogs name




pn2<-data.frame(species=pn[-which(duplicated(pn$species)),])
head(pn2); dim(pn2)

pn3<-merge(pn2, pnex, by="species", all.x=T, all.y=F)
head(pn3); dim(pn3)

which(pnex$species %in% pn2$species)

# write.table(pn3, file="plant_info_STJW.txt", row.names=F, quote=F, sep="\t")



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










