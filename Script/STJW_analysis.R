# updating script

dir()
dir("Data")


pn<-read.table("Data/Plant_info.txt",header=T)
s17<-read.table("Data/Survey_data_observed_sp_2017.txt",header=F)

head(pn); dim(pn)
head(s17[,1:10],8); dim(s17)



