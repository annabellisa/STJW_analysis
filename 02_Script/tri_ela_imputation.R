#update blanks in tri_ela 2017 count:
head(ct17)
rahead(ct17,3,7)
ct17$Tri_ela
Tri_eladat<-read.table(paste(data_dir,"tri_ela_2017.txt",sep=""),header=T)
head(Tri_eladat)
Tri_eladat$clump
Tri_eladat$reserve<-Tri_eladat$PLOT_ID #take the first character of PLOT_ID eg, just J and M
Tri_eladat$reserve<- substr(Tri_eladat$reserve,1,1)
Tri_eladat$reserve<-as.factor(Tri_eladat$reserve)
plot(Tri_eladat$clump,Tri_eladat$tuft)
cor.test(Tri_eladat$clump[Tri_eladat$clump<20],Tri_eladat$tuft[Tri_eladat$clump<20])
temod1<-lm(clump~tuft+reserve,data=Tri_eladat[Tri_eladat$clump<20,]) 
summary(temod1) #reserve not significant
temod2<-lm(clump~tuft,data=Tri_eladat[Tri_eladat$clump<20,])#remove reserve
summary(temod2)
head(Tri_eladat)
range(Tri_eladat$clump, na.rm = T)
range(Tri_eladat$tuft, na.rm = T)
te.nd<-data.frame(tuft=seq(min(Tri_eladat$tuft, na.rm = T),max(Tri_eladat$tuft, na.rm = T)))

te.pr<-data.frame(te.nd,clump=round(predict(temod2,newdata = te.nd),0))
Tri_eladat$clump
head(te.pr)

rahead(ct17,3,6)
ct17$Tri_ela
head(Tri_eladat)
tuft.topredict<-Tri_eladat$tuft[which(is.na(Tri_eladat$clump))]
head(te.pr)
which(te.pr$tuft %in% tuft.topredict)
relevanttufts<-te.pr[which(te.pr$tuft %in% tuft.topredict),]
relevanttufts$tuft[order(tuft.topredict)]
te.pr
Tri_eladat
which(is.na(Tri_eladat$clump))