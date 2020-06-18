# updating script
<<<<<<< HEAD

dir()


=======
# Lacey
>>>>>>> 455b0e3794d4b25501690abfb053f556dae98ea2
dir("Data")


pn<-read.table("Data/Plant_names.txt",header=T)
pnex<-read.table("Data/plant_info_example.txt",header=T)

head(pn); dim(pn)
head(pnex,2); dim(pnex)

pn2<-data.frame(species=pn[-which(duplicated(pn$species)),])
head(pn2); dim(pn2)

pn3<-merge(pn2, pnex, by="species", all.x=T, all.y=F)
head(pn3); dim(pn3)

which(pnex$species %in% pn2$species)

# write.table(pn3, file="plant_info_STJW.txt", row.names=F, quote=F, sep="\t")








