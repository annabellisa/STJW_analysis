

# These are the species with problems:
mismatch.sp<-sp.testout[sp.testout$lineup=="no",]

chk17<-read.table("00_Data/check17.txt",header=F)
chk18<-read.table("00_Data/check18.txt",header=F)
chk19<-read.table("00_Data/check19.txt",header=F)

# Update species code and check data for each:
sp.test<-"Lom_bra"

# Get count and cover data for mismatched species where they DO NOT line up:
xdat<-ct_dat[,c(1:4,which(colnames(ct_dat)==sp.test))]
colnames(xdat)[ncol(xdat)]<-"count_data"
xdat<-cbind(xdat,cv_dat[,which(colnames(cv_dat)==sp.test)])
colnames(xdat)[ncol(xdat)]<-"cover_data"
head(xdat)
mm.dat<-xdat[which(!xdat$count_data==xdat$cover_data),]
mm.dat$QUAD_ID<-paste(mm.dat$PLOT_ID, mm.dat$Treatment, sep="")
mm.quads<-unique(as.character(mm.dat$QUAD_ID))
mm.dat

# For each QUAD, the first row is cover, the second count
rahead(chk17,6,6)
full.name<-pinfo[pinfo$Sp==sp.test,]$Species
df.now<-data.frame(t(chk17[c(1:3,which(chk17$V1==full.name)),]))
colnames(df.now)<-df.now[1,]
df.now<-df.now[2:nrow(df.now),]

mm.dat
df.now[which(df.now$QUAD_ID %in% mm.quads),]
head(df.now)

# Conclusion # 1: cover is more wrong than count, but there are problems with the count data as well. 

# Gal_div, 2017:
# J3C count correct, cover wrong
# J4C count correct, cover wrong
# J4A count correct, cover wrong
# M3C count correct, cover wrong
# M3A count correct, cover wrong
# M3B count correct, cover wrong
# M4C count correct, cover wrong
# M7C count correct, cover wrong
# M7A count correct, cover wrong

# Gal_sp., 2017:
# J3C count correct (3), cover wrong (should be 5)
# "J4C" "J4A" "M3C" "M3A" "M3B" "M4C" "M7C" "M7A" count correct (all zeros), cover wrong (all 5, should be zero)

# Lom_fil and Lom_cor switched?
# What's happening with Lom_bra?
# Lom_bra = 2018, M2, M6; 2019, M6
# Lom_fil = almost everything mismatched
# Lom_cor = almost everything mismatched

# Ryt_sp2 and Ryt_sp4 switched?

# Vit_gre and Vit_gra switched?
# Vit_cun, 2017 and 2019, many rows
# Vit_gre 2017 J2, 2018 M1, J6
# Vit_gra 2017 J2, 2018 M1, J6

# Wah_com 2017, 2018 many rows
# Wah_lut, all years, many rows
# Wah_sp., all years, many rows
# Wah_sp2 2017 J5, J6, J7, J8

# Wur_dio, 2017, 2018, 2019, Mulangarri many rows
# Zor_dic 2017 J7, J8
