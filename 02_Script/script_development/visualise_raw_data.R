
#  VISUALISE raw data:    	# ----

head(gdf); dim(gdf)

rahead(rich,4,7); dim(rich)
rahead(shan,4,7); dim(shan)

# visualise raw data:

# year effects on species richness:

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(4,4,1,1), mgp=c(2.5,1,0))

for(i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  resp.dat<-rich[,which(colnames(rich)==resp.thisrun)]
  ylab.thisrun<-gdf$ylab[i]
  
  boxplot(resp.dat~rich$reserve+as.factor(rich$DATE), cex.axis=1, col=c("darkturquoise","darkolivegreen2"), ylab=ylab.thisrun, xlab="", las=1, xaxt="n")
  axis(side=1, at=c(1.5,3.5,5.5), labels=c(2017, 2018, 2019), cex.axis=0.8)
  arrows(c(2.5,4.5),0,c(2.5,4.5),40, length=0, col="grey70")
}
plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
legend(1,9,legend=c("Jerrabomberra","Mulangari"), col=c("darkturquoise","darkolivegreen2"), pch=15, bty="n", pt.cex = 3)

# treatement effects on species richness:

dev.new(width=11.69,height=8.27,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(5,6), mar=c(4,4,1,1), mgp=c(2.5,1,0))

for(i in 1:nrow(gdf)){
  
  resp.thisrun<-gdf$group[i]
  resp.dat<-rich[,which(colnames(rich)==resp.thisrun)]
  ylab.thisrun<-gdf$ylab[i]
  
  boxplot(resp.dat~rich$reserve+as.factor(rich$Treatment), cex.axis=1, col=c("darkturquoise","darkolivegreen2"), ylab=ylab.thisrun, xlab="", las=1, xaxt="n")
  axis(side=1, at=c(1.5,3.5,5.5), labels=c("C", "A", "B"), cex.axis=0.8)
  arrows(c(2.5,4.5),0,c(2.5,4.5),40, length=0, col="grey70")
}
plot(1:10, 1:10, type="n", bty="o", xaxt="n", yaxt="n", xlab="", ylab="")
legend(1,9,legend=c("Jerrabomberra","Mulangari"), col=c("darkturquoise","darkolivegreen2"), pch=15, bty="n", pt.cex = 3)

# close visualise ----
