
# exotic leg herbs richness: (repeat for exotic_legherb diversity, and for richness c4 and sigA)
# change data set and the response variable:
elh_mod<-lmer(exotic_legherb~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID), data=rich_sc)
summary(elh_mod)
elh_pr<-predictSE(elh_mod, nd1, se.fit = T)
elh_pr
elh_pr<-data.frame(nd1,fit=elh_pr$fit,se=elh_pr$se.fit)
elh_pr$lci<-elh_pr$fit-(elh_pr$se*1.96)
elh_pr$uci<-elh_pr$fit+(elh_pr$se*1.96)
head(elh_pr)

xofs<-0.2
arrowlgth<-0.02

dev.new(width=10,height=5,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(1,2), mar=c(2,4,4,1), mgp=c(2.5,1,0))

plot(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="C"]-xofs,elh_pr$fit[elh_pr$reserve=="J" & elh_pr$Treatment=="C"], pch=15, ylim=c(min(elh_pr$lci), max(elh_pr$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
arrows(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="C"]-xofs,elh_pr$lci[elh_pr$reserve=="J" & elh_pr$Treatment=="C"],elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="C"]-xofs,elh_pr$uci[elh_pr$reserve=="J" & elh_pr$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="A"],elh_pr$fit[elh_pr$reserve=="J" & elh_pr$Treatment=="A"], pch=15, col="red")
arrows(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="A"],elh_pr$lci[elh_pr$reserve=="J" & elh_pr$Treatment=="A"],elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="A"],elh_pr$uci[elh_pr$reserve=="J" & elh_pr$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="B"]+xofs,elh_pr$fit[elh_pr$reserve=="J" & elh_pr$Treatment=="B"], pch=15, col="blue")
arrows(elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="B"]+xofs,elh_pr$lci[elh_pr$reserve=="J" & elh_pr$Treatment=="B"],elh_pr$DATE[elh_pr$reserve=="J" & elh_pr$Treatment=="B"]+xofs,elh_pr$uci[elh_pr$reserve=="J" & elh_pr$Treatment=="B"], code=3, angle=90, length=arrowlgth)

legend(1.5,4,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)
mtext("Jerra",side=3, line=1)

# mulungarri:
plot(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="C"]-xofs,elh_pr$fit[elh_pr$reserve=="M" & elh_pr$Treatment=="C"], pch=15, ylim=c(min(elh_pr$lci), max(elh_pr$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
arrows(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="C"]-xofs,elh_pr$lci[elh_pr$reserve=="M" & elh_pr$Treatment=="C"],elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="C"]-xofs,elh_pr$uci[elh_pr$reserve=="M" & elh_pr$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="A"],elh_pr$fit[elh_pr$reserve=="M" & elh_pr$Treatment=="A"], pch=15, col="red")
arrows(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="A"],elh_pr$lci[elh_pr$reserve=="M" & elh_pr$Treatment=="A"],elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="A"],elh_pr$uci[elh_pr$reserve=="M" & elh_pr$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="B"]+xofs,elh_pr$fit[elh_pr$reserve=="M" & elh_pr$Treatment=="B"], pch=15, col="blue")
arrows(elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="B"]+xofs,elh_pr$lci[elh_pr$reserve=="M" & elh_pr$Treatment=="B"],elh_pr$DATE[elh_pr$reserve=="M" & elh_pr$Treatment=="B"]+xofs,elh_pr$uci[elh_pr$reserve=="M" & elh_pr$Treatment=="B"], code=3, angle=90, length=arrowlgth)
mtext("Mulangarri",side=3, line=1)
mtext("exotic_leg_forb richness", side=3, line=3)

#exotic leg_herb diversity:
# change data set and the response variable:
elhd_mod<-lmer(exotic_legherb~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID), data=shan_sc)
summary(elhd_mod)
elhd_pr<-predictSE(elhd_mod, nd1, se.fit = T)
elhd_pr
elhd_pr<-data.frame(nd1,fit=elhd_pr$fit,se=elhd_pr$se.fit)
elhd_pr$lci<-elhd_pr$fit-(elhd_pr$se*1.96)
elhd_pr$uci<-elhd_pr$fit+(elhd_pr$se*1.96)
head(elhd_pr)

xofs<-0.2
arrowlgth<-0.02

dev.new(width=10,height=5,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(1,2), mar=c(2,4,4,1), mgp=c(2.5,1,0))

plot(elhd_pr$DATE[elhd_pr$reserve=="J" & elhd_pr$Treatment=="C"]-xofs,elhd_pr$fit[elhd_pr$reserve=="J" & elhd_pr$Treatment=="C"], pch=15, ylim=c(min(elhd_pr$lci), max(elhd_pr$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
arrows(elhd_pr$DATE[elhd_pr$reserve=="J" & elhd_pr$Treatment=="C"]-xofs,elhd_pr$lci[elhd_pr$reserve=="J" & elhd_pr$Treatment=="C"],elhd_pr$DATE[elhd_pr$reserve=="J" & elhd_pr$Treatment=="C"]-xofs,elhd_pr$uci[elhd_pr$reserve=="J" & elh_pr$Treatment=="C"], code=3, angle=90, length=arrowlgth)
points(elhd_pr$DATE[elhd_pr$reserve=="J" & elhd_pr$Treatment=="A"],elhd_pr$fit[elhd_pr$reserve=="J" & elhd_pr$Treatment=="A"], pch=15, col="red")
arrows(elhd_pr$DATE[elhd_pr$reserve=="J" & elhd_pr$Treatment=="A"],elhd_pr$lci[elhd_pr$reserve=="J" & elhd_pr$Treatment=="A"],elhd_pr$DATE[elhd_pr$reserve=="J" & elhd_pr$Treatment=="A"],elhd_pr$uci[elhd_pr$reserve=="J" & elhd_pr$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(elhd_pr$DATE[elhd_pr$reserve=="J" & elhd_pr$Treatment=="B"]+xofs,elhd_pr$fit[elhd_pr$reserve=="J" & elhd_pr$Treatment=="B"], pch=15, col="blue")
arrows(elhd_pr$DATE[elhd_pr$reserve=="J" & elhd_pr$Treatment=="B"]+xofs,elhd_pr$lci[elhd_pr$reserve=="J" & elhd_pr$Treatment=="B"],elhd_pr$DATE[elhd_pr$reserve=="J" & elhd_pr$Treatment=="B"]+xofs,elhd_pr$uci[elhd_pr$reserve=="J" & elhd_pr$Treatment=="B"], code=3, angle=90, length=arrowlgth)

legend(1.5,4,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)
mtext("Jerra",side=3, line=1)

# mulungarri:
plot(elhd_pr$DATE[elhd_pr$reserve=="M" & elhd_pr$Treatment=="C"]-xofs,elhd_pr$fit[elhd_pr$reserve=="M" & elhd_pr$Treatment=="C"], pch=15, ylim=c(min(elhd_pr$lci), max(elhd_pr$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
arrows(elhd_pr$DATE[elhd_pr$reserve=="M" & elhd_pr$Treatment=="C"]-xofs,elhd_pr$lci[elhd_pr$reserve=="M" & elhd_pr$Treatment=="C"],elhd_pr$DATE[elhd_pr$reserve=="M" & elhd_pr$Treatment=="C"]-xofs,elhd_pr$uci[elhd_pr$reserve=="M" & elhd_pr$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(elhd_pr$DATE[elhd_pr$reserve=="M" & elhd_pr$Treatment=="A"],elhd_pr$fit[elhd_pr$reserve=="M" & elhd_pr$Treatment=="A"], pch=15, col="red")
arrows(elhd_pr$DATE[elhd_pr$reserve=="M" & elhd_pr$Treatment=="A"],elhd_pr$lci[elhd_pr$reserve=="M" & elhd_pr$Treatment=="A"],elhd_pr$DATE[elhd_pr$reserve=="M" & elhd_pr$Treatment=="A"],elhd_pr$uci[elhd_pr$reserve=="M" & elhd_pr$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(elhd_pr$DATE[elhd_pr$reserve=="M" & elhd_pr$Treatment=="B"]+xofs,elhd_pr$fit[elhd_pr$reserve=="M" & elhd_pr$Treatment=="B"], pch=15, col="blue")
arrows(elhd_pr$DATE[elhd_pr$reserve=="M" & elhd_pr$Treatment=="B"]+xofs,elhd_pr$lci[elhd_pr$reserve=="M" & elhd_pr$Treatment=="B"],elhd_pr$DATE[elhd_pr$reserve=="M" & elhd_pr$Treatment=="B"]+xofs,elhd_pr$uci[elhd_pr$reserve=="M" & elhd_pr$Treatment=="B"], code=3, angle=90, length=arrowlgth)
mtext("Mulangarri",side=3, line=1)
mtext("exotic_leg_forb diversity", side=3, line=3)


#native c4 grass richness:
nc4_mod<-lmer(native_c4~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID), data=rich_sc)
summary(nc4_mod)
nc4_pr<-predictSE(nc4_mod, nd1, se.fit = T)
nc4_pr
nc4_pr<-data.frame(nd1,fit=nc4_pr$fit,se=nc4_pr$se.fit)
nc4_pr$lci<-nc4_pr$fit-(nc4_pr$se*1.96)
nc4_pr$uci<-nc4_pr$fit+(nc4_pr$se*1.96)
head(nc4_pr)

xofs<-0.2
arrowlgth<-0.02

dev.new(width=10,height=5,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(1,2), mar=c(2,4,4,1), mgp=c(2.5,1,0))

plot(nc4_pr$DATE[nc4_pr$reserve=="J" & nc4_pr$Treatment=="C"]-xofs,nc4_pr$fit[nc4_pr$reserve=="J" & nc4_pr$Treatment=="C"], pch=15, ylim=c(min(nc4_pr$lci), max(nc4_pr$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
arrows(nc4_pr$DATE[nc4_pr$reserve=="J" & nc4_pr$Treatment=="C"]-xofs,nc4_pr$lci[nc4_pr$reserve=="J" & nc4_pr$Treatment=="C"],nc4_pr$DATE[nc4_pr$reserve=="J" & nc4_pr$Treatment=="C"]-xofs,nc4_pr$uci[nc4_pr$reserve=="J" & nc4_pr$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(nc4_pr$DATE[nc4_pr$reserve=="J" & nc4_pr$Treatment=="A"],nc4_pr$fit[nc4_pr$reserve=="J" & nc4_pr$Treatment=="A"], pch=15, col="red")
arrows(nc4_pr$DATE[nc4_pr$reserve=="J" & nc4_pr$Treatment=="A"],nc4_pr$lci[nc4_pr$reserve=="J" & nc4_pr$Treatment=="A"],nc4_pr$DATE[nc4_pr$reserve=="J" & nc4_pr$Treatment=="A"],nc4_pr$uci[nc4_pr$reserve=="J" & nc4_pr$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(nc4_pr$DATE[nc4_pr$reserve=="J" & nc4_pr$Treatment=="B"]+xofs,nc4_pr$fit[nc4_pr$reserve=="J" & nc4_pr$Treatment=="B"], pch=15, col="blue")
arrows(nc4_pr$DATE[nc4_pr$reserve=="J" & nc4_pr$Treatment=="B"]+xofs,nc4_pr$lci[nc4_pr$reserve=="J" & nc4_pr$Treatment=="B"],nc4_pr$DATE[nc4_pr$reserve=="J" & nc4_pr$Treatment=="B"]+xofs,nc4_pr$uci[nc4_pr$reserve=="J" & nc4_pr$Treatment=="B"], code=3, angle=90, length=arrowlgth)

legend(1.5,4,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)
mtext("Jerra",side=3, line=1)

# mulungarri:
plot(nc4_pr$DATE[nc4_pr$reserve=="M" & nc4_pr$Treatment=="C"]-xofs,nc4_pr$fit[nc4_pr$reserve=="M" & nc4_pr$Treatment=="C"], pch=15, ylim=c(min(nc4_pr$lci), max(nc4_pr$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
arrows(nc4_pr$DATE[nc4_pr$reserve=="M" & nc4_pr$Treatment=="C"]-xofs,nc4_pr$lci[nc4_pr$reserve=="M" & nc4_pr$Treatment=="C"],nc4_pr$DATE[nc4_pr$reserve=="M" & nc4_pr$Treatment=="C"]-xofs,nc4_pr$uci[nc4_pr$reserve=="M" & nc4_pr$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(nc4_pr$DATE[nc4_pr$reserve=="M" & nc4_pr$Treatment=="A"],nc4_pr$fit[nc4_pr$reserve=="M" & nc4_pr$Treatment=="A"], pch=15, col="red")
arrows(nc4_pr$DATE[nc4_pr$reserve=="M" & nc4_pr$Treatment=="A"],nc4_pr$lci[nc4_pr$reserve=="M" & nc4_pr$Treatment=="A"],nc4_pr$DATE[nc4_pr$reserve=="M" & nc4_pr$Treatment=="A"],nc4_pr$uci[nc4_pr$reserve=="M" & nc4_pr$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(nc4_pr$DATE[nc4_pr$reserve=="M" & nc4_pr$Treatment=="B"]+xofs,nc4_pr$fit[nc4_pr$reserve=="M" & nc4_pr$Treatment=="B"], pch=15, col="blue")
arrows(nc4_pr$DATE[nc4_pr$reserve=="M" & nc4_pr$Treatment=="B"]+xofs,nc4_pr$lci[nc4_pr$reserve=="M" & nc4_pr$Treatment=="B"],nc4_pr$DATE[nc4_pr$reserve=="M" & nc4_pr$Treatment=="B"]+xofs,nc4_pr$uci[nc4_pr$reserve=="M" & nc4_pr$Treatment=="B"], code=3, angle=90, length=arrowlgth)
mtext("Mulangarri",side=3, line=1)
mtext("native_c4 richness", side=3, line=3)


#significance A richness
siga_mod<-lmer(exotic_legherb~Treatment+DATE+reserve+Treatment:DATE+Treatment:DATE:reserve+(1|PLOT_ID), data=rich_sc)
summary(siga_mod)
siga_pr<-predictSE(siga_mod, nd1, se.fit = T)
siga_pr
siga_pr<-data.frame(nd1,fit=siga_pr$fit,se=siga_pr$se.fit)
siga_pr$lci<-siga_pr$fit-(siga_pr$se*1.96)
siga_pr$uci<-siga_pr$fit+(siga_pr$se*1.96)
head(siga_pr)

xofs<-0.2
arrowlgth<-0.02

dev.new(width=10,height=5,noRStudioGD = T,dpi=80, pointsize=12)
par(mfrow=c(1,2), mar=c(2,4,4,1), mgp=c(2.5,1,0))

plot(siga_pr$DATE[siga_pr$reserve=="J" & siga_pr$Treatment=="C"]-xofs,siga_pr$fit[siga_pr$reserve=="J" & siga_pr$Treatment=="C"], pch=15, ylim=c(min(siga_pr$lci), max(siga_pr$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
arrows(siga_pr$DATE[siga_pr$reserve=="J" & siga_pr$Treatment=="C"]-xofs,siga_pr$lci[siga_pr$reserve=="J" & siga_pr$Treatment=="C"],siga_pr$DATE[siga_pr$reserve=="J" & siga_pr$Treatment=="C"]-xofs,siga_pr$uci[siga_pr$reserve=="J" & siga_pr$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(siga_pr$DATE[siga_pr$reserve=="J" & siga_pr$Treatment=="A"],siga_pr$fit[siga_pr$reserve=="J" & siga_pr$Treatment=="A"], pch=15, col="red")
arrows(siga_pr$DATE[siga_pr$reserve=="J" & siga_pr$Treatment=="A"],siga_pr$lci[siga_pr$reserve=="J" & siga_pr$Treatment=="A"],siga_pr$DATE[siga_pr$reserve=="J" & siga_pr$Treatment=="A"],siga_pr$uci[siga_pr$reserve=="J" & siga_pr$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(siga_pr$DATE[siga_pr$reserve=="J" & siga_pr$Treatment=="B"]+xofs,siga_pr$fit[siga_pr$reserve=="J" & siga_pr$Treatment=="B"], pch=15, col="blue")
arrows(siga_pr$DATE[siga_pr$reserve=="J" & siga_pr$Treatment=="B"]+xofs,siga_pr$lci[siga_pr$reserve=="J" & siga_pr$Treatment=="B"],siga_pr$DATE[siga_pr$reserve=="J" & siga_pr$Treatment=="B"]+xofs,siga_pr$uci[siga_pr$reserve=="J" & siga_pr$Treatment=="B"], code=3, angle=90, length=arrowlgth)

legend(1.5,4,legend=c("Control","Spot spray","Boom spray"), col=c("black","red","blue"), pch=15, bty="n", pt.cex = 3)
mtext("Jerra",side=3, line=1)

# mulungarri:
plot(siga_pr$DATE[siga_pr$reserve=="M" & siga_pr$Treatment=="C"]-xofs,siga_pr$fit[siga_pr$reserve=="M" & siga_pr$Treatment=="C"], pch=15, ylim=c(min(siga_pr$lci), max(siga_pr$uci)), xlim=c(-0.3,2.3), xaxt="n", xlab="", ylab=ylab.thisrun, las=1)
axis(side = 1, at=c(0,1,2), labels=c(2017,2018,2019))
arrows(siga_pr$DATE[siga_pr$reserve=="M" & siga_pr$Treatment=="C"]-xofs,siga_pr$lci[siga_pr$reserve=="M" & siga_pr$Treatment=="C"],siga_pr$DATE[siga_pr$reserve=="M" & siga_pr$Treatment=="C"]-xofs,siga_pr$uci[siga_pr$reserve=="M" & siga_pr$Treatment=="C"], code=3, angle=90, length=arrowlgth)

points(siga_pr$DATE[siga_pr$reserve=="M" & siga_pr$Treatment=="A"],siga_pr$fit[siga_pr$reserve=="M" & siga_pr$Treatment=="A"], pch=15, col="red")
arrows(siga_pr$DATE[siga_pr$reserve=="M" & siga_pr$Treatment=="A"],siga_pr$lci[siga_pr$reserve=="M" & siga_pr$Treatment=="A"],siga_pr$DATE[siga_pr$reserve=="M" & siga_pr$Treatment=="A"],siga_pr$uci[siga_pr$reserve=="M" & siga_pr$Treatment=="A"], code=3, angle=90, length=arrowlgth, col="red")

points(siga_pr$DATE[siga_pr$reserve=="M" & siga_pr$Treatment=="B"]+xofs,siga_pr$fit[siga_pr$reserve=="M" & siga_pr$Treatment=="B"], pch=15, col="blue")
arrows(siga_pr$DATE[siga_pr$reserve=="M" & siga_pr$Treatment=="B"]+xofs,siga_pr$lci[siga_pr$reserve=="M" & siga_pr$Treatment=="B"],siga_pr$DATE[siga_pr$reserve=="M" & siga_pr$Treatment=="B"]+xofs,siga_pr$uci[siga_pr$reserve=="M" & siga_pr$Treatment=="B"], code=3, angle=90, length=arrowlgth)
mtext("Mulangarri",side=3, line=1)
mtext("significance A_ richness", side=3, line=3)




## THREE WAY INTERACTION:
head(gdf); dim(gdf)
rahead(rich_sc,4,7); dim(rich_sc)
rahead(shan_sc,4,7); dim(shan_sc)

# Native leg. forb, exotic leg forb

resp.thisrun<-gdf$group[26]
form.threeway<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+Treatment:reserve+reserve:DATE+(1|PLOT_ID)", sep="")

nfmod2<-lmer(form.threeway,data=rich_sc)
summary(nfmod2)
anova(nfmod2)

form.simple<-paste(resp.thisrun,"~Treatment*DATE+reserve+(1|PLOT_ID)", sep="")
mod.simp<-lmer(form.simple,data=rich_sc)
summary(mod.simp)
anova(mod.simp); AIC(mod.simp)
