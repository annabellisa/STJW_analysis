# SIMPLIFY MODEL: if the interaction between treatment and reserve is not significant, remove:

if(which(rownames(anova(m1))=="Treatment:reserve")>0.05){
  form.thisrun<-paste(resp.thisrun,"~Treatment+DATE+reserve+Treatment:DATE+(1|PLOT_ID)", sep="")
  m1<-lmer(formula = form.thisrun, data=data.set)
  summary(m1)
  anova(m1)
}

m1_coef<-coef.ext(m1)
m1_anova<-anova.ext(m1)
coef.out[[i]]<-m1_coef
anova.out[[i]]<-m1_anova

p1<-pred(model=m1, new.data = nd1,se.fit = T, type = "response")
preds.out[[i]]<-p1

# close models