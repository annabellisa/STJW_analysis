
#########################################
####  	       FUNCTIONS:    		 ####
#########################################

# author: Annabel Smith

# Drop levels and re-assign rownames to subsetted data frames:
tidy.df<-function(df){
df<-droplevels(df)
rownames(df)<-1:length(df[,1])
return(df)
}

# Randomly check rows of a data frame:
check.rows<-function(df,n=6){df[sample(1:length(df[,1]),n),]}

# head function for big genomic database:
ghead<-function(x) head(x[,1:10])

# head function for Raagini and Annabel:
rahead<-function(x,rows,cols) head(x[1:rows,1:cols])

# p value for lmer model:
lmer.pval<-function(x) 2*(1-pnorm(abs(x))) # x is the t-value from the model summary

# Variance inflation factor for multicollinearity analysis (from https://onlinecourses.science.psu.edu/stat501/node/347):
my_vif<-function(r2){1/(1-r2)}

# Blank plot. plot.new() does this but I don't understand the defaults; I defined this so I know what's where.
blankplot<-function()plot(1:10,1:10,bty="n",type="n",xaxt="n",yaxt="n",xlab="",ylab="")

# Polygon CI function
pg.ci<-function(x,data,x.subset=NULL,colour){

# DEFINITION:
# x, data and x.subset must be in quotes
# x: a vector of data on the x-axis
# data: a data frame which includes x 
# currently, the confidence intervals must be specified as $lci and $uci, need to fix this... 
# update 13th Feb 2017: I haven't fixed the $uci and $lci problem, but I've just added in an if to deal with two different cases that distinguish $uci.resp, etc.
# x.subset: a vector of data for subsetting x. If there no subsets, omit this argument or set it at NULL to plot a single polygon. If there are subsets, enter the colname of data frame x to be used as the subset

# For testing:

# i=3
# x<-as.character(step2spr.plot$x.variable[i])
# data<-as.character(step2spr.plot$data[i])
# x.subset<-as.character(step2spr.plot$x.subset[i])
# colour<-rgb(0,0,0,0.1)
# i=1

if(is.null(x.subset)==T){
xx<-paste(data,"$",x,sep="")
# lci<- paste(data,"$lci",sep="")
# uci<- paste(data,"$uci",sep="")
if(length(grep(".resp",colnames(get(data))))>0) lci<-paste(data,"$lci.resp",sep="") else lci<-paste(data,"$lci",sep="")
if(length(grep(".resp",colnames(get(data))))>0) uci<-paste(data,"$uci.resp",sep="") else uci<-paste(data,"$uci",sep="")

xvec <- c(eval(parse(text=xx)), tail(eval(parse(text=xx)), 1), rev(eval(parse(text=xx))), eval(parse(text=xx))[1])
yvec <- c(eval(parse(text=lci)), tail(eval(parse(text=uci)), 1), rev(eval(parse(text=uci))), eval(parse(text=lci))[1])
polygon(xvec, yvec, col=colour, border=NA)
} # close if no subsets

if(is.null(x.subset)==F){

# Get data and vector that is used for subsetting:
data.withsubset<-get(data)
subset.all<-data.withsubset[,x.subset]

# Specify subs.levs: levels for factors, unique numbers for binary variables, and first and third quartiles for continuous variables

if(is.factor(subset.all)) subs.levs<-levels(subset.all)

if(is.factor(subset.all)==F) {

if(length(unique(subset.all))==2) subs.levs<-unique(subset.all) 

} # close if subset is not a factor

for (i in 1:length(subs.levs)){

sub.thisrun<-subs.levs[i]
x.thisrun<-data.withsubset[which(subset.all==sub.thisrun),x]
# lci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"lci"]
# uci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"uci"]
if(length(grep(".resp",colnames(data.withsubset)))>0) lci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"lci.resp"] else lci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"lci"]
if(length(grep(".resp",colnames(data.withsubset)))>0) uci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"uci.resp"] else uci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"uci"]

xvec <- c(x.thisrun, tail(x.thisrun, 1), rev(x.thisrun), x.thisrun[1])
yvec <- c(lci.thisrun, tail(uci.thisrun, 1), rev(uci.thisrun), lci.thisrun[1])
polygon(xvec, yvec, col=colour, border=NA)

} # close for sub levels i

} # close if subset present

} # close polygon function

# Plot pairwise correlations between multiple terms in a data frame:
plot.pw<-function(data,cor.mat,pw.names){
# data: the data frame with variables to be correlated, IN QUOTES
# cor.mat: the correlation matrix, no quotes
# pw.names: pairs of names of variables to be correlated (result of a combn()), no quotes

# make matrix a data.frame so they can be sorted by correlation
cor.df<-data.frame(var1=pw.names[1,],var2=pw.names[2,],cor=cor.mat[lower.tri(cor.mat)])
cor.df<-cor.df[order(-abs(cor.df$cor)),]
cor.df<-tidy.df(cor.df)

for (i in 1:nrow(cor.df)){

data.thisrun<-cor.df[i,]

n1<-paste(data,"$",data.thisrun$var1,sep="")
n2<-paste(data,"$",data.thisrun$var2,sep="")
r.thisrun<-data.thisrun[,"cor"]

if (length(unique(eval(parse(text=n1))))>2) x.thisrun<-eval(parse(text=n1)) else x.thisrun<-as.factor(eval(parse(text=n1)))

plot(x.thisrun,eval(parse(text=n2)), xlab="", ylab=data.thisrun$var2, bty="l",pch=20, las=1,main="")
title(xlab=data.thisrun$var1,line=2)
par(xpd=NA)
mtext(bquote(bold("r = "~.(as.character(round(r.thisrun,3))))),adj=0, col="red",cex=1.5)

} # close i no. correlations

} # close plot.pw function

# Analyse multicollinearity among many terms in a data frame. This has been simplified from mcl_v2 in CATFORD_E93 to make it more general:
mcl_v3<-function(datasets){
# datasets: vector containing the datasets to be analysed, IN QUOTES

res.list<-list()

for (m in 1:length(datasets)){

data.thisrun<-get(datasets[m])
head(data.thisrun,3)

cnames<-colnames(data.thisrun)
out.mat1<-matrix(data=NA, nrow=length(cnames), ncol=3)

for(i in 1:length(cnames)){
resp.thisrun<-cnames[i]
preds.thisrun<-cnames[-which(cnames==resp.thisrun)]

formula.thisrun<-paste(resp.thisrun,"~",paste(preds.thisrun,collapse="+"), sep="")
r2<-summary(lm(eval(parse(text=formula.thisrun)),data=data.thisrun))$r.squared
VIF<-my_vif(r2)
out.mat1[i,1]<-resp.thisrun
out.mat1[i,2]<-r2
out.mat1[i,3]<-VIF
}
res.r2<-data.frame(out.mat1)
colnames(res.r2)<-c("response","r2","VIF")
res.r2[,2]<-as.numeric(as.character(res.r2[,2]))
res.r2[,3]<-as.numeric(as.character(res.r2[,3]))

res.list[[m]]<-res.r2
}
return(data.frame(do.call(rbind,res.list)))
} # close multicollinearity function

# PREDICT function
# extends predictSE to calculate CIs and dataframe-ise predictions with 'new data'
pred<-function(model,new.data,se.fit=T,type="response"){
  
  if(class(model)!="glmmadmb"){
    pr1<-predictSE(model,new.data,se.fit=se.fit, type=type)
    df1<-data.frame(new.data,fit=pr1$fit,se=pr1$se.fit, lci=pr1$fit-(1.96*pr1$se.fit), uci=pr1$fit+(1.96*pr1$se.fit))
  } # close lme4 models
  
  if(class(model)=="glmmadmb"){
    pr1<-suppressWarnings(predict(model,new.data,se.fit=se.fit, type="link"))
    df1<-data.frame(new.data,fit.link=pr1$fit,se.link=pr1$se.fit, lci.link=pr1$fit-(1.96*pr1$se.fit), uci.link=pr1$fit+(1.96*pr1$se.fit))	
    
    if(summary(model)$link=="log"){
      df1$fit.resp<-round(exp(df1$fit.link),4)
      df1$lci.resp<-round(exp(df1$lci.link),4)
      df1$uci.resp<-round(exp(df1$uci.link),4)
    } # close log
    
    if(summary(model)$link=="logit"){
      df1$fit.resp<-round(invlogit(df1$fit.link),6)
      df1$lci.resp<-round(invlogit(df1$lci.link),6)
      df1$uci.resp<-round(invlogit(df1$uci.link),6)
    } # close logit 
    
  } # close admb models
  
  return(df1)
  
} # close predict function

# Extract coefficient table from model summary
coef.ext<-function(model) {
  
  if(class(model)[1]=="lmerModLmerTest"){
    mod.store<-data.frame(term=dimnames(summary(model)$coefficients[,c(1,2,5)])[[1]],round(summary(model)$coefficients[,c(1,2,5)],4))
    mod.store<-tidy.df(mod.store)
    colnames(mod.store)<-c("term","est","se","P")
    return(mod.store)
  } # close lmerTest
  
  if(class(model)[1]=="lmerMod"){
    mod.store<-data.frame(term=dimnames(summary(model)$coefficients[,c(1,2,3)])[[1]],round(summary(model)$coefficients[,c(1,2,3)],4))
    mod.store<-tidy.df(mod.store)
    mod.store$P<-lmer.pval(mod.store[,4])
    mod.store<-mod.store[,c(1,2,3,5)]
    colnames(mod.store)<-c("term","est","se","P")
    return(mod.store)
  } # close lmer
  
  if(class(model)[1]=="glmmadmb"){
    mod.store<-data.frame(term=dimnames(summary(model)$coefficients[,c(1,2,4)])[[1]],round(summary(model)$coefficients[,c(1,2,4)],4))
    mod.store<-tidy.df(mod.store)
    colnames(mod.store)<-c("term","est","se","P")
    return(mod.store)
  } # close glmmadmb
  
  if(class(model)[1]=="glmerMod"){
    mod.store<-data.frame(term=dimnames(summary(model)$coefficients[,c(1,2,4)])[[1]],round(summary(model)$coefficients[,c(1,2,4)],4))
    mod.store<-tidy.df(mod.store)
    colnames(mod.store)<-c("term","est","se","P")
    return(mod.store)
  } # close glmerMod
  
} # close coef.ext

# Extract anova table from model summary
anova.ext<-function(model) {
  
  if(class(model)[1]=="lmerModLmerTest"){
    mod.store<-data.frame(term=rownames(data.frame(anova(m1))),round(data.frame(anova(m1)),4))
    mod.store<-tidy.df(mod.store)
    colnames(mod.store)<-c("term","sum.sq","mean.sq","NumDF","DenDF","F.value","p")
    return(mod.store)
  } # close lmerTest

} # close anova extract

# Effect size plots from coef table:
efs.plot<-function(coef.table,heading){
  # coef.table is the output from coef.sum
  # heading must be in quotes
  par(mar=c(4,8,2,1), mgp=c(2.6,1,0))
  plot(rev(coef.table$est),1:length(coef.table[,1]), pch=20, xlim=c(min(coef.table$lci.link),max(coef.table$uci.link)), xlab="Effect size", ylab="", yaxt="n", bty="l", main=heading, cex.main=1, font.main=1)
  arrows(rev(coef.table$lci.link),1:length(coef.table[,1]),rev(coef.table$uci.link),1:length(coef.table[,1]),code=0)
  arrows(0,0,0,50,code=0)
  axis(side=2, at=1:1:length(coef.table[,1]), labels=rev(as.character(coef.table$term)),las=1, cex.axis=0.8)
}

# From: http://glmm.wikidot.com/faq
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  if(model$family=="truncpoiss") rp <- summary(model)$residuals else rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# WALD TEST: to calculate p-vals from a factor with three or more levels:
# NOTE THAT q MUST BE THE LENGTH OF THE MATRIX IN QUESTION.
Wald<-function(object,R,q) {
  if (!is.matrix(R)) stop("Restrictions must be a matrix")
  b<-fixef(object) 
  vc<-vcov(object)
  w<-t(R%*%b-q)%*%solve(R%*%vc%*%t(R))%*%(R%*%b-q)
  pw<-1-pchisq(w[1],length(q))
  return(invisible(list(chisq=as.vector(w),pvalue=pw)))
} # close Wald test

# R MATRIX FOR WALD TEST. 
R.mat<-function(model,coef.name){
  if (class(model)=="glmmadmb") {
    coef.tab<-summary(model)$coefficients
    mod.fr<-"model$frame"}
  if (class(model)=="lmerMod") {
    coef.tab<-summary(model)$coefficients
    mod.fr<-"model@frame"}
  if (class(model)=="glmerMod") {
    coef.tab<-summary(model)$coefficients
    mod.fr<-"model@frame"}
  # STEP 1. Determine what type of variable (main or interaction) coef.name is. Get the names of the variables that make up the variable in question. t.names of main effects will = 1 and of interactions will = 2. t.names should never = 3 unless there are three way interactions.
  t.mat<-attr(terms(model),"factors")
  t.names<-names(which(t.mat[,colnames(t.mat)==coef.name]==1))
  
  # STEP 2. Get the levels of the factor so that they can be appended to the coef.name:
  if (length(t.names)==1) lev<-levels(eval(parse(text=paste(mod.fr,"$",coef.name,sep=""))))
  
  if (length(t.names)==2){
    # the classes of the interaction terms:
    int.terms<-data.frame(term=t.names,class=c(class(eval(parse(text=paste(mod.fr,"$",t.names[1],sep="")))),class(eval(parse(text=paste(mod.fr,"$",t.names[2],sep=""))))))
    factor.term<-int.terms$term[which(int.terms$class=="factor")]
    if(length(factor.term)==1) lev<-levels(eval(parse(text=paste(mod.fr,"$",factor.term,sep=""))))
    else stop("more than one factor in interaction term")
    # This line was added to make the function work on Alice's treatment:yr models. NOTE THAT IT MIGHT NOT WORK FOR ALL INTERACTION MODELS:
    lev<-paste(factor.term,lev,":",levels(factor.term)[which(levels(factor.term)!=factor.term)],sep="")
  }
  
  # STEP 3: generate the R matrix using the levels and the length of the levels:
  if (length(t.names)==1){
    R<-matrix(data=0,nrow=length(lev)-1,ncol=length(coef.tab[,1]))
    for(k in 1:(length(lev)-1)){
      R[k,which(rownames(coef.tab)==paste(coef.name,lev[k+1],sep=""))]<-1
    }
  }
  if(length(t.names)==2){
    R<-matrix(data=0,nrow=length(lev)-1,ncol=length(coef.tab[,1]))
    for(k in 1:(length(lev)-1)){
      R[k,which(rownames(coef.tab)==lev[k+1])]<-1
    }
  }
  R<-R
} # close R.mat function








