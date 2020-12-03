
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











