
# Exotic, for example, has a total of 11 zeros and 21 ones in the count data. For shannon's the ones become zeros, making 32 zeros. So data sets with lots of ones have many zeros. This is the point of a diversity estimate. A plot with a single species has no diversity
length(which(shan$exotic==0))-length(which(rich$exotic==0))
length(which(simp$exotic==0))-length(which(rich$exotic==0))

exotic_dat<-data.frame(rich=rich$exotic,shan=shan$exotic,simp=simp$exotic,invsimp=invsimp$exotic)
# shan and simp both have zero when richness==1
# invsimp has one when richness==1 (need to replace the Inf values with zeros for when richness is zero)

zero_one<-data.frame(rich=rich$exotic[c(which(rich$exotic==0),which(rich$exotic==1))],shan=shan$exotic[c(which(rich$exotic==0),which(rich$exotic==1))],simp=simp$exotic[c(which(rich$exotic==0),which(rich$exotic==1))],invsimp=invsimp$exotic[c(which(rich$exotic==0),which(rich$exotic==1))])

