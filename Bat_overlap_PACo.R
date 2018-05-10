overlap<-as.matrix(read.csv(file.choose(),header=T))
rownames(overlap)=colnames(overlap)
overlap=1-overlap
#overlap=as.matrix(dist(overlap, method = "maximum", diag = TRUE, upper = TRUE))

overlap <- overlap[rownames(host.D),rownames(host.D)]

mantel(overlap,host.D, method="pearson",permutations=10000)

HP <- as.matrix(read.table(file.choose(), header=TRUE))

overlap <- overlap[rownames(HP),rownames(HP)]
host.D <- host.D[rownames(HP), rownames(HP)]
para.D <- para.D[colnames(HP),colnames(HP)]

mantel.partial(host.D,overlap,para.D, method="pearson",permutations=10000)
partial.mantel.test(host.D,overlap,para.D, resamp = 1000, method = "pearson", quiet = FALSE)

# 
### 3. APPLY PACo FUNCTION  
PACo.fit <- PACo(host.D, overlap, HP)
HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo) #Procrustes Ordination 
NLinks = sum(HP) #Number of H-P links; needed for further computations
#
#3.1 Plot of host and parasite ordination:
#par(mfrow=c(2,2))
HostX <- HP.proc$X #host ordination matrix
ParY <- HP.proc$Yrot #parasite ordination matrix, scaled and rotated to fit HostX
#Plotting host and parasite ordinations
plot(HostX, asp=1, pch=1) 
points(ParY, pch=2)
arrows(ParY[,1], ParY[,2],HostX[,1], HostX[,2], length=0.12, angle=15, xpd=FALSE)
HostX <- unique(HP.proc$X) 
ParY <- unique(HP.proc$Yrot) #unique() removes duplicated points - convenient for labelling of points below
identify(ParY[,1], ParY[,2], rownames(ParY), offset=1, xpd=TRUE, cex=0.8) #interactive labelling
identify(HostX[,1], HostX[,2], rownames(HostX),offset=1, xpd=TRUE, cex= 0.8)
#write.csv(HostX,file="C:/Users/Clif McKee/Desktop/trees3/HostX.csv")
#write.csv(ParY,file="C:/Users/Clif McKee/Desktop/trees3/ParY.csv")
HostX<-HostX[rownames(HP),]
ParY<-ParY[colnames(HP),]
HostFam<-read.csv("C:/Users/Clif McKee/Desktop/trees4/HostFam.csv",header=T,sep=",")
ParFam<-read.csv("C:/Users/Clif McKee/Desktop/trees4/ParFam.csv",header=T,sep=",")
rownames(HostFam)=HostFam$X;HostFam<-HostFam[-1];HostFam<-HostFam[rownames(HP),]
rownames(ParFam)=ParFam$X;ParFam<-ParFam[-1];ParFam<-ParFam[colnames(HP),]
ordihull(HostX, HostFam$Family, label=T, lty = 1, col = "red")
ordihull(HostX, HostFam$Suborder, label=T, lty = 1, col = "blue")
ordihull(HostX, HostFam$Region, label=T, lty = 1, col = "darkgreen")
ordihull(ParY, ParFam$Family, label=T, lty = 1, col = "red")
ordihull(ParY, ParFam$Suborder, label=T, lty = 1, col = "blue")
ordihull(ParY, ParFam$Region, label=T, lty = 1, col = "darkgreen")
ordiellipse(HostX, HostFam$Family,kind="se",conf=0.95,col="red")

#
#3.2 Goodness-of-fit-test
m2.obs <- HP.proc$ss #observed sum of squares
N.perm = 10000 #set number of permutations for testing
P.value = 1
seed <-.Random.seed[trunc(runif(1,1,626))]
set.seed(seed)
#set.seed(5) ### use this option to obtain reproducible randomizations
for (n in c(1:N.perm))
{ 
  if (NLinks <= nrow(HP) | NLinks <= ncol(HP)) 	#control statement to avoid all parasites beig associated to a single host 
  {	flag2 <- TRUE 
    while (flag2 == TRUE)	{ 
      HP.perm <- t(apply(HP,1,sample))
      if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE else flag2 <- FALSE
    }  
  } else { HP.perm <- t(apply(HP,1,sample))} #permutes each HP row independently
  PACo.perm <- PACo(overlap, para.D, HP.perm)
  m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss #randomized sum of squares
  #write (m2.perm, file = "D:/m2_perm.txt", sep ="\t", append =TRUE) #option to save m2 from each permutation
  if (m2.perm <= m2.obs)
  {P.value = P.value + 1}
  print(P.value)
}
P.value <- P.value/N.perm
cat(" The observed m2 is ", m2.obs, "\n", "P-value = ", P.value, " based on ", N.perm," permutations.")

#3.3 Contribution of individual links
HP.ones <- which(HP > 0, arr.in=TRUE)
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)# empty matrix of jackknifed squared residuals
colnames (SQres.jackn) <- paste(rownames(HP.proc$X),rownames(HP.proc$Yrot), sep="-") #colnames identify the H-P link
t.critical = qt(0.975,NLinks-1) #Needed to compute 95% confidence intervals.
for(i in c(1:NLinks)) #PACo setting the ith link = 0
{HP.ind <- HP
 HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
 PACo.ind <- PACo(overlap, para.D, HP.ind)
 Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo) 
 res.Proc.ind <- c(residuals(Proc.ind))
 res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
 SQres.jackn [i, ] <- res.Proc.ind	#Append residuals to matrix of jackknifed squared residuals
} 
SQres.jackn <- SQres.jackn**2 #Jackknifed residuals are squared
SQres <- (residuals (HP.proc))**2 # Vector of original square residuals
#jackknife calculations:
SQres.jackn <- SQres.jackn*(-(NLinks-1))
SQres <- SQres*NLinks
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres)) #apply jackknife function to matrix
phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE) #mean jackknife estimate per link
phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE) #standard deviation of estimates
phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks) #upper 95% confidence interval
#barplot of squared jackknifed residuals
par(mfrow=c(1,1))
pat.bar <- barplot(names.arg = " ",phi.mean, space = 0.25, col="white",
                   xlab= "Host-parasite links", ylab= "Squared residuals",cex.lab=1.2)
#text(pat.bar, par("usr")[3] - 0.001, srt = 330, adj = 0, labels = colnames(SQres.jackn), xpd = TRUE, font = 1, cex=0.6)
#arrows(pat.bar, phi.mean, pat.bar, phi.UCI, length= 0.05, angle=90)
#abline(a=median(phi.mean), b=0, lty=2, xpd=FALSE) #draws a line across the median residual value
abline(a=median(phi.mean), b=0, lty=2, xpd=FALSE) #draws a line across the median residual value
IQR<-quantile(phi.mean)[4]-quantile(phi.mean)[2]
#abline(a=mean(phi.mean)+0.5*IQR, b=0, lty=1, col="darkgreen", xpd=FALSE)
abline(a=quantile(phi.mean)[4]+1*IQR, b=0, lty=1, xpd=FALSE)
abline(a=quantile(phi.mean)[4]+1.5*IQR, b=0, lty=1, col="red", xpd=FALSE)
phi.stat<-cbind(phi.mean,phi.UCI)
write.csv(phi.mean,file="C:/Users/Clif McKee/Desktop/trees4/phi.mean_overlap_alt_ML.csv")
    
parafit.out<-parafit(overlap, para.D, HP, nperm = 10000, test.links = T,
                      seed = NULL, correction = "cailliez", silent = FALSE)
cat(" The observed ParafitGlobal is ", parafit.out$ParaFitGlobal, "\n", "P-value = ", parafit.out$p.global, " based on ", parafit.out$nperm," permutations.")
write.csv(parafit.out[3],file="C:/Users/Clif McKee/Desktop/trees4/parafit.out_overlap_alt_ML.csv")

bias_Links<-as.matrix(read.csv(file.choose(),header=T))
rownames(bias_Links)=colnames(bias_Links)
Links <- bias_Links[rownames(HP),rownames(HP)]

superout<-NULL
z=1
for(i in seq(0,1,0.1)){
  for(j in seq(0,1,0.1)){
    super<-((i*host.D)+(1-i)*overlap)
    PACo.fit <- PACo(super, para.D, HP)
    HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit2$P.PCo)
    parafit.out<-parafit(super, para.D, HP, nperm = 1, test.links = FALSE,
                          seed = NULL, correction = "cailliez", silent = TRUE)
    
    superout<-rbind(superout,cbind(i,1-i,j,
                                   HP.proc$ss,
                                   parafit.out$ParaFitGlobal))
    print(z)
    z=z+1
  }
}
colnames(superout)<-c("i","1-i","j",
                      "m2.WoS","m2.SampleSize","m2.Links","m2.LinksSample",
                      "PFG.WoS","PFG.SampleSize","PFG.Links","PFG.LinksSample")
write.csv(superout,file="C:/Users/Clif McKee/Desktop/trees4/superout.csv")

library(ggplot2)
superout<-as.data.frame(superout)
plot1<-qplot(i,1-i,data=superout,colour=PACo.m2)
plot1+scale_colour_gradient(low="red",high="white")

plot2<-qplot(i,1-i,data=superout,colour=ParaFitGlobal)
plot2+scale_colour_gradient(low="red",high="white")