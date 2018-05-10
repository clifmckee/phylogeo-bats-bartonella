setwd("~/Bartonella/Developing manuscripts/Global bats and bartonella/Version4 files")
library (ape)
library(vegan)
							### 1. PACo FUNCTION: adjustemt prior to Procrustes analysis
PACo <- function (H.dist, P.dist, HP.bin)
{ 
HP.bin <- which(HP.bin > 0, arr.in=TRUE)
H.PCo <- pcoa(H.dist, correction="cailliez")$vectors #Performs PCo of Host distances 
P.PCo <- pcoa(P.dist, correction="cailliez")$vectors #Performs PCo of Parasite distances
H.PCo <- H.PCo[HP.bin[,1],] #adjust Host PCo vectors 
P.PCo <- P.PCo[HP.bin[,2],]  ##adjust Parasite PCo vectors
list (H.PCo = H.PCo, P.PCo = P.PCo)
}
							### 2. DATA INPUT
	#2.1 Host and parasite phylogenetic data (should be one of the following):
		#2.1.1 Phylogenetic trees:
#names<-as.matrix(read.csv("C:/Users/Clif McKee/Desktop/trees4/HostDnames.csv"), header=T)
TreeH <- read.tree(file.choose()) #this function reads Newick trees
TreeH <- read.nexus(file.choose()) #this function reads Nexus trees
#TreeH$tip.label<-colnames(names)
plot(TreeH)
TreeHdrop <- drop.tip(TreeH,c("Ornithorhynchus.anatinus","Rattus.rattus","Equus.caballus"))
plot(TreeHdrop)

TreeP <- read.tree(file.choose()) #this function reads Newick trees
TreeP <- read.nexus(file.choose()) #this function reads Nexus trees
plot(TreeP)
#TreePdrop <- drop.tip(TreeP,c("Brucella.melitensis","Rhizobium.leguminosarum","Ochrobactrum.anthropi"))
TreePdrop <- drop.tip(TreeP,c("'Brucella.melitensis'","'Rhizobium.leguminosarum'","'Ochrobactrum.anthropi'"))
plot(TreePdrop)
            #Compute patristic distances:
host.D <- cophenetic (TreeH)
host.D<-host.D/max(host.D)
para.D <- cophenetic (TreeP)
para.D<-para.D/max(para.D)

#
		#2.1.2 Aligned sequences
seqH <- read.dna(file.choose(),format="fasta")
seqP <- read.dna(file.choose(),format="fasta")
            # Compute distance matrices from sequence data
host.D <- dist.dna(seqH, model = "F84", as.matrix=TRUE) #Uses the Felsenstein (1984)

														#evolutionary model  
para.D <- dist.dna(seqP, model = "F84", as.matrix=TRUE) #(similar to HKY85).

#
		#2.1.3 Distance matrices
host.D <- as.matrix(read.table(file.choose(), header=TRUE))
para.D <- as.matrix(read.table(file.choose(), header=TRUE))
#			
	#2.2 ## Read HP: host-parasite association matrix
		#Hosts in rows, parasites in columns. Taxa names are included in the file and should match those in tree, sequence or distance files. 
HP <- as.matrix(read.table(file.choose(), header=TRUE))
        #Sort host and parasite taxa in distance matrices to match the HP matrix:
host.D <- host.D[rownames(HP),rownames(HP)]
para.D <- para.D[colnames(HP),colnames(HP)]
# 
						### 3. APPLY PACo FUNCTION  
PACo.fit <- PACo(host.D, para.D, HP)
HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo) #Procrustes Ordination 
NLinks = sum(HP) #Number of H-P links; needed for further computations
#
	#3.1 Plot of host and parasite ordination:
HostX <- HP.proc$X #host ordination matrix
ParY <- HP.proc$Yrot #parasite ordination matrix, scaled and rotated to fit HostX
        #Plotting host and parasite ordinations
#par(mfrow=c(2,2))
plot(HostX, asp=1, pch=1,xlim=c(-0.4,0.2),ylim=c(-0.5,0.3))
#plot(ParY, pch=2,xlim=c(-0.2,0.2),ylim=c(-0.25,0.15),xlab="Axis.1",ylab="Axis.2")
#points(ParY, pch=2)
#arrows(ParY[,1], ParY[,2],HostX[,1], HostX[,2], length=0.12, angle=15, xpd=FALSE)
HostX <- unique(HP.proc$X) 
#ParY <- unique(HP.proc$Yrot) #unique() removes duplicated points - convenient for labelling of points below
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

#Hulls
ordihull(HostX, HostFam$Family, lty = 1, label=T, col = "red")
ordihull(HostX, HostFam$Suborder, lty = 1, label=T, col = "blue")
ordihull(HostX, HostFam$Region, lty = 1, label=T, col = "darkgreen")
ordihull(ParY, ParFam$Family, lty = 1, label=T, col = "red")
ordihull(ParY, ParFam$Suborder, lty = 1, label=T, col = "blue")
ordihull(ParY, ParFam$Region, lty = 1, label=T, col = c("darkgreen"))

HostEig<-apply(HostX,2,sd)^2
HostPVE<-HostEig/sum(HostEig)
ParEig<-apply(ParY,2,sd)^2
ParPVE<-ParEig/sum(ParEig)

#
	#3.2 Goodness-of-fit-test
m2.obs <- HP.proc$ss #observed sum of squares
N.perm = 100 #set number of permutations for testing
P.value = 0.05
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
		PACo.perm <- PACo(host.D, para.D, HP.perm)
		m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss #randomized sum of squares
		if (m2.perm <= m2.obs)
		{P.value = P.value + 1} 
	print(P.value)
	}
P.value <- P.value/N.perm
cat(" The observed m2 is ", m2.obs, "\n", "P-value = ", P.value, " based on ", N.perm," permutations.")
#
	#3.3 Contribution of individual links
HP.ones <- which(HP > 0, arr.in=TRUE)
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)# empty matrix of jackknifed squared residuals
  colnames (SQres.jackn) <- paste(rownames(HP.proc$X),rownames(HP.proc$Yrot), sep="-") #colnames identify the H-P link
t.critical = qt(0.975,NLinks-1) #Needed to compute 95% confidence intervals.
for(i in c(1:NLinks)) #PACo setting the ith link = 0
{HP.ind <- HP
HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
PACo.ind <- PACo(host.D, para.D, HP.ind)
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
write.csv(phi.mean,file="C:/Users/Clif McKee/Desktop/trees4/phi.mean_alt_ML.csv")
#
## end of code ##

parafit.out<-parafit(host.D, para.D, HP, nperm = 10000, test.links = T,
        seed = NULL, correction = "cailliez", silent = FALSE)
cat(" The observed ParafitGlobal is ", parafit.out$ParaFitGlobal, "\n", "P-value = ", parafit.out$p.global, " based on ", parafit.out$nperm," permutations.")
write.csv(parafit.out[3],file="C:/Users/Clif McKee/Desktop/trees4/parafit.out_alt_ML.csv")


#draw host associations
par(mfrow=c(2,2))
assoc <- as.matrix(read.table(file.choose(), header=TRUE))
cophyloplot(TreeHdrop,TreePdrop,assoc,type = "phylogram",space=50,show.tip.label=F,rotate=T,col="red")

#test vector vs. host
library(beeswarm)
resid.vec<-read.csv(file.choose(),header=T,sep=",")
bats<-resid.vec[ which(resid.vec$Carrier=='bat'), ]
ectos<-resid.vec[ which(resid.vec$Carrier=='ectoparasite'),]

par(mfrow=c(2,2))
#ML
#PACo residuals
boxplot(ML.PACo.resid~Carrier,data=resid.vec,
        xlab="Carrier",ylab="PACo residuals",main="ML.PACo.resid",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$ML.PACo.resid,ectos$ML.PACo.resid,conf.int=T)
# ks.test(bats$ML.PACo.resid,ectos$ML.PACo.resid,conf.int=T)

#ParaFit F1
boxplot(ML.ParaFit.F1~Carrier,data=resid.vec,
        xlab="Carrier",ylab="ParaFit F1",main="ML.ParaFit.F1",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$ML.ParaFit.F1,ectos$ML.ParaFit.F1,conf.int=T)
# ks.test(bats$ML.ParaFit.F1,ectos$ML.ParaFit.F1,conf.int=T)

#BEAST
#PACo residuals
boxplot(BEAST.PACo.resid~Carrier,data=resid.vec,
        xlab="Carrier",ylab="PACo residuals",main="BEAST.PACo.resid",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$BEAST.PACo.resid,ectos$BEAST.PACo.resid)
# ks.test(bats$BEAST.PACo.resid,ectos$BEAST.PACo.resid)


#ParaFit F1
boxplot(BEAST.ParaFit.F1~Carrier,data=resid.vec,
        xlab="Carrier",ylab="ParaFit F1",main="BEAST.ParaFit.F1",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$BEAST.ParaFit.F1,ectos$BEAST.ParaFit.F1)
# ks.test(bats$BEAST.ParaFit.F1,ectos$BEAST.ParaFit.F1)

par(mfrow=c(2,2))
#ML overlap
#PACo residuals
boxplot(ML.overlap.PACo.resid~Carrier,data=resid.vec,
        xlab="Carrier",ylab="PACo residuals",main="ML.overlap.PACo.resid",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$ML.overlap.PACo.resid,ectos$ML.overlap.PACo.resid)
# ks.test(bats$ML.overlap.PACo.resid,ectos$ML.overlap.PACo.resid)

#ParaFit F1
boxplot(ML.overlap.ParaFit.F1~Carrier,data=resid.vec,
        xlab="Carrier",ylab="ParaFit F1",main="ML.overlap.ParaFit.F1",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$ML.overlap.ParaFit.F1,ectos$ML.overlap.ParaFit.F1)
# ks.test(bats$ML.overlap.ParaFit.F1,ectos$ML.overlap.ParaFit.F1)

#BEAST overlap
#PACo residuals
boxplot(BEAST.overlap.PACo.resid~Carrier,data=resid.vec,
        xlab="Carrier",ylab="PACo residuals",main="BEAST.overlap.PACo.resid",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$BEAST.overlap.PACo.resid,ectos$BEAST.overlap.PACo.resid)
# ks.test(bats$BEAST.overlap.PACo.resid,ectos$BEAST.overlap.PACo.resid)

#ParaFit F1
boxplot(BEAST.overlap.ParaFit.F1~Carrier,data=resid.vec,
        xlab="Carrier",ylab="ParaFit F1",main="BEAST.overlap.ParaFit.F1",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$BEAST.overlap.ParaFit.F1,ectos$BEAST.overlap.ParaFit.F1)
# ks.test(bats$BEAST.overlap.ParaFit.F1,ectos$BEAST.overlap.ParaFit.F1)

par(mfrow=c(2,2))
HostFam<-read.csv("C:/Users/Clif McKee/Desktop/trees4/HostFam.csv",header=T,sep=",")
fit1<-lm(log(Links)~log(WoS),data=HostFam);summary(fit1)
cor.test(log(HostFam$Links),log(HostFam$WoS))
plot(log(HostFam$WoS),log(HostFam$Links),
     xlab="log(Web of Science articles)",
     ylab="log(Bat-Bartonella links)",
     main="a")
abline(lm(fit1), col="red", lwd=2)

HostFam.noEido<-HostFam[ which(HostFam$X != "Eidolon.helvum"),]
fit2<-lm(log(Links)~log(WoS),data=HostFam.noEido);summary(fit2)
cor.test(log(HostFam.noEido$Links),log(HostFam.noEido$WoS))
plot(log(HostFam.noEido$WoS),log(HostFam.noEido$Links),
     xlab="log(Web of Science articles)",
     ylab="log(Bat-Bartonella links)",
     main="b")
abline(lm(fit2), col="red", lwd=2)

HostFam.noMyot<-HostFam[ which(HostFam$X != "Myotis.myotis"),]
fit3<-lm(log(Links)~log(Sample),data=HostFam.noMyot);summary(fit3)
cor.test(log(HostFam.noMyot$Links),log(HostFam.noMyot$Sample))
plot(log(HostFam.noMyot$Sample),log(HostFam.noMyot$Links),
     xlab="log(Study sample size)",
     ylab="log(Bat-Bartonella links)",
     main="c")
abline(lm(fit3), col="red", lwd=2)

HostFam.end<-HostFam.noMyot[ which(HostFam.noMyot$X != "Eidolon.helvum"),]
fit4<-lm(log(Links)~log(Sample),data=HostFam.end);summary(fit4)
cor.test(log(HostFam.end$Links),log(HostFam.end$Sample))
plot(log(HostFam.end$Sample),log(HostFam.end$Links),
     xlab="log(Study sample size)",
     ylab="log(Bat-Bartonella links)",
     main="d")
abline(lm(fit4), col="red", lwd=2)