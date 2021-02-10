# Set working directory
setwd("/Volumes/LaCie SSD/Documents/Completed projects/McKee CD et al [2016] Inf Gen Evol/Version4 files")

# Load packages
library(ape)
library(vegan)
library(paco)
library(beeswarm)

#################################
### CORE COPHYLOGENY ANALYSIS ###
#################################

### Can be run using different combinations of host trees and parasite trees

## Main combinations for cophylogeny tests
# ML trees: host = "./Bat/Bat tree (ML_GTR+G+I_4).nwk", parasite = "./Bartonella/Bartonella tree (ML_GTR+G+I_4).nwk"
# Bayesian trees: host = "./Bat/Bat alignment (MAFFT L)_combo-OUT.trees.txt", parasite = "./Bartonella/Bartonella alignment (MAFFT L)_combo-OUT.trees.txt"
# Host-parasite association matrix = "./Host associations.txt"

## Alternative combinations for cophylogeny tests
# ML trees: host = "./Bat alt/Bat alt tree (ML_GTR+G+I_4).nwk", parasite = "./Bartonella/Bartonella tree (ML_GTR+G+I_4).nwk"
# Bayesian trees: host = "./Bat alt/Bat alt alignment (MAFFT L)_combo-OUT.trees.txt", parasite = "./Bartonella/Bartonella alignment (MAFFT L)_combo-OUT.trees.txt"
# Host-parasite association matrix = "./Host associations alt.txt"

### 2. DATA INPUT

# 2.1 Host and parasite phylogenetic data (should be one of the following):
# 2.1.1 Phylogenetic trees:

## Main analysis
# ML host tree
TreeH <- read.tree("./Bat/Bat tree (ML_GTR+G+I_4).nwk") #this function reads Newick trees
# Bayesian host tree
# TreeH <- read.nexus("./Bat/Bat alignment (MAFFT L)_combo-OUT.trees.txt") #this function reads Nexus trees
## Alternative analysis
# Alternative ML host tree
# TreeH <- read.tree("./Bat alt/Bat alt tree (ML_GTR+G+I_4).nwk") #this function reads Newick trees
# Alternative Bayesian host tree
# TreeH <- read.nexus("./Bat alt/Bat alt alignment (MAFFT L)_combo-OUT.trees.txt") #this function reads Nexus trees
plot(TreeH)
TreeHdrop <- drop.tip(TreeH,c("Ornithorhynchus.anatinus","Rattus.rattus","Equus.caballus"))
plot(TreeHdrop)

# ML parasite tree
TreeP <- read.tree("./Bartonella/Bartonella tree (ML_GTR+G+I_4).nwk") #this function reads Newick trees
# Bayesian parasite tree
# TreeP <- read.nexus("./Bartonella/Bartonella alignment (MAFFT L)_combo-OUT.trees.txt") #this function reads Nexus trees
plot(TreeP)
TreePdrop <- drop.tip(TreeP,c("'Brucella.melitensis'","'Rhizobium.leguminosarum'","'Ochrobactrum.anthropi'"))
plot(TreePdrop)

# Compute patristic distances:
host.D <- cophenetic(TreeH)
host.D <- host.D/max(host.D)
para.D <- cophenetic(TreeP)
para.D <- para.D/max(para.D)

# 2.2 Read HP: host-parasite association matrix
# Hosts in rows, parasites in columns. Taxa names are included in the file and should match those in tree, sequence or distance files. 
## Main analysis
HP <- as.matrix(read.table("./Host associations.txt", header=TRUE))
## Alternative analysis
# HP <- as.matrix(read.table("./Host associations alt.txt", header=TRUE))
#Sort host and parasite taxa in distance matrices to match the HP matrix:
host.D <- host.D[rownames(HP),rownames(HP)]
para.D <- para.D[colnames(HP),colnames(HP)]

### 3. APPLY PACo FUNCTION
D <- prepare_paco_data(host.D, para.D, HP)
D <- add_pcoord(D, correction = "cailliez")
D <- PACo(D, seed = 1234, nperm = 100, symmetric = F, method = "r0", shuffled = T)
D <- paco_links(D)
# Plot null versus observed
top.D <- ceiling(max(D$shuffled))
hist(D$shuffled, xlim=c(0, top.D))
abline(v=D$gof$ss, col='red')
# Goodness of fit and p-value
D$gof
# Output residuals and jackknife values for individual links
D.resid = as.data.frame(residuals_paco(D$proc, type='interaction')); colnames(D.resid) = 'resid'
D.jk.median = median(D$jackknife)
D.jk.75th = quantile(D$jackknife, 0.75)
D.jk.IQR = IQR(D$jackknife)
D.jk.upper = D.jk.75th + 1.5*D.jk.IQR
sum(D$jackknife > D.jk.median)
sum(D$jackknife > D.jk.upper)

# Run ParaFit
parafit.out <- parafit(host.D, para.D, HP, nperm = 99, test.links = T,
        seed = 1234, correction = "cailliez", silent = F)
cat(" The observed ParafitGlobal is ", parafit.out$ParaFitGlobal, "\n", "P-value = ", parafit.out$p.global, " based on ", parafit.out$nperm," permutations.")

##############################
### Range overlap analysis ###
##############################

# Input range overlap matrix
# Main analysis
overlap <- as.matrix(read.csv("./overlaparea.csv",header=T))
# Alternative analysis
# overlap <- as.matrix(read.csv("./overlaparea_alt.csv",header=T))

# Modify matrix
rownames(overlap)=colnames(overlap)
overlap=1-overlap
overlap <- overlap[rownames(HP),rownames(HP)]

# Run Mantel test to look at correlation between phylogenetic distance and range overlap
mantel(overlap,host.D, method="pearson",permutations=99)

### 3. APPLY PACo FUNCTION
O <- prepare_paco_data(overlap, para.D, HP)
O <- add_pcoord(O, correction = "cailliez")
O <- PACo(O, seed = 1234, nperm = 100, symmetric = F, method = "r0", shuffled = T)
O <- paco_links(O)
# Plot null versus observed
top.O <- ceiling(max(O$shuffled))
hist(O$shuffled, xlim=c(0, top.O))
abline(v=O$gof$ss, col='red')
# Goodness of fit and p-value
O$gof
# Output residuals and jackknife values for individual links
O.resid = as.data.frame(residuals_paco(O$proc, type='interaction')); colnames(O.resid) = 'resid'
O.jk.median = median(O$jackknife)
O.jk.75th = quantile(O$jackknife, 0.75)
O.jk.IQR = IQR(O$jackknife)
O.jk.upper = O.jk.75th + 1.5*O.jk.IQR
sum(O$jackknife > O.jk.median)
sum(O$jackknife > O.jk.upper)

# Run ParaFit
parafit.out <- parafit(overlap, para.D, HP, nperm = 99, test.links = T,
                       seed = 1234, correction = "cailliez", silent = F)
cat(" The observed ParafitGlobal is ", parafit.out$ParaFitGlobal, "\n", "P-value = ", parafit.out$p.global, " based on ", parafit.out$nperm," permutations.")

# Test the effect of mixing the host phylogenetic matrix and host geographic overlap on cophylogeny fit
superout<-NULL
z=1
for(i in seq(0,1,0.1)){
  super<-((i*host.D)+(1-i)*overlap)
  M <- prepare_paco_data(super, para.D, HP)
  M <- add_pcoord(M, correction = "cailliez")
  M <- PACo(M, seed = 1234, nperm = 1, symmetric = F, method = "r0", shuffled = T)
  p.out <- parafit(super, para.D, HP, nperm = 1, test.links = F,
                   seed = 1234, correction = "cailliez", silent = T)
  
  superout <- rbind(superout, cbind(i,1-i,
                                    M$gof$ss,
                                    p.out$ParaFitGlobal))
  print(z)
  z=z+1
}
colnames(superout)<-c("i","1-i",
                      "m2.WoS",
                      "PFG.WoS")

#################
### Figure S9 ###
#################

#test vector vs. host
resid.vec<-read.csv("./resid.vec.csv",header=T,sep=",")
bats<-resid.vec[ which(resid.vec$Carrier=='bat'), ]
ectos<-resid.vec[ which(resid.vec$Carrier=='ectoparasite'),]

par(mfrow=c(2,2))
#ML
#PACo residuals
boxplot(ML.PACo.resid~Carrier,data=resid.vec,
        xlab="Carrier",ylab="PACo residuals",main="ML.PACo.resid",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$ML.PACo.resid,ectos$ML.PACo.resid,conf.int=T)

#ParaFit F1
boxplot(ML.ParaFit.F1~Carrier,data=resid.vec,
        xlab="Carrier",ylab="ParaFit F1",main="ML.ParaFit.F1",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$ML.ParaFit.F1,ectos$ML.ParaFit.F1,conf.int=T)

#BEAST
#PACo residuals
boxplot(BEAST.PACo.resid~Carrier,data=resid.vec,
        xlab="Carrier",ylab="PACo residuals",main="BEAST.PACo.resid",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$BEAST.PACo.resid,ectos$BEAST.PACo.resid)

#ParaFit F1
boxplot(BEAST.ParaFit.F1~Carrier,data=resid.vec,
        xlab="Carrier",ylab="ParaFit F1",main="BEAST.ParaFit.F1",
        col=c("#0000ff50","#ffff0050"))
wilcox.test(bats$BEAST.ParaFit.F1,ectos$BEAST.ParaFit.F1)

################################################################
### Extra figure - same as Figure S9 but using range overlap ###
################################################################

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

##################
### Figure S10 ###
##################

par(mfrow=c(2,2))
HostFam<-read.csv("./HostFam.csv",header=T,sep=",")
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
