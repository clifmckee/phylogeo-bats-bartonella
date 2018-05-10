##########################
# Code for combining shapefiles for each species
# Works: 2015.05.13
##########################

# Necessary packages
# gUnaryUnion() and gSimplify() are in the "rgeos" package
library(rgeos)

##########################
# Instructions:
# 1. Start out by using the gUnaryUnion() function on all species polygons.
#     This will collapse all polygons attributed to a species to just a single polygon.
# 2. Run loop A below to check all shapefiles for self-intersections.
#     These create holes in the polygons that interfere with calculating intersections.
# 3. Examine which shapefiles produce self-intersections.
# 4. Manually change the shapefiles with self-intersections with the gSimplify() function
#     Use tol=0.0001 which is the setting for the rebuffering
##########################

polygon.sp1<-gUnaryUnion(mammterr[13913:14114,])#Myot.daub
polygon.sp2<-gUnaryUnion(mammterr[29195:29256,])#Myot.myst
polygon.sp3<-gUnaryUnion(mammterr[24898:24918,])#Nyct.noct
polygon.sp4<-gUnaryUnion(mammterr[28254:28258,])#Myot.myot
polygon.sp5<-gUnaryUnion(mammterr[19587:19594,])#Epte.nils
polygon.sp6<-gUnaryUnion(mammterr[30014:30050,])#Pipi.pipi
polygon.sp7<-gUnaryUnion(gSimplify(mammterr[9198,],tol=0.0001))#Myot.nigr
polygon.sp8<-gUnaryUnion(mammterr[8500:8517,])#Myot.bran
polygon.sp9<-gUnaryUnion(gSimplify(mammterr[18058:18087,],tol=0.0001))#Tylo.pach
polygon.sp10<-gUnaryUnion(mammterr[33124:33148,])#Myot.natt
polygon.sp11<-gUnaryUnion(gSimplify(mammterr[25233:25235,],tol=0.0001))#Cory.town
polygon.sp12<-gUnaryUnion(gSimplify(mammterr[23284:23289,],tol=0.0001))#Rhin.pear
polygon.sp13<-gUnaryUnion(mammterr[18733:18735,])#Rhin.land
polygon.sp14<-gUnaryUnion(mammterr[17240:17273,])#Rous.aegy
polygon.sp15<-gUnaryUnion(gSimplify(mammterr[14186:14205,],tol=0.0001))#Eido.helv
polygon.sp16<-gUnaryUnion(gSimplify(mammterr[26509,],tol=0.0001))#Mcrp.pusi
polygon.sp17<-gUnaryUnion(mammterr[23545,])#Epom.gamb
polygon.sp18<-gUnaryUnion(gSimplify(mammterr[2030:2119,],tol=0.0001))#Pter.hypo
polygon.sp19<-gUnaryUnion(gSimplify(mammterr[4807:4835,],tol=0.0001))#Pten.jago
polygon.sp20<-gUnaryUnion(gSimplify(mammterr[4517:4536,],tol=0.0001))#Harp.whit
polygon.sp21<-gUnaryUnion(gSimplify(mammterr[31169:31173,],tol=0.0001))#Desm.rotu
polygon.sp22<-gUnaryUnion(gSimplify(mammterr[35157:35160,],tol=0.0001))#Caro.pers
polygon.sp23<-gUnaryUnion(gSimplify(mammterr[5808:5811,],tol=0.0001))#Glos.sori
polygon.sp24<-gUnaryUnion(gSimplify(mammterr[30252:30259,],tol=0.0001))#Stur.lili
polygon.sp25<-gUnaryUnion(gSimplify(mammterr[14261:14267,],tol=0.0001))#Mcny.micr
polygon.sp26<-gUnaryUnion(mammterr[8780:8783,])#Arti.tolt
polygon.sp27<-gUnaryUnion(gSimplify(mammterr[34406,],tol=0.0001))#Arti.plan
polygon.sp28<-gUnaryUnion(gSimplify(mammterr[25790:25791,],tol=0.0001))#Pyst.hast
polygon.sp29<-gUnaryUnion(gSimplify(mammterr[37069,],tol=0.0001))#Vamp.bide
polygon.sp30<-gUnaryUnion(gSimplify(mammterr[3172,],tol=0.0001))#Arti.obsc
polygon.sp31<-gUnaryUnion(gSimplify(mammterr[18798:18806,],tol=0.0001))#Mono.redm
polygon.sp32<-gUnaryUnion(mammterr[2246:2265,])#Brac.cave
polygon.sp33<-gUnaryUnion(gSimplify(mammterr[2138:2157,],tol=0.0001))#Arti.jama
polygon.sp34<-gUnaryUnion(gSimplify(mammterr[29413:29420,],tol=0.0001))#Mcts.wate
polygon.sp35<-gUnaryUnion(gSimplify(mammterr[34932:34935,],tol=0.0001))#Diae.youn
polygon.sp36<-gUnaryUnion(gSimplify(mammterr[17391:17393,],tol=0.0001))#Pyny.poey
polygon.sp37<-gUnaryUnion(gSimplify(mammterr[1536:1537,],tol=0.0001))#Loph.silv
polygon.sp38<-gUnaryUnion(gSimplify(mammterr[12884:12889,],tol=0.0001))#Arti.litu
polygon.sp39<-gUnaryUnion(gSimplify(mammterr[30296:30320,],tol=0.0001))#Noct.lepo
polygon.sp40<-gUnaryUnion(gSimplify(mammterr[28601:28608,],tol=0.0001))#Ptnt.davy
polygon.sp41<-gUnaryUnion(gSimplify(mammterr[25720,],tol=0.0001))#Ptnt.pers
polygon.sp42<-gUnaryUnion(gSimplify(mammterr[9310:9316,],tol=0.0001))#Mnpt.nata
polygon.sp43<-gUnaryUnion(mammterr[25553:25568,])#Mnpt.schr
polygon.sp44<-gUnaryUnion(mammterr[32978:32995,])#Tria.pers
polygon.sp45<-gUnaryUnion(gSimplify(mammterr[9887:9896,],tol=0.0001))#Hipp.giga
polygon.sp46<-gUnaryUnion(gSimplify(mammterr[27137:27226,],tol=0.0001))#Hipp.diad
polygon.sp47<-gUnaryUnion(gSimplify(mammterr[31206:31226,],tol=0.0001))#Hipp.larv
polygon.sp48<-gUnaryUnion(gSimplify(mammterr[12812:12816,],tol=0.0001))#Mgdm.lyra
polygon.sp49<-gUnaryUnion(gSimplify(mammterr[1941:1960,],tol=0.0001))#Rhin.born
polygon.sp50<-gUnaryUnion(gSimplify(mammterr[31456:31457,],tol=0.0001))#Mgps.niph
polygon.sp51<-gUnaryUnion(gSimplify(mammterr[36982:36993,],tol=0.0001))#Rhin.acum
polygon.sp52<-gUnaryUnion(gSimplify(mammterr[18909:18914,],tol=0.0001))#Rhin.sini
polygon.sp53<-gUnaryUnion(mammterr[32165:32176,])#Cole.afra
polygon.sp54<-gUnaryUnion(gSimplify(mammterr[1220:1225,],tol=0.0001))#Hipp.armi
polygon.sp55<-gUnaryUnion(gSimplify(mammterr[14254:14257,],tol=0.0001))#Hipp.fulv
polygon.sp56<-gUnaryUnion(gSimplify(mammterr[6247:6288,],tol=0.0001))#Taph.mela
polygon.sp57<-gUnaryUnion(gSimplify(mammterr[1430:1461,],tol=0.0001))#Chae.plic
polygon.sp58<-gUnaryUnion(gSimplify(mammterr[40977:40978,],tol=0.0001))#Eido.dupr
polygon.sp59<-gUnaryUnion(gSimplify(mammterr[14237:14240,],tol=0.0001))#Anou.geof
polygon.sp60<-gUnaryUnion(mammterr[19109,])#Stur.mord
polygon.sp61<-gUnaryUnion(gSimplify(mammterr[32192:32194,],tol=0.0001))#Vamp.thyo
polygon.sp62<-gUnaryUnion(gSimplify(mammterr[3092,],tol=0.0001))#Caro.cast
polygon.sp63<-gUnaryUnion(mammterr[3421,])#Plat.vitt
polygon.sp64<-gUnaryUnion(gSimplify(mammterr[9878,],tol=0.0001))#Caro.sowe
polygon.sp65<-gUnaryUnion(gSimplify(mammterr[1692:1695,],tol=0.0001))#Myot.keay
polygon.sp66<-gUnaryUnion(gSimplify(mammterr[19814:19815,],tol=0.0001))#Urod.bilo

##########################
# Error checking and centroid loops
##########################

#A. Loop for checking all species shapefiles for self-intersections
for (i in 1:length(allnames)){
  out<-gIsValid(get(paste("polygon.sp",i,sep="")),reason=TRUE,byid=TRUE)
  print(i)
  print(out)
}

centroids<-array(0,dim=c(66,3))
#B. Loop for calculating the centroid of all species shapefiles
for(i in 1:length(allnames)){
  out<-coordinates(get(paste("polygon.sp",i,sep="")))
  print(i)
  print(out)
  centroids[i,1:3]<-c(paste("polygon.sp",i,sep=""),out[1],out[2])
}
write.csv(centroids,"C:/Users/newUser/Documents/Thesis/Infection, Genetics, and Evolution draft/TERRESTRIAL_MAMMALS/centroids.csv")

##########################
# End of code
##########################