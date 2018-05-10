##########################
# Code for determining which species distributions overlap
# Works 2015.05.13
##########################

rm(list=ls()) # Clear workspace

##########################
# Import data
##########################

# Necessary packages
library(sp)
library(rgeos)
library(mapdata)
library(maptools)

# Read in the shape file
# Available from IUCN at http://www.iucnredlist.org/technical-documents/spatial-data/
# readShapeSpatial() is in the "sp" package
mammterr=readShapeSpatial("C:/Users/newUser/Documents/Thesis/Infection, Genetics, and Evolution draft/TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp") #shape files for the distributions of every terrestrial mammal in the IUCN database
summary(mammterr)
class(mammterr)

allnames<-c("Myotis daubentonii","Myotis mystacinus","Nyctalus noctula","Myotis myotis",
            "Eptesicus nilssonii","Pipistrellus pipistrellus","Myotis nigricans",
            "Myotis brandtii","Tylonycteris pachypus","Myotis nattereri","Corynorhinus townsendii",
            "Rhinolophus pearsonii","Rhinolophus landeri","Rousettus aegyptiacus",
            "Eidolon helvum","Micropteropus pusillus","Epomophorus gambianus","Pteropus hypomelanus",
            "Ptenochirus jagori","Harpyionycteris whiteheadi","Desmodus rotundus",
            "Carollia perspicillata","Glossophaga soricina","Sturnira lilium",
            "Micronycteris microtis","Artibeus toltecus","Artibeus planirostris","Phyllostomus hastatus",
            "Vampyressa bidens","Artibeus obscurus","Monophyllus redmani","Brachyphylla cavernarum",
            "Artibeus jamaicensis","Macrotus waterhousii","Diaemus youngi",
            "Phyllonycteris poeyi","Lophostoma silvicolum","Artibeus lituratus",
            "Noctilio leporinus","Pteronotus davyi","Pteronotus personatus","Miniopterus natalensis",
            "Miniopterus schreibersii","Triaenops persicus","Hipposideros gigas",
            "Hipposideros diadema","Hipposideros larvatus","Megaderma lyra","Rhinolophus borneensis",
            "Megaerops niphanae","Rhinolophus acuminatus","Rhinolophus sinicus","Coleura afra",
            "Hipposideros armiger","Hipposideros fulvus","Taphozous melanopogon","Chaerephon plicatus",
            "Eidolon dupreanum","Anoura geoffroyi","Sturnira mordax","Vampyressa thyone",
            "Carollia castanea","Platyrrhinus vittatus","Carollia sowelli","Myotis keaysi",
            "Uroderma bilobatum") #the species I care about
all.binomial<-mammterr$binomial

# Which rows of the data are the species I care about
keep=list()
for(i in 1:length(allnames)){
keep[[i]]=which(all.binomial==allnames[i])
}

x=keep[[1]]

for(i in 2:length(allnames)){
x=c(x,keep[[i]])
}

keep=x
myspecies.distr=mammterr[keep,]

# Helpful to output the "keep" numbers to Excel so you can manually inspect polygons and
#     see which species is assigned to each shapefile

##########################
# Checking data for errors
##########################

# See Bat_polygons.R for code to resolve self-intersections and other errors
source("C:/Users/newUser/Documents/Bartonella/Developing manuscripts/Global bats and bartonella/Version4 files/Bat_polygons.R")

##########################
# Species distribution plots
##########################

# Look at a couple of species distribution plots
plot(polygon.sp15,col="red") # Eidolon helvum
plot(polygon.sp17,add=TRUE,col="blue") # Epomophorus gambianus

# Plot all species distributions as transparent layers
greentrans<-rgb(0,100,0,50,maxColorValue=255)
# New World
map('worldHires',xlim=c(-180,-35),ylim=c(-90,90))
plot(myspecies.distr,add=TRUE,col=greentrans)
# Old World
map('worldHires',xlim=c(-35,180),ylim=c(-90,90))
plot(myspecies.distr,add=TRUE,col=greentrans)

##########################
# Overlap and area calculations
##########################

# gArea() and gIntersection() use the "rgeos" package
# over() uses the "sp" package

# Calculate area of polygon
gArea(polygon.sp15) # Eidolon helvum

# Calculate area of overlap between two polygons
gArea(gIntersection(polygon.sp15,polygon.sp17)) # Eidolon helvum and Epomophorus gambianus

# Calculate area of union between two polygons
gArea(gUnion(polygon.sp15,polygon.sp17)) # Eidolon helvum and Epomophorus gambianus

# Gives line segments that overlap
# All NAs means no overlaps
over(polygon.sp15,polygon.sp17) # Eidolon helvum and Epomophorus gambianus

# If they overlap, give a 1, if not 0
ifelse(sum(over(polygon.sp15,polygon.sp17),na.rm=TRUE)>0,1,0) # Eidolon helvum and Epomophorus gambianus

##########################
# Overlap and area loops
##########################

### Overlap loop

# Set up an empty matrix
overlapmatrix=matrix(NA,length(allnames),length(allnames),dimnames=list(allnames, allnames))

# Create an index for the loop
index=expand.grid(1:length(allnames),1:length(allnames))

# Loops through index and calculates whether two shapefiles have lines that overlap
for(i in 1:dim(index)[1]){
if(exists(paste("polygon.sp",index[i,1],sep=""))&exists(paste("polygon.sp",index[i,2],sep=""))){
overlapmatrix[index[i,1],index[i,2]]=
ifelse(sum(over(get(paste("polygon.sp",index[i,1],sep="")),
    get(paste("polygon.sp",index[i,2],sep=""))),na.rm=TRUE)>0,1,0)
    }
}

# Output matrix is sent to Excel sheet
write.csv(overlapmatrix,file="C:/Users/Clif McKee/Desktop/trees3/overlap.csv")

### Area loop

# Start with matrix of 0s and 1s from overlap matrix
areamatrix=overlapmatrix

# Loops through overlap matrix and calculates area of intersection between overlapping shapefiles
for(i in 1:dim(index)[1]){
if(exists(paste("polygon.sp",index[i,1],sep=""))&exists(paste("polygon.sp",index[i,2],sep=""))
   &overlapmatrix[index[i,1],index[i,2]]!=0){
  areamatrix[index[i,1],index[i,2]]=
    gArea(
      gIntersection(get(paste("polygon.sp",index[i,1],sep="")),
                    get(paste("polygon.sp",index[i,2],sep=""))))/
    gArea(
      gUnion(get(paste("polygon.sp",index[i,1],sep="")),
             get(paste("polygon.sp",index[i,2],sep=""))))
  print(areamatrix[index[i,1],index[i,2]]) # Print as loop goes just to make sure it is working
    }
}

# Output matrix is sent to Excel sheet
write.csv(areamatrix,file="C:/Users/Clif McKee/Desktop/trees3/overlaparea.csv")

##########################
# End of code
##########################

# Further modifications to area overlap matrix:
#     1. Calculation of percent overlap - do in Excel manually
#         Creates an asymmetrical matrix that can be imported into R
#     2. Calculate distance matrix - use maximum distance
#         (i.e., minimum percent overlap between two shapefiles)
# Example code
pct.overlap<-as.matrix(read.csv(file.choose(),header=T))
rownames(pct.overlap)=colnames(pct.overlap)
pct.distance=dist(pct.overlap, method = "maximum", diag = TRUE, upper = TR