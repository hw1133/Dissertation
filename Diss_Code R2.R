#Graphs and Stats for "Complexity in Origins: A Study of Migratory Individuals in Ketton Quarry through d34S Analysis and Statistical Modelling" by Hannah Wilson
#Code below by H. Wilson, Uni of Edinburgh over course of MSc Dissertation
# May 2023 - August 2023 

#set working directory
setwd("~/")

#### IQR rule for outliers
## find Quantile 1 (Q1) and Quantile 3 (Q3) for the data
## 1.5*Q1 = x, 1.5*Q3=y
## any values below x or above y are outliers

#libraries for descriptive and inferential stats
library(car)
library(ggplot2)
library(lattice)
library(FSA)
library(stats)

####read in data: all samples of d34S in database
Data1 <- read.csv("~/Everybody 34S csv.csv") #Appendix G

#create factor levels and models for d34S data
data1.f <- factor(Data1$Site)
data1.f2 <- factor(Data1$Group1)
data1.f3 <- factor(Data1$Country)
data1.f4 <- factor(Data1$Region)
mod1 <- Data1$X34S~data1.f
mod2 <- Data1$X34S~data1.f2
mod3 <- Data1$X34S~data1.f3
mod11 <- Data1$X34S~data1.f4
sub1 <- subset(Data1,Group1 %in% c("1")) #KQ
sub2 <- subset(Data1,Group1 %in% c("2")) #All other sites

#variances of d34and boxplot for KQ versus compilation of all other sites
var(sub1$X34S)
var(sub2$X34S)
boxplot(sub1$X34S, xlab="KQ", ylab="d34S Values", main="d34S Values for KQ", col=c("red"))

#variances of d34S for each site (I know there are faster ways, sorry)
sub10 <- subset(Data1,Site %in% c("Ketton Quarry"))
sub11 <- subset(Data1,Site %in% c("Auldhame"))
sub12 <- subset(Data1,Site %in% c("Stavanger"))
sub13 <- subset(Data1,Site %in% c("Sigtuna"))
sub14 <- subset(Data1,Site %in% c("Bjorned"))
sub15 <- subset(Data1,Site %in% c("Vivallen"))
sub16 <- subset(Data1,Site %in% c("St Olofsholm"))
sub17 <- subset(Data1,Site %in% c("Varnhem"))
sub18 <- subset(Data1,Site %in% c("Hofstadir"))
sub19 <- subset(Data1,Site %in% c("Ostriv"))
sub20 <- subset(Data1,Site %in% c("Baltics"))
sub21 <- subset(Data1,Site %in% c("Eura"))
var(sub10$X34S)
var(sub11$X34S)
var(sub12$X34S)
var(sub13$X34S)
var(sub14$X34S)
var(sub15$X34S)
var(sub16$X34S)
var(sub17$X34S)
var(sub18$X34S)
var(sub19$X34S)
var(sub20$X34S)
var(sub21$X34S)

## descriptive stats of d34S by region
# All regions together not including KQ
sub22 <- subset(Data1,Site %in% c("Eura","Auldhame","Stavanger","Sigtuna","Bjorned","Vivallen","St Olofsholm","Varnhem","Hofstadir","Ostriv","Baltics"))
mean(sub22$X34S)
median(sub22$X34S)
min(sub22$X34S)
max(sub22$X34S)
sd(sub22$X34S)
var(sub22$X34S)
# Britain not including KQ
mean(sub11$X34S) 
median(sub11$X34S)
min(sub11$X34S)
max(sub11$X34S)
sd(sub11$X34S)
var(sub11$X34S)
# Norway
mean(sub12$X34S)
median(sub12$X34S)
min(sub12$X34S)
max(sub12$X34S)
sd(sub12$X34S)
var(sub12$X34S)
# Sweden
sub23 <- subset(Data1,Site %in% c("Sigtuna","Bjorned","Vivallen","St Olofsholm","Varnhem"))
mean(sub23$X34S)
median(sub23$X34S)
min(sub23$X34S)
max(sub23$X34S)
sd(sub23$X34S)
var(sub23$X34S)
# Iceland
mean(sub18$X34S)
median(sub18$X34S)
min(sub18$X34S)
max(sub18$X34S)
sd(sub18$X34S)
var(sub18$X34S)
# Ukraine
mean(sub19$X34S)
median(sub19$X34S)
min(sub19$X34S)
max(sub19$X34S)
sd(sub19$X34S)
var(sub19$X34S)
# Baltics
mean(sub20$X34S)
median(sub20$X34S)
min(sub20$X34S)
max(sub20$X34S)
sd(sub20$X34S)
var(sub20$X34S)
# Finland
mean(sub21$X34S)
median(sub21$X34S)
min(sub21$X34S)
max(sub21$X34S)
sd(sub21$X34S)
var(sub21$X34S)

###plot
plot(mod1, xlab="Site", ylab="34S", cex.axis=0.7)
plot(mod2, cex.axis=0.7, col=c("red","snow3"), xlab="d34S: KQ vs Other Sites", ylab="d34S")
plot(mod3, xlab="KQ and Regions", ylab="d34S", main="Comparison of d34S Values in KQ to Regions of Interest", cex.axis=0.7)
plot(mod11, xlab="KQ and Regions", ylab="d34S", main="Comparison of d34S Values in KQ \nto Regions of Interest \n     ", cex.axis=0.7, col=c("red", "snow3", "snow3", "snow3","snow3","snow3","snow3","snow3"))
histogram(~X34S|Group1, data=Data1, xlab="d34S", main="Distribution of d34S Values \nKetton Quarry vs Other Sites", breaks=seq(from=-20, to=25, by=5))
densityPlot(Data1$X34S, xlab="d34S", main = "Distribution of d34S in all Sites")
qqPlot(mod2, ylab="dd34S")

## inferential stats
kruskal.test(X34S~Site, Data1) #By site
wilcox.test(X34S~Group1, Data1) #KQ vs not KQ
kruskal.test(X34S~Country, Data1) #By country/region
levels(data1.f3)
dunnTest(X34S~data1.f3, Data1, method="bonferroni")

####read in 2nd set of data: d34S only within KQ
Data2 <- read.csv("~/34s KQ only.csv") #Appendix G
Data2

##create factor levels and models for d34S data
data2.f <- factor(Data2$Group) #bone or dentinal Collagen
data2.f2 <- factor(Data2$Outliers) #High outliers, central cluster, low outliers
mod4 <- Data2$d34S~data2.f
mod5 <- Data2$d34S~data2.f2
sub3 <- subset(Data2,Group %in% c("1")) #dentinal collagen
sub4 <- subset(Data2,Group %in% c("2")) #human bone collagen
sub5 <- subset(Data2,Group %in% c("3")) #equid bone collagen
sub6 <- subset(Data2,Group %in% c("2", "3")) #all bone collagen
sub7 <- subset(Data2,Outliers %in% c("1")) # central cluster
sub8 <- subset(Data2,Outliers %in% c("2")) # high outliers
sub9 <- subset(Data2,Outliers %in% c("3")) # low outliers

#descriptive stats
var(sub3$d34S)
mean(sub3$d34S)
median(sub3$d34S)
sd(sub3$d34S)
var(sub4$d34S)
var(sub5$d34S)
var(sub6$d34S)
mean(sub6$d34S)
median(sub6$d34S)
sd(sub6$d34S)
min(sub6$d34S)
max(sub6$d34S)
mean(sub7$d34S)
min(sub7$d34S)
sd(sub7$d34S)
mean(sub9$d34S)
sd(sub9$d34S)

##plot
car::Boxplot(Data2$d34S~data2.f2, main="Majority Range d34S Values vs Outlier Groups", xlab="Groups", ylab="34S Value", col=c("royalblue3","steelblue2", "lightsteelblue1"))
plot(mod4, xlab="Group", ylab="34S Value", main="d34S Values vs Collagen Type", col=c("darkgreen","palegreen3", "darkseagreen1"))
histogram(~d34S|Group, data=Data2)

##inferential stats
kruskal.test(d34S~data2.f, data=Data2)
dunnTest(d34S~data2.f, data=Data2, method="bonferroni")
kruskal.test(d34S~data2.f2, data=Data2)
dunnTest(d34S~data2.f2, data=Data2, method="bonferroni")

#read in data for scatterplot
KQData <- read.csv("~/KQ Scatterplot csv.csv") #Appendix G
KQData

#plot
plot(KQData$FormAge, KQData$X34S, main="KQ d34S Value by Probable Collagen Formation Age", xlab="Probable Formation Age", ylab="d34S (‰)")
abline(a=15.0897, b=-0.3434)
abline(a=-20.5303, b=0.5107)
abline(a=12.9489, b=-0.2412)
abline(a=15.4446, b=-0.3393)
abline(a=-12.7168, b=0.1954)
abline(a=1.577778, b=0, col="red")

####read in data: oxygen comparison data 
dataA <- read.csv("Everybody 18O.csv") #Appendix G
dataA

#create factor levels and models for d18Odw data
dataA.f <- factor(dataA$Site) #by site
dataA.f2 <- factor(dataA$Group) #KQ vs not KQ
dataA.f3 <- factor(dataA$Region) #by region
modA1 <- dataA$d18O~dataA.f
modA2 <- dataA$d18O~dataA.f2
modA3 <- dataA$d18O~dataA.f3
subA1 <- subset(dataA,Group %in% c("1")) #KQ
subA2 <- subset(dataA,Group %in% c("2")) #not KQ

#descriptive stats
var(subA1$d18O)
var(subA2$d18O)

#plot
plot(modA3, xlab="Region", ylab="d18O", main="d18O by Region vs KQ" , cex.axis=0.7, col=c("red", "snow3", "snow3", "snow3","snow3","snow3"))
hist(dataA$d18O, xlab="d18O", main = "Distribution of d18O")
histogram(~d18O|Group, data=dataA, xlab="d18O", main="Distribution of d18O Values \nKetton Quarry vs Other Sites", breaks=seq(from=-18, to=3, by=3))
histogram(~d18O|Region, data=dataA, xlab="d18O", main="Distribution of d18O Values \nby Region", breaks=seq(from=-18, to=3, by=3), )

#inferential stats
leveneTest(dataA$d18O, dataA.f2)
leveneTest(dataA$d18O, dataA.f3)
kruskal.test(d18O~dataA.f2, dataA) # O by KQ/not
kruskal.test(d18O~dataA.f3, dataA) # O by region
dunnTest(d18O~dataA.f3, dataA)

####read in data: strontium comparison data
dataB <- read.csv("Everybody 87Sr.csv") #Appendix G
dataB

#create factor levels and models ofr 87/86Sr data
dataB.f <- factor(dataB$Site)
dataB.f2 <- factor(dataB$Group)
dataB.f3 <- factor(dataB$Region)
modB1 <- dataB$Sr~dataB.f
modB2 <- dataB$Sr~dataB.f2
modB3 <- dataB$Sr~dataB.f3
subB1 <- subset(dataB,Group %in% c("1")) #KQ
subB2 <- subset(dataB,Group %in% c("2")) #not KQ
subB3 <- subset(dataB,Region %in% c("1")) #Britain without KQ
subB4 <- subset(dataB,Region %in% c("2")) #Norway
subB5 <- subset(dataB,Region %in% c("3")) #Sweden
subB6 <- subset(dataB,Region %in% c("4")) #Denmark
subB7 <- subset(dataB,Region %in% c("5")) #Russia

#descriptive stats
var(subB1$Sr)*100000 #easier to see variation when not in sci notation
var(subB2$Sr)*100000 #easier to see variation when not in sci notation
var(subB3$Sr)*100000 #easier to see variation when not in sci notation
var(subB4$Sr)*100000 #easier to see variation when not in sci notation
var(subB5$Sr)*100000 #easier to see variation when not in sci notation
var(subB6$Sr)*100000 #easier to see variation when not in sci notation
var(subB7$Sr)*100000 #easier to see variation when not in sci notation

#plot
plot(modB3, xlab="Region", ylab="86/87Sr", cex.axis=0.7, col=c("red", "snow3", "snow3", "snow3"))
subB11 <- subset(dataB,Region %in% c("0","1","4"))
dataB.f4 <- factor(subB11$Region)
modB12 <- subB11$Sr~dataB.f4
plot(modB12, xlab="Region", ylab="86/87Sr", cex.axis=0.7, col=c("red", "snow3", "snow3"))
hist(dataB$Sr, xlab="86/87Sr", main = "Distribution of 86/87Sr")
histogram(~Sr|Group, data=dataB, xlab="86/87Sr", main="Distribution of 86/87Sr Values \nKetton Quarry vs Other Sites", breaks=seq(from=0.69, to=0.75, by=0.01))
histogram(~Sr|Region, data=dataB, xlab="86/87Sr", main="Distribution of 86/87Sr Values \nby Region", breaks=seq(from=0.68, to=0.75, by=0.01))

#inferential stats
kruskal.test(Sr~dataB.f3, data=dataB) #Sr by region
dunnTest(Sr~dataB.f3, data=dataB)

#####libraries for isoscape mapping: some may not be used, but better to have too many libraries than make something break
library(readr)
library(readxl)
library("rgdal")
library("raster")
library("maps")
library("maptools")
library("rasterVis")
library("mvnmle")
library("mixtools")
library("fossil")
library(sp)
library(maptools)
library(rgdal)
library(classInt)
library(assignR)
library(colorRamps)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(cellWise)
library(ggsci)
library(ggpmisc)
library(ggpubr)
library(ggridges)
library(magrittr)
library(latex2exp)
library(viridis)
library(ggrepel)
library(cowplot)
library(lavaan)
library(plyr)
library(caTools)
library(bitops)
library(Hmisc)
library(ggdist)
library(tidyquant)
library(see)
library(assignR)

####################################################################################
# Isoscape probability mapping from Colleter et al 2021 and Bataille et al 2021
####################################################################################

####################################################################################
# Define new functions
####################################################################################

## function returns raster of posterior probabilities for univariate normal data
calcCellProb <- function(x,isoscape,std){
  m <- getValues(isoscape)
  s <- getValues(std) # use this is you have raster of variances
  m <- dnorm(x,mean=m,sd=s)
  cell.dens <- setValues(isoscape,m)
  return(cell.dens)
}

## function returns raster of posterior probabilities for bivariate normal data
## x is the unknown tissue of interest, will have two values, one for each isotope
## m is a 2-D vector, all the values in the raster for each isotope
## v is the same as m, but for variances
## r is a single number - the covariance. Can be vector if estimated as non-stationary
## ras is a raster that will serve as a template for the final product
calcCellProb2D <- function(x,m,v,r,ras) {
  pd <- 1/(2*pi*sqrt(v[,1])*sqrt(v[,2])*sqrt(1-r^2))*exp(-(1/(2*(1-r^2)))*
                                                           ((x[1]-m[,1])^2/v[,1]+(x[2]-m[,2])^2/v[,2]-(2*r*(x[1]-m[,1])*
                                                                                                         (x[2]-m[,2]))/(sqrt(v[,1])*sqrt(v[,2]))))
  pdras <- setValues(ras,pd)
  return(pdras)
}

## function returns raster of posterior probability distribution
calcPostProb <- function(x){
  pp <- x/cellStats(x,sum)
  return(pp)
}


## function to get normalized cell probabilites
calcNormProb <- function(x){
  np <- x/cellStats(x,max)
  return(np)
}

####################################################################################
## Data input: read in tissue isotope data, GIS raster models - have to split data
## KQ individuals with only d34S values first 
####################################################################################
iso.data <- read.csv("~/KQ S csv.csv") # tissue isotope data, Appendix G
view(iso.data)

data(wrld_simpl)

### Download d18O isoscape rasters at RCWIP http://www-naweb.iaea.org/napc/ih/IHS_resources_rcwip.html
### input the annual amount-weighted mean and sd from RCWIP

rain_d18O_aaw <- raster("~/RCWIP_grid_5_14_100_20130502.tif")  # this is baseline isoscape
rain_d18O_err_aaw <- raster("~/RCWIP_grid_err_5_14_100_20130502.tif")  # this is baseline isoscape

### Download 87Sr/86Sr isoscape rasters and uncertainty at https://drive.google.com/drive/folders/1g9rCGo3Kd3hz2o5JKkSbgNsGJclvsuQm?usp=sharing
### See Bataille et al. 2020 for details
rf_sr<-raster("~/rf_plantsoilmammal1.tif")
rf_sr_err <-raster("~/srse.tif")

#Clip data to study area Europe
#europe<-as.vector(c(-500000,1000000,5000000,6500000))
europe<-as.vector(c(-860000,2700000,4300000,7900000))
rangemap<-crop(rf_sr/rf_sr, europe, snap='near')
writeRaster(rangemap, "~/range.tif",overwrite=TRUE)
sr<-rf_sr*rangemap
sr.se<-rf_sr_err*rangemap

###Reproject d18O isoscape to Sr projection and resolution
d18O <- projectRaster(rain_d18O_aaw,sr, method="bilinear")
d18O.se<-projectRaster(rain_d18O_err_aaw,sr, method="bilinear")+1###Add 1 per mile uncertainty due to regression equation

###The d34S isoscape is available at https://drive.google.com/drive/folders/1g9rCGo3Kd3hz2o5JKkSbgNsGJclvsuQm?usp=sharing
d34S<-raster("~/rf_d34S.tif")
d34S<-d34S*rangemap
d34S.se<-d34S/d34S*3

d34S_proj<-project(as.matrix(iso.data[,c("Lat","Long")]), "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
iso.data$X<-d34S_proj[,1]
iso.data$Y<-d34S_proj[,2]
#writeRaster(d18O, "d18O.tif", overwrite=TRUE)
#writeRaster(sr, "sr.tif",overwrite=TRUE)

########################### single isotope assignments for sulphur #########################################################

distance.list <- vector("list", length(iso.data[,1]))
origins <- stack()
for (i in seq(along=iso.data[,1])){
  pp <- calcCellProb(iso.data$d34S[i],d34S,d34S.se) # compute probs
  np <- calcNormProb(pp)
  origins <- raster::stack(origins,np)
  
}

#re-project the site coordinates to get the site to plot (thanks to P. Schauer)
myproj <- proj4string(origins)#get the CRS string
iso.data.sp <- iso.data #make a copy
coordinates(iso.data.sp) <- ~ Long + Lat #specify the coordinates
proj4string(iso.data.sp) <- CRS("+proj=longlat +datum=WGS84") #project first
iso.data.sp <- spTransform(iso.data.sp, myproj) #then re-project using the projection from the base map

wrld_simpl2 <- wrld_simpl #make a copy
wrld_simpl2 <- spTransform(wrld_simpl2, myproj) #this is to add in country borders etc. if you want to as they have the same projection issues as the site coords

## summary map plots
pdf("maps_S_KQ.pdf")
opar<-par()
par(mfrow=c(2,2),mar=c(2,3,2,2)+0.01)
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  #plot(wrld_simpl,add=T)
  # plot known origin and individual label info
  #points(s.data$Long,s.data$Lat,pch=21, col="black", bg=colcode, lwd=0.4, cex=symbol.size)
  points(iso.data.sp[i,],col="black",cex=0.2, pch=16)
  #text(-5,47,paste("??34S=",iso.data$d34S[i],"???"),cex=0.70)
  text(-55000, 7650000,paste(iso.data$Grave_ID[i]),cex=0.70)
  #legend(-13, 65, legend=names(attr(colcode, "table")),
  #      fill=attr(colcode, "palette"), cex=0.6, bty="n")
  #hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)

###Store posterior probability maps in raster stack
S34<-origins

####################################################################################
## Data input: read in tissue isotope data, GIS raster models - have to split data 
## Triple-isotope individuals (KQ and Auldhame)
####################################################################################
iso.data <- read.csv("~/KQ S_Sr_O csv final.csv") # tissue isotope data, Appendix G
view(iso.data)

data(wrld_simpl)

### d18O single isotope assignments for chenery dw values 

distance.list <- vector("list", length(iso.data[,1]))
origins <- stack()
for (i in seq(along=iso.data[,1])){
  pp <- calcCellProb(iso.data$d18O_C[i],d18O,d18O.se) # compute probs
  np <- calcNormProb(pp)
  origins <- raster::stack(origins,np)
}

#re-project the site coordinates to get the site to plot (thanks to P. Schauer)
myproj <- proj4string(origins)#get the CRS string
iso.data.sp <- iso.data #make a copy
coordinates(iso.data.sp) <- ~ Long + Lat #specify the coordinates
proj4string(iso.data.sp) <- CRS("+proj=longlat +datum=WGS84") #project first
iso.data.sp <- spTransform(iso.data.sp, myproj) #then re-project using the projection from the base map

wrld_simpl2 <- wrld_simpl #make a copy
wrld_simpl2 <- spTransform(wrld_simpl2, myproj) #this is to add in country borders etc. if you want to as they have the same projection issues as the site coords


## summary map plots
pdf("maps_d18OKQ_Final.pdf")
opar<-par()
par(mfrow=c(2,2),mar=c(2,3,2,2)+0.01)
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  #plot(wrld_simpl,add=T)
  # plot known origin and individual label info
  #points(s.data$Long,s.data$Lat,pch=21, col="black", bg=colcode, lwd=0.4, cex=symbol.size)
  points(iso.data.sp[i,],col="black",cex=0.2, pch=16)
  #text(-5,47,paste("??34S=",iso.data$d34S[i],"???"),cex=0.70)
  text(-55000, 7650000,paste(iso.data$Grave_ID[i]),cex=0.70)
  #legend(-13, 65, legend=names(attr(colcode, "table")),
  #      fill=attr(colcode, "palette"), cex=0.6, bty="n")
  #hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)

###Store posterior probability maps in raster stack
O<-origins

### SR single isotope assignments

distance.list <- vector("list", length(iso.data[,1]))
origins <- stack()
for (i in seq(along=iso.data[,1])){
  pp <- calcCellProb(iso.data$Sr[i],sr,sr.se) # compute probs
  np <- calcNormProb(pp)
  origins <- raster::stack(origins,np)
  
}

#re-project the site coordinates to get the site to plot (thanks to P. Schauer)
myproj <- proj4string(origins)#get the CRS string
iso.data.sp <- iso.data #make a copy
coordinates(iso.data.sp) <- ~ Long + Lat #specify the coordinates
proj4string(iso.data.sp) <- CRS("+proj=longlat +datum=WGS84") #project first
iso.data.sp <- spTransform(iso.data.sp, myproj) #then re-project using the projection from the base map

wrld_simpl2 <- wrld_simpl #make a copy
wrld_simpl2 <- spTransform(wrld_simpl2, myproj) #this is to add in country borders etc. if you want to as they have the same projection issues as the site coords

## summary map plots
pdf("maps_srKQ_Final.pdf")
opar<-par()
par(mfrow=c(2,2),mar=c(2,3,2,2)+0.01)
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  #plot(wrld_simpl,add=T)
  # plot known origin and individual label info
  #points(s.data$Long,s.data$Lat,pch=21, col="black", bg=colcode, lwd=0.4, cex=symbol.size)
  points(iso.data.sp[i,],col="black",cex=0.2, pch=16)
  #text(-5,47,paste("??34S=",iso.data$d34S[i],"???"),cex=0.70)
  text(-55000, 7650000,paste(iso.data$Grave_ID[i]),cex=0.70)
  #legend(-13, 65, legend=names(attr(colcode, "table")),
  #      fill=attr(colcode, "palette"), cex=0.6, bty="n")
  #hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)

###Store posterior probability maps in raster stack
Sr87_86Sr<-origins

###########################d18O and 87Sr/86Sr dual assignments#########################################################
###Combine d18O and 87Sr/86Sr assuming independence, using the chenery outputs
OSr<-O*Sr87_86Sr

pdf("maps_OSrKQ_final.pdf")
opar<-par()
par(mfrow=c(2,2),mar=c(2,3,2,2)+0.01)
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(OSr[[i]],axes = FALSE)
  #plot(wrld_simpl,add=T)
  # plot known origin and individual label info
  #points(s.data$Long,s.data$Lat,pch=21, col="black", bg=colcode, lwd=0.4, cex=symbol.size)
  points(iso.data.sp[i,],col="black",cex=0.2, pch=16)
  #text(-5,47,paste("??34S=",iso.data$d34S[i],"???"),cex=0.70)
  text(-55000, 7650000,paste(iso.data$Grave_ID[i]),cex=0.70)
  #legend(-13, 65, legend=names(attr(colcode, "table")),
  #      fill=attr(colcode, "palette"), cex=0.6, bty="n")
  #hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)

### single isotope assignments for sulphur 

distance.list <- vector("list", length(iso.data[,1]))
origins <- stack()
for (i in seq(along=iso.data[,1])){
  pp <- calcCellProb(iso.data$d34S[i],d34S,d34S.se) # compute probs
  np <- calcNormProb(pp)
  origins <- raster::stack(origins,np)
  
}

#re-project the site coordinates to get the site to plot (thanks to P. Schauer)
myproj <- proj4string(origins)#get the CRS string
iso.data.sp <- iso.data #make a copy
coordinates(iso.data.sp) <- ~ Long + Lat #specify the coordinates
proj4string(iso.data.sp) <- CRS("+proj=longlat +datum=WGS84") #project first
iso.data.sp <- spTransform(iso.data.sp, myproj) #then re-project using the projection from the base map

wrld_simpl2 <- wrld_simpl #make a copy
wrld_simpl2 <- spTransform(wrld_simpl2, myproj) #this is to add in country borders etc. if you want to as they have the same projection issues as the site coords

## summary map plots
pdf("maps_SKQ_Final.pdf")
opar<-par()
par(mfrow=c(2,2),mar=c(2,3,2,2)+0.01)
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(origins[[i]],axes = FALSE)
  #plot(wrld_simpl,add=T)
  # plot known origin and individual label info
  #points(s.data$Long,s.data$Lat,pch=21, col="black", bg=colcode, lwd=0.4, cex=symbol.size)
  points(iso.data.sp[i,],col="black",cex=0.2, pch=16)
  #text(-5,47,paste("??34S=",iso.data$d34S[i],"???"),cex=0.70)
  text(-55000, 7650000,paste(iso.data$Grave_ID[i]),cex=0.70)
  #legend(-13, 65, legend=names(attr(colcode, "table")),
  #      fill=attr(colcode, "palette"), cex=0.6, bty="n")
  #hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)

###Store posterior probability maps in raster stack
S34<-origins

###########################d18O, 87Sr/86Sr, 34S triple assignments#########################################################
###Combine d18O and 87Sr/86Sr assuming independence, using the chenery outputs
SOSr<-O*Sr87_86Sr*S34

pdf("maps_SOSrKQ_Final.pdf")
opar<-par()
par(mfrow=c(2,2),mar=c(2,3,2,2)+0.01)
xl<-length(iso.data[,1])
for (i in seq(along=iso.data[,1])){
  plot(SOSr[[i]],axes = FALSE)
  #plot(wrld_simpl,add=T)
  # plot known origin and individual label info
  #points(s.data$Long,s.data$Lat,pch=21, col="black", bg=colcode, lwd=0.4, cex=symbol.size)
  points(iso.data.sp[i,],col="black",cex=0.2, pch=16)
  #text(-5,47,paste("??34S=",iso.data$d34S[i],"???"),cex=0.70)
  text(-55000, 7650000,paste(iso.data$Grave_ID[i]),cex=0.70)
  #legend(-13, 65, legend=names(attr(colcode, "table")),
  #      fill=attr(colcode, "palette"), cex=0.6, bty="n")
  #hist(distance.list[[i]],xlim=c(0,xl),xlab="",ylab="",main="")
}
par(opar)
dev.off()
rm(i)
