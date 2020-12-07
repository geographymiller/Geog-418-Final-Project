##################################################################################### 
##################################################################################### 
##################################################################################### 
                               #DATA PREP AND CLEANING

#Install Packages
install.packages("spgwr")

#Load Libraries
library("spgwr")
library("spatstat")
library("tmap")
library("gstat")
library("sf")
library("raster")
library("rgdal")
library("e1071")
library("spdep")
library("maptools") 
library("sp")
library("gtable")
library("rgeos")
library("gridExtra")
library("grid")
library("ggplot2")

#Set working directory
dir <- "/Users/Anna/Desktop/Geography 418 Lab Assignments/Geog 418 Final Project"
setwd(dir)

#Reading in elevation dataset
elev <- readOGR(dsn = "./DataFinal", layer = "ElevSample") #Read in data
elev <- spTransform(elev, CRS("+init=epsg:26910"))

#Reading in VRI data
VRI <- readOGR(dsn = "./DataFinal", layer = "WatershedVRI") #Read in shapefile
VRI <- spTransform(VRI, CRS("+init=epsg:26910"))
head(VRI@data)

vriCleanCols <- c("FID_VEG_CO", "POLYGON_ID", "PROJ_AGE_1",
                  "SITE_INDEX", "SPECIES__4", "SPECIES__5",
                  "PROJ_HEI_1", "SPECIES_PC", "SPECIES__6",
                  "VRI_LIVE_S", "BASAL_AREA", "WHOLE_STEM",
                  "CROWN_CL_1")

vriClean <- VRI[,vriCleanCols]

newNames <- c("FID", "PolyID", "Stand_Age", "Site_Index",
              "CoDom_Sp", "Dom_Sp", "Stand_HT", "DomSP_Perc", 
              "CDomSP_Perc", "Stand_Dens", "Stand_BA", "Stand_StemBio", "Stand_CrownCl")

colnames(vriClean@data) <- newNames
head(vriClean@data)

#Choose a Variable and get rid of NA values
vriClean <- vriClean[!is.na(vriClean@data$Stand_HT), ]



#CREATE MAP OF STAND HEIGHT AND STUDY AREA
map_StdHT <- tm_shape(vriClean) +
  tm_polygons(col = "Stand_HT",
              title = "Forest Stand Height (meters)",
              style = "jenks",
              palette = "viridis", n = 6, border.alpha = 0.5) +
  tm_legend(legend.outside=TRUE) +
  tm_scale_bar(position=c("left", "bottom")) +
  tm_compass(type="arrow", position=c("right", "top"), show.labels = 1) +
  tm_layout(main.title = "Forest Stand Height in the Greater Victoria Water Supply Area", main.title.position = "center")

map_StdHT

tmaptools::palette_explorer()

#####################################################################################
#####################################################################################
#####################################################################################
                        #DESCRIPTIVE STATISTICS OBJECTIVE #1

#FOREST STAND MEASURES OF CENTRAL TENDENCY
#Calculate the mean, mode and median for forest stand
meanHeight <- mean(vriClean$Stand_HT, na.rm = TRUE)

modeHeight <- as.numeric(names(sort(table(vriClean$Stand_HT), decreasing = TRUE))[1])

medianHeight <- median(vriClean$Stand_HT, na.rm = TRUE)

#FOREST STAND MEASURES OF VARIABILITY
#Calculate SD, kurtosis, skewness, Co-variance and range for forest stand
sdHeight <- sd(vriClean$Stand_HT, na.rm = TRUE)

kurtHeight <- kurtosis(vriClean$Stand_HT, na.rm = TRUE)

skewHeight <- skewness(vriClean$Stand_HT, na.rm = TRUE)

cvHeight <- (sdHeight / meanHeight) * 100

rangeHeight <- range(vriClean$Stand_HT)

#ELEVATION MEASURES OF CENTRAL TENDENCY
#Calculate the mean, mode and median for elevation
meanElevation <- mean(elev$grid_code, na.rm = TRUE)

modeElevation <- as.numeric(names(sort(table(elev$grid_code), decreasing = TRUE))[1])

medianElevation <- median(elev$grid_code, na.rm = TRUE)

#ELEVATION MEASURES OF VARIABILITY
#Calculate SD, kurtosis, skewness c0-variance and range for elevation
sdElevation <- sd(elev$grid_code, na.rm = TRUE)

kurtElevation <- kurtosis(elev$grid_code, na.rm = TRUE)

skewElevation <- skewness(elev$grid_code, na.rm = TRUE)

cvElevation <- (sdElevation / meanElevation) * 100

rangeElevation <- range(elev$grid_code)

#TABLE FOR FOREST HEIGHT AND ELEVATION DESCRIPTIVE STATS
Variables = c("Forest Stand Height", "Elevation") 
Mean = c(meanHeight, meanElevation) 
Mode <- c(modeHeight, modeElevation)
Median = c(medianHeight, medianElevation)
Standard_Deviation = c(sdHeight, sdElevation) 
Kurtosis <- c(kurtHeight, kurtElevation)
Skewness <- c(skewHeight, skewElevation) 
Coefficient_of_Variation <- c(cvHeight, cvElevation) 
Range <- c(round(rangeHeight, 3), round(rangeElevation, 3))
temp1 <- c(Range[1], "-", Range[2])
temp1
temp1 <- as.character(temp1)
temp1 <- toString(temp1)
temp1 <- sub(",","",temp1)
temp1 <- sub(",","",temp1)

temp2 <- c(Range[3], "-", Range[4])
temp2
temp2 <- as.character(temp2)
temp2 <- toString(temp2)
temp2 <- sub(",","",temp2)
temp2 <- sub(",","",temp2)
Range <- c(temp1, temp2)
Range

#Changing SigFigs...

Mean
Mean <- round(Mean,3)
Mean

Median
Median <- round(Median,3)
Median

Mode
Mode <- round(Mode,3)
Mode

Skewness
Skewness <-round(Skewness,3)
Skewness

Coefficient_of_Variation
Coefficient_of_Variation <-round(Coefficient_of_Variation,3)
Coefficient_of_Variation

Kurtosis
Kurtosis <-round(Kurtosis,3)
Kurtosis


#Make Table #1
data.for.table1 = data.frame(Variables, Mean, Median, Mode, Skewness, Kurtosis, Coefficient_of_Variation, Range)

table1 <- tableGrob(data.for.table1, rows = c("",""))  
t1Caption <- textGrob("Table 1: Statistical Anaylsis for Stand Height and Elevation in the Greater Victoria Water Supply Area", gp = gpar(fontsize = 10))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1)

grid.arrange(table1, newpage = TRUE)


#Create Histogram
hist(vriClean$Stand_HT, breaks = 30, main = "Frequency of Forest Stand Height (in meters) in the GVWSA", xlab = "Tree Height (meters)") 

######################################################################################
######################################################################################
######################################################################################
                    #SPATIAL AUTOCORRELATION OF FOREST HEIGHT OBJECTIVE #2

#Queens Weight Nieghbourhood Defining
vri.nb <- poly2nb(vriClean)
vri.net <- nb2lines(vri.nb, coords=coordinates(vriClean))
crs(vri.net) <- crs(vriClean)

tm_shape(vriClean) + tm_borders(col='lightgrey') + 
  tm_shape(vri.net) + tm_lines(col='red') +
  tm_layout(main.title = "Queen Weights for Neighbourhood Definition of Stand Height", main.title.position = "TOP")

#Weights Matrix
vri.lw <- nb2listw(vri.nb, zero.policy = TRUE, style = "W")
print.listw(vri.lw, zero.policy = TRUE)

#Global Moran's I TEST
mi <- moran.test(vriClean$Stand_HT, vri.lw, zero.policy = TRUE)
mi

mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

z <- (mI-eI)/(sqrt(var))

#Global Moran's I values in Table
Forest_Characteristic = c("Stand Height") #Create an object for the labels
Global_Morans_I = c(mI)
Expected_Morans_I = c(eI) #Create an object for the means
Variance = c(var) #Create an object for the standard deviations
Z_Test = c(z) #Create an object for the medians

data.for.table1 = data.frame(Forest_Characteristic, Global_Morans_I, Expected_Morans_I, Variance, Z_Test)

table1 <- tableGrob(data.for.table1, rows = c(""))
t1Caption <- textGrob("Table 2: Global Moran's I Statistical Analysis for Stand Height", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1)

grid.arrange(table1, newpage = TRUE)
#### TABLE TESTER END

#Local Moran's I TEST
lisa.test <- localmoran(vriClean$Stand_HT, vri.lw, zero.policy = TRUE)

vriClean$Ii <- lisa.test[,1]
vriClean$E.Ii<- lisa.test[,2]
vriClean$Var.Ii<- lisa.test[,3]
vriClean$Z.Ii<- lisa.test[,4]
vriClean$P<- lisa.test[,5]

########################
tmaptools::palette_explorer()

map_LISA <- tm_shape(vriClean) + 
  tm_polygons(col = "Z.Ii", 
              title = "Z-Test Score", 
              style = "fixed", breaks = c(-Inf, -1.96, 1.96, +Inf), 
              palette = "PiYG", n = 3) +
  tm_layout(main.title = "Local Spatial Autocorrelation for Stand Height in the GVWSA", main.title.position = "center")


map_LISA

######
moran.plot(vriClean$Stand_HT, vri.lw, zero.policy=TRUE, spChk=NULL, labels=NULL, xlab="Stand Height", 
           ylab="Spatially Lagged Stand Height", quiet=NULL, main = "Scatterplot of Local Moran's I for Stand Height")


#######################################################################################
#######################################################################################
#######################################################################################
                    #SPATIAL INTERPOLATION OF ELEVATION OBJECTIVE #3

## INVERSE DISTANCE WEIGHTED
#Create a grid called grd to use in your interpolation

grd <- as.data.frame(spsample(elev, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(elev)

elev = na.omit(elev)

#IDW Interpolation *running the interpolation*
P.idw <- gstat::idw(grid_code ~ 1, elev, newdata=grd, idp=4)
r       <- raster(P.idw)
r.m     <- mask(r, vriClean)

IDW_map <- tm_shape(r.m) + 
  tm_raster(n=10,palette = "-PiYG",
            title="Predicted Elevation \n(in meters)") + 
  tm_shape(elev) + tm_dots(size=0.2) +
  tm_layout(main.title = "Elevation Interpolation using IDW", main.title.position = "center") +
  tm_legend(legend.outside=TRUE) 

IDW_map

#################################################
# Leave-one-out validation routine
IDW.out <- vector(length = length(elev))
for (i in 1:length(elev)) {
  IDW.out[i] <- idw(grid_code ~ 1, elev[-i,], elev[i,], idp=4)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ elev$grid_code, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ elev$grid_code), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
sqrt(sum((IDW.out - elev$grid_code)^2) / length(elev))


#################################################
# Implementation of a jackknife technique to estimate a confidence interval at each unsampled point.
# Create the interpolated surface
img <- gstat::idw(grid_code~1, elev, newdata=grd, idp=4)
n   <- length(elev)
Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)

# Remove a point then interpolate (do this n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(grid_code~1, elev[-i,], newdata=grd, idp=4)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
}

# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- img
img.sig$v <- CI /img$var1.pred 

# Clip the confidence raster to Southern California
r.conf <- raster(img.sig, layer="v")
r.m.conf <- mask(r, vriClean)

# Plot the map
#### WHAT DO YOU PUT IN FOR POINTS??????
tm_shape(r.m.conf) + tm_raster(n=7,title="95% Confidence Interval \n(in ppm)") +
  tm_shape(elev) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE) +
  tm_layout(main.title = "Confidence Interval for IDW interpolation technique for Elevation \nin the GVWSA", main.title.position = "center")


#######################################################################################
#######################################################################################
#######################################################################################
            #REGRESSION ANALYSIS OF ELEVATION AND FOREST STAND HEIGHT #Objective 3

#Combining Elevation and Forest Height

#These steps will help you combine the outputs 
#from your spatial interpolation with your income data.
#Convert your interpolation into a raster and map it:
r <- raster(P.idw)
sufaceMap <- tm_shape(r.m) + 
  tm_raster(n=5,palette = "viridis",
            title="Elevation (m)") +
  tm_shape(elev) + tm_dots(size=0.2) +
  tm_layout(legend.outside = TRUE) +
  tm_layout(main.title = "Interpolation of Elevation over the GVWSA Study Area")
sufaceMap

#If you have too many cells, 
#you can reduce the number by aggregating values
#agg <- aggregate(yourRasterFromKriging, fact=??, fun=mean)

#Extract average elev for each polygon
vriClean$Elev <- extract(r, vriClean, fun = mean)[,1]


map_Bio <- tm_shape(vriClean) +
  tm_polygons(col = "Elev", title = "Elevation (meters)", style = "jenks", palette = "viridis", n = 6, border.alpha = 0.5) +
  tm_legend(legend.position = c("LEFT", "BOTTOM")) +
  tm_layout(main.title = "Average Elevation for Each Polygon in the \n Greater Victoria Water Supply Area", main.title.position = "center")
map_Bio

#####################################################################################
#####################################################################################
#####################################################################################
                        #Linear Regression OBJECTIVE #3 CONTINUED....

#Let's say your dataset with both Elev and Height are stored in a dataset called VRI.
#Plot Height and Elev from the VRI dataset you created
plot(vriClean$Stand_HT ~ vriClean$Elev)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
VRI.no0 <-  vriClean[which(vriClean$Stand_HT > 0), ]
VRI.no0 <-  VRI.no0[which(VRI.no0$Elev > 0), ]

#Now plot the data again
plot(VRI.no0$Stand_HT ~ VRI.no0$Elev)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(VRI.no0$Stand_HT ~ VRI.no0$Elev)

#Add the regression model to the plot you created
plot(VRI.no0$Stand_HT ~ VRI.no0$Elev)
abline(lm.model, col = "red")

#Get the summary of the results
summary(lm.model)

#add the fitted values to your spatialpolygon dataframe
VRI.no0$predictlm <- lm.model$fitted.values

#You want to determine if the model residuals are spatially clustered. 
#add the residuals to your spatialpolygon dataframe
VRI.no0$residuals <- residuals.lm(lm.model)

#Observe the result to make sure it looks correct
head(VRI.no0@data)

#Now, create choropleth map of residuals
map_resid <- tm_shape(VRI.no0) +
  tm_polygons(col = "residuals",
              title = "Stand Height Residuals",
              style = "jenks",
              palette = "viridis", n = 6) +
  tm_layout(legend.outside = TRUE)

map_resid
##################################################
# GLOBAL MORANS I OF THE RESIDUALS
#Queens Weight Nieghbourhood Defining
vri.nb <- poly2nb(VRI.no0)
vri.net <- nb2lines(vri.nb, coords=coordinates(VRI.no0))
crs(vri.net) <- crs(VRI.no0)

#Weights Matrix
vri.lw <- nb2listw(vri.nb, zero.policy = TRUE, style = "W")
print.listw(vri.lw, zero.policy = TRUE)


#Global Moran's I TEST
mi.res <- moran.test(VRI.no0$residuals, vri.lw, zero.policy = TRUE)
mi.res

mI.res <- mi$estimate[[1]]
eI.res <- mi$estimate[[2]]
var.res <- mi$estimate[[3]]

z.res <- (mI-eI)/(sqrt(var))

#Global Moran's I values in Table
Variable = c("Residuals") #Create an object for the labels
Global_Morans_I = c(mI.res)
Expected_Morans_I = c(eI.res) #Create an object for the means
Variance = c(var.res) #Create an object for the standard deviations
Z_Test = c(z.res) #Create an object for the medians

data.for.table1 = data.frame(Variable, Global_Morans_I, Expected_Morans_I, Variance, Z_Test)

table1 <- tableGrob(data.for.table1, rows = c(""))
t1Caption <- textGrob("Table 3: Global Moran's I Analysis for Residuals of OLS Regression", gp = gpar(fontsize = 09))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1)

grid.arrange(table1, newpage = TRUE)

#######################################################################################
#######################################################################################
#######################################################################################
#GEOGRAPHICALLY WEIGHTED REGRESSION



#Let's say you are continuing with your data from the regression analysis. 
#The first thing you need to do is to add the polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the "coordinates" function from the sp library
VRI.no0.coords <- sp::coordinates(VRI.no0)
#Observe the result:
head(VRI.no0.coords)
#Now add the coordinates back to the spatialpolygondataframe
VRI.no0$X <- VRI.no0.coords[,1]
VRI.no0$Y <- VRI.no0.coords[,2]

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(VRI.no0$Stand_HT ~ VRI.no0$Elev, 
                        data=VRI.no0, coords=cbind(VRI.no0$X,VRI.no0$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(VRI.no0$Stand_HT ~ VRI.no0$Elev, 
                data=VRI.no0, coords=cbind(VRI.no0$X,VRI.no0$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 


#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
VRI.no0$localr <- results$localR2

#Create choropleth map of r-square values
map_r2 <- tm_shape(VRI.no0) +
  tm_polygons(col = "localr",
              title = "R2 values",
              style = "jenks",
              palette = "viridis", n = 6) +
  tm_layout(main.title = "R-squared Values of Geographically Weighted Regression", main.title.position = "center")
map_r2

#Time for more magic. Let's map the coefficients
VRI.no0$coeff <- results$VRI.no0.Elev
#Create choropleth map of the coefficients
map_coef <- tm_shape(VRI.no0) +
  tm_polygons(col = "coeff",
              title = "Coefficients",
              style = "jenks",
              palette = "viridis", n = 6) +
  tm_layout(main.title = "Coefficients Values of Geographically Weighted Regression", main.title.position = "center")
map_coef





#######################################################################################
#######################################################################################
#######################################################################################
                            #POINT PATTERN ANALYSIS 

kma <- elev
kma$x <- coordinates(kma)[,1]
kma$y <- coordinates(kma)[,2]

#check for and remove duplicated points
#first, finds zero distance among points to see if there are any duplicates
zd <- zerodist(kma)
zd

#if there are duplicates, remove them
kma <- remove.duplicates(kma)

#create an "extent" object which can be used to create the observation window for spatstat
kma.ext <- as.matrix(extent(vriClean)) 

#observation window
window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,]))

#create ppp oject from spatstat
kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)

####################################################################
####################################################################
#NEAREST NEIGHBOUR ANALYSIS

nearestNeighbour <- nndist(kma.ppp)

##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"


##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour
nnd = (sum(nearestNeighbour$Distance))/(353)

#STUDY AREA CALCULATION
wholearea <- gArea(vriClean, byid = FALSE)

#mean nearest neighbour for random spatial distribution

studyArea <- wholearea

pointDensity <- (353)/(wholearea)

r.nnd = ((1)/((2)*(sqrt(pointDensity))))

d.nnd = ((1.07453)/(sqrt(pointDensity)))

R = ((nnd)/(r.nnd))

SE.NND <- ((0.26136)/(sqrt(353*pointDensity)))

z = ((nnd)-(r.nnd))/(SE.NND)


Topographic_Feature = c("Elevation") #Create an object for the labels
Number_Of_Points = c(353)
Nearest_Neighbour = c(nnd) #Create an object for the means
Random_Pattern = c(r.nnd) #Create an object for the standard deviations
Standardized_Mean = c(R) #Create an object for the medians
Z_Test <- c(z) #Create an object for the modes

Nearest_Neighbour
Nearest_Neighbour <- round(Nearest_Neighbour,3)
Nearest_Neighbour

Random_Pattern
Random_Pattern <- round(Random_Pattern,3)
Random_Pattern

Standardized_Mean
Standardized_Mean <-round(Standardized_Mean,3)
Standardized_Mean

Z_Test
Z_Test <-round(Z_Test,3)
Z_Test

data.for.table1 = data.frame(Topographic_Feature, Number_Of_Points, Nearest_Neighbour, Random_Pattern, Standardized_Mean, Z_Test)

table1 <- tableGrob(data.for.table1, rows = c("")) #make a table "Graphical Object" (GrOb) 
t1Caption <- textGrob("Table 4: Nearest Neighbour Analysis for Elevation in the GVWSA", gp = gpar(fontsize = 10))
padding <- unit(5, "mm")

table1 <- gtable_add_rows(table1, 
                          heights = grobHeight(t1Caption) + padding, 
                          pos = 0)

table1 <- gtable_add_grob(table1,
                          t1Caption, t = 1, l = 2, r = ncol(data.for.table1) + 1)

grid.arrange(table1, newpage = TRUE)


###KERNEL DENSITY ESTIMATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#2D (gaussian) kernel, compare how bandwidth (sigma) selection influences the point density estimates
#since data are projected, sigma is represented in metres
#eps is the width and height of the pixels (1000m X 1000m)
#coerce to a SpatialGridDataFrame for plotting
kde.100 <- density(kma.ppp, sigma = 100, at = "pixels", eps = c(400, 400))
kde.SG <- as(kde.100, "SpatialGridDataFrame")
kde.500 <- density(kma.ppp, sigma = 500, at = "pixels", eps = c(400, 400))
kde.SG <- cbind(kde.SG, as(kde.500, "SpatialGridDataFrame"))
kde.250 <- density(kma.ppp, sigma = 250, at = "pixels", eps = c(400, 400))
kde.SG <- cbind(kde.SG, as(kde.250, "SpatialGridDataFrame"))
kde.50 <- density(kma.ppp, sigma = 50, at = "pixels", eps = c(400, 400))
kde.SG <- cbind(kde.SG, as (kde.50, "SpatialGridDataFrame"))


names(kde.SG) <- c("Size100", "Size500", "Size250", "Size50")
#plot
x11() #opens a new plot window
spplot(kde.SG, main = "Kernel Density Estimation for Elevation in the GVWSA")

#can see how the bandwidth selection influences the density estimates
summary(kde.SG)

#use cross-validation to get the bandwidth that minimizes MSE
bw.d <- bw.diggle(kma.ppp)
#plot the "optimal" bandwidth
plot(bw.d, ylim=c(-10, 10), main= "Kernel Density Estimation for Elevation in the GVWSA")

#density using the cross-validation bandwidth
kde.bwo <- density(kma.ppp, sigma = bw.d, at = "pixels", eps = c(100, 100))
plot(kde.bwo, main = "Kernel Density Estimation for Elevation in the GVWSA")


