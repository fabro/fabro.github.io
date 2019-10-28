#Script to create a virtual geographic gradient of species richness
#with a spatially random knowledge shortfall and apply Regression-Kriging 
#and Ordinary Kriging model to recover this virtual gradient

#--------------------------------------------------------------------------------#
#Summary
#1- Simulated gradient with a spatially random shortfall scenario
#2- Regression-Kriging to recover the simulated gradient
  #with 60% of degradation
  #with 75% of degradation
  #with 80% of degradation
  #with 85% of degradation
  #with 90% of degradation
  #with 95% of degradation

#--------------------------------------------------------------------------------#

#1- Simulated gradient with a spatially random shortfall scenario

#Function  

#Description: all raster cells (e.g. of species richness) will be multiplied by a value
#sampled from a uniforme distribution (i.e. a proxy of local scientific knowledge)

#Argument: 
#rich: species richness raster

#Value:
#"result"[[1]]: "degraded" species richness gradient
#"result"[[2]]: raster with the "degraded" values (sample from the uniform distribution)

rich_deg1 <- function(rich){
  
  deg <- raster(extent(rich), res = res(rich))#raster that will receive the degraded values
  
  size <-  dim(rich)
  
  for (i in 1:size[1]){#"run" across raster rows
    
    for(j in 1:size[2]){#"run" across raster columns
      
      um <- rich[i,j]
      
      dois <- is.na(um)
      
      if(dois!=TRUE){#if the "rich" cell is not NA, insert a degraded value in the same 
                     #position in "deg"
        
        deg[i,j] <- runif(1)#sampling a degradation value from a uniform distribution
        
      }
      
    }
    
  }
  
  rich.deg <- rich*deg#degrading the gradient
  
  resu <- list(rich.deg, deg)
  
  return(resu)
}  

#Analyzis
res.rich_deg1 <- rich_deg1 (mapaRiq)

#Figures
par(mfrow=c(1,3))
plot(mapaRiq)
plot(res.rich_deg1 [[1]],) 
plot(res.rich_deg1 [[2]])  
dev.off()

#---------------------------------------------------------------------------------------------#

#2- Regression-Kriging and Ordinary Kriging

#Function to search for the best variogram model

#Description: select cells from the richness raster based on the degradation percentage;use
#them to create the semi-variogram, fit different variogram models and select the best one

#Arguments: 
#rich: richness raster
#deg: degraded raster
#prop: degradation proportion

#Values:
#"results"[[1]]: results from "autofitVariogram()"
#"results"[[2]]: SpatialPointsDataFrame (cells with values < "prop" = 0)
#results[[3]]: SpatialPointsDataFrame (only with cells with values > "prop")

fit.variog <- function(rich, deg, prop = 0.5){
  
  samp.rast <- rich * (deg > prop)
  
  samp.points1 <- rasterToPoints(samp.rast)
  
  samp.points2 <- as.data.frame(samp.points1[-c(which(samp.points1[,3]==0)),])
  
  coordinates(samp.points2) <- ~x+y
  
  variog <- autofitVariogram(layer~1, samp.points2)
  
  resu <- list(variog, samp.points1, samp.points2)
  
  return (resu)
  
}

#Input data
res.rich_deg1 [[1]]#degraded gradient raster
res.rich_deg1 [[2]]#raster with degraded values


#---60% degradation---#

#semi-variogram
fit.ran.var.06 <- fit.variog(res.rich_deg1 [[1]], res.rich_deg1 [[2]], 0.6)
plot(fit.ran.var.06[[1]])
plot(rasterFromXYZ(fit.ran.var.06[[2]])) 
plot(rasterFromXYZ(fit.ran.var.06[[3]]))
dev.off()

#Sample points
samp.ran.06 <- as.data.frame(fit.ran.var.06[[3]])

#Co-variable sample points
temp.ran.06.rs <- mask(crop(temp,rasterFromXYZ(fit.ran.var.06[[3]])), rasterFromXYZ(fit.ran.var.06[[3]]))
temp.ran.06.df <- rasterToPoints(temp.ran.06.rs)

prec.ran.06.rs <- mask(crop(prec,rasterFromXYZ(fit.ran.var.06[[3]])) , rasterFromXYZ(fit.ran.var.06[[3]]))
prec.ran.06.df <- rasterToPoints(prec.ran.06.rs)


#Full sample points 
samp.ran.06 <- data.frame(samp.ran.06, temp.ran.06.df[,3], 
                          prec.ran.06.df[,3])
colnames(samp.ran.06) <- c("x","y","rich","temp","prec")
samp.ran.06.sp <- samp.ran.06
coordinates(samp.ran.06.sp) <- ~x+y
samp.ran.06.sp2 <- samp.ran.06.sp 
proj4string(samp.ran.06.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.ran.06 <- autofitVariogram(rich~temp+prec, input_data = samp.ran.06.sp2)

#Regression kriging
r.v.ran.06#insert the parameters from this variogram model in the function "krige()"

r.krig.ran.06 <- krige(formula=rich~temp+prec, locations=samp.ran.06.sp2, newdata=full.sp2,
                       model=vgm(psill=245.20836,model="Ste",range=1426.541,kappa=0.3,nugget=61.78844))

#Determing 0 for cells with "rich < 0"
r.krig.ran.06.df <- as.data.frame(r.krig.ran.06)
r.krig.ran.06.df.rich <- ifelse(r.krig.ran.06.df[,3]<0,0,r.krig.ran.06.df[,3])
r.krig.ran.06.df.rich <- data.frame(r.krig.ran.06.df[,1:2],r.krig.ran.06.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.ran.06 <- r.krig.ran.06.df.rich[,3] - rich.df[,3]
rmse.r.ran.06 <- sqrt((sum((r.krig.diff.ran.06)^2))/length(r.krig.diff.ran.06))#accuracy (RMSE)
r.krig.ran.06.df.rich <- data.frame(r.krig.ran.06.df.rich,r.krig.diff.ran.06) 
colnames(r.krig.ran.06.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.ran.06 <- autofitVariogram(rich~1, input_data = samp.ran.06.sp2)

#Kriging
v.ran.06#insert the parameters from this variogram model in the function "krige()"

krig.ran.06 <- krige(formula=rich~1, locations=samp.ran.06.sp2, newdata=full.sp2,
                       model=vgm(psill=10199.5640,model="Ste",range=309655.5,kappa=0.3,nugget=45.7688))

#Determing 0 for cells with "rich < 0"
krig.ran.06.df <- as.data.frame(krig.ran.06)
krig.ran.06.df.rich <- ifelse(krig.ran.06.df[,3]<0,0,krig.ran.06.df[,3])
krig.ran.06.df.rich <- data.frame(krig.ran.06.df[,1:2],krig.ran.06.df.rich)

#Root mean squate error (RMSE)
krig.diff.ran.06 <- krig.ran.06.df.rich[,3] - rich.df[,3]
rmse.ran.06 <- sqrt((sum((krig.diff.ran.06)^2))/length(krig.diff.ran.06))#accuracy (RMSE)
krig.ran.06.df.rich <- data.frame(krig.ran.06.df.rich,krig.diff.ran.06) 
colnames(krig.ran.06.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg1[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.ran.var.06[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(r.krig.ran.06.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.ran.06.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()

#---65% degradation---#

#semi-variogram
fit.ran.var.065 <- fit.variog(res.rich_deg1 [[1]], res.rich_deg1 [[2]], 0.65)
plot(fit.ran.var.065[[1]])
plot(rasterFromXYZ(fit.ran.var.065[[2]])) 
plot(rasterFromXYZ(fit.ran.var.065[[3]]))
dev.off()

#Sample points
samp.ran.065 <- as.data.frame(fit.ran.var.065[[3]])

#Co-variable sample points
temp.ran.065.rs <- mask(crop(temp,rasterFromXYZ(fit.ran.var.065[[3]])), rasterFromXYZ(fit.ran.var.065[[3]]))
temp.ran.065.df <- rasterToPoints(temp.ran.065.rs)

prec.ran.065.rs <- mask(crop(prec,rasterFromXYZ(fit.ran.var.065[[3]])) , rasterFromXYZ(fit.ran.var.065[[3]]))
prec.ran.065.df <- rasterToPoints(prec.ran.065.rs)

#Full sample points 
samp.ran.065 <- data.frame(samp.ran.065, temp.ran.065.df[,3], 
                          prec.ran.065.df[,3])
colnames(samp.ran.065) <- c("x","y","rich","temp","prec")
samp.ran.065.sp <- samp.ran.065
coordinates(samp.ran.065.sp) <- ~x+y
samp.ran.065.sp2 <- samp.ran.065.sp 
proj4string(samp.ran.065.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.ran.065 <- autofitVariogram(rich~temp+prec, input_data = samp.ran.065.sp2)

#Regression kriging
r.v.ran.065#insert the parameters from this variogram model in the function "krige()"

r.krig.ran.065 <- krige(formula=rich~temp+prec, locations=samp.ran.065.sp2, newdata=full.sp2,
                       model=vgm(psill=265.89140,model="Ste",range=1578.583,kappa=0.3,nugget=42.30319))

#Determing 0 for cells with "rich < 0"
r.krig.ran.065.df <- as.data.frame(r.krig.ran.065)
r.krig.ran.065.df.rich <- ifelse(r.krig.ran.065.df[,3]<0,0,r.krig.ran.065.df[,3])
r.krig.ran.065.df.rich <- data.frame(r.krig.ran.065.df[,1:2],r.krig.ran.065.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.ran.065 <- r.krig.ran.065.df.rich[,3] - rich.df[,3]
rmse.r.ran.065 <- sqrt((sum((r.krig.diff.ran.065)^2))/length(r.krig.diff.ran.065))#accuracy (RMSE)
r.krig.ran.065.df.rich <- data.frame(r.krig.ran.065.df.rich,r.krig.diff.ran.065) 
colnames(r.krig.ran.065.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.ran.065 <- autofitVariogram(rich~1, input_data = samp.ran.065.sp2)

#Kriging
v.ran.065#insert the parameters from this variogram model in the function "krige()"

krig.ran.065 <- krige(formula=rich~1, locations=samp.ran.065.sp2, newdata=full.sp2,
                        model=vgm(psill=1162.17427,model="Ste",range=5940.474,kappa=0.4,nugget=55.57243))

#Determing 0 for cells with "rich < 0"
krig.ran.065.df <- as.data.frame(krig.ran.065)
krig.ran.065.df.rich <- ifelse(krig.ran.065.df[,3]<0,0,krig.ran.065.df[,3])
krig.ran.065.df.rich <- data.frame(krig.ran.065.df[,1:2],krig.ran.065.df.rich)

#Root mean squate error (RMSE)
krig.diff.ran.065 <- krig.ran.065.df.rich[,3] - rich.df[,3]
rmse.ran.065 <- sqrt((sum((krig.diff.ran.065)^2))/length(krig.diff.ran.065))#accuracy (RMSE)
krig.ran.065.df.rich <- data.frame(krig.ran.065.df.rich,krig.diff.ran.065) 
colnames(krig.ran.065.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg1[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.ran.var.065[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(r.krig.ran.065.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.ran.065.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()

#---70% degradation---#

#semi-variogram
fit.ran.var.07 <- fit.variog(res.rich_deg1 [[1]], res.rich_deg1 [[2]], 0.7)
plot(fit.ran.var.07[[1]])
plot(rasterFromXYZ(fit.ran.var.07[[2]]))
plot(rasterFromXYZ(fit.ran.var.07[[3]]))

#Sample points
samp.ran.07<-as.data.frame(fit.ran.var.07[[3]])

#Co-variable sample points
temp.ran.07.rs <- mask(crop(temp,rasterFromXYZ(fit.ran.var.07[[3]])), rasterFromXYZ(fit.ran.var.07[[3]]))
temp.ran.07.df <- rasterToPoints(temp.ran.07.rs)

prec.ran.07.rs <- mask(crop(prec,rasterFromXYZ(fit.ran.var.07[[3]])) , rasterFromXYZ(fit.ran.var.07[[3]]))
prec.ran.07.df <- rasterToPoints(prec.ran.07.rs)

#Full sample points
samp.ran.07 <- data.frame(samp.ran.07, temp.ran.07.df[,3], 
                           prec.ran.07.df[,3])
colnames(samp.ran.07) <- c("x","y","rich","temp","prec")
samp.ran.07.sp <- samp.ran.07
coordinates(samp.ran.07.sp) <- ~x+y
samp.ran.07.sp2 <- samp.ran.07.sp 
proj4string(samp.ran.07.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.ran.07 <- autofitVariogram(rich~temp+prec, input_data = samp.ran.07.sp2)

#r-Kriging
r.v.ran.07#insert the parameters from this variogram model in the function "krige()"

r.krig.ran.07 <- krige(formula=rich~temp+prec, locations=samp.ran.07.sp2, newdata=full.sp2,
                        model=vgm(psill=273.64079,model="Ste",range=1474.363,kappa=0.3,nugget=26.94171))

#Determing 0 for cells with "rich < 0"
r.krig.ran.07.df <- as.data.frame(r.krig.ran.07)
r.krig.ran.07.df.rich <- ifelse(r.krig.ran.07.df[,3]<0,0,r.krig.ran.07.df[,3])
r.krig.ran.07.df.rich <- data.frame(r.krig.ran.07.df[,1:2],r.krig.ran.07.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.ran.07 <- r.krig.ran.07.df.rich[,3] - rich.df[,3]
rmse.r.ran.07 <- sqrt((sum((r.krig.diff.ran.07)^2))/length(r.krig.diff.ran.07))
r.krig.ran.07.df.rich <- data.frame(r.krig.ran.07.df.rich,r.krig.diff.ran.07) 
colnames(r.krig.ran.07.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.ran.07 <- autofitVariogram(rich~1, input_data = samp.ran.07.sp2)

#Kriging
v.ran.07#insert the parameters from this variogram model in the function "krige()"

krig.ran.07 <- krige(formula=rich~1, locations=samp.ran.07.sp2, newdata=full.sp2,
                       model=vgm(psill=1400.12280,model="Ste",range=7483.598,kappa=0.4,nugget=46.25454))

#Determing 0 for cells with "rich < 0"
krig.ran.07.df <- as.data.frame(krig.ran.07)
krig.ran.07.df.rich <- ifelse(krig.ran.07.df[,3]<0,0,krig.ran.07.df[,3])
krig.ran.07.df.rich <- data.frame(krig.ran.07.df[,1:2],krig.ran.07.df.rich)

#Root mean squate error (RMSE)
krig.diff.ran.07 <- krig.ran.07.df.rich[,3] - rich.df[,3]
rmse.ran.07 <- sqrt((sum((krig.diff.ran.07)^2))/length(krig.diff.ran.07))
krig.ran.07.df.rich <- data.frame(krig.ran.07.df.rich,krig.diff.ran.07) 
colnames(krig.ran.07.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg1[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.ran.var.07[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(r.krig.ran.07.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.ran.07.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()


#---75% degradation---#

#semi-variogram
fit.ran.var.075 <- fit.variog(res.rich_deg1 [[1]], res.rich_deg1 [[2]], 0.75)
plot(fit.ran.var.075[[1]])
plot(rasterFromXYZ(fit.ran.var.075[[2]]))
plot(rasterFromXYZ(fit.ran.var.075[[3]]))

#Sample points
samp.ran.075<-as.data.frame(fit.ran.var.075[[3]])

#Co-variable sample points
temp.ran.075.rs <- mask(crop(temp,rasterFromXYZ(fit.ran.var.075[[3]])), rasterFromXYZ(fit.ran.var.075[[3]]))
temp.ran.075.df <- rasterToPoints(temp.ran.075.rs)

prec.ran.075.rs <- mask(crop(prec,rasterFromXYZ(fit.ran.var.075[[3]])) , rasterFromXYZ(fit.ran.var.075[[3]]))
prec.ran.075.df <- rasterToPoints(prec.ran.075.rs)

#Full sample points
samp.ran.075 <- data.frame(samp.ran.075, temp.ran.075.df[,3], 
                           prec.ran.075.df[,3])
colnames(samp.ran.075) <- c("x","y","rich","temp","prec")
samp.ran.075.sp <- samp.ran.075
coordinates(samp.ran.075.sp) <- ~x+y
samp.ran.075.sp2 <- samp.ran.075.sp 
proj4string(samp.ran.075.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.ran.075 <- autofitVariogram(rich~temp+prec, input_data = samp.ran.075.sp2)

#r-Kriging
r.v.ran.075#insert the parameters from this variogram model in the function "krige()"

r.krig.ran.075 <- krige(formula=rich~temp+prec, locations=samp.ran.075.sp2, newdata=full.sp2,
                        model=vgm(psill=273.17720,model="Ste",range=1632.248,kappa=0.3,nugget=23.98099))

#Determing 0 for cells with "rich < 0"
r.krig.ran.075.df <- as.data.frame(r.krig.ran.075)
r.krig.ran.075.df.rich <- ifelse(r.krig.ran.075.df[,3]<0,0,r.krig.ran.075.df[,3])
r.krig.ran.075.df.rich <- data.frame(r.krig.ran.075.df[,1:2],r.krig.ran.075.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.ran.075 <- r.krig.ran.075.df.rich[,3] - rich.df[,3]
rmse.r.ran.075 <- sqrt((sum((r.krig.diff.ran.075)^2))/length(r.krig.diff.ran.075))
r.krig.ran.075.df.rich <- data.frame(r.krig.ran.075.df.rich,r.krig.diff.ran.075) 
colnames(r.krig.ran.075.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.ran.075 <- autofitVariogram(rich~1, input_data = samp.ran.075.sp2)

#Kriging
v.ran.075#insert the parameters from this variogram model in the function "krige()"

krig.ran.075 <- krige(formula=rich~1, locations=samp.ran.075.sp2, newdata=full.sp2,
                        model=vgm(psill=1357.99663,model="Ste",range=6530.193,kappa=0.4,nugget=38.48108))

#Determing 0 for cells with "rich < 0"
krig.ran.075.df <- as.data.frame(krig.ran.075)
krig.ran.075.df.rich <- ifelse(krig.ran.075.df[,3]<0,0,krig.ran.075.df[,3])
krig.ran.075.df.rich <- data.frame(krig.ran.075.df[,1:2],krig.ran.075.df.rich)

#Root mean squate error (RMSE)
krig.diff.ran.075 <- krig.ran.075.df.rich[,3] - rich.df[,3]
rmse.ran.075 <- sqrt((sum((krig.diff.ran.075)^2))/length(krig.diff.ran.075))
krig.ran.075.df.rich <- data.frame(krig.ran.075.df.rich,krig.diff.ran.075) 
colnames(krig.ran.075.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg1[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.ran.var.075[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(r.krig.ran.075.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.ran.075.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()


#---80% degradation---#

#semi-variogram
fit.ran.var.080 <- fit.variog(res.rich_deg1 [[1]], res.rich_deg1 [[2]], 0.80)
plot(fit.ran.var.080[[1]])
plot(rasterFromXYZ(fit.ran.var.080[[2]]))
plot(rasterFromXYZ(fit.ran.var.080[[3]]))

#Sample points
samp.ran.080<-as.data.frame(fit.ran.var.080[[3]])

#Co-variable sample points
temp.ran.080.rs <- mask(crop(temp,rasterFromXYZ(fit.ran.var.080[[3]])), rasterFromXYZ(fit.ran.var.080[[3]]))
temp.ran.080.df <- rasterToPoints(temp.ran.080.rs)

prec.ran.080.rs <- mask(crop(prec,rasterFromXYZ(fit.ran.var.080[[3]])) , rasterFromXYZ(fit.ran.var.080[[3]]))
prec.ran.080.df <- rasterToPoints(prec.ran.080.rs)

#Full sample points
samp.ran.080 <- data.frame(samp.ran.080, temp.ran.080.df[,3], 
                           prec.ran.080.df[,3])
colnames(samp.ran.080) <- c("x","y","rich","temp","prec")
samp.ran.080.sp <- samp.ran.080
coordinates(samp.ran.080.sp) <- ~x+y
samp.ran.080.sp2 <- samp.ran.080.sp 
proj4string(samp.ran.080.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.ran.080 <- autofitVariogram(rich~temp+prec, input_data = samp.ran.080.sp2)

#r-Kriging
r.v.ran.080#insert the parameters from this variogram model in the function "krige()"

r.krig.ran.080 <- krige(formula=rich~temp+prec, locations=samp.ran.080.sp2, newdata=full.sp2,
                        model=vgm(psill=287.03907,model="Ste",range=1652.689,kappa=0.3,nugget=16.13323))

#Determing 0 for cells with "rich < 0"
r.krig.ran.080.df <- as.data.frame(r.krig.ran.080)
r.krig.ran.080.df.rich <- ifelse(r.krig.ran.080.df[,3]<0,0,r.krig.ran.080.df[,3])
r.krig.ran.080.df.rich <- data.frame(r.krig.ran.080.df[,1:2],r.krig.ran.080.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.ran.080 <- r.krig.ran.080.df.rich[,3] - rich.df[,3]
rmse.r.ran.080 <- sqrt((sum((r.krig.diff.ran.080)^2))/length(r.krig.diff.ran.080))
r.krig.ran.080.df.rich <- data.frame(r.krig.ran.080.df.rich,r.krig.diff.ran.080) 
colnames(r.krig.ran.080.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.ran.080 <- autofitVariogram(rich~1, input_data = samp.ran.080.sp2)

#Kriging
v.ran.080#insert the parameters from this variogram model in the function "krige()"

krig.ran.080 <- krige(formula=rich~1, locations=samp.ran.080.sp2, newdata=full.sp2,
                        model=vgm(psill=1355.92072,model="Ste",range=6151.23,kappa=0.4,nugget=31.84657))

#Determing 0 for cells with "rich < 0"
krig.ran.080.df <- as.data.frame(krig.ran.080)
krig.ran.080.df.rich <- ifelse(krig.ran.080.df[,3]<0,0,krig.ran.080.df[,3])
krig.ran.080.df.rich <- data.frame(krig.ran.080.df[,1:2],krig.ran.080.df.rich)

#Root mean squate error (RMSE)
krig.diff.ran.080 <- krig.ran.080.df.rich[,3] - rich.df[,3]
rmse.ran.080 <- sqrt((sum((krig.diff.ran.080)^2))/length(krig.diff.ran.080))
krig.ran.080.df.rich <- data.frame(krig.ran.080.df.rich,krig.diff.ran.080) 
colnames(krig.ran.080.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg1[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.ran.var.080[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(r.krig.ran.080.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.ran.080.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()


#---85% degradation---#

#semi-variogram
fit.ran.var.085 <- fit.variog(res.rich_deg1 [[1]], res.rich_deg1 [[2]], 0.85)
plot(fit.ran.var.085[[1]])
plot(rasterFromXYZ(fit.ran.var.085[[2]]))
plot(rasterFromXYZ(fit.ran.var.085[[3]]))

#Sample points
samp.ran.085<-as.data.frame(fit.ran.var.085[[3]])

#Co-variable sample points
temp.ran.085.rs <- mask(crop(temp,rasterFromXYZ(fit.ran.var.085[[3]])), rasterFromXYZ(fit.ran.var.085[[3]]))
temp.ran.085.df <- rasterToPoints(temp.ran.085.rs)

prec.ran.085.rs <- mask(crop(prec,rasterFromXYZ(fit.ran.var.085[[3]])) , rasterFromXYZ(fit.ran.var.085[[3]]))
prec.ran.085.df <- rasterToPoints(prec.ran.085.rs)

#Full sample points
samp.ran.085 <- data.frame(samp.ran.085, temp.ran.085.df[,3], 
                           prec.ran.085.df[,3])
colnames(samp.ran.085) <- c("x","y","rich","temp","prec")
samp.ran.085.sp <- samp.ran.085
coordinates(samp.ran.085.sp) <- ~x+y
samp.ran.085.sp2 <- samp.ran.085.sp 
proj4string(samp.ran.085.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.ran.085 <- autofitVariogram(rich~temp+prec, input_data = samp.ran.085.sp2)

#r-Kriging
r.v.ran.085#insert the parameters from this variogram model in the function "krige()"

r.krig.ran.085 <- krige(formula=rich~temp+prec, locations=samp.ran.085.sp2, newdata=full.sp2,
                        model=vgm(psill=287.530726,model="Ste",range=1348.192,kappa=0.3,nugget=5.996012))

#Determing 0 for cells with "rich < 0"
r.krig.ran.085.df <- as.data.frame(r.krig.ran.085)
r.krig.ran.085.df.rich <- ifelse(r.krig.ran.085.df[,3]<0,0,r.krig.ran.085.df[,3])
r.krig.ran.085.df.rich <- data.frame(r.krig.ran.085.df[,1:2],r.krig.ran.085.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.ran.085 <- r.krig.ran.085.df.rich[,3] - rich.df[,3]
rmse.r.ran.085 <- sqrt((sum((r.krig.diff.ran.085)^2))/length(r.krig.diff.ran.085))
r.krig.ran.085.df.rich <- data.frame(r.krig.ran.085.df.rich,r.krig.diff.ran.085) 
colnames(r.krig.ran.085.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.ran.085 <- autofitVariogram(rich~1, input_data = samp.ran.085.sp2)

#Kriging
v.ran.085#insert the parameters from this variogram model in the function "krige()"

krig.ran.085 <- krige(formula=rich~1, locations=samp.ran.085.sp2, newdata=full.sp2,
                        model=vgm(psill=988.55030,model="Ste",range=3519.208,kappa=0.4,nugget=24.44609))

#Determing 0 for cells with "rich < 0"
krig.ran.085.df <- as.data.frame(krig.ran.085)
krig.ran.085.df.rich <- ifelse(krig.ran.085.df[,3]<0,0,krig.ran.085.df[,3])
krig.ran.085.df.rich <- data.frame(krig.ran.085.df[,1:2],krig.ran.085.df.rich)

#Root mean squate error (RMSE)
krig.diff.ran.085 <- krig.ran.085.df.rich[,3] - rich.df[,3]
rmse.ran.085 <- sqrt((sum((krig.diff.ran.085)^2))/length(krig.diff.ran.085))
krig.ran.085.df.rich <- data.frame(krig.ran.085.df.rich,krig.diff.ran.085) 
colnames(krig.ran.085.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg1[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.ran.var.085[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(r.krig.ran.085.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.ran.085.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()


#---90% degradation---#

#semi-variogram
fit.ran.var.09 <- fit.variog(res.rich_deg1 [[1]], res.rich_deg1 [[2]], 0.9)
plot(fit.ran.var.09[[1]])
plot(rasterFromXYZ(fit.ran.var.09[[2]]))
plot(rasterFromXYZ(fit.ran.var.09[[3]]))

#Sample points
samp.ran.09<-as.data.frame(fit.ran.var.09[[3]])

#Co-variable sample points

temp.ran.09.rs <- mask(crop(temp,rasterFromXYZ(fit.ran.var.09[[3]])), rasterFromXYZ(fit.ran.var.09[[3]]))
temp.ran.09.df <- rasterToPoints(temp.ran.09.rs)

prec.ran.09.rs <- mask(crop(prec,rasterFromXYZ(fit.ran.var.09[[3]])) , rasterFromXYZ(fit.ran.var.09[[3]]))
prec.ran.09.df <- rasterToPoints(prec.ran.09.rs)

#Full sample points
samp.ran.09 <- data.frame(samp.ran.09, temp.ran.09.df[,3], 
                          prec.ran.09.df[,3])
colnames(samp.ran.09) <- c("x","y","rich","temp","prec")
samp.ran.09.sp <- samp.ran.09
coordinates(samp.ran.09.sp) <- ~x+y
samp.ran.09.sp2 <- samp.ran.09.sp 
proj4string(samp.ran.09.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig) 
r.v.ran.09 <- autofitVariogram(rich~temp+prec, input_data = samp.ran.09.sp2)

#Regression kriging
r.v.ran.09#insert the parameters from this variogram model in the function "krige()"

r.krig.ran.09 <- krige(formula=rich~temp+prec, locations=samp.ran.09.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                       model=vgm(psill=290.9848,model="Ste",range=1307.071,nugget=0,kappa=0.3))

#Determing 0 for cells with rich<0
r.krig.ran.09.df <- as.data.frame(r.krig.ran.09)
r.krig.ran.09.df.rich <- ifelse(r.krig.ran.09.df[,3]<0,0,r.krig.ran.09.df[,3])
r.krig.ran.09.df.rich <- data.frame(r.krig.ran.09.df[,1:2],r.krig.ran.09.df.rich)

#Root mean square error (RMSE)
r.krig.diff.ran.09 <- r.krig.ran.09.df.rich[,3] - rich.df[,3]
rmse.r.ran.09 <- sqrt((sum((r.krig.diff.ran.09)^2))/length(r.krig.diff.ran.09))
r.krig.ran.09.df.rich <- data.frame(r.krig.ran.09.df.rich,r.krig.diff.ran.09) 
colnames(r.krig.ran.09.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.ran.09 <- autofitVariogram(rich~1, input_data = samp.ran.09.sp2)

#Kriging
v.ran.09#insert the parameters from this variogram model in the function "krige()"

krig.ran.09 <- krige(formula=rich~1, locations=samp.ran.09.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                       model=vgm(psill=1055.62698,model="Ste",range=3755.748,nugget=24.27605,kappa=0.4))

#Determing 0 for cells with rich<0
krig.ran.09.df <- as.data.frame(krig.ran.09)
krig.ran.09.df.rich <- ifelse(krig.ran.09.df[,3]<0,0,krig.ran.09.df[,3])
krig.ran.09.df.rich <- data.frame(krig.ran.09.df[,1:2],krig.ran.09.df.rich)

#Root mean square error (RMSE)
krig.diff.ran.09 <- krig.ran.09.df.rich[,3] - rich.df[,3]
rmse.ran.09 <- sqrt((sum((krig.diff.ran.09)^2))/length(krig.diff.ran.09))
krig.ran.09.df.rich <- data.frame(krig.ran.09.df.rich,krig.diff.ran.09) 
colnames(krig.ran.09.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg1[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.ran.var.09[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))#sample raster
plot(rasterFromXYZ(r.krig.ran.09.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.ran.09.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()



#---95% degradation---#

#semi-variogram
fit.ran.var.095 <- fit.variog(res.rich_deg1 [[1]], res.rich_deg1 [[2]], 0.95)
plot(fit.ran.var.095[[1]])
plot(rasterFromXYZ(fit.ran.var.095[[2]]))
plot(rasterFromXYZ(fit.ran.var.095[[3]]))

#Sample points
samp.ran.095<-as.data.frame(fit.ran.var.095[[3]])

#Co-variable sample points
temp.ran.095.rs <- mask(crop(temp,rasterFromXYZ(fit.ran.var.095[[3]])), rasterFromXYZ(fit.ran.var.095[[3]]))
temp.ran.095.df <- rasterToPoints(temp.ran.095.rs)

prec.ran.095.rs <- mask(crop(prec,rasterFromXYZ(fit.ran.var.095[[3]])) , rasterFromXYZ(fit.ran.var.095[[3]]))
prec.ran.095.df <- rasterToPoints(prec.ran.095.rs)

#Full sample points
samp.ran.095 <- data.frame(samp.ran.095, temp.ran.095.df[,3], 
                          prec.ran.095.df[,3])
colnames(samp.ran.095) <- c("x","y","rich","temp","prec")
samp.ran.095.sp <- samp.ran.095
coordinates(samp.ran.095.sp) <- ~x+y
samp.ran.095.sp2 <- samp.ran.095.sp 
proj4string(samp.ran.095.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.ran.095 <- autofitVariogram(rich~temp+prec, input_data = samp.ran.095.sp2)

#Regression kriging
r.v.ran.095#insert the parameters from this variogram model in the function "krige()"

r.krig.ran.095 <- krige(formula=rich~temp+prec, locations=samp.ran.095.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                       model=vgm(psill=264.22266,model="Ste",kappa=0.3,range=1319.741,nugget=12.91133))

#Determing 0 for cells with rich<0
r.krig.ran.095.df <- as.data.frame(r.krig.ran.095)
r.krig.ran.095.df.rich <- ifelse(r.krig.ran.095.df[,3]<0,0,r.krig.ran.095.df[,3])
r.krig.ran.095.df.rich <- data.frame(r.krig.ran.095.df[,1:2],r.krig.ran.095.df.rich)


#Root mean square error (RMSE)
r.krig.diff.ran.095 <- r.krig.ran.095.df.rich[,3] - rich.df[,3]
rmse.r.ran.095 <- sqrt((sum((r.krig.diff.ran.095)^2))/length(r.krig.diff.ran.095))
r.krig.ran.095.df.rich <- data.frame(r.krig.ran.095.df.rich,r.krig.diff.ran.095) 
colnames(r.krig.ran.095.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.ran.095 <- autofitVariogram(rich~1, input_data = samp.ran.095.sp2)

#Regression kriging
v.ran.095#insert the parameters from this variogram model in the function "krige()"

krig.ran.095 <- krige(formula=rich~1, locations=samp.ran.095.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                        model=vgm(psill=2068.80010,model="Ste",kappa=0.4,range=9615.904,nugget=47.26203))

#Determing 0 for cells with rich<0
krig.ran.095.df <- as.data.frame(krig.ran.095)
krig.ran.095.df.rich <- ifelse(krig.ran.095.df[,3]<0,0,krig.ran.095.df[,3])
krig.ran.095.df.rich <- data.frame(krig.ran.095.df[,1:2],krig.ran.095.df.rich)


#Root mean square error (RMSE)
krig.diff.ran.095 <- krig.ran.095.df.rich[,3] - rich.df[,3]
rmse.ran.095 <- sqrt((sum((krig.diff.ran.095)^2))/length(krig.diff.ran.095))
krig.ran.095.df.rich <- data.frame(krig.ran.095.df.rich,krig.diff.ran.095) 
colnames(krig.ran.095.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg1[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.ran.var.095[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))#sample raster
plot(rasterFromXYZ(r.krig.ran.095.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.ran.095.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()


#Plot with all RMSE's

#r-kriging
ran_rmse <- c(rmse.r.ran.06,rmse.r.ran.065,rmse.r.ran.07,rmse.r.ran.075,rmse.r.ran.080,
              rmse.r.ran.085,rmse.r.ran.09,rmse.r.ran.095)
comp <- c(60, 65, 70, 75, 80, 85, 90, 95)
plot(comp,ran_rmse)

#kriging
k.ran_rmse <- c(rmse.ran.06,rmse.ran.065,rmse.ran.07,rmse.ran.075,rmse.ran.080,
              rmse.ran.085,rmse.ran.09,rmse.ran.095)
comp <- c(60, 65, 70, 75, 80, 85, 90, 95)
plot(comp,k.ran_rmse)



