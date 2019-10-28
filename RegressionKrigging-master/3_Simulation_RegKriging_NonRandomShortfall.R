#Script to create a virtual geographic gradient of species richness
#with a spatially non-random knowledge shortfall and apply 
#the Regression-Kriging and Ordinary Kriging model to recover 
#this virtual gradient

#-------------------------------------------------------------------------------#
#Summary
#1- Simulated gradient with a spatially non-random shortfall scenario
#2- Regression-Kriging and Ordinary Kriging to recover the simulated gradient
  #with 60% of degradation
  #with 75% of degradation
  #with 80% of degradation
  #with 85% of degradation
  #with 90% of degradation
  #with 95% of degradation




#--------------------------------------------------------------------------------#

#1- Simulated gradient with a spatially non-random shortfall scenariO

#Function  

#Description: a focal raster (e.g. species richness) is multiplied by a raster representing
             #the scientific knowledge about that information (e.g. human density)

#Argument: 
#rich: raster representing species richness
#fact: raster representing the scientific knowledge about species richness (obs: 
       #the function transforms this raster)

#Results:
#"results"[[1]]: "degraded" species richness gradient
#resultado[[2]]:  raster with the "degraded" values (transformed "fact" raster)

rich_deg2 <- function(rich, fact){
  
  #log-transform the values of the "factor" raster into relative values (i_value/max_value)
  log.factor <- log(round(fact,0))
  max.v <- maxValue(log.factor)
  rel.factor <- log.factor/max.v
  
  size <-  dim(rich)#number of rows and columns of the raster
  
  for (i in 1:size[1]){#run accross the raster rows
    
    for(j in 1:size[2]){#run accross the raster columns
      
      um <- rich[i,j]#select one cell
      dois <- is.na(um)#is NA?
      
      if(dois!=TRUE){ #if not, do this...
        
        rich[i,j] <- rich[i,j] * rel.factor[i,j]#multiply by the factor value
      }
    }
  }
  
  resu <- list(rich, rel.factor)
  
  return(resu)
}  

#Human density (e.g. http://sedac.ciesin.columbia.edu/data)
h.pop <- raster("HumanPopulation2015.asc")#read the raster used as proxy of "scientific knowledge"
h.pop <- disaggregate(h.pop, fact=2)#setting the right resolution
h.pop <- crop(h.pop, temp)


log.factort <- log2(round(h.pop,0))
max.vt <- maxValue(log.factort)
rel.factort <- log.factort/max.vt
rel.factort2 <- (round(h.pop,0)-cellStats(h.pop,median))/max.vt
plot(rel.factort2)



#"Degrading" the gradient
res.rich_deg2 <- rich_deg2 (mapaRiq, h.pop)
res.rich_deg2 [[1]]
res.rich_deg2 [[2]]
plot(res.rich_deg2 [[1]])

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Richness Pattern")
plot(h.pop,main="Human Population Pattern")
plot(res.rich_deg2 [[2]],main="Log_relative Human Population")  
plot(res.rich_deg2 [[1]],main="Observed Richness Pattern") 


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

#---60% degradation---#

#semi-variogram
fit.hum.var.06 <- fit.variog(res.rich_deg2 [[1]], res.rich_deg2 [[2]], 0.6)
plot(fit.hum.var.06[[1]])
plot(rasterFromXYZ(fit.hum.var.06[[2]]))
plot(rasterFromXYZ(fit.hum.var.06[[3]]))

#Sample points
samp.hum.06<-as.data.frame(fit.hum.var.06[[3]])

#Co-variable sample points
temp.hum.06.rs <- mask(crop(temp,rasterFromXYZ(fit.hum.var.06[[3]])), rasterFromXYZ(fit.hum.var.06[[3]]))
temp.hum.06.df <- rasterToPoints(temp.hum.06.rs)

prec.hum.06.rs <- mask(crop(prec,rasterFromXYZ(fit.hum.var.06[[3]])) , rasterFromXYZ(fit.hum.var.06[[3]]))
prec.hum.06.df <- rasterToPoints(prec.hum.06.rs)

#Full sample points
samp.hum.06 <- data.frame(samp.hum.06, temp.hum.06.df[,3], 
                          prec.hum.06.df[,3])
colnames(samp.hum.06) <- c("x","y","rich","temp","prec")
samp.hum.06.sp <- samp.hum.06
coordinates(samp.hum.06.sp) <- ~x+y
samp.hum.06.sp2 <- samp.hum.06.sp 
proj4string(samp.hum.06.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.hum.06 <- autofitVariogram(rich~temp+prec, input_data = samp.hum.06.sp2)

#Regression-Kriging
r.v.hum.06#insert the parameters from this variogram model in the function "krige()"

r.krig.hum.06 <- krige(formula=rich~temp+prec, locations=samp.hum.06.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                       model=vgm(psill=288.849698,model="Ste",range=1970.132,kappa=0.3,nugget=8.686873))


#Determining 0 for cells with "rich < 0"
r.krig.hum.06.df <- as.data.frame(r.krig.hum.06)
r.krig.hum.06.df.rich <- ifelse(r.krig.hum.06.df[,3]<0,0,r.krig.hum.06.df[,3]) 
r.krig.hum.06.df.rich <- data.frame(r.krig.hum.06.df[,1:2],r.krig.hum.06.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.hum.06 <- r.krig.hum.06.df.rich[,3] - rich.df[,3]
rmse.r.hum.06 <- sqrt((sum((r.krig.diff.hum.06)^2))/length(r.krig.diff.hum.06))
r.krig.hum.06.df.rich <- data.frame(r.krig.hum.06.df.rich,r.krig.diff.hum.06) 
colnames(r.krig.hum.06.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.hum.06 <- autofitVariogram(rich~1, input_data = samp.hum.06.sp2)

#Kriging
v.hum.06#insert the parameters from this variogram model in the function "krige()"

krig.hum.06 <- krige(formula=rich~1, locations=samp.hum.06.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                       model=vgm(psill=1308.584,model="Ste",range=6483.441,kappa=0.4,nugget=26.58955))

#Determining 0 for cells with "rich < 0"
krig.hum.06.df <- as.data.frame(krig.hum.06)
krig.hum.06.df.rich <- ifelse(krig.hum.06.df[,3]<0,0,krig.hum.06.df[,3]) 
krig.hum.06.df.rich <- data.frame(krig.hum.06.df[,1:2],krig.hum.06.df.rich)

#Root mean squate error (RMSE)
krig.diff.hum.06 <- krig.hum.06.df.rich[,3] - rich.df[,3]
rmse.hum.06 <- sqrt((sum((krig.diff.hum.06)^2))/length(krig.diff.hum.06))
krig.hum.06.df.rich <- data.frame(krig.hum.06.df.rich,krig.diff.hum.06) 
colnames(krig.hum.06.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg2[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.hum.var.06[[2]]),zlim=c(0,maxValue(mapaRiq$layer)),main="Sample Points")
plot(rasterFromXYZ(r.krig.hum.06.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.hum.06.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))

#---65% degradation---#

#semi-variogram
fit.hum.var.065 <- fit.variog(res.rich_deg2 [[1]], res.rich_deg2 [[2]], 0.65)
plot(fit.hum.var.065[[1]])
plot(rasterFromXYZ(fit.hum.var.065[[2]]))
plot(rasterFromXYZ(fit.hum.var.065[[3]]))

#Sample points
samp.hum.065<-as.data.frame(fit.hum.var.065[[3]])

#Co-variable sample points
temp.hum.065.rs <- mask(crop(temp,rasterFromXYZ(fit.hum.var.065[[3]])), rasterFromXYZ(fit.hum.var.065[[3]]))
temp.hum.065.df <- rasterToPoints(temp.hum.065.rs)

prec.hum.065.rs <- mask(crop(prec,rasterFromXYZ(fit.hum.var.065[[3]])) , rasterFromXYZ(fit.hum.var.065[[3]]))
prec.hum.065.df <- rasterToPoints(prec.hum.065.rs)

#Full sample points
samp.hum.065 <- data.frame(samp.hum.065, temp.hum.065.df[,3], 
                          prec.hum.065.df[,3])
colnames(samp.hum.065) <- c("x","y","rich","temp","prec")
samp.hum.065.sp <- samp.hum.065
coordinates(samp.hum.065.sp) <- ~x+y
samp.hum.065.sp2 <- samp.hum.065.sp 
proj4string(samp.hum.065.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.hum.065 <- autofitVariogram(rich~temp+prec, input_data = samp.hum.065.sp2)

#Regression-Kriging
r.v.hum.065#insert the parameters from this variogram model in the function "krige()"

r.krig.hum.065 <- krige(formula=rich~temp+prec, locations=samp.hum.065.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                       model=vgm(psill=1191.421519,model="Ste",range=89087.75,kappa=0.2,nugget=1.685316))


#Determining 0 for cells with "rich < 0"
r.krig.hum.065.df <- as.data.frame(r.krig.hum.065)
r.krig.hum.065.df.rich <- ifelse(r.krig.hum.065.df[,3]<0,0,r.krig.hum.065.df[,3]) 
r.krig.hum.065.df.rich <- data.frame(r.krig.hum.065.df[,1:2],r.krig.hum.065.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.hum.065 <- r.krig.hum.065.df.rich[,3] - rich.df[,3]
rmse.r.hum.065 <- sqrt((sum((r.krig.diff.hum.065)^2))/length(r.krig.diff.hum.065))
r.krig.hum.065.df.rich <- data.frame(r.krig.hum.065.df.rich,r.krig.diff.hum.065) 
colnames(r.krig.hum.065.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.hum.065 <- autofitVariogram(rich~1, input_data = samp.hum.065.sp2)

#Kriging
v.hum.065#insert the parameters from this variogram model in the function "krige()"

krig.hum.065 <- krige(formula=rich~1, locations=samp.hum.065.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                        model=vgm(psill=18311.125,model="Ste",range=267594.5,kappa=0.4,nugget=49.601))


#Determining 0 for cells with "rich < 0"
krig.hum.065.df <- as.data.frame(krig.hum.065)
krig.hum.065.df.rich <- ifelse(krig.hum.065.df[,3]<0,0,krig.hum.065.df[,3]) 
krig.hum.065.df.rich <- data.frame(krig.hum.065.df[,1:2],krig.hum.065.df.rich)

#Root mean squate error (RMSE)
krig.diff.hum.065 <- krig.hum.065.df.rich[,3] - rich.df[,3]
rmse.hum.065 <- sqrt((sum((krig.diff.hum.065)^2))/length(krig.diff.hum.065))
krig.hum.065.df.rich <- data.frame(krig.hum.065.df.rich,krig.diff.hum.065) 
colnames(krig.hum.065.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg2[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.hum.var.065[[2]]),zlim=c(0,maxValue(mapaRiq$layer)),main="Sample Points")
plot(rasterFromXYZ(r.krig.hum.065.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.hum.065.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))


#---70% degradation---#

#semi-variogram
fit.hum.var.07 <- fit.variog(res.rich_deg2 [[1]], res.rich_deg2 [[2]], 0.7)
plot(fit.hum.var.07[[1]])
plot(rasterFromXYZ(fit.hum.var.07[[2]]))
plot(rasterFromXYZ(fit.hum.var.07[[3]]))

#Sample points
samp.hum.07 <- as.data.frame(fit.hum.var.07[[3]])

#Co-variable sample points
temp.hum.07.rs <- mask(crop(temp,rasterFromXYZ(fit.hum.var.07[[3]])), rasterFromXYZ(fit.hum.var.07[[3]]))
temp.hum.07.df <- rasterToPoints(temp.hum.07.rs)

prec.hum.07.rs <- mask(crop(prec,rasterFromXYZ(fit.hum.var.07[[3]])) , rasterFromXYZ(fit.hum.var.07[[3]]))
prec.hum.07.df <- rasterToPoints(prec.hum.07.rs)

#Full sample points
samp.hum.07 <- data.frame(samp.hum.07, temp.hum.07.df[,3], 
                           prec.hum.07.df[,3])
colnames(samp.hum.07) <- c("x","y","rich","temp","prec")
samp.hum.07.sp <- samp.hum.07
coordinates(samp.hum.07.sp) <- ~x+y
samp.hum.07.sp2 <- samp.hum.07.sp 
proj4string(samp.hum.07.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#Variogram model (r-krig)
r.v.hum.07 <- autofitVariogram(rich~temp+prec, input_data = samp.hum.07.sp2)

#Regression-Kriging
r.v.hum.07#insert the parameters from this variogram model in the function "krige()"

r.krig.hum.07 <- krige(formula=rich~temp+prec, locations=samp.hum.07.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                        model=vgm(psill=978.896800,model="Ste",range=47969.12,kappa=0.2,nugget=6.304754))

#Determining 0 for cells with "rich < 0"
r.krig.hum.07.df <- as.data.frame(r.krig.hum.07)
r.krig.hum.07.df.rich <- ifelse(r.krig.hum.07.df[,3]<0,0,r.krig.hum.07.df[,3])
r.krig.hum.07.df.rich <- data.frame(r.krig.hum.07.df[,1:2],r.krig.hum.07.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.hum.07 <- r.krig.hum.07.df.rich[,3] - rich.df[,3]
rmse.r.hum.07 <- sqrt((sum((r.krig.diff.hum.07)^2))/length(r.krig.diff.hum.07))#accuracy (RMSE)
r.krig.hum.07.df.rich <- data.frame(r.krig.hum.07.df.rich,r.krig.diff.hum.07) 
colnames(r.krig.hum.07.df.rich ) <- c("x","y","rich","error")

#Variogram model (krig)
v.hum.07 <- autofitVariogram(rich~1, input_data = samp.hum.07.sp2)

#Kriging
v.hum.07#insert the parameters from this variogram model in the function "krige()"

krig.hum.07 <- krige(formula=rich~1, locations=samp.hum.07.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                       model=vgm(psill=894068.54627,model="Ste",range=1923709,kappa=0.6,nugget=86.01934))

#Determining 0 for cells with "rich < 0"
krig.hum.07.df <- as.data.frame(krig.hum.07)
krig.hum.07.df.rich <- ifelse(krig.hum.07.df[,3]<0,0,krig.hum.07.df[,3])
krig.hum.07.df.rich <- data.frame(krig.hum.07.df[,1:2],krig.hum.07.df.rich)

#Root mean squate error (RMSE)
krig.diff.hum.07 <- krig.hum.07.df.rich[,3] - rich.df[,3]
rmse.hum.07 <- sqrt((sum((krig.diff.hum.07)^2))/length(krig.diff.hum.07))#accuracy (RMSE)
krig.hum.07.df.rich <- data.frame(krig.hum.07.df.rich,krig.diff.hum.07) 
colnames(krig.hum.07.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg2[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.hum.var.07[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(r.krig.hum.07.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.hum.07.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()
x11()

#---75% degradation---#

#semi-variogram
fit.hum.var.075 <- fit.variog(res.rich_deg2 [[1]], res.rich_deg2 [[2]], 0.75)
plot(fit.hum.var.075[[1]])
plot(rasterFromXYZ(fit.hum.var.075[[2]]))
plot(rasterFromXYZ(fit.hum.var.075[[3]]))

#Sample points
samp.hum.075<-as.data.frame(fit.hum.var.075[[3]])

#Co-variable sample points
temp.hum.075.rs <- mask(crop(temp,rasterFromXYZ(fit.hum.var.075[[3]])), rasterFromXYZ(fit.hum.var.075[[3]]))
temp.hum.075.df <- rasterToPoints(temp.hum.075.rs)

prec.hum.075.rs <- mask(crop(prec,rasterFromXYZ(fit.hum.var.075[[3]])) , rasterFromXYZ(fit.hum.var.075[[3]]))
prec.hum.075.df <- rasterToPoints(prec.hum.075.rs)

#Full sample points
samp.hum.075 <- data.frame(samp.hum.075, temp.hum.075.df[,3], 
                           prec.hum.075.df[,3])
colnames(samp.hum.075) <- c("x","y","rich","temp","prec")
samp.hum.075.sp <- samp.hum.075
coordinates(samp.hum.075.sp) <- ~x+y
samp.hum.075.sp2 <- samp.hum.075.sp 
proj4string(samp.hum.075.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#Variogram model (r-krig)
r.v.hum.075 <- autofitVariogram(rich~temp+prec, input_data = samp.hum.075.sp2)

#Regression-Kriging
r.v.hum.075#insert the parameters from this variogram model in the function "krige()"

r.krig.hum.075 <- krige(formula=rich~temp+prec, locations=samp.hum.075.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                        model=vgm(psill=181.1703,model="Ste",range=138.3405,kappa=0.6,nugget=0))

#Determining 0 for cells with "rich < 0"
r.krig.hum.075.df <- as.data.frame(r.krig.hum.075)
r.krig.hum.075.df.rich <- ifelse(r.krig.hum.075.df[,3]<0,0,r.krig.hum.075.df[,3])
r.krig.hum.075.df.rich <- data.frame(r.krig.hum.075.df[,1:2],r.krig.hum.075.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.hum.075 <- r.krig.hum.075.df.rich[,3] - rich.df[,3]
rmse.r.hum.075 <- sqrt((sum((r.krig.diff.hum.075)^2))/length(r.krig.diff.hum.075))#accuracy (RMSE)
r.krig.hum.075.df.rich <- data.frame(r.krig.hum.075.df.rich,r.krig.diff.hum.075) 
colnames(r.krig.hum.075.df.rich ) <- c("x","y","rich","error")

#Variogram model (r-krig)
v.hum.075 <- autofitVariogram(rich~1, input_data = samp.hum.075.sp2)

#Regression-Kriging
v.hum.075#insert the parameters from this variogram model in the function "krige()"

krig.hum.075 <- krige(formula=rich~1, locations=samp.hum.075.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                        model=vgm(psill=170494.922,model="Ste",range=483254.7,kappa=0.6,nugget=114.425))

#Determining 0 for cells with "rich < 0"
krig.hum.075.df <- as.data.frame(krig.hum.075)
krig.hum.075.df.rich <- ifelse(krig.hum.075.df[,3]<0,0,krig.hum.075.df[,3])
krig.hum.075.df.rich <- data.frame(krig.hum.075.df[,1:2],krig.hum.075.df.rich)

#Root mean squate error (RMSE)
krig.diff.hum.075 <- krig.hum.075.df.rich[,3] - rich.df[,3]
rmse.hum.075 <- sqrt((sum((krig.diff.hum.075)^2))/length(krig.diff.hum.075))#accuracy (RMSE)
krig.hum.075.df.rich <- data.frame(krig.hum.075.df.rich,krig.diff.hum.075) 
colnames(krig.hum.075.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg2[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.hum.var.075[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(r.krig.hum.075.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.hum.075.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()


#---80% degradation---#

#semi-variogram
fit.hum.var.080 <- fit.variog(res.rich_deg2 [[1]], res.rich_deg2 [[2]], 0.80)
plot(fit.hum.var.080[[1]])
plot(rasterFromXYZ(fit.hum.var.080[[2]]))
plot(rasterFromXYZ(fit.hum.var.080[[3]]))

#Sample points
samp.hum.080<-as.data.frame(fit.hum.var.080[[3]])

#Co-variable sample points
temp.hum.080.rs <- mask(crop(temp,rasterFromXYZ(fit.hum.var.080[[3]])), rasterFromXYZ(fit.hum.var.080[[3]]))
temp.hum.080.df <- rasterToPoints(temp.hum.080.rs)

prec.hum.080.rs <- mask(crop(prec,rasterFromXYZ(fit.hum.var.080[[3]])) , rasterFromXYZ(fit.hum.var.080[[3]]))
prec.hum.080.df <- rasterToPoints(prec.hum.080.rs)

#Full sample points
samp.hum.080 <- data.frame(samp.hum.080, temp.hum.080.df[,3], 
                           prec.hum.080.df[,3])
colnames(samp.hum.080) <- c("x","y","rich","temp","prec")
samp.hum.080.sp <- samp.hum.080
coordinates(samp.hum.080.sp) <- ~x+y
samp.hum.080.sp2 <- samp.hum.080.sp 
proj4string(samp.hum.080.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#Variogram model (r-krig)
r.v.hum.080 <- autofitVariogram(rich~temp+prec, input_data = samp.hum.080.sp2)

#Regression-Kriging
r.v.hum.080#insert the parameters from this variogram model in the function "krige()"

r.krig.hum.080 <- krige(formula=rich~temp+prec, locations=samp.hum.080.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                        model=vgm(psill=159.1654,model="Gau",range=65.83032,nugget=0))

#Determining 0 for cells with "rich < 0"
r.krig.hum.080.df <- as.data.frame(r.krig.hum.080)
r.krig.hum.080.df.rich <- ifelse(r.krig.hum.080.df[,3]<0,0,r.krig.hum.080.df[,3])
r.krig.hum.080.df.rich <- data.frame(r.krig.hum.080.df[,1:2],r.krig.hum.080.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.hum.080 <- r.krig.hum.080.df.rich[,3] - rich.df[,3]
rmse.r.hum.080 <- sqrt((sum((r.krig.diff.hum.080)^2))/length(r.krig.diff.hum.080))#accuracy (RMSE)
r.krig.hum.080.df.rich <- data.frame(r.krig.hum.080.df.rich,r.krig.diff.hum.080) 
colnames(r.krig.hum.080.df.rich ) <- c("x","y","rich","error")

#Variogram model (krig)
v.hum.080 <- autofitVariogram(rich~1, input_data = samp.hum.080.sp2)

#Kriging
v.hum.080#insert the parameters from this variogram model in the function "krige()"

krig.hum.080 <- krige(formula=rich~1, locations=samp.hum.080.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                        model=vgm(psill=294999.1603,model="Ste",range=99467.25,kappa=1.4,nugget=169.6302))

#Determining 0 for cells with "rich < 0"
krig.hum.080.df <- as.data.frame(krig.hum.080)
krig.hum.080.df.rich <- ifelse(krig.hum.080.df[,3]<0,0,krig.hum.080.df[,3])
krig.hum.080.df.rich <- data.frame(krig.hum.080.df[,1:2],krig.hum.080.df.rich)

#Root mean squate error (RMSE)
krig.diff.hum.080 <- krig.hum.080.df.rich[,3] - rich.df[,3]
rmse.hum.080 <- sqrt((sum((krig.diff.hum.080)^2))/length(krig.diff.hum.080))#accuracy (RMSE)
krig.hum.080.df.rich <- data.frame(krig.hum.080.df.rich,krig.diff.hum.080) 
colnames(krig.hum.080.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg2[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.hum.var.080[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(r.krig.hum.080.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.hum.080.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()


#---85% degradation---#

#semi-variogram
fit.hum.var.085 <- fit.variog(res.rich_deg2 [[1]], res.rich_deg2 [[2]], 0.85)
plot(fit.hum.var.085[[1]])
plot(rasterFromXYZ(fit.hum.var.085[[2]]))
plot(rasterFromXYZ(fit.hum.var.085[[3]]))

#Sample points
samp.hum.085<-as.data.frame(fit.hum.var.085[[3]])

#Co-variable sample points
temp.hum.085.rs <- mask(crop(temp,rasterFromXYZ(fit.hum.var.085[[3]])), rasterFromXYZ(fit.hum.var.085[[3]]))
temp.hum.085.df <- rasterToPoints(temp.hum.085.rs)

prec.hum.085.rs <- mask(crop(prec,rasterFromXYZ(fit.hum.var.085[[3]])) , rasterFromXYZ(fit.hum.var.085[[3]]))
prec.hum.085.df <- rasterToPoints(prec.hum.085.rs)

#Full sample points
samp.hum.085 <- data.frame(samp.hum.085, temp.hum.085.df[,3], 
                           prec.hum.085.df[,3])
colnames(samp.hum.085) <- c("x","y","rich","temp","prec")
samp.hum.085.sp <- samp.hum.085
coordinates(samp.hum.085.sp) <- ~x+y
samp.hum.085.sp2 <- samp.hum.085.sp 
proj4string(samp.hum.085.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.hum.085 <- autofitVariogram(rich~temp+prec, input_data = samp.hum.085.sp2)

#Regression-Kriging
r.v.hum.085#insert the parameters from this variogram model in the function "krige()"

r.krig.hum.085 <- krige(formula=rich~temp+prec, locations=samp.hum.085.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                        model=vgm(psill=106.22520,model="Ste",range=116.4211,kappa=2,nugget=48.97641))

#Determining 0 for cells with "rich < 0"
r.krig.hum.085.df <- as.data.frame(r.krig.hum.085)
r.krig.hum.085.df.rich <- ifelse(r.krig.hum.085.df[,3]<0,0,r.krig.hum.085.df[,3])
r.krig.hum.085.df.rich <- data.frame(r.krig.hum.085.df[,1:2],r.krig.hum.085.df.rich)

#Root mean squate error (RMSE)
r.krig.diff.hum.085 <- r.krig.hum.085.df.rich[,3] - rich.df[,3]
rmse.r.hum.085 <- sqrt((sum((r.krig.diff.hum.085)^2))/length(r.krig.diff.hum.085))#accuracy (RMSE)
r.krig.hum.085.df.rich <- data.frame(r.krig.hum.085.df.rich,r.krig.diff.hum.085) 
colnames(r.krig.hum.085.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.hum.085 <- autofitVariogram(rich~1, input_data = samp.hum.085.sp2)

#Kriging
v.hum.085#insert the parameters from this variogram model in the function "krige()"

krig.hum.085 <- krige(formula=rich~1, locations=samp.hum.085.sp2, newdata=full.sp2,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                        model=vgm(psill=161.98401,model="Sph",range=257.5932,nugget=88.18795))

#Determining 0 for cells with "rich < 0"
krig.hum.085.df <- as.data.frame(krig.hum.085)
krig.hum.085.df.rich <- ifelse(krig.hum.085.df[,3]<0,0,krig.hum.085.df[,3])
krig.hum.085.df.rich <- data.frame(krig.hum.085.df[,1:2],krig.hum.085.df.rich)

#Root mean squate error (RMSE)
krig.diff.hum.085 <- krig.hum.085.df.rich[,3] - rich.df[,3]
rmse.hum.085 <- sqrt((sum((krig.diff.hum.085)^2))/length(krig.diff.hum.085))#accuracy (RMSE)
krig.hum.085.df.rich <- data.frame(krig.hum.085.df.rich,krig.diff.hum.085) 
colnames(krig.hum.085.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg2[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.hum.var.085[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(r.krig.hum.085.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.hum.085.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()



#---90% degradation---#

#semi-variogram
fit.hum.var.09 <- fit.variog(res.rich_deg2 [[1]], res.rich_deg2 [[2]], 0.9)
plot(fit.hum.var.09[[1]])
plot(rasterFromXYZ(fit.hum.var.09[[2]]))
plot(rasterFromXYZ(fit.hum.var.09[[3]]))

#Sample points
samp.hum.09<-as.data.frame(fit.hum.var.09[[3]])

#Co-variable sample points
temp.hum.09.rs <- mask(crop(temp,rasterFromXYZ(fit.hum.var.09[[3]])), rasterFromXYZ(fit.hum.var.09[[3]]))
temp.hum.09.df <- rasterToPoints(temp.hum.09.rs)

prec.hum.09.rs <- mask(crop(prec,rasterFromXYZ(fit.hum.var.09[[3]])) , rasterFromXYZ(fit.hum.var.09[[3]]))
prec.hum.09.df <- rasterToPoints(prec.hum.09.rs)

#Full sample points
samp.hum.09 <- data.frame(samp.hum.09, temp.hum.09.df[,3], 
                          prec.hum.09.df[,3])
colnames(samp.hum.09) <- c("x","y","rich","temp","prec")
samp.hum.09.sp <- samp.hum.09
coordinates(samp.hum.09.sp) <- ~x+y
samp.hum.09.sp2 <- samp.hum.09.sp 
proj4string(samp.hum.09.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig) 
r.v.hum.09 <- autofitVariogram(rich~temp+prec, input_data = samp.hum.09.sp2)

#Regression-Kriging

r.v.hum.09#insert the parameters from this variogram model in the function "krige()"

r.krig.hum.09 <- krige(formula=rich~temp+prec, locations=samp.hum.09.sp2, newdata=full.sp2,
                       model=vgm(psill=131.1698,model="Ste",range=57.8766, kappa=0.8,nugget=0))

#Determining 0 for cells with rich<0
r.krig.hum.09.df <- as.data.frame(r.krig.hum.09)
r.krig.hum.09.df.rich <- ifelse(r.krig.hum.09.df[,3]<0,0,r.krig.hum.09.df[,3])
r.krig.hum.09.df.rich <- data.frame(r.krig.hum.09.df[,1:2],r.krig.hum.09.df.rich)

#Root mean square error (RMSE)
r.krig.diff.hum.09 <- r.krig.hum.09.df.rich[,3] - rich.df[,3]
rmse.r.hum.09 <- sqrt((sum((r.krig.diff.hum.09)^2))/length(r.krig.diff.hum.09))
r.krig.hum.09.df.rich <- data.frame(r.krig.hum.09.df.rich,r.krig.diff.hum.09) 
colnames(r.krig.hum.09.df.rich ) <- c("x","y","rich","error")

#variogram model (krig) 
v.hum.09 <- autofitVariogram(rich~1, input_data = samp.hum.09.sp2)

#Kriging
v.hum.09#insert the parameters from this variogram model in the function "krige()"

krig.hum.09 <- krige(formula=rich~1, locations=samp.hum.09.sp2, newdata=full.sp2,
                       model=vgm(psill=139.2843,model="Ste",range=353.506, kappa=1.6,nugget=133.1158))

#Determining 0 for cells with rich<0
krig.hum.09.df <- as.data.frame(krig.hum.09)
krig.hum.09.df.rich <- ifelse(krig.hum.09.df[,3]<0,0,krig.hum.09.df[,3])
krig.hum.09.df.rich <- data.frame(krig.hum.09.df[,1:2],krig.hum.09.df.rich)

#Root mean square error (RMSE)
krig.diff.hum.09 <- krig.hum.09.df.rich[,3] - rich.df[,3]
rmse.hum.09 <- sqrt((sum((krig.diff.hum.09)^2))/length(krig.diff.hum.09))
krig.hum.09.df.rich <- data.frame(krig.hum.09.df.rich,krig.diff.hum.09) 
colnames(krig.hum.09.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg2[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.hum.var.09[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))#SpatialPointsDataFrame (cels. com valores > "x" de "prop" sao = 0) 
plot(rasterFromXYZ(r.krig.hum.09.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.hum.09.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()
x11()


#---95% degradation---#

#semi-variogram
fit.hum.var.095 <- fit.variog(res.rich_deg2 [[1]], res.rich_deg2 [[2]], 0.95)
plot(fit.hum.var.095[[1]])
plot(rasterFromXYZ(fit.hum.var.095[[2]]))
plot(rasterFromXYZ(fit.hum.var.095[[3]]))

#Sample points
samp.hum.095 <- as.data.frame(fit.hum.var.095[[3]])

#Co-variable sample points
temp.hum.095.rs <- mask(crop(temp,rasterFromXYZ(fit.hum.var.095[[3]])), rasterFromXYZ(fit.hum.var.095[[3]]))
temp.hum.095.df <- rasterToPoints(temp.hum.095.rs)

prec.hum.095.rs <- mask(crop(prec,rasterFromXYZ(fit.hum.var.095[[3]])) , rasterFromXYZ(fit.hum.var.095[[3]]))
prec.hum.095.df <- rasterToPoints(prec.hum.095.rs)

#Full sample points
samp.hum.095 <- data.frame(samp.hum.095, temp.hum.095.df[,3], 
                          prec.hum.095.df[,3])
colnames(samp.hum.095) <- c("x","y","rich","temp","prec")
samp.hum.095.sp <- samp.hum.095
coordinates(samp.hum.095.sp) <- ~x+y
samp.hum.095.sp2 <- samp.hum.095.sp 
proj4string(samp.hum.095.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#variogram model (r-krig)
r.v.hum.095 <- autofitVariogram(rich~temp+prec, input_data = samp.hum.095.sp2)

#Regression-Kriging
r.v.hum.095#insert the parameters from this variogram model in the function "krige()"

r.krig.hum.095 <- krige(formula=rich~temp+prec, locations=samp.hum.095.sp2, newdata=full.sp2,
                       model=vgm(psill=3764.52483,model="Ste",range=12778.07, kappa=10,nugget=93.59895))

#Determining 0 for cells with rich<0
r.krig.hum.095.df <- as.data.frame(r.krig.hum.095)
r.krig.hum.095.df.rich <- ifelse(r.krig.hum.095.df[,3]<0,0,r.krig.hum.095.df[,3])
r.krig.hum.095.df.rich <- data.frame(r.krig.hum.095.df[,1:2],r.krig.hum.095.df.rich)

#Root mean square error (RMSE)
r.krig.diff.hum.095 <- r.krig.hum.095.df.rich[,3] - rich.df[,3]
rmse.r.hum.095 <- sqrt((sum((r.krig.diff.hum.095)^2))/length(r.krig.diff.hum.095))
r.krig.hum.095.df.rich <- data.frame(r.krig.hum.095.df.rich,r.krig.diff.hum.095) 
colnames(r.krig.hum.095.df.rich ) <- c("x","y","rich","error")

#variogram model (krig)
v.hum.095 <- autofitVariogram(rich~1, input_data = samp.hum.095.sp2)

#Regression-Kriging
v.hum.095#insert the parameters from this variogram model in the function "krige()"

krig.hum.095 <- krige(formula=rich~1, locations=samp.hum.095.sp2, newdata=full.sp2,
                        model=vgm(psill=374.07604,model="Ste",range=3427.833, kappa=0.2,nugget=83.17597))

#Determining 0 for cells with rich<0
krig.hum.095.df <- as.data.frame(krig.hum.095)
krig.hum.095.df.rich <- ifelse(krig.hum.095.df[,3]<0,0,krig.hum.095.df[,3])
krig.hum.095.df.rich <- data.frame(krig.hum.095.df[,1:2],krig.hum.095.df.rich)

#Root mean square error (RMSE)
krig.diff.hum.095 <- krig.hum.095.df.rich[,3] - rich.df[,3]
rmse.hum.095 <- sqrt((sum((krig.diff.hum.095)^2))/length(krig.diff.hum.095))
krig.hum.095.df.rich <- data.frame(krig.hum.095.df.rich,krig.diff.hum.095) 
colnames(krig.hum.095.df.rich ) <- c("x","y","rich","error")

#Figures
par(mfrow=c(2,2))
plot(mapaRiq,main="Simulated Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(res.rich_deg2[[1]],main="Observed Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(fit.hum.var.095[[2]]),main="Sample Points",zlim=c(0,maxValue(mapaRiq$layer)))#SpatialPointsDataFrame (cels. com valores > "x" de "prop" sao = 0) 
plot(rasterFromXYZ(r.krig.hum.095.df.rich[,-4]),main="R-Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
plot(rasterFromXYZ(krig.hum.095.df.rich[,-4]),main="Krig Pattern",zlim=c(0,maxValue(mapaRiq$layer)))
dev.off()

#Plot with all RMSE's 

#r-kriging
hum_rmse <- c(rmse.r.hum.06,rmse.r.hum.065,rmse.r.hum.07,rmse.r.hum.075,rmse.r.hum.080,
              rmse.r.hum.085,rmse.r.hum.09,rmse.r.hum.095)
comp <- c(60, 65, 70, 75, 80, 85, 90, 95)
plot(comp,hum_rmse)
dev.off()

#kriging
k.hum_rmse <- c(rmse.hum.06,rmse.hum.065,rmse.hum.07,rmse.hum.075,rmse.hum.080,
              rmse.hum.085,rmse.hum.09,rmse.hum.095)
comp <- c(60, 65, 70, 75, 80, 85, 90, 95)
plot(comp,k.hum_rmse)
dev.off()
