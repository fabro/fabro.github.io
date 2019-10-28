#Application of the Regression-Kriging model to unveil the geographic
#gradient of Furnariidae species richness for South America


#1- Data
#2- Selecting well-sampled cells
#3- Selecting the co-variates for Regression-Kriging
#4- Regression Kriging 



#--- 1- Data ---#

#Fabaceae occurrence data
fur <- read.table("Furnariidae.txt", h=T)

#Temperature
temp.raw <- raster("Temperature")
proj4string(temp.raw) <- "+proj=aitoff +ellps=WGS84"
temp.raw <- crop(temp.raw,c(-111, -34.5, -56, 28))
temp <- aggregate(temp.raw,fact=0.5/res(temp.raw)[1],fun=mean)
plot(temp)

#Temperature seasonality (standard deviation)
temp.sd.raw <- raster("Temperature sd")
proj4string(temp.sd.raw) <- "+proj=aitoff +ellps=WGS84"
temp.sd.raw <- crop(temp.sd.raw,c(-111, -34.5, -56, 28))
temp.sd <- aggregate(temp.sd.raw,fact=0.5/res(temp.sd.raw)[1],fun=mean)
temp.sd <- mask(temp.sd,temp)
plot(temp.sd)

#canopy height
can.raw <- raster("Canopy height")
proj4string(can.raw) <- "+proj=aitoff +ellps=WGS84"
can.raw <- crop(can.raw,c(-111, -34.5, -56, 28))
can <- aggregate(can.raw, fact=0.5/res(can.raw)[1], fun=mean )
can <- mask(can,temp)
plot(can)

#AET
aet.raw <- raster("AET")
proj4string(aet.raw) <- "+proj=aitoff +ellps=WGS84"
aet.raw <- crop(aet.raw,c(-111, -34.5, -56, 28))
aet <- aggregate(aet.raw, fact=0.5/res(aet.raw)[1], fun=mean )
aet <- mask(aet,temp)
plot(aet)

#Elevation 
elev.raw <- raster("Elevation")
proj4string(elev.raw) <- "+proj=aitoff +ellps=WGS84"
elev.raw <- crop(elev.raw,c(-111, -34.5, -56, 28))
elev<-aggregate(elev.raw,factor=0.5/res(elev.raw)[1],fun=mean)
elev <- mask(elev,temp)
plot(elev)

#Elevation heterogeneity (standard deviation)
elev.sd <- aggregate(elev.raw,fact=0.5/res(elev.raw)[1], fun=sd )
elev.sd <- mask(elev.sd, temp)
plot(elev.sd)

#Precipitation
prec.raw <- raster("Precipitation")
proj4string(prec.raw) <- "+proj=aitoff +ellps=WGS84"
prec.raw <- crop(prec.raw,c(-111, -34.5, -56, 28))
prec <- aggregate(prec.raw,fact=0.5/res(prec.raw)[1],fun=mean)
prec <- mask(prec,temp)
plot(prec)

#Precipitation seasonality (standard deviation)
prec.sd.raw <- raster("Precipitation coefficient variation")
proj4string(prec.sd.raw) <- "+proj=aitoff +ellps=WGS84"
prec.sd.raw <- crop(prec.sd.raw,c(-111, -34.5, -56, 28))
prec.sd <- aggregate(prec.sd.raw,fact=0.5/res(prec.sd.raw)[1],fun=mean)
prec.sd <- mask(prec.sd,temp)
plot(prec.sd)

#Historical temperature (LGM minimum temperature - current minimum temp)
temp.hmin.raw <- raster("LGM minimum temperature")
proj4string(temp.hmin) <- "+proj=aitoff +ellps=WGS84"
temp.hmin.raw <- crop(temp.hmin.raw,c(-111, -34.5, -56, 28))#lgm min temp
temp.hmin <- aggregate(temp.hmin.raw,fact=0.5/res(temp.hmin.raw)[1],fun=mean)
temp.hmin <- mask(temp.hmin,temp)

temp.min.raw <- raster("Minimum temperature")
temp.min.raw <- crop(temp.min.raw,c(-111, -34.5, -56, 28))
temp.min <- aggregate(temp.min.raw,fact=0.5/res(temp.min.raw)[1],fun=mean)
temp.min <- mask(temp.min,temp)

temp.hmin = temp.hmin - temp.min
plot(temp.hmin)
plot(temp.min)

#Historical precipitation (LGM driest quarter - current driest quarter)
prec.hmin.raw <- raster("LGM precipitation from the driest quarter")
proj4string(prec.hmin) <- "+proj=aitoff +ellps=WGS84"
prec.hmin.raw <- crop(prec.hmin.raw,c(-111, -34.5, -56, 28))
prec.hmin <- aggregate(prec.hmin.raw,fact=0.5/res(prec.hmin.raw)[1],fun=mean)
prec.hmin <- mask(prec.hmin,temp)

prec.min.raw <- raster("Current precipitation from the driest quarter")
proj4string(prec.min) <- "+proj=aitoff +ellps=WGS84"
prec.min.raw <- crop(prec.min.raw,c(-111, -34.5, -56, 28))
prec.min <- aggregate(prec.min.raw,fact=0.5/res(prec.min.raw)[1],fun=mean)
prec.min <- mask(prec.min,temp)

prec.hmin = prec.hmin - prec.min
plot(prec.hmin)
plot(prec.min)

#Stack
env.clim <- stack(aet,temp,temp.sd,prec,prec.sd,elev,
                  elev.sd,can,temp.hmin,prec.hmin)

#--- 2- Selecting well-sampled cells ---# 

#KnowB: Completeness > 85%, a Ratio > 8, and a Slope < 0.02
data(adworld)
KnowB(data=fur, format="A", cell=30,  cutoff = 8,
      cutoffCompleteness = 85, cutoffSlope = 0.02, save="RData")

#resultant matrix
load("Estimators.RData")
head(values)

#selecting cells with: Numer of Records > 500, Completeness > 85%, 
#Ratio > 8, and Slope < 0.02
c <- which(values[,3]>500 &values[,7] > 95 & values[,8] > 15 & values[,6]<0.02 )#thresholds 
estim <- values[c,]
dim(estim)

#raster to determine the well-sampled cells
rst1 <- temp 
res(rst1) <- 0.5 

#well-sampled cells
completeness <- rasterize(estim[,1:2],rst1,estim[,7])
proj4string(completeness) <- "+proj=aitoff +ellps=WGS84"
completeness <- crop(completeness,temp)

x11()
plot(completeness>95,axes=F,box=F)
points(estim[,1:2],col="red", pch=18,cex=1)
plot(la,add=T)

#number of records
records <- rasterize(values[,1:2],rst1,values[,3])
proj4string(records) <- "+proj=aitoff +ellps=WGS84"
records <- crop(records,temp)

x11()
plot(log10(records),axes=F,box=F, col=rev(gray.colors(50)))
plot(la,add=T)

#observed richness
richness <- rasterize(values[,1:2],rst1,values[,4])
proj4string(richness) <- "+proj=aitoff +ellps=WGS84"
richness <- crop(richness,c(-111, -34.5, -56, 28))
richness <- mask(richness,temp)

x11()
plot(richness,axes=F,box=F)

#--- 3- Selecting the co-variates for Regression-Kriging ---#

#Data

#sampled points: geographic coordinates and co-variates from the well-sampled cells
aet.var <- raster::extract(x=env.clim[[1]], y=estim[,1:2])
temp.var <- raster::extract(x=env.clim[[2]], y=estim[,1:2])
temp.sd.var <- raster::extract(x=env.clim[[3]], y=estim[,1:2])
prec.var <- raster::extract(x=env.clim[[4]], y=estim[,1:2])
prec.sd.var <- raster::extract(x=env.clim[[5]], y=estim[,1:2])
elev.var <- raster::extract(x=env.clim[[6]], y=estim[,1:2])
elev.sd.var <- raster::extract(x=env.clim[[7]], y=estim[,1:2])
can.var <- raster::extract(x=env.clim[[8]], y=estim[,1:2])
temp.hmin.var <- raster::extract(x=env.clim[[9]], y=estim[,1:2])
prec.hmin.var <- raster::extract(x=,env.clim[[10]], y=estim[,1:2])


samp.points.rat <- data.frame(y=estim[,1], x=estim[,2], rich=estim[,4],
                              aet=aet.var,temp=temp.var,temp.sd=temp.sd.var,
                              prec=prec.var,prec.sd=prec.sd.var,
                              elev=elev.var,elev.sd=elev.sd.var,
                              can=can.var,temp.hmin=temp.hmin.var,
                              prec.hmin=prec.hmin.var)

#excluding NAs
samp.points.rat<- samp.points.rat[-c(which(is.na(temp.sd.var))),]

#Decreasing collinearity between covariates (criterion: r.pearson < 0.7)
head(samp.points.rat)
cor(samp.points.rat[,4:13])#
#exclude: aet, elev, prec.sd,prec.hmin, elev.sd, temp.hmin
cor(samp.points.rat[,c(5,6,7,11)])
#maintain:temp+temp.sd+prec+can
samp.points.rat2=samp.points.rat[,c(1,2,3,5,6,7,11)]

#Spatializing "sampled points" 
samp.points.rat.sp <- samp.points.rat3 
sp::coordinates(samp.points.rat.sp) <- ~x+y
proj4string(samp.points.rat.sp) <- "+proj=aitoff +ellps=WGS84"

#full points: geographic coordinates and and co-variates from all cells
full.points <- data.frame(aet=data.frame(rasterToPoints(env.clim[[1]],spatial=T))[,c(2,3,1)],
                          temp=data.frame(rasterToPoints(env.clim[[2]],spatial=T))[,1],
                          temp.sd=data.frame(rasterToPoints(env.clim[[3]],spatial=T))[,1],
                          prec=data.frame(rasterToPoints(env.clim[[4]],spatial=T))[,1],
                          prec.sd=data.frame(rasterToPoints(env.clim[[5]],spatial=T))[,1],
                          elev=data.frame(rasterToPoints(env.clim[[6]],spatial=T))[,1],
                          elev.sd=data.frame(rasterToPoints(env.clim[[7]],spatial=T))[,1],
                          can=data.frame(rasterToPoints(env.clim[[8]],spatial=T))[,1],
                          temp.hmin=data.frame(rasterToPoints(env.clim[[9]],spatial=T))[,1],
                          prec.hmin=data.frame(rasterToPoints(env.clim[[10]],spatial=T))[,1])

#spatializing "full points"
colnames(full.points) <- c("x", "y", "aet", "temp","temp.sd", "prec","prec.sd",
                           "elev","elev.sd","can","temp.hmin","prec.hmin")

head(full.points)
full.points.sp <- full.points 
coordinates(full.points.sp) <- ~x+y
proj4string(full.points.sp) <- "+proj=aitoff +ellps=WGS84"


#--- 4- Regression-Kriging ---# 

#selecting the best sample variogram model
(sample.var.rat <- autofitVariogram(rich~prec+temp.sd+can, input_data = samp.points.rat.sp))

#obs: insert r-k model's parameters (e.g., psill, range...) from "sample.var.rat"
rk.rat <- krige(formula=rich~prec+temp.sd+can, locations=samp.points.rat.sp, newdata=full.points.sp,#inserir os parametros do melhor modelo do semi-variograma dos residuos aqui: veja em "r.v.ran.09"
                model=vgm(psill=208.7099,model="Ste",range= 0.4592364   , 
                          nugget=0,kappa=1.7))
proj4string(rk.rat) <- "+proj=aitoff +ellps=WGS84" 

#estimated richness
est.rich <- rasterFromXYZ(rk.rat)#estimated richness
est.rich2 <- est.rich#estimated richness
est.rich2[is.na(est.rich2)] <- 0#transforming NA's in 0
est.rich2[est.rich2<0] <- 0#estimated richness without negative values
est.rich2=mask(est.rich2,temp)

x11()
plot(est.rich2, add=T,xaxt="n", yaxt="n", box=F, axes=F,zlim=c(0, maxValue(richness)), main="Estimated Richness")

#observed richness 
richness[is.na(richness)] <- 0
richness[richness<0] <- 0#estimated richness without negative values
richness<-mask(richness,temp)

plot(richness, add=T,main="Observed Richness",xaxt="n", yaxt="n", box=F, axes=0)

#residual
richness2=richness
richness2[is.na(richness2)]=0
resid <- richness - est.rich2#residuals
proj4string(resid) <- "+proj=aitoff +ellps=WGS84"
resid <- mask(resid,temp)

plot(resid,add=T, col=heat.colors(50, alpha=1),xaxt="n", yaxt="n", box=F, axes=0,main="Residual")

#Correlation between estimated richness and observed richness
rst.ob.rich <- getValues(richness)#observed rich vector
rst.est.rich <- getValues(est.rich2)#estimated rich vector
rst.resid <- getValues(resid)#residual raster vector
variab <- data.frame(ob=rst.ob.rich, est=rst.est.rich, res=rst.resid)
variab <- variab[-which(is.na(variab[,2])),]
cuts <- classIntervals(variab[,3], style="fixed", fixedBreaks=variab[,3] )
plotclr <- heat.colors(5,alpha=0.7)
colcode <- findColours(cuts, plotclr)

plot(x=variab[,2], y=variab[,1], col=colcode, xlim=c(0,80),
     ylim=c(0,120), pch=18, xlab="Estimated Richness", 
     ylab="Observed Richness", cex=1.2, cex.axis=1.3, cex.lab=1.3)
abline(reg=c(0,1), lty=3)





