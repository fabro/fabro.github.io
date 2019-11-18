#Script to simulate a virtual geographic gradient of species richness



#Summary:
#1- Co-variables: temperature and precipitation
#2- Simulating 100 species' geographic distributions 
  #(according to Roberts et al. 2017) 
#3- Generating the geographic gradient of species richness



#--------------------------------------------------------------------------------#

#1- Co-variables: temperature and precipitation

#Reading the environmental rasters:
temp <- raster("Temp_SA.asc") 
prec <- raster("Prec_SA.asc")

#Stack
predictors <- stack(temp, prec)#Joining both rasters 

#Data.frames:
temp.sp <- rasterToPoints(temp)
temp.df <- as.data.frame(temp.sp)

prec.sp <- rasterToPoints(prec)
prec.df <- as.data.frame(prec.sp)

grd <- temp.sp[,1:2] 

#------------------------------------------------------------------------#

#2- Simulating species' geographic ranges (according to the scripts available in
#Roberts et al. 2017)

n.sims=100

for(i in 1:n.sims){
  
  ## FIRST, simulate the spatially structured environment
  
  # Some other env var is just standard (could be soil or topography or whatever)
  expCov <- RMexp(var = 0.3, scale = 5)
  x.1 <- as.vector(t(RFsimulate(expCov, x = grd[,1], y = grd[,2], spConform = FALSE)))
  x.1 <- scale(x.1)
  
  # "Precip" has a higher variance
  #expCov <- RMexp(var = 0.3, scale = 0.1)
  #x.2 <- as.vector(t(RFsimulate(expCov, x = grd[,1], y = grd[,2], spConform = FALSE)))
  #x.2 <- scale(x.2)
  x.2 <- scale(prec.df[,3])
  
  # "Temp" has a longer range
  #expCov <- RMgauss(var = 0.1, scale = 0.3) #0.4 was a bit much
  #x.3 <- as.vector(t(RFsimulate(expCov, x = grd[,1], y = grd[,2], spConform = FALSE)))
  #x.3 <- scale(x.3)
  x.3 <- scale(temp.df[,3])
  
  ## SECOND, simulate the biotic interaction and the disease
  
  # Prevalence of disease
  #p.d <- 0.05
  
  # Add constant to x.2 and x.3 to keep them positive
  #t.ps <- -floor(min(x.2, x.3))
  
  # Set up disease with 1 where it occurs
  #x.4 <- rep(1, prod(gridsize))
  
  # Change locations to zero where the ratio of precip to temp isn't large enough to be 
  # in the percentile specified by prevalence (the top wet and warm locations have disease)
  #x.4[((x.2 + t.ps)/(x.3 + t.ps)) < quantile(((x.2 + t.ps)/(x.3 + t.ps)), (1 - p.d))] <- 0
  
  # Putting in "food" as a linear combination of previous variables
  x.5 <- (x.1 + x.2 + x.3 + x.2*x.3)
  x.5 <- scale(x.5)
  
  #  Add more "normal" variables as additional material for overfitting
  #expCov <- RMexp(var = 0.1, scale = 0.1)
  #x.6 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  #x.6 <- scale(x.6)
  #x.7 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  #x.7 <- scale(x.7)
  #x.8 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  #x.8 <- scale(x.8)
  #expCov <- RMgauss(var = 0.1, scale = 0.3)
  #x.9 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  #x.9 <- scale(x.9)
  #expCov <- RMgauss(var = 0.1, scale = 0.3)
  #x.10 <- as.vector(t(RFsimulate(expCov, x = Xvec, y = Yvec, spConform = FALSE)))
  #x.10 <- scale(x.10)
  
  
  ## THIRD, smulate the species
  
  # Species should depend linearly on food (x.5).
  # Species should be limited by an interaction between precip and temp.
  # In reality, water potential depends linearly on temp, but temp is in K, so the
  # real variance is small and relatively unimportant compared to precip. I will simulate that
  # here by turning the simulated temp into somewhat realistic temp and then doing calcs (precip/temp)
  # as if it were in K (so add 273 first and in the end standardize the variable again)
  # dependence on temp is unimodal f(x) = 1/(sqrt(2*pi)*sigma)*e^(-((x - mean)^2/(2*sigma^2)))
  # I'll have sigma = 1 here and mean = 0
  
  # "water availability"
  #x.11 <- x.2/(x.3 + 273)
  #x.11 <- scale(x.11)
  
  # Gaussian dependence on temperature
  x.12 <- 1/(sqrt(2*pi))*exp(-(x.3^2/4))
  x.12 <- scale(x.12)
  
  # Gaussian dependence on water
  x.13 <- 1/(sqrt(2*pi))*exp(-(x.2^2/4))
  x.13 <- scale(x.13)
  
  # x.1 = some standard var such as soil - unknown to model
  # x.5 = "food" a combination of x.1, x.2, x.3 and x.2*x.3
  # x.12 and x.13, Gaussian response of temp and precip
  y <- x.1 + x.5 + x.12 + x.13 #+ x.6
  #y <- x.12 + x.13
  y <- scale(y)
  
  # Then use the water potential as limiting factor by reducing y to water potential where 
  # water potential is lower (x.11, "water potential" derived from x.2 and x.3)
  #y[y>x.11] <- x.11[y>x.11]
  #y[x.4 == 1] <- min(y)
  
  # Empty list for simulated landscapes (created on first loop)
  if(i==1){sim.data <- list()}
  sim.data[[i]] <- data.frame(y, mget(paste0("x.", c(1,2,3,5,12,13))))
  
} 

#Binarizing species' geographic ranges 

sim.data2<-list()
for(i in 1:length(sim.data)){
  sim.data2[[i]] <- sim.data[[i]][,1] > -0.25#range threshold
}


#-----------------------------------------------------------------------#

#3- Generating the geographic gradient of species richness

sim.data3<-rowSums(do.call(cbind,sim.data2))#spp richness
mapaRiq=rasterize(grd,raster(x=temp),as.numeric(sim.data3))#rasterizing spp richness
rich.df <- as.data.frame(rasterToPoints(mapaRiq))
plot(mapaRiq)

#Objects: richness and co-variables for all cells
full.df <- data.frame(rich.df, temp.df[,3], prec.df[,3])#data frame
colnames(full.df) <- c("x","y","rich","temp","prec")
full.sp <- full.df#spatial data.frame
coordinates(full.sp) <- ~x+y
full.sp2 <- full.sp#object with spatial projection
proj4string(full.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#Geographic coordinates for all cells
pred.loc.df <- full.df[,c(1,2)]#data frame
pred.loc.sp <- as.data.frame(pred.loc.df)
coordinates(pred.loc.sp) <- ~x+y#spatial data.frame
pred.loc.sp2 <- pred.loc.sp#object with spatial projection
proj4string(pred.loc.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs" 

