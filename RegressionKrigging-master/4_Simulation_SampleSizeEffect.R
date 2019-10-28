#Script to evaluate the effect of sample size (number of "sampled" cells) on 
#the accuracy of the Regression-Kriging

#----------------------------------------------------#
#Summary:
#1-Data
#2-Function
#3-Results:
  #for the non-random shortfall scenario
  #for the random shortfall scenario

#-------------------------------------------------------

#1-Data:

#full spatial data.frame
full.sp2
#full observed richness matrix
rich.df
#sample matrix (x,y,rich,temp,prec)

#non-random scenario
samp.hum.06
samp.hum.065
samp.hum.07
samp.hum.075
samp.hum.080
samp.hum.085
samp.hum.09
samp.hum.095

#random scenario
samp.ran.06
samp.ran.065
samp.ran.07
samp.ran.075
samp.ran.080
samp.ran.085
samp.ran.09
samp.ran.095

#------------------------------------------------------------------------#

#2-Function

#Description: evaluate the sample size effect upon RMSE (of Regression-Kriging) 
#with n iterations

#arguments: 
#data: sample matrix (x,y,rich,temp,prec)
#n1: sample size
#n2: number of iterations

#Value: RMSE for each iteration

sam.size.ef <- function(data, n1, n2){

resu <- numeric()#receive the RMSE's from loop  
  
for(i in 1:n2){#iterations  

print(i)
  
#curating data
data.n <- data[sample(dim(data)[1],n1),]#sample table
data.n.sp <- data.n
coordinates(data.n.sp) <- ~x+y
data.n.sp2 <- data.n.sp 
proj4string(data.n.sp2) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"

#fitting variogram models
r.v.data <- autofitVariogram(rich~temp+prec, input_data = data.n.sp2)

#r-kriging
if(any(colnames(r.v.data$var_model)=="kappa")){

r.krig.data <- krige(formula=rich~temp+prec, locations=data.n.sp2, newdata=full.sp2,
                     model=vgm(psill=r.v.data$var_model$psill[2],
                     model=paste(r.v.data$var_model$model[2],"",sep=""),
                     range=r.v.data$var_model$range[2],nugget=r.v.data$var_model$psill[1],
                     kappa=r.v.data$var_model$kappa[2]))

}else{
  
r.krig.data <- krige(formula=rich~temp+prec, locations=data.n.sp2, newdata=full.sp2,
                       model=vgm(psill=r.v.data$var_model$psill[2],
                                 model=paste(r.v.data$var_model$model[2],"",sep=""),
                                 range=r.v.data$var_model$range[2],nugget=r.v.data$var_model$psill[1]))
  
}

#RMSE
r.krig.data.df <- as.data.frame(r.krig.data)
r.krig.data.df.rich <- ifelse(r.krig.data.df[,3]<0,0,r.krig.data.df[,3])
r.krig.data.df.rich <- data.frame(r.krig.data.df[,1:2],r.krig.data.df.rich)
r.krig.data.df.rich.diff <- r.krig.data.df.rich[,3] - rich.df[,3]
resu[i] <- sqrt((sum((r.krig.data.df.rich.diff)^2))/length(r.krig.data.df.rich.diff))#accuracy (RMSE)

}

return(resu)

}


#-----------------------------------------------------------------------#

#3-Results:

#---Non-random shortfall scenario---#

#size: 10
res.60.10 <- sam.size.ef(samp.hum.06,10,100)
res.65.10 <- sam.size.ef(samp.hum.065,10,100)
res.70.10 <- sam.size.ef(samp.hum.07,10,100)
res.75.10 <- sam.size.ef(samp.hum.075,10,100)
res.80.10 <- sam.size.ef(samp.hum.080,10,100)
res.85.10 <- sam.size.ef(samp.hum.085,10,100)
res.90.10 <- sam.size.ef(samp.hum.09,10,100)
res.95.10 <- sam.size.ef(samp.hum.095,10,100)

#size: 50
res.60.50 <- sam.size.ef(samp.hum.06,50,100)
res.65.50 <- sam.size.ef(samp.hum.065,50,100)
res.70.50 <- sam.size.ef(samp.hum.07,50,100)
res.75.50 <- sam.size.ef(samp.hum.075,50,100)
res.80.50 <- sam.size.ef(samp.hum.080,50,100)
res.85.50 <- sam.size.ef(samp.hum.085,50,100)

#size: 100
res.60.100 <- sam.size.ef(samp.hum.06,100,100)
res.65.100 <- sam.size.ef(samp.hum.065,100,100)
res.70.100 <- sam.size.ef(samp.hum.07,100,100)
res.75.100 <- sam.size.ef(samp.hum.075,100,100)
res.80.100 <- sam.size.ef(samp.hum.080,100,100)
res.85.100 <- sam.size.ef(samp.hum.085,100,100)

#Joining the results
ss10 <- cbind(res.60.10,res.65.10,res.70.10,res.75.10,res.80.10,res.85.10,res.90.10,res.95.10) 
ss100 <- cbind(res.60.100,res.65.100,res.70.100,res.75.100,res.80.100,res.85.100) 
ss50 <- cbind(res.60.50,res.65.50,res.70.50,res.75.50,res.80.50,res.85.50) 

#plots

x11()

#100 cells
mean100 <- colMeans(ss100)
se100 <- apply(ss100, 2, function(x) sd(x)/sqrt(10))
ci100 <- rbind(mean100+(1.96*(se100)), mean100-(1.96*(se100)))  
comp <- c(60, 65, 70, 75, 80, 85, 90, 95)
comp2 <- c(60, 65, 70, 75, 80, 85) 
plotCI(comp2, mean100, ui=ci100[1,], ylim=c(17,36), li=ci100[2,], pch=16, xlab="Completeness", ylab="RMSE")

#50 cells
mean50 <- colMeans(ss50)
se50 <- apply(ss50, 2, function(x) sd(x)/sqrt(10))
ci50 <- rbind(mean50+(1.96*(se50)), mean50-(1.96*(se50)))  
plotCI(comp2, mean50, ui=ci50[1,], li=ci50[2,], ylim=c(17,36), pch=16, xlab="Completeness", ylab="RMSE")

#10 cells
mean10 <- colMeans(ss10)
se10 <- apply(ss10, 2, function(x) sd(x)/sqrt(10))
ci10 <- rbind(mean10+(1.96*(se10)), mean10-(1.96*(se10)))  
comp <- c(60, 65, 70, 75, 80, 85, 90, 95)
plotCI(comp, mean10, ui=ci10[1,], li=ci10[2,],ylim=c(17,36), pch=16, xlab="Completeness", ylab="RMSE")

#United plots
x11()
par(mfrow=c(1,3))
plotCI(comp, mean10, ui=ci10[1,], li=ci10[2,],ylim=c(17,36), pch=16,
       cex.axis = 1.2,cex.lab=1.2, xlab="Completeness", ylab="RMSE")
plotCI(comp, c(mean50,0,0), ui=c(ci50[1,],0,0), li=c(ci50[2,],0,0), ylim=c(17,36), pch=16,
        cex.axis = 1.2,cex.lab=1.2, xlab="Completeness", ylab="RMSE")
plotCI(comp, c(mean100,0,0), ui=c(ci100[1,],0,0), li=c(ci100[2,],0,0),ylim=c(17,36), pch=16,
       cex.axis = 1.2,cex.lab=1.2, xlab="Completeness", ylab="RMSE")
dev.off()

#---Random shortfall scenario---#

#size: 10
res2.60.10 <- sam.size.ef(samp.ran.06,10,100)
res2.65.10 <- sam.size.ef(samp.ran.065,10,100)
res2.70.10 <- sam.size.ef(samp.ran.07,10,100)
res2.75.10 <- sam.size.ef(samp.ran.075,10,100)
res2.80.10 <- sam.size.ef(samp.ran.080,10,100)
res2.85.10 <- sam.size.ef(samp.ran.085,10,100)
res2.90.10 <- sam.size.ef(samp.ran.09,10,100)
res2.95.10 <- sam.size.ef(samp.ran.095,10,100)

#size: 50
res2.60.50 <- sam.size.ef(samp.ran.06,50,100)
res2.65.50 <- sam.size.ef(samp.ran.065,50,100)
res2.70.50 <- sam.size.ef(samp.ran.07,50,100)
res2.75.50 <- sam.size.ef(samp.ran.075,50,100)
res2.80.50 <- sam.size.ef(samp.ran.080,50,100)
res2.85.50 <- sam.size.ef(samp.ran.085,50,100)
res2.90.50 <- sam.size.ef(samp.ran.09,50,100)
res2.95.50 <- sam.size.ef(samp.ran.095,50,100)

#size: 100
res2.60.100 <- sam.size.ef(samp.ran.06,100,100)
res2.65.100 <- sam.size.ef(samp.ran.065,100,100)
res2.70.100 <- sam.size.ef(samp.ran.07,100,100)
res2.75.100 <- sam.size.ef(samp.ran.075,100,100)
res2.80.100 <- sam.size.ef(samp.ran.080,100,100)
res2.85.100 <- sam.size.ef(samp.ran.085,100,100)
res2.90.100 <- sam.size.ef(samp.ran.09,100,100)
res2.95.100 <- sam.size.ef(samp.ran.095,100,100)

#Joining the results
ss2.10 <- cbind(res2.60.10,res2.65.10,res2.70.10,res2.75.10,res2.80.10,res2.85.10,res2.90.10,res2.95.10) 
ss2.100 <- cbind(res2.60.100,res2.65.100,res2.70.100,res2.75.100,res2.80.100,res2.85.100,res2.90.100,res2.95.100) 
ss2.50 <- cbind(res2.60.50,res2.65.50,res2.70.50,res2.75.50,res2.80.50,res2.85.50,res2.90.50,res2.95.50) 

#plots

x11()

#100 cells
mean2.100 <- colMeans(ss2.100)
se2.100 <- apply(ss2.100, 2, function(x) sd(x)/sqrt(10))
ci2.100 <- rbind(mean2.100+(1.96*(se2.100)), mean2.100-(1.96*(se2.100)))  
plotCI(comp, mean2.100, ui=ci2.100[1,], li=ci2.100[2,], pch=16, xlab="Completeness", ylab="RMSE")

#50 cells
mean2.50 <- colMeans(ss2.50)
se2.50 <- apply(ss2.50, 2, function(x) sd(x)/sqrt(10))
ci2.50 <- rbind(mean2.50+(1.96*(se2.50)), mean2.50-(1.96*(se2.50)))  
plotCI(comp, mean2.50, ui=ci2.50[1,], li=ci2.50[2,],  pch=16, xlab="Completeness", ylab="RMSE")

#10 cells
mean2.10 <- colMeans(ss2.10)
se2.10 <- apply(ss2.10, 2, function(x) sd(x)/sqrt(10))
ci2.10 <- rbind(mean2.10+(1.96*(se2.10)), mean2.10-(1.96*(se2.10)))  
plotCI(comp, mean2.10, ui=ci2.10[1,], li=ci2.10[2,], pch=16, xlab="Completeness", ylab="RMSE")

#United plots
x11()
par(mfrow=c(1,3))
plotCI(comp, mean2.10, ui=ci2.10[1,], li=ci2.10[2,],ylim=c(11,32), 
       pch=16, xlab="Completeness", ylab="RMSE", cex.axis = 1.2,cex.lab=1.2)
plotCI(comp, mean2.50, ui=ci2.50[1,], li=ci2.50[2,],ylim=c(11,32), 
       pch=16, xlab="Completeness", ylab="RMSE", cex.axis = 1.2,cex.lab=1.2)
plotCI(comp, mean2.100, ui=ci2.100[1,], li=ci2.100[2,],ylim=c(11,32),
       pch=16, xlab="Completeness", ylab="RMSE", cex.axis = 1.2,cex.lab=1.2)
dev.off()

