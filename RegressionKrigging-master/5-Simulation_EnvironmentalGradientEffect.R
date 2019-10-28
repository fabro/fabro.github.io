#Script to evaluate the effect of selecting sampled cells at different
#environmnetal distances in the accuracy of the Regression-Kriging model

#Summary:

  #1- Environmental filtering function (according to Varela et al. 2014)
  #2- "Environmental gradient effect" function
  #3- Application in the Non-random shortfall scenario
  #4- Application in the random shortfall scenario

#-----------------------------------------------------------------------#

#1- Environmental filtering function (according to Varela et al. 2014)

#reference: https://github.com/SaraVarela/envSample/blob/master/R/envSample.R

#Description: select points in a environmental space at a given resolution.

#Arguments:
#coord: geographic coordinates of the points
#filters: two environmental vectors to create the environmental space
#res: resolution of the environmental grid

#Value: matrix with the geographic coordinates of the selected points
envSample<- function (coord, filters, res, do.plot=TRUE){
  
  n<- length (filters)
  pot_points<- list ()
  for (i in 1:n){
    k<- filters [[i]] [!is.na(filters[[i]])]
    ext1<- range (k)
    ext1 [1]<- ext1[1]- 1
    x<- seq(ext1[1],ext1[2], by=res[[i]])
    pot_points[[i]]<- x
  }
  pot_p<- expand.grid(pot_points)
  
  ends<- NULL
  for (i in 1:n){
    fin<- pot_p [,i] + res[[i]]
    ends<- cbind (ends, fin)
  }
  
  pot_pp<- data.frame (pot_p, ends)
  pot_pp<- data.frame (pot_pp, groupID=c(1:nrow (pot_pp)))
  rows<- length (filters[[1]])
  filter<- data.frame(matrix(unlist(filters), nrow=rows))
  real_p<- data.frame (coord, filter)
  
  names_real<- c("lon", "lat")
  names_pot_st<- NULL
  names_pot_end<- NULL
  sql1<- NULL
  for (i in 1:n){
    names_real<- c(names_real, paste ("filter", i, sep=""))
    names_pot_st<- c(names_pot_st, paste ("start_f", i, sep=""))
    names_pot_end<- c(names_pot_end, paste ("end_f", i, sep=""))
    sql1<- paste (sql1, paste ("real_p.filter", i, sep=""), sep=", ")   
  }
  
  names (real_p)<- names_real
  names (pot_pp)<- c(names_pot_st, names_pot_end, "groupID")
  
  conditions<- paste ("(real_p.filter", 1, "<= pot_pp.end_f", 1,") and (real_p.filter", 1, "> pot_pp.start_f", 1, ")", sep="")
  for (i in 2:n){
    conditions<- paste (conditions, 
                        paste ("(real_p.filter", i, "<= pot_pp.end_f", i,") and (real_p.filter", i, "> pot_pp.start_f", i, ")", sep=""), 
                        sep="and")
  }
  
  selection_NA<- sqldf(paste ("select real_p.lon, real_p.lat, pot_pp.groupID",   
                              sql1, "from pot_pp left join real_p on", conditions, sep=" "))
  
  selection<- selection_NA [complete.cases(selection_NA),]
  
  final_points<- selection[!duplicated(selection$groupID), ]
  coord_filter<- data.frame (final_points$lon, final_points$lat) 
  names (coord_filter)<- c("lon", "lat")
  
  if (do.plot==TRUE){
    par (mfrow=c(1,2), mar=c(4,4,0,0.5))
    plot (filters[[1]], filters[[2]], pch=19, 
          col="grey50", xlab="Filter 1", ylab="Filter 2")
    points (final_points$filter1, final_points$filter2, 
            pch=19, col="#88000090")
    plot (coord, pch=19, col="grey50")
    map(add=T)
    points (coord_filter, pch=19, col="#88000090")
    
  }
  coord_filter
}

#-------------------------------------------------------------------------#

#2- "Environmental gradient effect" function

#Description: evaluate the environmental gradient effect upon RMSE 
#(of Regression-Kriging) with n1 iterations and n2 samples (based on envSamp;
#Varela et al. 2014)

#arguments: 
#samp.ran: sample matrix (x,y,rich,temp,prec)
#n1: number of iterations
#n2: number of samples
#res1: resolution of temperature
#resolution of precipitation

#Value: RMSE for each iteration


env.gr.ef <- function(samp.ran, n1, n2, res1, res2){
  
  resu <- numeric()#receive the RMSE's from loop  
  
  for(i in 1:n1){#iterations  
    
    print(i)
    
    #geographic coordinates
    envsamp1 = samp.ran[,-c(3:5)]
    #environmental coordinates
    envsamp2 = samp.ran[,4:5]
    #envSample
    envsamp3 = envSample(envsamp1,filters=list(envsamp2[,1],envsamp2[,2]),
                         res=list(res1,res2),do.plot=F)
    #n2 samples
    samp.temp = extract(temp, envsamp3)
    samp.prec = extract(prec, envsamp3)
    samp.rich = extract(mapaRiq, envsamp3)
    envsamp4 <- data.frame(envsamp3, samp.rich, samp.temp, samp.prec)
    envsamp4 <- envsamp4[sample(1:dim(envsamp4)[1],n2),]
    #spatial points
    colnames(envsamp4) <- c("x","y","rich","temp","prec")
    envsamp5 <- envsamp4
    coordinates(envsamp5) <- ~x+y
    envsamp6 <- envsamp5 
    proj4string(envsamp6) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
    
    #fitting variogram models
    r.v.data <- autofitVariogram(rich~temp+prec, input_data = envsamp6)
    
    #r-kriging
    if(any(colnames(r.v.data$var_model)=="kappa")){
      
      r.krig.data <- krige(formula=rich~temp+prec, locations=envsamp6, newdata=full.sp2,
                           model=vgm(psill=r.v.data$var_model$psill[2],
                                     model=paste(r.v.data$var_model$model[2],"",sep=""),
                                     range=r.v.data$var_model$range[2],nugget=r.v.data$var_model$psill[1],
                                     kappa=r.v.data$var_model$kappa[2]))
      
    }else{
      
      r.krig.data <- krige(formula=rich~temp+prec, locations=envsamp6, newdata=full.sp2,
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

#--------------------------------------------------------------------------#

#3- Application to the Random shortfall scenario

#60%
resu.06 = env.gr.ef(samp.ran.06,100,50,5,50)
resu1.06 = env.gr.ef(samp.ran.06,100,50,10,100)
resu2.06 = env.gr.ef(samp.ran.06,100,50,15,150)
resu3.06 = env.gr.ef(samp.ran.06,100,50,20,200)
resu4.06 = env.gr.ef(samp.ran.06,100,50,25,250)
resu5.06 = env.gr.ef(samp.ran.06,100,50,30,300)
resu6.06 = env.gr.ef(samp.ran.06,100,50,35,350)
resu7.06 = env.gr.ef(samp.ran.06,100,50,40,400)

resu.06.mat = cbind(resu.06,resu1.06,resu2.06,
                  resu3.06, resu4.06, resu5.06, 
                  resu6.06, resu7.06)
mean.06 <- colMeans(resu.06.mat)
se.06 <- apply(resu.06.mat, 2, function(x) sd(x)/sqrt(10))
ci.06 <- rbind(mean.06+(1.96*(se.06)), mean.06-(1.96*(se.06)))  
env.grad <- c(1:8)
plotCI(env.grad, mean.06, ui=ci.06[1,], li=ci.06[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
dev.off()

#65%
resu.065 = env.gr.ef(samp.ran.065,100,50,5,50)
resu1.065 = env.gr.ef(samp.ran.065,100,50,10,100)
resu2.065 = env.gr.ef(samp.ran.065,100,50,15,150)
resu3.065 = env.gr.ef(samp.ran.065,100,50,20,200)
resu4.065 = env.gr.ef(samp.ran.065,100,50,25,250)
resu5.065 = env.gr.ef(samp.ran.065,100,50,30,300)
resu6.065 = env.gr.ef(samp.ran.065,100,50,35,350)
resu7.065 = env.gr.ef(samp.ran.065,100,50,40,400)

resu.065.mat = cbind(resu.065,resu1.065,resu2.065,
                    resu3.065, resu4.065, resu5.065, 
                    resu6.065, resu7.065)
mean.065 <- colMeans(resu.065.mat)
se.065 <- apply(resu.065.mat, 2, function(x) sd(x)/sqrt(10))
ci.065 <- rbind(mean.065+(1.96*(se.065)), mean.065-(1.96*(se.065)))  
env.grad <- c(1:8)
plotCI(env.grad, mean.065, ui=ci.065[1,], li=ci.065[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#70%
resu.07 = env.gr.ef(samp.ran.07,100,50,5,50)
resu1.07 = env.gr.ef(samp.ran.07,100,50,10,100)
resu2.07 = env.gr.ef(samp.ran.07,100,50,15,150)
resu3.07 = env.gr.ef(samp.ran.07,100,50,20,200)
resu4.07 = env.gr.ef(samp.ran.07,100,50,25,250)
resu5.07 = env.gr.ef(samp.ran.07,100,50,30,300)
resu6.07 = env.gr.ef(samp.ran.07,100,50,35,350)
resu7.07 = env.gr.ef(samp.ran.07,100,50,40,400)

resu.07.mat = cbind(resu.07,resu1.07,resu2.07,
                     resu3.07, resu4.07, resu5.07, 
                     resu6.07, resu7.07)
mean.07 <- colMeans(resu.07.mat)
se.07 <- apply(resu.07.mat, 2, function(x) sd(x)/sqrt(10))
ci.07 <- rbind(mean.07+(1.96*(se.07)), mean.07-(1.96*(se.07)))  
env.grad <- c(1:8)
plotCI(env.grad, mean.07, ui=ci.07[1,], li=ci.07[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#75%
resu.075 = env.gr.ef(samp.ran.075,100,50,5,50)
resu1.075 = env.gr.ef(samp.ran.075,100,50,10,100)
resu2.075 = env.gr.ef(samp.ran.075,100,50,15,150)
resu3.075 = env.gr.ef(samp.ran.075,100,50,20,200)
resu4.075 = env.gr.ef(samp.ran.075,100,50,25,250)
resu5.075 = env.gr.ef(samp.ran.075,100,50,30,300)
resu6.075 = env.gr.ef(samp.ran.075,100,50,35,350)
resu7.075 = env.gr.ef(samp.ran.075,100,50,40,400)

resu.075.mat = cbind(resu.075,resu1.075,resu2.075,
                    resu3.075, resu4.075, resu5.075, 
                    resu6.075, resu7.075)
mean.075 <- colMeans(resu.075.mat)
se.075 <- apply(resu.075.mat, 2, function(x) sd(x)/sqrt(10))
ci.075 <- rbind(mean.075+(1.96*(se.075)), mean.075-(1.96*(se.075)))  
env.grad <- c(1:8)
plotCI(env.grad, mean.075, ui=ci.075[1,], li=ci.075[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#80%
resu.080 = env.gr.ef(samp.ran.080,100,50,5,50)
resu1.080 = env.gr.ef(samp.ran.080,100,50,10,100)
resu2.080 = env.gr.ef(samp.ran.080,100,50,15,150)
resu3.080 = env.gr.ef(samp.ran.080,100,50,20,200)
resu4.080 = env.gr.ef(samp.ran.080,100,50,25,250)
resu5.080 = env.gr.ef(samp.ran.080,100,50,30,300)
resu6.080 = env.gr.ef(samp.ran.080,100,50,35,350)
resu7.080 = env.gr.ef(samp.ran.080,100,50,40,400)

resu.080.mat = cbind(resu.080,resu1.080,resu2.080,
                     resu3.080, resu4.080, resu5.080, 
                     resu6.080, resu7.080)
mean.080 <- colMeans(resu.080.mat)
se.080 <- apply(resu.080.mat, 2, function(x) sd(x)/sqrt(10))
ci.080 <- rbind(mean.080+(1.96*(se.080)), mean.080-(1.96*(se.080)))  
env.grad <- c(1:8)
plotCI(env.grad, mean.080, ui=ci.080[1,], li=ci.080[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#85%
resu.085 = env.gr.ef(samp.ran.085,100,50,5,50)
resu1.085 = env.gr.ef(samp.ran.085,100,50,10,100)
resu2.085 = env.gr.ef(samp.ran.085,100,50,15,150)
resu3.085 = env.gr.ef(samp.ran.085,100,50,20,200)
resu4.085 = env.gr.ef(samp.ran.085,100,50,25,250)
resu5.085 = env.gr.ef(samp.ran.085,100,50,30,300)
resu6.085 = env.gr.ef(samp.ran.085,100,50,35,350)
resu7.085 = env.gr.ef(samp.ran.085,100,50,40,400)

resu.085.mat = cbind(resu.085,resu1.085,resu2.085,
                     resu3.085, resu4.085, resu5.085, 
                     resu6.085, resu7.085)
mean.085 <- colMeans(resu.085.mat)
se.085 <- apply(resu.085.mat, 2, function(x) sd(x)/sqrt(10))
ci.085 <- rbind(mean.085+(1.96*(se.085)), mean.085-(1.96*(se.085)))  
env.grad <- c(1:8)
plotCI(env.grad, mean.085, ui=ci.085[1,], li=ci.085[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#90%
resu.09 = env.gr.ef(samp.ran.09,100,50,5,50)
resu1.09 = env.gr.ef(samp.ran.09,100,50,10,100)
resu2.09 = env.gr.ef(samp.ran.09,100,50,15,150)
resu3.09 = env.gr.ef(samp.ran.09,100,50,20,200)
resu4.09 = env.gr.ef(samp.ran.09,100,50,25,250)
resu5.09 = env.gr.ef(samp.ran.09,100,50,30,300)
resu6.09 = env.gr.ef(samp.ran.09,100,50,35,350)
resu7.09 = env.gr.ef(samp.ran.09,100,50,40,400)

resu.09.mat = cbind(resu.09,resu1.09,resu2.09,
                     resu3.09, resu4.09, resu5.09, 
                     resu6.09, resu7.09)
mean.09 <- colMeans(resu.09.mat)
se.09 <- apply(resu.09.mat, 2, function(x) sd(x)/sqrt(10))
ci.09 <- rbind(mean.09+(1.96*(se.09)), mean.09-(1.96*(se.09)))  
env.grad <- c(1:8)
plotCI(env.grad, mean.09, ui=ci.09[1,], li=ci.09[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#95%
resu.095 = env.gr.ef(samp.ran.095,100,50,5,50)
resu1.095 = env.gr.ef(samp.ran.095,100,50,10,100)
resu2.095 = env.gr.ef(samp.ran.095,100,50,15,150)
resu3.095 = env.gr.ef(samp.ran.095,100,50,20,200)
resu4.095 = env.gr.ef(samp.ran.095,100,50,25,250)
resu5.095 = env.gr.ef(samp.ran.095,100,50,30,300)
resu6.095 = env.gr.ef(samp.ran.095,100,50,35,350)
resu7.095 = env.gr.ef(samp.ran.095,100,50,40,400)

resu.095.mat = cbind(resu.095,resu1.095,resu2.095,
                    resu3.095, resu4.095, resu5.095, 
                    resu6.095,NA) #resu7.095)
mean.095 <- colMeans(resu.095.mat)
se.095 <- apply(resu.095.mat, 2, function(x) sd(x)/sqrt(10))
ci.095 <- rbind(mean.095+(1.96*(se.095)), mean.095-(1.96*(se.095)))  
env.grad <- c(1:8)
plotCI(env.grad, mean.095, ui=ci.095[1,], li=ci.095[2,], pch=16, xlab="Environmental distance", ylab="RMSE")


#all plots

x11()
par(mfrow=c(2,4))
plotCI(env.grad, mean.06, ui=ci.06[1,], li=ci.06[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, mean.065, ui=ci.065[1,], li=ci.065[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, mean.07, ui=ci.07[1,], li=ci.07[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, mean.075, ui=ci.075[1,], li=ci.075[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, mean.080, ui=ci.080[1,], li=ci.080[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, mean.085, ui=ci.085[1,], li=ci.085[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, mean.09, ui=ci.09[1,], li=ci.09[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, mean.095, ui=ci.095[1,], li=ci.095[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
dev.off()


#--------------------------------------------------------------------------#

#4- Application to the Non-random shortfall scenario

#60%
s.resu.06 = env.gr.ef(samp.hum.06,100,50,5,50)
s.resu1.06 = env.gr.ef(samp.hum.06,100,50,10,100)
s.resu2.06 = env.gr.ef(samp.hum.06,100,50,15,150)
s.resu3.06 = env.gr.ef(samp.hum.06,100,50,20,200)
s.resu4.06 = env.gr.ef(samp.hum.06,100,50,25,250)
s.resu5.06 = env.gr.ef(samp.hum.06,100,50,30,300)
s.resu6.06 = env.gr.ef(samp.hum.06,100,50,35,350)
s.resu7.06 = env.gr.ef(samp.hum.06,100,50,40,400)

s.resu.06.mat = cbind(s.resu.06,s.resu1.06,s.resu2.06,
                      s.resu3.06, s.resu4.06, s.resu5.06, 
                      s.resu6.06, s.resu7.06)
s.mean.06 <- colMeans(s.resu.06.mat)
s.se.06 <- apply(s.resu.06.mat, 2, function(x) sd(x)/sqrt(10))
s.ci.06 <- rbind(s.mean.06+(1.96*(s.se.06)), s.mean.06-(1.96*(s.se.06)))  
env.grad <- c(1:8)
plotCI(env.grad, s.mean.06, ui=s.ci.06[1,], li=s.ci.06[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
dev.off()

#65%
s.resu.065 = env.gr.ef(samp.hum.065,100,50,5,50)
s.resu1.065 = env.gr.ef(samp.hum.065,100,50,10,100)
s.resu2.065 = env.gr.ef(samp.hum.065,100,50,15,150)
s.resu3.065 = env.gr.ef(samp.hum.065,100,50,20,200)
s.resu4.065 = env.gr.ef(samp.hum.065,100,50,25,250)
s.resu5.065 = env.gr.ef(samp.hum.065,100,50,30,300)
s.resu6.065 = env.gr.ef(samp.hum.065,100,50,35,350)
s.resu7.065 = env.gr.ef(samp.hum.065,100,50,40,400)

s.resu.065.mat = cbind(s.resu.065,s.resu1.065,s.resu2.065,
                       s.resu3.065, s.resu4.065, s.resu5.065, 
                       s.resu6.065, s.resu7.065)
s.mean.065 <- colMeans(s.resu.065.mat)
s.se.065 <- apply(s.resu.065.mat, 2, function(x) sd(x)/sqrt(10))
s.ci.065 <- rbind(s.mean.065+(1.96*(s.se.065)), s.mean.065-(1.96*(s.se.065)))  
env.grad <- c(1:8)
plotCI(env.grad, s.mean.065, ui=s.ci.065[1,], li=s.ci.065[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#70%
s.resu.07 = env.gr.ef(samp.hum.07,100,50,5,50)
s.resu1.07 = env.gr.ef(samp.hum.07,100,50,10,100)
s.resu2.07 = env.gr.ef(samp.hum.07,100,50,15,150)
s.resu3.07 = env.gr.ef(samp.hum.07,100,50,20,200)
s.resu4.07 = env.gr.ef(samp.hum.07,100,50,25,250)
s.resu5.07 = env.gr.ef(samp.hum.07,100,50,30,300)
s.resu6.07 = env.gr.ef(samp.hum.07,100,50,35,350)
s.resu7.07 = env.gr.ef(samp.hum.07,100,50,40,400)

s.resu.07.mat = cbind(s.resu.07,s.resu1.07,s.resu2.07,
                      s.resu3.07, s.resu4.07, s.resu5.07, 
                      s.resu6.07, s.resu7.07)
s.mean.07 <- colMeans(s.resu.07.mat)
s.se.07 <- apply(s.resu.07.mat, 2, function(x) sd(x)/sqrt(10))
s.ci.07 <- rbind(s.mean.07+(1.96*(s.se.07)), s.mean.07-(1.96*(s.se.07)))  
env.grad <- c(1:8)
plotCI(env.grad, s.mean.07, ui=s.ci.07[1,], li=s.ci.07[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#75%
s.resu.075 = env.gr.ef(samp.hum.075,100,50,5,50)
s.resu1.075 = env.gr.ef(samp.hum.075,100,50,10,100)
s.resu2.075 = env.gr.ef(samp.hum.075,100,50,15,150)
s.resu3.075 = env.gr.ef(samp.hum.075,100,50,20,200)
s.resu4.075 = env.gr.ef(samp.hum.075,100,50,25,250)
s.resu5.075 = env.gr.ef(samp.hum.075,100,50,30,300)
s.resu6.075 = env.gr.ef(samp.hum.075,100,50,35,350)
s.resu7.075 = env.gr.ef(samp.hum.075,100,50,40,400)


s.resu.075.mat = cbind(s.resu.075,s.resu1.075,s.resu2.075,
                       s.resu3.075, s.resu4.075, s.resu5.075, 
                     s.resu6.075,s.resu7.075)
s.mean.075 <- colMeans(s.resu.075.mat)
s.se.075 <- apply(s.resu.075.mat, 2, function(x) sd(x)/sqrt(10))
s.ci.075 <- rbind(s.mean.075+(1.96*(s.se.075)), s.mean.075-(1.96*(s.se.075)))  
env.grad <- c(1:8)
plotCI(env.grad, s.mean.075, ui=s.ci.075[1,], li=s.ci.075[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#80%
s.resu.080 = env.gr.ef(samp.hum.080,100,45,5,50)
s.resu1.080 = env.gr.ef(samp.hum.080,100,45,10,100)
s.resu2.080 = env.gr.ef(samp.hum.080,100,45,15,150)
s.resu3.080 = env.gr.ef(samp.hum.080,100,45,20,200)
s.resu4.080 = env.gr.ef(samp.hum.080,100,45,25,250)
s.resu5.080 = env.gr.ef(samp.hum.080,100,45,30,300)
s.resu6.080 = env.gr.ef(samp.hum.080,100,45,35,350)
s.resu7.080 = env.gr.ef(samp.hum.080,100,45,40,400)

s.resu.080.mat = cbind(s.resu.080,s.resu1.080,s.resu2.080,
                     s.resu3.080, s.resu4.080, s.resu5.080, 
                     s.resu6.080,s.resu7.080)
s.mean.080 <- colMeans(s.resu.080.mat)
s.se.080 <- apply(s.resu.080.mat, 2, function(x) sd(x)/sqrt(10))
s.ci.080 <- rbind(s.mean.080+(1.96*(s.se.080)), s.mean.080-(1.96*(s.se.080)))  
env.grad3 <- c(1:8)
plotCI(env.grad, s.mean.080, ui=s.ci.080[1,], li=s.ci.080[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#85%
s.resu.085 = env.gr.ef(samp.hum.085,100,30,5,50)
s.resu1.085 = env.gr.ef(samp.hum.085,100,30,10,100)
s.resu2.085 = env.gr.ef(samp.hum.085,100,30,15,150)
s.resu3.085 = env.gr.ef(samp.hum.085,100,30,20,200)
s.resu4.085 = env.gr.ef(samp.hum.085,100,30,25,250)
s.resu5.085 = env.gr.ef(samp.hum.085,100,30,30,300)
s.resu6.085 = env.gr.ef(samp.hum.085,100,30,35,350)
s.resu7.085 = env.gr.ef(samp.hum.085,100,30,40,400)

s.resu.085.mat = cbind(s.resu.085,s.resu1.085,s.resu2.085,
                     s.resu3.085, s.resu4.085,s.resu5.085, 
                     s.resu6.085, s.resu7.085)
s.mean.085 <- colMeans(s.resu.085.mat)
s.se.085 <- apply(s.resu.085.mat, 2, function(x) sd(x)/sqrt(10))
s.ci.085 <- rbind(s.mean.085+(1.96*(s.se.085)), s.mean.085-(1.96*(s.se.085)))  
env.grad <- c(1:8)
plotCI(env.grad2, s.mean.085, ui=s.ci.085[1,], li=s.ci.085[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
dev.off()
#90%
s.resu.09 = env.gr.ef(samp.hum.09,100,15,5,50)
s.resu1.09 = env.gr.ef(samp.hum.09,100,15,10,100)
s.resu2.09 = env.gr.ef(samp.hum.09,100,15,15,150)
s.resu3.09 = env.gr.ef(samp.hum.09,100,15,20,200)
s.resu4.09 = env.gr.ef(samp.hum.09,100,15,25,250)
s.resu5.09 = env.gr.ef(samp.hum.09,100,15,30,300)
s.resu6.09 = env.gr.ef(samp.hum.09,100,15,35,350)
s.resu7.09 = env.gr.ef(samp.hum.09,100,15,40,400)

s.resu.09.mat = cbind(s.resu.09,s.resu1.09,s.resu2.09,
                    s.resu3.09, s.resu4.09, s.resu5.09, 
                    s.resu6.09, s.resu7.09)
s.mean.09 <- colMeans(s.resu.09.mat)
s.se.09 <- apply(s.resu.09.mat, 2, function(x) sd(x)/sqrt(10))
s.ci.09 <- rbind(s.mean.09+(1.96*(s.se.09)), s.mean.09-(1.96*(s.se.09)))  
env.grad <- c(1:8)
plotCI(env.grad, s.mean.09, ui=s.ci.09[1,], li=s.ci.09[2,], pch=16, xlab="Environmental distance", ylab="RMSE")

#95%
s.resu.095 = env.gr.ef(samp.hum.095,100,10,5,50)
s.resu1.095 = env.gr.ef(samp.hum.095,100,10,10,100)
s.resu2.095 = env.gr.ef(samp.hum.095,100,10,10,150)
s.resu3.095 = env.gr.ef(samp.hum.095,100,10,10,200)
s.resu4.095 = env.gr.ef(samp.hum.095,100,10,10,250)
s.resu5.095 = env.gr.ef(samp.hum.095,100,10,10,300)
s.resu6.095 = env.gr.ef(samp.hum.095,100,10,10,350)
s.resu7.095 = env.gr.ef(samp.hum.095,100,10,10,400)

s.resu.095.mat = cbind(s.resu.095,s.resu1.095,s.resu2.095,
                       s.resu3.095, s.resu4.095, s.resu5.095, 
                       s.resu6.095,s.resu7.095) 
s.mean.095 <- colMeans(s.resu.095.mat)
s.se.095 <- apply(s.resu.095.mat, 2, function(x) sd(x)/sqrt(10))
s.ci.095 <- rbind(s.mean.095+(1.96*(s.se.095)), s.mean.095-(1.96*(s.se.095)))  
env.grad <- c(1:8)
plotCI(env.grad, s.mean.095, ui=s.ci.095[1,], li=s.ci.095[2,], pch=16, xlab="Environmental distance", ylab="RMSE")


#all plots

x11()
par(mfrow=c(2,4))
plotCI(env.grad, s.mean.06, ui=s.ci.06[1,], li=s.ci.06[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, s.mean.065, ui=s.ci.065[1,], li=s.ci.065[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, s.mean.07, ui=s.ci.07[1,], li=s.ci.07[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, s.mean.075, ui=s.ci.075[1,], li=s.ci.075[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, s.mean.080, ui=s.ci.080[1,], li=s.ci.080[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, s.mean.085, ui=s.ci.085[1,], li=s.ci.085[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, s.mean.09, ui=s.ci.09[1,], li=s.ci.09[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
plotCI(env.grad, s.mean.095, ui=s.ci.095[1,], li=s.ci.095[2,], pch=16, xlab="Environmental distance", ylab="RMSE")
dev.off()

#obs: some datasets presented different sample sizes 
#80% completeness: 45 points
#85% completeness: 30 points 
#90% completeness: 15 points
#95% completeness: 10 points



