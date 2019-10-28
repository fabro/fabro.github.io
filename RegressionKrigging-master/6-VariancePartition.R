#Script to create a Variance partition of the simulated "real" richness 
#between the estimated richness from the Regrssion-Kriging and Ordinary
#Kriging model

#We made multiple regressions between the simulated "real' richness (response variable)
#with the estimated richness from r-kriging and kriging for each shortfall
#scenario and completeness percentage.

#Summary: 
  
  #1- simulated "real" richness
  #2- Random shortfall scenario
  #3- Non-random shortfall scenario
  #4- Barplots
  
#------------------------------------------------------------#

#1- simulated "real" richness
rich.vt = rasterToPoints(mapaRiq)


#--------------------------------------------------------------#
#2- Random shortfall scenario 

#60% completeness
#full matrix
rich.06 = data.frame(y=rich.vt[,3],y.k=krig.ran.06.df.rich[,3],y.r=r.krig.ran.06.df.rich[,3])
#regressions
part.06 = varpart(rich.06[,1], ~ y.k, ~ y.r,data=rich.06)

#65% completeness
#full matrix
rich.065 = data.frame(y=rich.vt[,3],y.k=krig.ran.065.df.rich[,3],y.r=r.krig.ran.065.df.rich[,3])
#regressions
part.065 = varpart(rich.065[,1], ~ y.k, ~ y.r,data=rich.065)

#70% completeness
#full matrix
rich.07 = data.frame(y=rich.vt[,3],y.k=krig.ran.07.df.rich[,3],y.r=r.krig.ran.07.df.rich[,3])
#regressions
part.07 = varpart(rich.07[,1], ~ y.k, ~ y.r,data=rich.07)

#75% completeness
#full matrix
rich.075 = data.frame(y=rich.vt[,3],y.k=krig.ran.075.df.rich[,3],y.r=r.krig.ran.075.df.rich[,3])
#regressions
part.075 = varpart(rich.075[,1], ~ y.k, ~ y.r,data=rich.075)

#80% completeness
#full matrix
rich.08 = data.frame(y=rich.vt[,3],y.k=krig.ran.080.df.rich[,3],y.r=r.krig.ran.080.df.rich[,3])
#regressions
part.08 = varpart(rich.08[,1], ~ y.k, ~ y.r,data=rich.08)

#85% completeness
#full matrix
rich.085 = data.frame(y=rich.vt[,3],y.k=krig.ran.085.df.rich[,3],y.r=r.krig.ran.085.df.rich[,3])
#regressions
part.085 = varpart(rich.085[,1], ~ y.k, ~ y.r,data=rich.085)

#90% completeness
#full matrix
rich.09 = data.frame(y=rich.vt[,3],y.k=krig.ran.09.df.rich[,3],y.r=r.krig.ran.09.df.rich[,3])
#regressions
part.09 = varpart(rich.09[,1], ~ y.k, ~ y.r,data=rich.09)

#95% completeness
#full matrix
rich.095 = data.frame(y=rich.vt[,3],y.k=krig.ran.095.df.rich[,3],y.r=r.krig.ran.095.df.rich[,3])
#regressions
part.095 = varpart(rich.095[,1], ~ y.k, ~ y.r,data=rich.095)

#uniting all partitions
ran.part <- cbind(part.06$part$indfract$Adj.R.squared,
                  part.065$part$indfract$Adj.R.squared,
                  part.07$part$indfract$Adj.R.squared,
                  part.075$part$indfract$Adj.R.squared,
                  part.08$part$indfract$Adj.R.squared,
                  part.085$part$indfract$Adj.R.squared,
                  part.09$part$indfract$Adj.R.squared,
                  part.095$part$indfract$Adj.R.squared)
colnames(ran.part) <- c("60","65","70","75","80","85","90","95")
row.names(ran.part) <- c("a","b","c","d") 



#---------------------------------------------------------------------#

#3- Non-random shortfall scenario

#60% completeness
#full matrix
s.rich.06 = data.frame(y=rich.vt[,3],y.k=krig.hum.06.df.rich[,3],y.r=r.krig.hum.06.df.rich[,3])
#regressions
s.part.06 = varpart(s.rich.06[,1], ~ y.k, ~ y.r,data=s.rich.06)

#65% completeness
#full matrix
s.rich.065 = data.frame(y=rich.vt[,3],y.k=krig.hum.065.df.rich[,3],y.r=r.krig.hum.065.df.rich[,3])
#regressions
s.part.065 = varpart(s.rich.065[,1], ~ y.k, ~ y.r,data=s.rich.065)

#70% completeness
#full matrix
s.rich.07 = data.frame(y=rich.vt[,3],y.k=krig.hum.07.df.rich[,3],y.r=r.krig.hum.07.df.rich[,3])
#regressions
s.part.07 = varpart(s.rich.07[,1], ~ y.k, ~ y.r,data=s.rich.07)

#75% completeness
#full matrix
s.rich.075 = data.frame(y=rich.vt[,3],y.k=krig.hum.075.df.rich[,3],y.r=r.krig.hum.075.df.rich[,3])
#regressions
s.part.075 = varpart(s.rich.075[,1], ~ y.k, ~ y.r,data=s.rich.075)

#80% completeness
#full matrix
s.rich.08 = data.frame(y=rich.vt[,3],y.k=krig.hum.080.df.rich[,3],y.r=r.krig.hum.080.df.rich[,3])
#regressions
s.part.08 = varpart(s.rich.08[,1], ~ y.k, ~ y.r,data=s.rich.08)

# 85%
#full matrix
s.rich.085 = data.frame(y=rich.vt[,3],y.k=krig.hum.085.df.rich[,3],y.r=r.krig.hum.085.df.rich[,3])
#regressions

s.part.085 = varpart(s.rich.085[,1], ~ y.k, ~ y.r,data=s.rich.085)

#90% completeness
#full matrix
s.rich.09 = data.frame(y=rich.vt[,3],y.k=krig.hum.09.df.rich[,3],y.r=r.krig.hum.09.df.rich[,3])
#regressions
s.part.09 = varpart(s.rich.09[,1], ~ y.k, ~ y.r,data=s.rich.09)

#95% completeness
#full matrix
s.rich.095 = data.frame(y=rich.vt[,3],y.k=krig.hum.095.df.rich[,3],y.r=r.krig.hum.095.df.rich[,3])
#regressions
s.part.095 = varpart(s.rich.095[,1], ~ y.k, ~ y.r,data=s.rich.095)

#uniting the partitions
str.part <- cbind(s.part.06$part$indfract$Adj.R.squared,
                  s.part.065$part$indfract$Adj.R.squared,
                  s.part.07$part$indfract$Adj.R.squared,
                  s.part.075$part$indfract$Adj.R.squared,
                  s.part.08$part$indfract$Adj.R.squared,
                  s.part.085$part$indfract$Adj.R.squared,
                  s.part.09$part$indfract$Adj.R.squared,
                  s.part.095$part$indfract$Adj.R.squared)
colnames(str.part) <- c("60","65","70","75","80","85","90","95")
row.names(str.part) <- c("a","b","c","d") 


#---------------------------------------------------------------------#

#4- barplots
x11()
par(mfrow=c(1,2))

color=c("black","red","blue","green")
barplot(str.part,col=color)
legend(0.2,0.5, legend=row.names(str.part),fill=color, cex = 1,bg="white")

barplot(ran.part,col=color)
legend(0.2,0.5, row.names(ran.part),fill=color, cex = 1,bg="white")
