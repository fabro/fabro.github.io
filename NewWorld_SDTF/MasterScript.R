library(picante)
library(phytools)
library(phytools)
library(picante)
library(fossil)

######## DATA preparation #########

setwd("~/Your_PATH")
biosp<-read.csv('Biogeo_sp.csv',header=T)
smithtree<-read.tree('ALLMBMOD.TXT')
matriz<-create.matrix(biosp,tax.name="A_N",locality='codigo')#creates a community matrix
comat2<-t(matriz)
sum(comat2[1,])

# Calculate diversity measures
match<-match.phylo.comm(smithtree,comat2)
pd1<-pd(match$comm,match$phy,include.root=F)#Calculates phylogenetic Diversity
SESMPD<-ses.mpd(match$comm,cophenetic(match$phy),null.model='taxa.label',runs=999,iterations=1000)#calculates standarized size MPD (to be transformed later to NRI)
loessmod<-loess(PD~SR,data=pd1)#LOESS model to calculate rPD
rPD<-loessmod$residuals
NRI<-SESMPD$mpd.obs.z*-1#NRI= ses.mpd*-1
Phylomeasures<-cbind(NRI,rPD,pd1$SR,pd1$PD)
Phylom<-as.data.frame(Phylomeasures);names(Phylom)<-c('NRI','rPD','SR','PD')
bio_codes<-read.csv("~/Desktop/MasterScript/biogeo_codes.csv",header=T)#getting the site list with codes
NRTI1<-cbind(Phylom,bio_codes)
NRTI1<-NRTI1[,-5]
head(NRTI1)
setwd('~/Documents')
write.csv(NRTI1,"Phylom.csv")#creates biodiversity measurementes data.frame
#a quick view of the data:
tapply(NRTI1$NRI,NRTI1$biosp.biogeo,mean)
tapply(NRTI1$rPD,NRTI1$biosp.biogeo,mean)
####GAM
library(mgcv)
ABSLATP<-abs(NRTI1$lat)
MPEN<-gam(NRTI1$PD~s(ABSLATP))#GAM where PD is predicted by absolute latitude
summary(MPEN)
SEQ<-seq(min(ABSLATP),max(ABSLATP),length=50)
pred<-predict(MPEN,data.frame(ABSLATP=SEQ))
plot(abs(NRTI1$lat),NRTI1$PD)
lines(SEQ,pred,col='red')

####Spatial analyses
library(raster)
library(SpatialPack)
library(spData)
library(spdep)
library(spatialreg)
NRI_coords<-data.frame(NRTI1$lon,NRTI1$lat);names(NRI_coords)<-c('lon','lat')#
W01<-raster('/Volumes/HAZEDRIVE/zproyecciones/WorldClim_Actual/wc2.0_2.5m_bio/wc2-5/bio1.bil')
W12<-raster('/Volumes/HAZEDRIVE/zproyecciones/WorldClim_Actual/wc2.0_2.5m_bio/wc2-5/bio12.bil')
W01CC<-raster('/Volumes/HAZEDRIVE/zproyecciones/WorldClim_Paleo/cclgmbi_2-5m/cclgmbi1.tif')
W12CC<-raster('/Volumes/HAZEDRIVE/zproyecciones/WorldClim_Paleo/cclgmbi_2-5m/cclgmbi12.tif')
W01MR<-raster('/Volumes/HAZEDRIVE/zproyecciones/WorldClim_Paleo/mrlgmbi_2-5m/mrlgmbi1.tif')
W12MR<-raster('/Volumes/HAZEDRIVE/zproyecciones/WorldClim_Paleo/mrlgmbi_2-5m/mrlgmbi12.tif') 
####climatic extractions and calculations
w01<-extract(W01,NRI_coords)
w12<-extract(W12,NRI_coords)
w01cc<-extract(W01CC,NRI_coords)
w12cc<-extract(W12CC,NRI_coords)
w01mr<-extract(W01MR,NRI_coords)
w12mr<-extract(W12MR,NRI_coords)
w01mean<-(w01cc+w01mr)/2
w12mean<-(w12cc+w12mr)/2
AN01<-(w01-w01mean)
AN12<-(w12-w12mean)
IN01<-abs(AN01)
IN12<-abs(AN12)
####creating the spatial object for the SARs
long<-NRTI1$lon
lat<-NRTI1$lat
pd.data2<-as.data.frame(cbind(NRTI1$PD,NRTI1$rPD,NRTI1$NRI,w01,w12,w01mean,w12mean,AN01,AN12,IN01,IN12,long,lat))
dim(NRI_coords)
names(pd.data2)<-c('PD','rPD','NRI','W01','W12','W01P','W12P','AN01mean','AN12mean','IN01mean','IN12mean','longitude','latitude')
pd.data.comp<-na.omit(pd.data2)
pd.sp<-pd.data.comp
coordinates(pd.sp)=~longitude+latitude
crs(pd.sp)<-crs(W01)
crs_pd<-crs(pd.sp)
pd.spatial<-spTransform(pd.sp,crs_pd)
# to create neighbor matrix:
near<-knearneigh(pd.spatial,k=1)
nbr<-knn2nb(near)
tw<-nb2listw(nbr)
#####formulas for the SARs
presente01<-rPD~W01
presente12<-rPD~W12

pasado01<-rPD~W01P
pasado12<-rPD~W12P


InTemp<-rPD~IN01mean
InPrep<-rPD~IN12mean
ANTemp<-rPD~AN01mean
ANPrep<-rPD~AN12mean
####SAR models

sarpresente01<-errorsarlm(formula=presente01,data=pd.spatial,listw=tw)
sarpresente12<-errorsarlm(formula=presente12,data=pd.spatial,listw=tw)

sarpasado01<-errorsarlm(formula=pasado01,data=pd.spatial,listw=tw)
sarpasado12<-errorsarlm(formula=pasado12,data=pd.spatial,listw=tw)


SarAn01<-errorsarlm(formula=ANTemp,data=pd.spatial,listw=tw)
SarAn12<-errorsarlm(formula=ANPrep,data=pd.spatial,listw=tw)

SarIn01<-errorsarlm(formula=InTemp,data=pd.spatial,listw=tw)
SarIn12<-errorsarlm(formula=InPrep,data=pd.spatial,listw=tw)

#summary of the models


summary(sarpresente01,Nagelkerke=T)#rPD~current mean annual temperature
summary(sarpresente12,Nagelkerke=T)#rPD~current mean annual precipitation
summary(sarpasado01,Nagelkerke=T)#rPD~LGM mean annual temperature
summary(sarpasado12,Nagelkerke=T)#rPD~LGM mean annual precipitation
summary(SarAn01,Nagelkerke=T)#rPD~Anomaly of the mean annual temperature
summary(SarAn12,Nagelkerke=T)#rPD~Anomaly of the mean annual precipitation
summary(SarIn01,Nagelkerke=T)#rPD~Instability of the mean annual temperature
summary(SarIn12,Nagelkerke=T)#rPD~Instability of the mean annual precipitation


####Biogeographic pools
rm(list=ls())
setwd("~/Desktop/MasterScript")
smithtree<-read.tree('ALLMBMOD.TXT')
biosp<-read.csv("biogeo_sp.csv",header=T)
N1<-biosp[which(biosp$biogeo=="Mesoamerica"),]
N2<-biosp[which(biosp$biogeo=="North southamerica"),]
N3<-biosp[which(biosp$biogeo=="Caribean"),]
N4<-rbind(N1,N2,N3)
northpool<-rbind(N1,N2);write.csv(northpool,"northpool.csv");northpool<-read.csv("northpool.csv",header=T)#creates the Mesoamerica North Southamerica pool
southpool<-biosp[-which(!is.na(match(biosp$biogeo,N4$biogeo))),];write.csv(southpool,"southpool.csv");southpool<-read.csv("southpool.csv",header=T) #creates the southamerican pool
northsites<-unique(data.frame(northpool$codigo,northpool$site,northpool$lon,northpool$lat,northpool$biogeo))
southsites<-unique(data.frame(southpool$codigo,southpool$site,southpool$lon,southpool$lat,southpool$biogeo))##creates a site list for every pool



northmatrix<-t(create.matrix(northpool,tax.name ="A_N",locality = "codigo" ))
southmatrix<-t(create.matrix(southpool,tax.name ="A_N",locality = "codigo" ))
northmatch<-match.phylo.comm(smithtree,northmatrix)
southmatch<-match.phylo.comm(smithtree,southmatrix)
sesnorth<-ses.mpd(northmatch$comm,cophenetic(northmatch$phy),null.model = "taxa.labels",runs=999,iterations=1000)
sessouth<-ses.mpd(southmatch$comm,cophenetic(southmatch$phy),null.model = "taxa.labels",runs=999,iterations=1000)
NRInorth<-data.frame(row.names(northmatrix),-1*sesnorth$mpd.obs.z);names(NRInorth)<-c("codigo","NRI")
NRIsouth<-data.frame(row.names(southmatrix),-1*sessouth$mpd.obs.z );names(NRIsouth)<-c("codigo","NRI")


#creating a data frame per pool
namez<-c("codigo","sitio","lon","lat","biogeo")
names(northsites)<-namez
names(southsites)<-namez
nrin<-merge(NRInorth,northsites,by="codigo")
nris<-merge(NRIsouth,southsites,by="codigo")
#a quick view of the data
tapply(nrin$NRI,nrin$biogeo,mean)
tapply(nris$NRI,nris$biogeo,mean)
mean(nrin$NRI)
mean(nris$NRI)

#save.image("finalscript.Rdata")
##rm(list=ls())