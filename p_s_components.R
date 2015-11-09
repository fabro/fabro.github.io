library(ape)
library(picante)
library(vegan)


p.s.comps <- function(phy,sppnames,data,colnum){
  comps.list <- vector("list",3)
  #prune the tree, retaining only spp of interest
  phy.spp <- drop.tip(phy,tip = which(!phy$tip.label%in%sppnames)) 
  #create the phylogenetic distance and variance-covariance matrices
  phydis <-cophenetic(phy.spp)
  phycor <-vcv(phy.spp,corr=TRUE)
  diag(phycor)<-0
  #match the data and phylogeny
  phy.data <- match.phylo.data(phy.spp, data[,colnum])
  y1 <- phy.data$data
  #PVR (PCOA - from ape package. Matrix has to be a squared one)
  phy.pcord <- pcoa(phydis,correction="none", rn=NULL)
  phy.vecold <- phy.pcord$vectors
  #sequential eigenvector selection minimizing Moran's I
  phy.my1 <-Moran.I(y1,phycor)
  n <- length(phy.spp$tip.label)
  
  for(i in 1:n-1){
    m1 <-lm(y1~phy.vecold[,1:i])
    res <-m1$residuals		
    mor <-Moran.I(res,phycor)
    PMPVR <-mor$p.value
    if (PMPVR>0.05){phy.sel <-(i+2)}
    if (PMPVR>0.05){break}
  }
  phy.vecsel <-phy.vecold[,1:phy.sel]
  #linear modeling
  pvr1 <- lm(y1~phy.vecsel)
  r2pvr1 <- RsquareAdj(pvr1)$r.squared
  ppvr1 <- as.numeric(format.pval(pf(summary(pvr1)$fstatistic[1L], summary(pvr1)$fstatistic[2L], summary(pvr1)$fstatistic[3L], lower.tail =FALSE), digits = 4))
  fitpvr1 <- pvr1$fitted.values
  respvr1 <- pvr1$residuals
  phy.m <- Moran.I(respvr1,phycor)
  PMPVR <- phy.m$p.value
  #obtain the phylogenetic (fitted) and specific (residuals) components
  fitted <- fitpvr1[order(names(fitpvr1))]
  resids <- respvr1[order(names(respvr1))]
  #join the P and S components in one matrix
  comps <- cbind(fitted,resids)
  colnames(comps) <- c("Phy_component", "Spec_component")
  #save the results, including the main parameters of the linear model and autocorrelation test (after adding the phylogenetic vectors)
  comps.list[[1]] <- c(r2pvr1,ppvr1,phy.m$observed,PMPVR)
  names(comps.list[[1]]) <- c("R_squared","Rsq_pvalue","Moran_I","MoranI_pvalue")
  comps.list[[2]] <- comps
  comps.list[[3]] <- paste(ncol(phy.vecsel),'selected_vectors',sep=" ")
  return(comps.list)
}