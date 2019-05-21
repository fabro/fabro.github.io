### CLASE 5: ORDENACIÓN EN UN ESPACIO REDUCIDO

# Cargar paquertes, funciones y datos
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(missMDA)
library(FactoMineR)

# Funciones adicionales
source("funciones/cleanplot.pca.R")
source("funciones/PCA.newr.R")
source("funciones/CA.newr.R")

# Cargar los datos
load("Doubs.RData")

# Quitar el sitio vacio
spe <- spe[-8, ]
env <- env[-8, ]
spa <- spa[-8, ]

# Cargar los datos de oribatidos
load("mite.RData")

#### Análisis de componentes principales (PCA) ####
## PCA con todos los datos ambientales

# Recordatorio del contenido de la base de datos env
summary(env)

# PCA basado en una matriz de correlación
# Argumento scale=TRUES estandariza las variables
env.pca <- rda(env, scale = TRUE)
env.pca
summary(env.pca) # Default scaling 2
summary(env.pca, scaling = 1)

# Revisar PDF con explicación

###############============= Diapos ===###

#### Examinar y graficar los resultados parciales de salida del PCA ####
?cca.object   # Explains how an ordination object
              # produced by vegan is structured and how to
              # extract its results.

# Eigenvalues (valores propios)
(ev <- env.pca$CA$eig)

# Gráfica del modelo del palo roto
screeplot(env.pca, bstick = TRUE, npcs = length(env.pca$CA$eig))

## dos biplots pca: scaling 1 and scaling 2

# Plots using biplot.rda
par(mfrow = c(1, 2))
biplot(env.pca, scaling = 1, main = "PCA - scaling 1")
biplot(env.pca, main = "PCA - scaling 2")  # Default scaling 2

# Plots using cleanplot.pca

par(mfrow = c(1, 2))
cleanplot.pca(env.pca, scaling = 1, mar.percent = 0.08)
cleanplot.pca(env.pca, scaling = 2, mar.percent = 0.04)

#########============ Diapos ===###

#### Combinando resultados de clasificación y ordenación ####

# Agrupando los objetos usando datos ambientales: distancia euclideana despues de
# estandarizar las variables, seguido por una agrupación Ward
env.w <- hclust(dist(scale(env)), "ward.D")

# Cortar el dendrograma a 4 grupos
gr <- cutree(env.w, k = 4)
grl <- levels(factor(gr))

# Extraer los pesos de los sitios, scaling 1
sit.sc1 <- scores(env.pca, display = "wa", scaling = 1)

# Graficar los sitios con simbolos de grupos y colores (scaling 1)
par(mfrow = c(1, 1))
p <- plot(
  env.pca,
  display = "wa",
  scaling = 1,
  type = "n",
  main = "PCA correlation + clusters"
)
abline(v = 0, lty = "dotted")
abline(h = 0, lty = "dotted")
for (i in 1:length(grl)) {
  points(sit.sc1[gr == i, ],
         pch = (14 + i),
         cex = 2,
         col = i + 1)
}
text(sit.sc1, row.names(env), cex = 0.7, pos = 3)
# Agregar el dendrograma
ordicluster(p, env.w, col = "dark grey")

# Agregar la leyenda interactivamente
legend(
  locator(1),
  paste("Cluster", c(1:length(grl))),
  pch = 14 + c(1:length(grl)),
  col = 1 + c(1:length(grl)),
  pt.cex = 2
)

#############========== Diapos ===###

#### PCA en los datos de abundancias de peces ####

# Hellinger pre-transformation of the species data
spe.h <- decostand(spe, "hellinger")
(spe.h.pca <- rda(spe.h))

# Scree plot and broken stick model
screeplot(spe.h.pca,
  bstick = TRUE, 
  npcs = length(spe.h.pca$CA$eig)
)

# PCA biplots
spe.pca.sc1 <- scores(spe.h.pca, display = "species", scaling = 1)
spe.pca.sc2 <- scores(spe.h.pca, display = "species", scaling = 2)

par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.06)
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.06)

# Proyección a posteriori de variables ambientales en un PCA
par(mfrow = c(1, 1))

biplot(spe.h.pca, main = "PCA fish abundances - scaling 2")
(spe.h.pca.env <-
    envfit(spe.h.pca, env, scaling = 2)) # Scaling 2 is default
# Graficar variables significativas con un color seleccionado
plot(spe.h.pca.env, p.max = 0.05, col = 3)

# Esto ha agregado las variables ambientales significativas al último biplot creado
# Atención: envfit se debe dar en la misma escala que el gráfico para que se
# se agregue el resultado


# PCA on the environmental data using PCA.newr() and 
#    biplot.PCA.newr()
par(mfrow = c(1, 2))
# PCA; scaling 1 is the default for biplots in this function
env.PCA.PL <- PCA.newr(env, stand = TRUE)
biplot.PCA.newr(env.PCA.PL)

# PCA; scaling 2 in the biplot
biplot.PCA.newr(env.PCA.PL, scaling = 2)

#################=========== Diapos ===###

#### Análisis de correspondencia (CA) ####

# CA de la base de datos crudos
# Calcula CA
(spe.ca <- cca(spe))
summary(spe.ca)		# default scaling 2
summary(spe.ca, scaling = 1)

# Modelo del palo roto  screeplot.cca()
screeplot(spe.ca, bstick = TRUE, npcs = length(spe.ca$CA$eig))

# CA biplots
par(mfrow = c(1, 2))
# Scaling 1: sites are centroids of species
plot(spe.ca, 
     scaling = 1, 
     main = "CA fish abundances - biplot scaling 1"
)
# Scaling 2 (default): species are centroids of sites
plot(spe.ca, main = "CA fish abundances - biplot scaling 2")

###############============ Diapos === ###

#### Ajuste de curvas en un biplot CA ####
par(mfrow = c(1, 1))

plot(spe.ca, main = "CA fish abundances - scaling 2", 
     sub = "Fitted curves: discharge (red), ammonium (green)")
spe.ca.env <- envfit(spe.ca ~ dis + amm, env)
plot(spe.ca.env)  # Two arrows
ordisurf(spe.ca, env$dis, add = TRUE)
ordisurf(spe.ca, env$amm, add = TRUE, col = "green")

###############============ Diapos === ###

#### Cuadro de especies ordenados despues del resultado CA ####
vegemite(spe, spe.ca)
# CA-ordered species table illustrated as a heat map
tabasco(spe, spe.ca)

# CA using CA.newr() function

spe.CA.PL <- CA.newr(spe)
par(mfrow = c(1, 2))
biplot.CA(spe.CA.PL, scaling = 1, cex = 1)
biplot.CA(spe.CA.PL, scaling = 2, cex = 1)

# Ordering of the data table following the first CA axis
# The table is transposed, as in the vegemite() output
summary(spe.CA.PL)
t(spe[order(as.vector(spe.CA.PL$scaling1$sites[, 1])), 
      order(as.vector(spe.CA.PL$scaling1$species[, 1]))])

# Cumulative fit of species
spe.CA.PL$fit$cumulfit.spe
# Cumulative fit of sites
spe.CA.PL$fit$cumulfit.obj

###############============ Diapos === ###

#### Análisis de correspondencia múltiple (MCA) ####

## MCA on the oribatid mite environmental data set

# Preparation of supplementary data: mite classification 
# into 4 groups
# Hellinger transformation of oribatid mite species data
mite.h <- decostand(mite, "hel")
# Ward clustering of mite data
mite.h.ward <- hclust(dist(mite.h), "ward.D2")
# Cut the dendrogram to 4 groups
mite.h.w.g <- cutree(mite.h.ward, 4)
# Assembly of the data set
mite.envplus <- data.frame(mite.env, mite.h.w.g)

# MCA of qualitative environmental data plus supplementary
# variables: 
#   (1) quantitative environmental data, 
#   (2) 4-group mite classification.
#   Default: graph=TRUE.

mite.env.MCA <- MCA(mite.envplus, quanti.sup = 1:2, quali.sup = 6)
mite.env.MCA

# Contingency table crossing variables "Shrub" and "Topo"
table(mite.env$Shrub, mite.env$Topo)

###############============ Diapos === ###

# Principal coordinate analysis (PCoA) ============================

## PCoA on a percentage difference dissimilarity matrix of
## fish species

spe.bray <- vegdist(spe)
spe.b.pcoa <- cmdscale(spe.bray, k = (nrow(spe) - 1), eig = TRUE)
# Plot of the sites
dev.new(
   title = "PCoA on fish species - Percentage difference",
   noRStudioGD = TRUE
)
ordiplot(scores(spe.b.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Add weighted average projection of species
spe.wa <- wascores(spe.b.pcoa$points[, 1:2], spe)
text(spe.wa, rownames(spe.wa), cex = 0.7, col = "red")
# A posteriori projection of environmental variables
(spe.b.pcoa.env <- envfit(spe.b.pcoa, env))
# Plot significant variables with a user-selected colour
plot(spe.b.pcoa.env, p.max = 0.05, col = 3)


# PCoA and projection of species vectors using function pcoa()
spe.h.pcoa <- pcoa(dist(spe.h))
# Biplots
dev.new(title = "PCoA with species vectors",
        width = 12,
        height = 7,
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
# First biplot: Hellinger-transformed species data
biplot.pcoa(spe.h.pcoa, spe.h, dir.axis1 = -1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
text(-0.5, 0.45, "a", cex = 2)
# Second biplot: standardized Hellinger-transformed species data
spe.std <- scale(spe.h)
biplot.pcoa(spe.h.pcoa, spe.std, dir.axis1 = -1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
text(-2.7, 2.45, "b", cex = 2)
# Third biplot: standardized Hellinger-transformed species data;
# only four species in the plot (figure not shown in the book)
dev.new(title = "PCoA with species vectors - 4 species",
        noRStudioGD = TRUE
)
spe.std <- scale(spe.h)
biplot.pcoa(spe.h.pcoa, spe.h[, c(2, 5, 11, 21)], dir.axis1 = -1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)


# Comparison of PCoA results with Euclidean and non-Euclidean
# dissimilarity matrices

# PCoA on a Hellinger distance matrix
is.euclid(dist(spe.h))
summary(spe.h.pcoa)
spe.h.pcoa$values

# PCoA on a percentage difference dissimilarity matrix
is.euclid(spe.bray)
spe.bray.pcoa <- pcoa(spe.bray)
spe.bray.pcoa$values        # Observe eigenvalues 18 and following

# PCoA on the square root of a percentage difference
# dissimilarity matrix
is.euclid(sqrt(spe.bray))
spe.braysq.pcoa <- pcoa(sqrt(spe.bray))
spe.braysq.pcoa$values	# Observe the eigenvalues

# PCoA on a percentage difference dissimilarity matrix with
# Lingoes correction
spe.brayl.pcoa <- pcoa(spe.bray, correction = "lingoes")
spe.brayl.pcoa$values	   # Observe the eigenvalues, col. 1 and 2

# PCoA on a percentage difference dissimilarity matrix with
# Cailliez correction
spe.brayc.pcoa <- pcoa(spe.bray, correction = "cailliez")
spe.brayc.pcoa$values	   # Observe the eigenvalues, col. 1 and 2



# Nonmetric multidimensional scaling (NMDS) =======================

## NMDS applied to the Doubs fish species - percentage difference
## dissimilarity matrix

spe.nmds <- metaMDS(spe, distance = "bray")
spe.nmds
spe.nmds$stress
dev.new(title = "NMDS on fish species - Percentage difference",
   noRStudioGD = TRUE
)
plot(
  spe.nmds,
  type = "t",
  main = paste(
    "NMDS/Percentage difference - Stress =",
    round(spe.nmds$stress, 3)
  )
)

# Shepard plot and goodness of fit
dev.new(title = "NMDS - Shepard plot",
        width = 12,
        height = 6,
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
stressplot(spe.nmds, main = "Shepard plot")
gof <- goodness(spe.nmds)
plot(spe.nmds, type = "t", main = "Goodness of fit")
points(spe.nmds, display = "sites", cex = gof * 300)


# Add colours from a clustering result to an NMDS plot

# Ward clustering of percentage difference dissimilarity matrix
# and extraction of four groups
spe.bray.ward <-
  hclust(spe.bray, "ward.D") # Here better than ward.D2 for 4 groups
spe.bw.groups <- cutree(spe.bray.ward, k = 4)
grp.lev <- levels(factor(spe.bw.groups))

# Combination with NMDS result
sit.sc <- scores(spe.nmds)
dev.new(title = "NMDS plot with cluster colors",
   noRStudioGD = TRUE
)
p <-
  ordiplot(sit.sc, type = "n", 
           main = "NMDS/% difference + clusters Ward/% difference")
for (i in 1:length(grp.lev))
{
  points(sit.sc[spe.bw.groups == i, ],
         pch = (14 + i),
         cex = 2,
         col = i + 1)
}
text(sit.sc, row.names(spe), pos = 4, cex = 0.7)
# Add the dendrogram
ordicluster(p, spe.bray.ward, col = "dark grey")
# Add a legend interactively
legend(
  locator(1),
  paste("Group", c(1:length(grp.lev))),
  pch = 14 + c(1:length(grp.lev)),
  col = 1 + c(1:length(grp.lev)),
  pt.cex = 2
)


# -----------------------------------------------------------------
# The Code It Yourself Corner #2

# A simple function to perform PCA

myPCA <- function(Y) {
  Y.mat <- as.matrix(Y)
  object.names <- rownames(Y)
  var.names <- colnames(Y)
  
  # Centre the data (needed to compute matrix F)
  Y.cent <- scale(Y.mat, center = TRUE, scale = FALSE)
  
  # Covariance matrix S
  Y.cov <- cov(Y.cent)
  
  # Eigenvectors and eigenvalues of S (Legendre and Legendre 2012,
  # eq. 9.1 and 9.2)
  Y.eig <- eigen(Y.cov)
  
  # Copy the eigenvectors to matrix U (used to represent variables
  # in scaling 1 biplots)
  U <- Y.eig$vectors
  rownames(U) <- var.names
  
  # Compute matrix F (used to represent objects in scaling 1 plots)
  F <- Y.cent %*% U			# eq. 9.4
  rownames(F) <- object.names
  
  # Compute matrix U2 (to represent variables in scaling 2 plots)
  # eq. 9.8
  U2 <- U %*% diag(Y.eig$values ^ 0.5)
  rownames(U2) <- var.names
  
  # Compute matrix G (to represent objects in scaling 2 plots)
  # eq. 9.14
  G <- F %*% diag(Y.eig$values ^ 0.5)
  rownames(G) <- object.names
  
  # Output of a list containing all the results
  result <- list(Y.eig$values, U, F, U2, G)
  names(result) <- c("eigenvalues", "U", "F", "U2", "G")
  result
}

# -----------------------------------------------------------------

# PCA on fish species using hand-written function
fish.PCA <- myPCA(spe.h)
summary(fish.PCA)
# Eigenvalues
fish.PCA$eigenvalues
# Eigenvalues expressed as percentages
(pv <- 
  round(100 * fish.PCA$eigenvalues / sum(fish.PCA$eigenvalues), 
  2))
# Alternate computation of total variation (denominator)
round(100 * fish.PCA$eigenvalues / sum(diag(cov(spe.h))), 2)
# Cumulative eigenvalues expressed as percentages
round(
  cumsum(100 * fish.PCA$eigenvalues / sum(fish.PCA$eigenvalues)), 
  2)

# Biplots
dev.new(title = "PCA using homemade function",
        width = 12,
        height = 8)
par(mfrow = c(1, 2))
# Scaling 1 biplot
biplot(fish.PCA$F, fish.PCA$U)
# Scaling 2 biplot
biplot(fish.PCA$G, fish.PCA$U2)

# Plots using generic R plot() function
dev.new(
   title = "PCA using homemade function - generic plot function",
   width = 12,
   height = 8,
   noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
# Scaling 1
# Plot objects
plot(
  fish.PCA$F[, 1],
  fish.PCA$F[, 2],
  asp = 1,
  main = "PCA scaling 1",
  xlab = paste("Axis 1 (", pv[1], "%)", sep = ""),
  ylab = paste("Axis 2 (", pv[2], "%)", sep = "")
)
# Plot variables
arrows(
  x0 = 0,
  y0 = 0,
  fish.PCA$U[, 1],
  fish.PCA$U[, 2],
  length = 0.1,
  col = "red"
)
# Add object numbers
text(
  fish.PCA$F[, 1],
  fish.PCA$F[, 2],
  labels = row.names(spe),
  pos = 3,
  cex = 0.8
)
# Add variable names
text(
  fish.PCA$U[, 1],
  fish.PCA$U[, 2],
  labels = colnames(spe),
  adj = c(-0.2, 0.2),
  col = "red",
  cex = 0.8
)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Scaling 2
plot(
  fish.PCA$G[, 1],
  fish.PCA$G[, 2],
  asp = 1,
  main = "PCA scaling 2",
  xlab = paste("Axis 1 (", pv[1], "%)", sep = ""),
  ylab = paste("Axis 2 (", pv[2], "%)", sep = "")
)
arrows(
  x0 = 0,
  y0 = 0,
  fish.PCA$U2[, 1],
  fish.PCA$U2[, 2],
  length = 0.1,
  col = "red"
)
text(
  fish.PCA$G[, 1],
  fish.PCA$G[, 2],
  labels = row.names(spe),
  pos = 3,
  cex = 0.8
)
text(
  fish.PCA$U2[, 1],
  fish.PCA$U2[, 2],
  labels = colnames(spe),
  col = "red",
  adj = c(-0.2, 0.2),
  cex = 0.8
)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

   ```
