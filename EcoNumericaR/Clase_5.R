### CLASE 5: ORDENACIÓN CANÓNICA

# Load packages, functions and data ===============================
library(ade4)
library(adegraphics)
library(adespatial)
library(cocorresp)
library(vegan)
library(vegan3d)
library(ape)   # For PCoA with Lingoes correction (not in the book)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)

# Funciones adicionales
source("funciones/hcoplot.R")
source("funciones/triplot.rda.R")
source("funciones/plot.lda.R")
source("funciones/polyvars.R")
source("funciones/screestick.R")

# Datos
load("Doubs.RData")  

# Quitar el sitio vacio 8
spe <- spe[-8, ]
env <- env[-8, ]
spa <- spa[-8, ]

# Tomar la variable "dfs" (distancia al origen) para un uso posterior
dfs <- env[, 1]

# Quitar la variable "dfs" de la base de datos
env2 <- env[, -1]

# Recodificar la variable pendiente (slo) en un factor (cualitativo) para mostrar como son
# manejadas en las ordenaciones
slo2 <- rep(".very_steep", nrow(env))
slo2[env$slo <= quantile(env$slo)[4]] <- ".steep"
slo2[env$slo <= quantile(env$slo)[3]] <- ".moderate"
slo2[env$slo <= quantile(env$slo)[2]] <- ".low"
slo2 <- factor(slo2, 
  levels = c(".low", ".moderate", ".steep", ".very_steep"))
table(slo2)

# Crear una base de datos env3 con pendiente como una variable cualitativa
env3 <- env2
env3$slo <- slo2

# Crear dos subconjuntos de variables explicativas
# Fisiografia ()
# Physiography (gradiente corriente arriba-corriente abajo)
envtopo <- env2[, c(1 : 3)]
names(envtopo)

# Calidad del agua
envchem <- env2[, c(4 : 10)]
names(envchem)

# Transformación Hellinger de las especies
spe.hel <- decostand(spe, "hellinger")

########################### ============ Diapos ===###

#### Análisis de redundancia ####
# RDA de los datos de peces transformados a Hellinger, restringidos por todas las variables
# ambientales contenidas en env3
(spe.rda <- rda(spe.hel ~ ., env3)) # Observe the shortcut formula
summary(spe.rda)	# Scaling 2 (default)
# Leer pdf RDA

# Coeficientes canónicos del objeto rda
coef(spe.rda)

########################### ============ Diapos ===###
# R^2 sin ajustar tomado del objeto rda
(R2 <- RsquareAdj(spe.rda)$r.squared)

# R^2 ajustado tomado del objeto rda
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)

# Triplots de los resultados rda (puntuaciones lc)
# Puntuaciones de los sitios como combinaciones lineares de las variables ambientales
par(mfrow = c(2, 1))
plot(spe.rda,
  scaling = 1,
  display = c("sp", "lc", "cn"),
  main = "Triplot RDA spe.hel ~ env3 - scaling 1 - lc scores"
)
spe.sc1 <- 
  scores(spe.rda, 
         choices = 1:2, 
         scaling = 1, 
         display = "sp"
)
arrows(0, 0, 
  spe.sc1[, 1] * 0.92,
  spe.sc1[, 2] * 0.92,
  length = 0, 
  lty = 1, 
  col = "red"
)
text(-0.75, 0.7, "a", cex = 1.5)

# Scaling 2
plot(spe.rda, 
  display = c("sp", "lc", "cn"), 
  main = "Triplot RDA spe.hel ~ env3 - scaling 2 - lc scores"
)
spe.sc2 <- 
  scores(spe.rda, 
         choices = 1:2, 
         display = "sp"
)
arrows(0, 0, 
  spe.sc2[, 1] * 0.92, 
  spe.sc2[, 2] * 0.92,
  length = 0,
  lty = 1,
  col = "red"
)
text(-0.82, 0.55, "b", cex = 1.5)

# Triplots de los restultados rda (puntuaciones wa)
# Puntuaciones de los sitios como promedios ponderados (default vegan)
# Scaling 1 :  distance triplot
plot(spe.rda, 
  scaling = 1, 
  main = "Triplot RDA spe.hel ~ env3 - scaling 1 - wa scores"
)
arrows(0, 0, 
  spe.sc1[, 1] * 0.92, 
  spe.sc1[, 2] * 0.92, 
  length = 0, 
  lty = 1, 
  col = "red"
)

# Scaling 2 (default) :  correlation triplot
plot(spe.rda, 
  main = "Triplot RDA spe.hel ~ env3 - scaling 2 - wa scores")
arrows(0, 0, 
  spe.sc2[, 1] * 0.92, 
  spe.sc2[, 2] * 0.92, 
  length = 0, 
  lty = 1, 
  col = "red"
)

# Seleccionar especies con bondad de ajuste de al menos 0.6 en el plano de ordenación
# formado por los ejes 1 y 2
spe.good <- goodness(spe.rda)
sel.sp <- which(spe.good[, 2] >= 0.6)

# Función de triplots
par(mfrow = c(2, 1))
triplot.rda(spe.rda, 
  site.sc = "lc", 
  scaling = 1, 
  cex.char2 = 0.7, 
  pos.env = 3, 
  pos.centr = 1, 
  mult.arrow = 1.1, 
  mar.percent = 0.05, 
  select.spe = sel.sp
)
text(-0.92, 0.72, "a", cex = 2)
triplot.rda(spe.rda, 
  site.sc = "lc", 
  scaling = 2, 
  cex.char2 = 0.7, 
  pos.env = 3, 
  pos.centr = 1, 
  mult.arrow = 1.1, 
  mar.percent = 0.05, 
  select.spe = sel.sp
)
text(-2.82, 2, "b", cex = 2)

################### ================ Diapos === ###

#### Prueba global de los resultados RDA ####
anova(spe.rda, permutations = how(nperm = 999))

# Prueba de todos los ejes canónicos
anova(spe.rda, by = "axis", permutations = how(nperm = 999))

################### ================ Diapos === ###
#### RDA parcial: efecto de la química del agua, manteniendo fisiografia constante ####

# Sintaxis simple, X y W pueden estar en bases separadas de variables cuantitativas
(spechem.physio <- rda(spe.hel, envchem, envtopo))
summary(spechem.physio)

# Interface formula, X y W deben estar en la misma base de datos
(spechem.physio2 <- 
  rda(spe.hel ~ pH + har + pho + nit + amm + oxy + bod 
	    + Condition(ele + slo + dis), data = env2))


# Test of the partial RDA, using the results with the formula 
# interface to allow the tests of the axes to be run
anova(spechem.physio2, permutations = how(nperm = 999))
anova(spechem.physio2, permutations = how(nperm = 999), by = "axis")

# Partial RDA triplots (with fitted site scores) 
# with function triplot.rda
# Scaling 1
par(mfrow = c(2, 1))
triplot.rda(spechem.physio, 
  site.sc = "lc", 
  scaling = 1, 
  cex.char2 = 0.8, 
  pos.env = 3, 
  mar.percent = 0
)
text(-0.58, 0.64, "a", cex = 2)
# Scaling 2
triplot.rda(spechem.physio, 
  site.sc = "lc", 
  scaling = 2, 
  cex.char2 = 0.8, 
  pos.env = 3, 
  mult.spe = 1.1, 
  mar.percent = 0.04
)
text(-3.34, 3.64, "b", cex = 2)

# Alternative plot using plot.cca ---------------------------------
# Partial RDA triplots (with fitted site scores) 
# with function triplot.rda
# Scaling 1
plot(spechem.physio, 
     scaling = 1, 
     display = c("sp", "lc", "cn"), 
     main = "Triplot RDA spe.hel ~ chem | Topo - scaling 1 - lc scores")
spe3.sc <- 
  scores(spechem.physio, 
         choices = 1:2, 
         scaling = 1, 
         display = "sp"
)
arrows(0, 0, 
  spe3.sc[, 1] * 0.92, 
  spe3.sc[, 2] * 0.92, 
  length = 0, 
  lty = 1, 
  col = "red"
)

# Scaling 2
plot(spechem.physio, 
     display = c("sp", "lc", "cn"), 
     main = "Triplot RDA spe.hel ~ chem | Topo - scaling 2 - lc scores")
spe4.sc <- 
  scores(spechem.physio, 
  choices = 1:2, 
  display = "sp"
)
arrows(0, 0, 
  spe4.sc[, 1] * 0.88, 
  spe4.sc[, 2] * 0.88, 
  length = 0, 
  lty = 1, 
  col = "red"
)
# End Alternative plot using plot.cca -----------------------------
################### ================ Diapos === ###
#### Factores de inflación varianza (VIF) en dos RDAs ####
vif.cca(spe.rda)
# Partial RDA â€“ physiographic variables only
vif.cca(spechem.physio) 
################### ================ Diapos === ###

#### Selección hacia adelante de variables explicativas ####

# RDA with all explanatory variables except dfs
spe.rda.all <- rda(spe.hel ~ ., data = env2)
# Global adjusted R^2
(R2a.all <- RsquareAdj(spe.rda.all)$adj.r.squared)

# Forward selection using forward.sel()
forward.sel(spe.hel, env2, adjR2thresh = R2a.all)

# Forward selection using vegan's ordistep()
# This function allows the use of factors. 
mod0 <- rda(spe.hel ~ 1, data = env2)
step.forward <- 
  ordistep(mod0, 
           scope = formula(spe.rda.all), 
           direction = "forward", 
           permutations = how(nperm = 499)
)
RsquareAdj(step.forward)

############# ========== Diapos === ###

# Backward elimination using vegan's ordistep()
step.backward <-
  ordistep(spe.rda.all,
           permutations = how(nperm = 499))
# With redundant argument direction = "backward":
# step.backward <-
#  ordistep(spe.rda.all,
#           direction = "backward",
#           permutations = how(nperm = 499)
RsquareAdj(step.backward)

# Forward selection using vegan's ordiR2step()
# using a double stopping criterion (Blanchet et al. 2008a)
# and object env containing only quantitative variables.
step2.forward <- 
  ordiR2step(mod0, 
             scope = formula(spe.rda.all), 
             direction = "forward", 
             R2scope = TRUE,
             permutations = how(nperm = 199)
  )
RsquareAdj(step2.forward)

# Forward selection using vegan's ordiR2step()
# using a double stopping criterion (Blanchet et al. 2008a)
# and object env3 containing a factor.
mod00 <- rda(spe.hel ~ 1, data = env3)
spe.rda2.all <- rda(spe.hel ~ ., data = env3)
step3.forward <- 
  ordiR2step(mod00, 
             scope = formula(spe.rda2.all), 
             direction = "forward", 
             permutations = how(nperm = 199)
  )
RsquareAdj(step3.forward)
# Note that the adjusted R^2 of the complete model is smaller
# than that of the complete RDA with only quantitative
# variables.
# Some information has been lost when transforming the
# quantitative slo variable into a factor with 4 levels.

# Partial forward selection with variable slo held constant
mod0p <- rda(spe.hel ~ Condition(slo), data = env2)
mod1p <- rda(spe.hel ~ . + Condition(slo), data = env2)
step.p.forward <- 
  ordiR2step(mod0p, 
             scope = formula(mod1p), 
             direction = "forward", 
             permutations = how(nperm = 199)
)


## Parsimonious RDA

(spe.rda.pars <- rda(spe.hel ~ ele + oxy + bod, data = env2))
anova(spe.rda.pars, permutations = how(nperm = 999))
anova(spe.rda.pars, permutations = how(nperm = 999), by = "axis")
(R2a.pars <- RsquareAdj(spe.rda.pars)$adj.r.squared)
# Compare the variance inflation factors
vif.cca(spe.rda.all)
vif.cca(spe.rda.pars)

# Triplots of the parsimonious RDA (with fitted site scores)
par(mfrow = c(2, 1))
# Scaling 1
triplot.rda(spe.rda.pars, 
  site.sc = "lc", 
  scaling = 1, 
  cex.char2 = 0.8, 
  pos.env = 2, 
  mult.spe = 0.9, 
  mult.arrow = 0.92, 
  mar.percent = 0.01
)
# Scaling 2
triplot.rda(spe.rda.pars, 
  site.sc = "lc", 
  scaling = 2, 
  cex.char2 = 0.8, 
  pos.env = 2, 
  mult.spe = 1.1, 
  mar.percent = -0.02
)


# Alternate code using plot.cca -----------------------------------
# Triplots of the parsimonious RDA (with fitted site scores)
par(mfrow = c(1, 2))
# Scaling 1
plot(spe.rda.pars, 
     scaling = 1, 
     display = c("sp", "lc", "cn"), 
     main = "Triplot RDA spe.hel ~ ele+oxy+bod - scaling 1 - lc scores")
spe4.sc <- 
  scores(spe.rda.pars, 
         choices = 1:2, 
         scaling = 1, 
         display = "sp"
)
arrows(0, 0, 
  spe4.sc[, 1] * 0.92, 
  spe4.sc[, 2] * 0.92, 
  length = 0, 
  lty = 1, 
  col = "red"
)
# Scaling 2
plot(spe.rda.pars, 
     display = c("sp", "lc", "cn"), 
     main = "Triplot RDA spe.hel ~ ele+oxy+bod - scaling 2 - lc scores")
spe5.sc <- 
  scores(spe.rda.pars, 
         choices = 1:2, 
         display = "sp"
)
arrows(0, 0, 
  spe5.sc[, 1] * 0.9, 
  spe5.sc[, 2] * 0.9, 
  length = 0, 
  lty = 1, 
  col = "red"
)
# End Alternate code using plot.cca --------------------------------


## Environmental reconstruction (calibration) with RDA

# New (fictitious) objects with fish abundances
# Variables(species) must match those in the original data set in 
# name, number and order
# New site 1 is made from rounded means of species in sites 1 to 15
site1.new <- round(apply(spe[1:15, ], 2, mean))
# New site 2 is made from rounded means of species in sites 16 - 29
site2.new <- round(apply(spe[16:29, ], 2, mean))
obj.new <- t(cbind(site1.new, site2.new))

# Hellinger transformation of the new sites
obj.new.hel <- decostand(obj.new, "hel")

# Calibration
calibrate(spe.rda.pars, obj.new.hel)

# Compare with real values at sites 7 to 9 and 22 to 24: 
env2[7:9, c(1, 9, 10)]
env2[22:24, c(1, 9, 10)]


## Variation partitioning with two sets of explanatory variables

# Explanation of fraction labels (two, three and four explanatory 
# matrices) with optional colours
dev.new(
  title = "Symbols of variation partitioning fractions", 
  width = 6, 
  height = 2.3, noRStudioGD = TRUE
)
par(mfrow = c(1, 3), mar = c(1, 1, 1, 1))
showvarparts(2, bg = c("red", "blue"))
showvarparts(3, bg = c("red", "blue", "yellow"))
showvarparts(4, bg = c("red", "blue", "yellow", "green"))

## 1. Variation partitioning with all explanatory variables
##    (except dfs)
(spe.part.all <- varpart(spe.hel, envchem, envtopo))

# Plot of the partitioning results
dev.new(title = "Variation partitioning - all variables", 
        noRStudioGD = TRUE)
plot(spe.part.all, digits = 2, bg = c("red", "blue"))


## 2. Variation partitioning after forward selection of explanatory 
##    variables
# Separate forward selection in each subset of environmental 
# variables
spe.chem <- rda(spe.hel, envchem)
R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
forward.sel(spe.hel, 
  envchem, 
  adjR2thresh = R2a.all.chem, 
  nperm = 9999
)

spe.topo <- rda(spe.hel, envtopo)
R2a.all.topo <- RsquareAdj(spe.topo)$adj.r.squared
forward.sel(spe.hel, 
  envtopo, 
  adjR2thresh = R2a.all.topo, 
  nperm = 9999
)

# Parsimonious subsets of explanatory variables, based on forward 
# selections
names(envchem)
envchem.pars <- envchem[, c(4, 6, 7)]
names(envtopo)
envtopo.pars <- envtopo[, c(1, 2)]

# Variation partitioning
(spe.part <- varpart(spe.hel, envchem.pars, envtopo.pars))
dev.new(
  title = "Variation partitioning - parsimonious subsets", 
  noRStudioGD = TRUE
)
plot(spe.part, 
  digits = 2, 
  bg = c("red", "blue"), 
  Xnames = c("Chemistry", "Physiography"), 
  id.size = 0.7
)

# Tests of all testable fractions
# Test of fraction [a+b]
anova(rda(spe.hel, envchem.pars), permutations = how(nperm = 999))
# Test of fraction [b+c]
anova(rda(spe.hel, envtopo.pars), permutations = how(nperm = 999))
# Test of fraction [a+b+c]
env.pars <- cbind(envchem.pars, envtopo.pars)
anova(rda(spe.hel, env.pars), permutations = how(nperm = 999))
# Test of fraction [a]
anova(rda(spe.hel, envchem.pars, envtopo.pars), 
      permutations = how(nperm = 999)
)
# Test of fraction [c]
anova(rda(spe.hel, envtopo.pars, envchem.pars), 
      permutations = how(nperm = 999)
)


## 3. Variation partitioning without the 'nit' variable
envchem.pars2 <- envchem[, c(6, 7)]
(spe.part2 <- varpart(spe.hel, envchem.pars2, envtopo.pars))
dev.new(
  title = "Variation partitioning - parsimonious subset 2", 
  noRStudioGD = TRUE
)
plot(spe.part2, digits = 2)


## Two-way MANOVA by RDA

# Creation of a factor 'elevation' (3 levels, 9 sites each)
ele.fac <- gl(3, 9, labels = c("high", "mid", "low"))
# Creation of a factor mimicking 'pH'
pH.fac <- 
  as.factor(c(1, 2, 3, 2, 3, 1, 3, 2, 1, 2, 1, 3, 3, 2, 
              1, 1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 3))
# Is the two-way factorial design balanced?
table(ele.fac, pH.fac)

# Creation of Helmert contrasts for the factors and the interaction
ele.pH.helm <- 
  model.matrix(~ ele.fac * pH.fac, 
               contrasts = list(ele.fac = "contr.helmert", 
                                pH.fac = "contr.helmert"))[, -1]
ele.pH.helm
ele.pH.helm2 <- 
  model.matrix(~ ele.fac + pH.fac, 
	   contrasts = list(ele.fac = "contr.helmert", 
	                    pH.fac = "contr.helmert"))[, -1]
colnames(ele.pH.helm2)

# Check property 1 of Helmert contrasts : all variables sum to 0
apply(ele.pH.helm, 2, sum)
# Check property 2 of Helmert contrasts: their crossproducts 
# must be 0 within and between groups (factors and interaction)
crossprod(ele.pH.helm)

# Verify multivariate homogeneity of within-group covariance
# matrices using the betadisper() function (vegan package)
# implementing Marti Anderson's testing method (Anderson 2006)

# To avoid the rist of heterogeneity of variances with respect to 
# one factor because of the dispersion in the other (in case of 
# interaction), creation of a factor crossing the two factors, i.e. 
# defining the cell-by-cell attribution of the data
cell.fac <- gl(9, 3) 

spe.hel.d1 <- dist(spe.hel[1:27, ])

# Test of homogeneity of within-cell dispersions
(spe.hel.cell.MHV <- betadisper(spe.hel.d1, cell.fac))
anova(spe.hel.cell.MHV)     # Parametric test (not recommended here)
permutest(spe.hel.cell.MHV)

# Alternatively, test homogeneity of dispersions within each
# factor. 
# These tests ore more robust with this small example because 
# there are now 9 observations per group instead of 3. 
# Factor "elevation"
(spe.hel.ele.MHV <- betadisper(spe.hel.d1, ele.fac))
anova(spe.hel.ele.MHV)     # Parametric test (not recommended here)
permutest(spe.hel.ele.MHV) # Permutation test
# Factor "pH"
(spe.hel.pH.MHV <- betadisper(spe.hel.d1, pH.fac))
anova(spe.hel.pH.MHV)
permutest(spe.hel.pH.MHV) # Permutation test


## Step-by-step procedure using function rda()

# Test the interaction first. Factors ele and pH (columns 1-4)  
# are assembled to form the matrix of covariables for the test.
interaction.rda <- 
  rda(spe.hel[1:27, ], 
      ele.pH.helm[, 5:8], 
      ele.pH.helm[, 1:4])
anova(interaction.rda, permutations = how(nperm = 999))

# Test the main factor ele. The factor pH and the interaction
# are assembled to form the matrix of covariables. 
factor.ele.rda <- 
  rda(spe.hel[1:27, ], 
      ele.pH.helm[, 1:2], 
      ele.pH.helm[, 3:8])
anova(factor.ele.rda, 
      permutations = how(nperm = 999), 
      strata = pH.fac
)

# Test the main factor pH. The factor ele and the interaction
# are assembled to form the matrix of covariables. 
factor.pH.rda <- 
  rda(spe.hel[1:27, ], 
      ele.pH.helm[, 3:4], 
      ele.pH.helm[, c(1:2, 5:8)]) 
anova(factor.pH.rda, 
  permutations = how(nperm = 999), 
  strata = ele.fac
)

# RDA with the significant factor ele
ele.rda.out <- rda(spe.hel[1:27, ]~ ., as.data.frame(ele.fac))
# Triplot with "wa" sites related to factor centroids, and species 
# arrows
dev.new(
  title = "Multivariate ANOVA - elevation", 
  noRStudioGD = TRUE
)
plot(ele.rda.out, 
  scaling = 1, 
  display = "wa", 
  main = "Multivariate ANOVA, factor elevation - scaling 1 - 
          wa scores")
ordispider(ele.rda.out, ele.fac, 
  scaling = 1, 
  label = TRUE, 
  col = "blue"
)
spe.sc1 <- 
  scores(ele.rda.out, 
  scaling = 1, 
  display = "species")
arrows(0, 0, 
  spe.sc1[, 1] * 0.3, 
  spe.sc1[, 2] * 0.3, 
  length = 0.1, 
  angle = 10, 
  col = "red"
)
text(
  spe.sc1[, 1] * 0.3, 
  spe.sc1[, 2] * 0.3, 
  labels = rownames(spe.sc1), 
  pos = 4, 
  cex = 0.8, 
  col = "red"
)

## Permutational MANOVA using adonis2()
adonis2(spe.hel[1:27, ] ~ ele.fac * pH.fac, 
  method = "euc", 
  by = "term"
)


## RDA with a single second degree explanatory variable

# Create a matrix of dfs and its orthogonal second degree term 
# using function poly()
dfs.df <- poly(dfs, 2)
colnames(dfs.df) <- c("dfs", "dfs2")
# Verify that the polynomial terms are orthogonal
cor(dfs.df)
# Find out if both variables are significant
forward.sel(spe.hel, dfs.df)

# RDA and test
spe.dfs.rda <- rda(spe.hel ~ ., as.data.frame(dfs.df))
anova(spe.dfs.rda)

# Triplot using "lc" (model) site scores and scaling 2
dev.new(
  title = "RDA w. 2nd-degree variable - scaling 2", 
  noRStudioGD = TRUE
)
triplot.rda(spe.dfs.rda, 
  site.sc = "lc", 
  scaling = 2, 
  plot.sites = FALSE, 
  pos.env = 1, 
  mult.arrow = 0.9, 
  move.origin = c(-0.25, 0), 
  mar.percent = 0
)

# Alternate code using plot.cca -----------------------------
# Triplot using lc (model) site scores and scaling 2
dev.new(
  title = "RDA w. one 2nd-degree variable - scaling 2", 
  noRStudioGD = TRUE
)
plot(spe.dfs.rda, 
  scaling = 2, 
  display = c("sp", "lc", "cn"), 
  main = "Triplot RDA spe ~ dfs+dfs2 - scaling 2 - lc scores")
spe6.sc <- 
  scores(spe.dfs.rda, 
         choices = 1:2, 
         scaling = 2, 
         display = "sp")
arrows(0, 0, 
  spe6.sc[, 1] * 0.9, 
  spe6.sc[, 2] * 0.9, 
  length = 0, 
  lty = 1, 
  col = "red"
)
# End Alternate code using plot.cca --------------------------

# Maps of four fish species
dev.new(title = "Four fish species", noRStudioGD = TRUE)
par(mfrow = c(2, 2))
plot(spa, 
  asp = 1, 
  col = "brown", 
  cex = spe$Satr, 
  xlab = "x (km)", 
  ylab = "y (km)", 
  main = "Brown trout"
)
lines(spa, col = "light blue")
plot(spa, 
  asp = 1, 
  col = "brown", 
  cex = spe$Thth, 
  xlab = "x (km)", 
  ylab = "y (km)", 
  main = "Grayling"
)
lines(spa, col = "light blue")
plot(spa, 
  asp = 1, 
  col = "brown", 
  cex = spe$Alal, 
  xlab = "x (km)", 
  ylab = "y (km)", 
  main = "Bleak"
)
lines(spa, col = "light blue")
plot(spa, 
  asp = 1, 
  col = "brown", 
  cex = spe$Titi, 
  xlab = "x (km)", 
  ylab = "y (km)", 
  main = "Tench"
)
lines(spa, col = "light blue")


## Polynomial RDA, second degree, with forward selection
## among all environmental variables
 
env.square <- polyvars(env2, degr = 2)
names(env.square)
spe.envsq.rda <- rda(spe.hel ~ ., env.square)
R2ad <- RsquareAdj(spe.envsq.rda)$adj.r.squared
spe.envsq.fwd <- 
   forward.sel(spe.hel, 
               env.square, 
               adjR2thresh = R2ad)
spe.envsq.fwd
envsquare.red <- env.square[, sort(spe.envsq.fwd$order)]
(spe.envsq.fwd.rda <- rda(spe.hel ~., envsquare.red))
RsquareAdj(spe.envsq.fwd.rda)
summary(spe.envsq.fwd.rda)

# Triplot using lc (model) site scores and scaling 2
dev.new(
  title = "RDA w. 2nd-order variables - scaling 2", 
  noRStudioGD = TRUE
)
triplot.rda(spe.envsq.fwd.rda, 
  site.sc = "lc", 
  scaling = 2, 
  plot.sites = FALSE, 
  pos.env = 1, 
  mult.arrow = 0.9, 
  mult.spe = 0.9,
  mar.percent = 0
)


## Distance-based redundancy analysis (db-RDA) :  square-rooted
## percentage difference response matrix and 'ele' factor with 
## factor pH and interaction as covariables

# Rename columns of matrix of Helmert contrasts (for convenience)
colnames(ele.pH.helm) <- 
  c("ele1", "ele2", "pH1", "pH2", "ele1pH1", "ele1pH2", 
    "ele2pH1", "ele2pH2" )

# Create the matrix of covariables. MUST be of class matrix, 
# NOT data.frame
covariables <- ele.pH.helm[, 3:8]

# Compute the dissimilarity response matrix with veganâ€™s vegdist()
spe.bray27 <- vegdist(spe[1:27, ], "bray")
# â€¦ or with function dist.ldc() of adespatial
spe.bray27 <- dist.ldc(spe[1:27, ], "percentdiff")

# 1. dbrda() on the square-rooted dissimilarity matrix
bray.ele.dbrda <-  dbrda(
    sqrt(spe.bray27) ~ ele.pH.helm[, 1:2] + Condition(covariables))
anova(bray.ele.dbrda, permutations = how(nperm = 999)) 

# 2. capscale() with raw (site by species) data
# Rename factor (cosmetic for plot)
ele.fac. <- ele.fac
bray.env.cap <- 
  capscale(spe[1:27, ] ~ ele.fac. + Condition(covariables), 
           data = as.data.frame(ele.pH.helm), 
           distance = "bray", 
           add = "lingoes", 
           comm = spe[1:27, ])
anova(bray.env.cap, permutations = how(nperm = 999))
# Plot with "wa" scores to see dispersion of sites around the 
# factor levels
dev.new(
  title = "db-RDA - wa - scaling 1", 
  noRStudioGD = TRUE
)
triplot.rda(bray.env.cap, site.sc = "wa", scaling = 1)

# The results of the two analyses are slightly different because
# (1) the test is not performed in the same manner and (2) the 
# correction to make the response matrix Euclidean is not the same.


### NOT IN THE BOOK ###

## Alternative ways of computing db-RDA

# 1. capscale() with raw (site by species) data
# Alternate coding with explicit covariables coming from same
# object as the constraining variables : 
bray.env.capscale <- 
  capscale(spe[1:27, ] ~ ele1 + ele2 + 
       Condition(pH1 + pH2 + ele1pH1 + ele1pH2 + ele2pH1 + ele2pH2), 
       data = as.data.frame(ele.pH.helm), 
       distance = "bray", 
       add = "cailliez", 
       comm = spe[1:27, ])
anova(bray.env.capscale, permutations = how(nperm = 999))

# 2. PCoA with Lingoes (1971) correction
#    Explicit steps
spe.bray27.lin <- pcoa(spe.bray27, correction = "lingoes") 
spe.bray27.lingoes <- spe.bray27.lin$vectors.cor 
# Test of the factor ele. Factor pH and interaction, Helmert-coded,
# form the matrix of covariables
spe.L.ele.dbrda <- 
  rda(spe.bray27.lingoes, 
      ele.pH.helm[, 1:2], 
      covariables) 
anova(spe.L.ele.dbrda, permutations = how(nperm = 999))

# Same by staying in {vegan} and using wcmdscale() : 

spe.lingoes2 <- wcmdscale(spe.bray27, add = "lingoes") 
anova(
  rda(spe.lingoes2 ~ ele.pH.helm[, 1:2] + Condition(covariables))
)

### END NOT IN THE BOOK ###



# -----------------------------------------------------------------
# The Code It Yourself corner #3

myRDA <- function(Y, X)
{

	## 1. Preparation of the data

	Y.mat <- as.matrix(Y)
	Yc <- scale(Y.mat, scale = FALSE)

	X.mat <- as.matrix(X)
	Xcr <- scale(X.mat)

	# Dimensions
	n <- nrow(Y)
	p <- ncol(Y)
	m <- ncol(X)

	## 2. Computation of the multivariate linear regression

	# Matrix of regression coefficients (eq. 11.9)
	B <- solve(t(Xcr) %*% Xcr) %*% t(Xcr) %*% Yc

	# Matrix of fitted values (eq. 11.10)
	Yhat <- Xcr %*% B

	# Matrix of residuals
	Yres <- Yc - Yhat


	## 3. PCA on fitted values

	# Covariance matrix (eq. 11.12)
	S <- cov(Yhat)

	# Eigenvalue decomposition
	eigenS <- eigen(S)

	# How many canonical axes?
	kc <- length(which(eigenS$values > 0.00000001))

	# Eigenvalues of canonical axes
	ev <- eigenS$values[1 : kc]
	# Total variance (inertia) of the centred matrix Yc
	trace = sum(diag(cov(Yc)))
	
	# Orthonormal eigenvectors (contributions of response 
	# variables, scaling 1)
	U <- eigenS$vectors[, 1 : kc]
	row.names(U) <- colnames(Y)

	# Site scores (vegan's wa scores, scaling 1; eq.11.17)
	F <- Yc %*% U
	row.names(F) <- row.names(Y)

	# Site constraints (vegan's 'lc' scores, scaling 1; 
	# eq. 11.18)
	Z <- Yhat %*% U
	row.names(Z) <- row.names(Y)

	# Canonical coefficients (eq. 11.19)
	CC <- B %*% U
	row.names(CC) <- colnames(X)

	# Explanatory variables
	# Species-environment correlations
	corXZ <- cor(X, Z)

	# Diagonal matrix of weights
	D <- diag(sqrt(ev / trace))

	# Biplot scores of explanatory variables
	coordX <- corXZ %*% D    # Scaling 1
	coordX2 <- corXZ         # Scaling 2
	row.names(coordX) <- colnames(X)
	row.names(coordX2) <- colnames(X)

	# Scaling to sqrt of the relative eigenvalue
	# (for scaling 2)
	U2 <- U %*% diag(sqrt(ev))
	row.names(U2) <- colnames(Y)
	F2 <- F %*% diag(1/sqrt(ev))
	row.names(F2) <- row.names(Y)
	Z2 <- Z %*% diag(1/sqrt(ev))
	row.names(Z2) <- row.names(Y)

	# Unadjusted R2
	R2 <- sum(ev/trace)
	# Adjusted R2
	R2a <- 1 - ((n - 1)/(n - m - 1)) * (1 - R2)


	# 4. PCA on residuals
	# Write your own code as in Chapter 5. It could begin 
	# with : 
	#     eigenSres <- eigen(cov(Yres))
	#     evr <- eigenSres$values
	

	# 5. Output

	result <- 
	  list(trace, R2, R2a, ev, CC, U, F, Z, coordX, 
	       U2, F2, Z2, coordX2)
	names(result) <- 
	  c("Total_variance", "R2", "R2adj", "Can_ev", 
	    "Can_coeff", "Species_sc1", "wa_sc1", "lc_sc1", 
	    "Biplot_sc1", "Species_sc2", "wa_sc2", "lc_sc2", 
	    "Biplot_sc2") 

	result
}


doubs.myRDA <- myRDA(spe.hel, env2)
summary(doubs.myRDA)
# Retrieve adjusted R-square
doubs.myRDA$R2adj

###################################################################



#### Análisis de correspondencia canónica (CCA) ####
# CCA de los datos de especies crudos, restringidos por todas las variables ambientales
(spe.cca <- cca(spe ~ ., env3))
summary(spe.cca)	# Scaling 2 (default)

# Unadjusted and adjusted R^2 - like statistics
RsquareAdj(spe.cca)

## CCA triplots (using lc site scores)
par(mfrow = c(2, 1))
# Scaling 1: species scores scaled to the relative eigenvalues, 
# sites are weighted averages of the species
plot(spe.cca, 
  scaling = 1, 
  display = c("sp", "lc", "cn"), 
  main = "Triplot CCA spe ~ env3 - scaling 1"
)
text(-2.3, 4.1, "a", cex = 1.5)
# Default scaling 2: site scores scaled to the relative 
# eigenvalues, species are weighted averages of the sites
plot(spe.cca, 
  display = c("sp", "lc", "cn"), 
  main = "Triplot CCA spe ~ env3 - scaling 2")
text(-2.3, 2.6, "b", cex = 1.5)

par(mfrow = c(1, 2))
# CCA scaling 1 biplot without species (using lc site scores)
plot(spe.cca, 
  scaling = 1, 
  display = c("lc", "cn"), 
  main = "Biplot CCA spe ~ env3 - scaling 1"
)
# CCA scaling 2 biplot with species but without sites
plot(spe.cca, 
  scaling = 2, 
  display = c("sp", "cn"), 
  main = "Biplot CCA spe ~ env3 - scaling 2"
)

# Permutation test of the overall analysis
anova(spe.cca, permutations = how(nperm = 999))
# Permutation test of each axis
anova(spe.cca, by = "axis", permutations = how(nperm = 999))

###########==================== Diapos ===###

#### Selección hacia adelante CCA ####
# This function allows the use of factors like 'slo' in env3
cca.step.forward <- 
  ordistep(cca(spe ~ 1, data = env3), 
           scope = formula(spe.cca), 
           direction = "forward", 
           permutations = how(nperm = 199))


## Parsimonious CCA using ele, oxy and bod
spe.cca.pars <- cca(spe ~ ele + oxy + bod, data = env3)
anova(spe.cca.pars, permutations = how(nperm = 999))
anova(spe.cca.pars, permutations = how(nperm = 999), by = "axis")
# R-square â€“ like statistics
RsquareAdj(spe.cca.pars)
# Compare variance inflation factors
vif.cca(spe.cca)
vif.cca(spe.cca.pars)

# -----------------------------------------------------------------

## A simple function to compute CCA-based variation partitioning
## with bootstrap adjusted R-square

# This function is limited to two explanatory matrices

varpart.cca <- function(Y, X, W)
{
  # Computation of the three CCA
    yx.cca <- cca(Y, X)
    yw.cca <- cca(Y, W)
    yxw.cca <- cca(Y, cbind(X, W))

  # Computation of adjusted inertia
    fract.ab <- RsquareAdj(yx.cca)$adj.r.squared
    fract.bc <- RsquareAdj(yw.cca)$adj.r.squared
    fract.abc <- RsquareAdj(yxw.cca)$adj.r.squared
    fract.a <- fract.abc-fract.bc
    fract.b <- fract.ab-fract.a
    fract.c <- fract.abc-fract.ab
    fract.d <- 1-fract.abc

  # Output of results
    res <- matrix(0, 7, 1)
    rownames(res) <- 
        c("[ab]", "[bc]", "[abc]", "[a]", "[b]", "[c]", "[d]")
    colnames(res) <- "Value"
    res[1, 1] <- round(fract.ab, 4)
    res[2, 1] <- round(fract.bc, 4)
    res[3, 1] <- round(fract.abc, 4)
    res[4, 1] <- round(fract.a, 4)
    res[5, 1] <- round(fract.b, 4)
    res[6, 1] <- round(fract.c, 4)
    res[7, 1] <- round(fract.d, 4)

  res
}

# -----------------------------------------------------------------


## Three-dimensional interactive ordination plots
## (requires vegan3d package)

# Plot of the sites only (wa scores)
ordirgl(spe.cca.pars, type = "t", scaling = 1)

# Connect weighted average scores to linear combination scores
orglspider(spe.cca.pars, scaling = 1, col = "purple")

# Plot the sites (wa scores) with a clustering result
# Colour sites according to cluster membership
gr <- cutree(hclust(vegdist(spe.hel, "euc"), "ward.D2"), 4)
ordirgl(spe.cca.pars, 
  type = "t", 
  scaling = 1, 
  ax.col = "black", 
  col = gr + 1
)
# Connect sites to cluster centroids
orglspider(spe.cca.pars, gr, scaling = 1)

# Complete CCA 3D triplot
ordirgl(spe.cca.pars, type = "t", scaling = 2)
orgltext(spe.cca.pars, 
  display = "species", 
  type = "t", 
  scaling = 2, 
  col = "cyan"
)

# Plot species groups (Jaccard dissimilarity, useable in R mode)
gs <- 
  cutree(
      hclust(vegdist(t(spe), method = "jaccard"), "ward.D2"), 
      k = 4)
ordirgl(spe.cca.pars, 
         display = "species", 
         type = "t", 
         col = gs + 1)


################# ============== Diapos === ###

#### Análisis linear discriminante (LDA) ####

# Ward clustering result of Hellinger-transformed species data, 
# cut into 4 groups
gr <- cutree(hclust(vegdist(spe.hel, "euc"), "ward.D2"), k = 4)

# Environmental matrix with only 3 variables (ele, oxy and bod)
env.pars2 <- as.matrix(env2[, c(1, 9, 10)])

# Verify multivariate homogeneity of within-group covariance
# matrices using the betadisper() function {vegan}
env.pars2.d1 <- dist(env.pars2)
(env.MHV <- betadisper(env.pars2.d1, gr))
anova(env.MHV)
permutest(env.MHV)	# Permutational test

# Log transform ele and bod
env.pars3 <- cbind(log(env2$ele), env2$oxy, log(env2$bod))
colnames(env.pars3) <- c("ele.ln", "oxy", "bod.ln") 
rownames(env.pars3) <- rownames(env2)
env.pars3.d1 <- dist(env.pars3)
(env.MHV2 <- betadisper(env.pars3.d1, gr))
permutest(env.MHV2)

# Preliminary test :  do the means of the explanatory variable 
# differ among groups?
# Compute Wilk'S lambda test
# First way: with function Wilks.test() of package rrcov, Ï‡2 test
Wilks.test(env.pars3, gr)
# Second way: with function manova() of stats, which uses
#             an F-test approximation
lw <-  manova(env.pars3 ~ as.factor(gr))
summary(lw, test = "Wilks")


## Computation of LDA - identification functions (on unstandardized 
## variables)

env.pars3.df <- as.data.frame(env.pars3)
(spe.lda <- lda(gr ~ ele.ln + oxy + bod.ln, data = env.pars3.df))
# Alternate coding without formula interface :  
#    spe.lda <- lda(env.pars3.df, gr)
# The result object contains the information necessary to interpret 
# the LDA
summary(spe.lda)

# Display the group means for the 3 variables
spe.lda$means

# Extract the unstandardized identification functions (matrix C, 
# eq. 11.33 in Legendre and Legendre 2012)
(C <- spe.lda$scaling)

# Classification of two new objects (identification)
# A new object is created with two sites: 
#     (1) ln(ele) = 6.8, oxygen = 9 and ln(bod) = 0.8 
# and (2) ln(ele) = 5.5, oxygen = 10 and ln(bod) = 1.0
newo <- data.frame(c(6.8, 5.5), c(9, 10), c(0.8, 1))
colnames(newo) <- colnames(env.pars3)
newo
(predict.new <- predict(spe.lda, newdata = newo))


## Computation of LDA - discrimination functions (on standardized 
## variables)
 
env.pars3.sc <- as.data.frame(scale(env.pars3.df))
spe.lda2 <- lda(gr ~ ., data = env.pars3.sc)

# Display the group means for the 3 variables
spe.lda2$means

# Extract the classification functions
(C2 <- spe.lda2$scaling)

# Compute the canonical eigenvalues
spe.lda2$svd^2

# Position the objects in the space of the canonical variates
(Fp2 <- predict(spe.lda2)$x)
# alternative way :  Fp2 <- as.matrix(env.pars3.sc) %*% C2

# Classification of the objects
(spe.class2 <- predict(spe.lda2)$class)

# Posterior probabilities of the objects to belong to the groups
# (rounded for easier interpretation)
(spe.post2 <- round(predict(spe.lda2)$posterior, 2))

# Contingency table of prior versus predicted classifications
(spe.table2 <- table(gr, spe.class2))

# Proportion of correct classification (classification success)
diag(prop.table(spe.table2, 1))

# Plot the LDA results using the homemade function plot.lda()
par(mfrow=c(1,1))
plot.lda(lda.out = spe.lda2, 
  groups = gr, 
  plot.sites = 2, 
  plot.centroids = 1, 
  mul.coef = 2.35
)
################# ============== Diapos === ###

#### LDA con una clasificación basada en jacknife ####
(spe.lda.jac <- 
  lda(gr ~ ele.ln + oxy + bod.ln, 
      data = env.pars3.sc, 
      CV = TRUE))
summary(spe.lda.jac)

# Numbers and proportions of correct classification
spe.jac.class <- spe.lda.jac$class
spe.jac.table <- table(gr, spe.jac.class)
# Classification success
diag(prop.table(spe.jac.table, 1))

# Example of Legendre and Legendre (2012, p. 683)

grY <- c(1, 1, 1, 2, 2, 3, 3)
x1 <- c(1, 2, 2, 8, 8, 8, 9)
x2 <- c(2, 2, 1, 7, 6, 3, 3)
X <- as.data.frame(cbind(x1, x2))

# Computation of unstandardized identification functions
unstand.lda <- lda(grY ~ ., data = X)

# Computation of standardized discriminant functions
X.sc <- as.data.frame(scale(X))
stand.lda <- lda(grY ~ ., data = X.sc)


################ ===============Diapos  === ###


# Principal response curves (PRC) =================================

# Code from the prc() help file, with additional comments
# Chlorpyrifos experiment and experimental design :  Pesticide
# treatment in ditches (replicated) and followed over from 4 weeks
# before to 24 weeks after exposure 

# Extract the data (available in vegan)
data(pyrifos)

# Create factors for time (week) and treatment (dose). Create an
# additional factor "ditch" representing the mesocosm, for testing 
# purposes
week <-
  gl(11, 12, 
     labels = c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose <- 
  factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 
         11))
ditch <- gl(12, 1, length = 132)

# PRC
mod <- prc(pyrifos, dose, week)
mod            # Modified RDA
summary(mod)   # Results formatted as PRC

# PRC plot; at the right of it, only species with large total
# (log-transformed) abundances are reported
logabu <- colSums(pyrifos)
dev.new(title = "PRC", noRStudioGD = TRUE)
plot(mod, select = logabu > 200)

# Statistical test
# Ditches are randomized, we have a time series, and are only
# interested in the first axis
ctrl <- 
  how(plots = Plots(strata = ditch, type = "free"), 
      within = Within(type = "series"), nperm = 999)
anova(mod, permutations = ctrl, first = TRUE)



# Predictive co-correspondence analysis (CoCA) ====================

data(bryophyte)
data(vascular)

# Co-correspondence analysis is computed using the function coca()
# The default option is method = "predictive" 
(carp.pred <- coca(bryophyte ~ ., data = vascular))

# Leave-one-out cross-validation
crossval(bryophyte, vascular)

# Permutation test
(carp.perm <- permutest(carp.pred, permutations = 99))

# Only two significant axes: refit
(carp.pred <- coca(bryophyte ~ ., data = vascular, n.axes = 2))

# Extract the site scores and the species loadings used in
# the biplots
carp.scores <- scores(carp.pred)
load.bryo <- carp.pred$loadings$Y
load.plant <- carp.pred$loadings$X

 # We have generated two plots. As in ter Braak and Schaffers
 # (2004, Fig. 3), in both plots the site scores are derived
 # from the vascular plants (carp.scores$sites$X) and the
 # species scores are the "loadings with respect to
 # normalized site scores"

# Printing options:
?plot.predcoca  

dev.new(
  title = "Predictive co-correspondence analysis (CoCA)", 
  width = 12, 
  height = 6, 
  noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
plot(carp.pred, 
  type = "none", 
  main = "Bryophytes", 
  xlim = c(-2, 3), 
  ylim = c(-3, 2)
)
points(carp.scores$sites$X, pch = 16, cex = 0.5)
text(load.bryo, 
  labels = rownames(load.bryo), 
  cex = 0.7, 
  col = "red"
) 
plot(carp.pred, 
  type = "none", 
  main = "Vascular plants", 
  xlim = c(-2, 3), 
  ylim = c(-3, 2)
)
points(carp.scores$sites$X, pch = 16, cex = 0.5)
text(load.plant, 
  labels = rownames(load.plant), 
  cex = 0.7, 
  col = "blue"
) 

# Detach package cocorresp to avoid conflicts with ade4:
detach("package:cocorresp", unload = TRUE)
# If not sufficient:
unloadNamespace("cocorresp")


# Canonical correlation analysis (CCorA) ==========================

# Preparation of data (transformations to make variable 
# distributions approximately symmetrical)
envchem2 <- envchem
envchem2$pho <- log(envchem$pho)
envchem2$nit <- sqrt(envchem$nit)
envchem2$amm <- log1p(envchem$amm)
envchem2$bod <- log(envchem$bod)
envtopo2 <- envtopo
envtopo2$ele <- log(envtopo$ele)
envtopo2$slo <- log(envtopo$slo)
envtopo2$dis <- sqrt(envtopo$dis)

# CCorA (on standardized variables)
chem.topo.ccora <- 
  CCorA(envchem2, envtopo2, 
        stand.Y = TRUE, 
        stand.X = TRUE, 
        permutations = how(nperm = 999))
chem.topo.ccora

dev.new(
  title = "Canonical correlation analysis (CCorA)", 
  width = 9, 
  height = 6, 
  noRStudioGD = TRUE
)
biplot(chem.topo.ccora, plot.type = "biplot")



# Co-inertia analysis =============================================

# PCA on both matrices using ade4 functions
dudi.chem <- dudi.pca(envchem2, 
         scale = TRUE, 
         scannf = FALSE)
dudi.topo <- dudi.pca(envtopo2, 
         scale = TRUE, 
         scannf = FALSE)
# Cumulated relative variation of eigenvalues
cumsum(dudi.chem$eig / sum(dudi.chem$eig))
# Cumulated relative variation of eigenvalues
cumsum(dudi.topo$eig / sum(dudi.topo$eig))

# Are the row weights equal in the 2 analyses?
all.equal(dudi.chem$lw, dudi.topo$lw)

# Co-inertia analysis
coia.chem.topo <- coinertia(dudi.chem, dudi.topo, 
            scannf = FALSE, 
            nf = 2)
summary(coia.chem.topo)

# Relative variation on first eigenvalue
coia.chem.topo$eig[1] / sum(coia.chem.topo$eig)
# Permutation test
randtest(coia.chem.topo, nrepet = 999)

# Plot results
dev.new(title = "Co-inertia analysis", noRStudioGD = TRUE)
plot(coia.chem.topo)


# Multiple factor analysis (MFA) ==================================

# MFA on 3 groups of variables : 
# Regroup the 3 tables (Hellinger-transformed species, 
# physiographic variables, chemical variables) 
tab3 <- data.frame(spe.hel, envtopo, envchem)
dim(tab3)
# Number of variables in each group
(grn <- c(ncol(spe), ncol(envtopo), ncol(envchem)))

# Close the previous graphic windows
graphics.off()
# Compute the MFA without multiple plots
t3.mfa <- MFA(
 tab3,
 group = grn,
 type = c("c", "s", "s"),
 ncp = 2,
 name.group = c("Fish community", "Physiography", "Water quality"),
 graph = FALSE
)
t3.mfa
summary(t3.mfa)
t3.mfa$ind

# Plot the results
dev.new(title = "Partial axes", noRStudioGD = TRUE)
plot(t3.mfa,
     choix = "axes",
     habillage = "group",
     shadowtext = TRUE)
dev.new(title = "Sites", noRStudioGD = TRUE)
plot(
  t3.mfa,
  choix = "ind",
  partial = "all",
  habillage = "group")
dev.new(title = "Quantitative variables", noRStudioGD = TRUE)
plot(t3.mfa,
     choix = "var",
     habillage = "group",
     shadowtext = TRUE)
dev.new(title = "Groups", noRStudioGD = TRUE)
plot(t3.mfa, choix = "group")


# RV coefficients with tests (p-values above the diagonal of 
# the matrix)
rvp <- t3.mfa$group$RV
rvp[1, 2] <- coeffRV(spe.hel, scale(envtopo))$p.value
rvp[1, 3] <- coeffRV(spe.hel, scale(envchem))$p.value
rvp[2, 3] <- coeffRV(scale(envtopo), scale(envchem))$p.value
round(rvp[-4, -4], 6)

# Eigenvalues, scree plot and broken stick model
ev <- t3.mfa$eig[, 1]
names(ev) <- paste("MFA", 1 : length(ev))
dev.new(
  title = "MFA eigenvalues and broken stick model", 
  noRStudioGD = TRUE
)
screestick(ev, las = 2)

# Alternative to the standard, automatic MFA plots : 
# Plot only the significant variables (correlations)

# Select the most characteristic variables
aa <- dimdesc(t3.mfa, axes = 1:2, proba = 0.0001)

# Plot
dev.new(
  title = "MFA :  correlations among significant variables", 
  noRStudioGD = TRUE
)
varsig <- 
  t3.mfa$quanti.var$cor[unique(c(rownames(aa$Dim.1$quanti), 
	            rownames(aa$Dim.2$quanti))), ]
plot(varsig[, 1:2], 
  asp = 1, 
  type = "n", 
  xlim = c(-1, 1),
  ylim = c(-1, 1)
)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
symbols(0, 0, circles = 1, inches = FALSE, add = TRUE)
arrows(0, 0, varsig[, 1], varsig[, 2], length = 0.08, angle = 20)
for (v in 1 : nrow(varsig))
{
	if (abs(varsig[v, 1]) > abs(varsig[v, 2]))
	{
		if (varsig[v, 1] >=  0) pos <- 4
		else pos <- 2
	}
	else
	{
		if (varsig[v, 2] >=  0) pos <- 3
		else pos <- 1
	}
	text(varsig[v, 1], varsig[v, 2], 
	     labels = rownames(varsig)[v], 
	     pos = pos)
}



# RLQ and fourth-corner analyses ==================================

data(aravo)
dim(aravo$spe)
dim(aravo$traits)
dim(aravo$env)

# Preliminary analyses: CA, Hill-Smith and PCA
afcL.aravo <- dudi.coa(aravo$spe, scannf = FALSE)
acpR.aravo <- dudi.hillsmith(aravo$env, 
                 row.w = afcL.aravo$lw,
                 scannf = FALSE)
acpQ.aravo <- dudi.pca(aravo$traits, 
                 row.w = afcL.aravo$cw,
                 scannf = FALSE)

# RLQ analysis
rlq.aravo <- rlq(
                 dudiR = acpR.aravo, 
                 dudiL = afcL.aravo, 
                 dudiQ = acpQ.aravo,
                 scannf = FALSE)
dev.new(
  title = "RLQ", 
  noRStudioGD = TRUE
)
plot(rlq.aravo)
# Traits by environment crossed table
rlq.aravo$tab


# Since the plots are crowded, one can plot them one by one 
# in large graphical windows.

dev.new(title = "RLQ - site (L) scores", noRStudioGD = TRUE)
s.label(rlq.aravo$lR, 
  plabels.boxes.draw = FALSE, 
  ppoints.alpha = 0,
  psub.text = "a",
  psub.cex = 2, 
  psub.position = "topleft"
)
dev.new(title = "RLQ - species abundances", noRStudioGD = TRUE)
s.label(rlq.aravo$lQ, 
  plabels.boxes.draw = FALSE, 
  ppoints.alpha = 0,
  psub.text = "b",
  psub.cex = 2, 
  psub.position = "topleft"
)
dev.new(title = "RLQ - environmental variables", 
  noRStudioGD = TRUE)
s.arrow(rlq.aravo$l1,
  psub.text = "c",
  psub.cex = 2, 
  psub.position = "topleft"
)
dev.new(title = "RLQ - species traits", noRStudioGD = TRUE)
s.arrow(rlq.aravo$c1,
  psub.text = "d",
  psub.cex = 2, 
  psub.position = "topleft"
)

# Global test
randtest(rlq.aravo, nrepet = 999, modeltype = 6)

## Fourth-corner analysis (takes time with 49999 permutations!)
fourth.aravo <- fourthcorner(
                      tabR = aravo$env, 
                      tabL = aravo$spe, 
                      tabQ = aravo$traits,
                      modeltype = 6,
                      p.adjust.method.G = "none", 
                      p.adjust.method.D = "none", 
                      nrepet = 49999)
# Correction for multiple testing, here using FDR
fourth.aravo.adj <- p.adjust.4thcorner(
                      fourth.aravo,
                      p.adjust.method.G = "fdr", 
                      p.adjust.method.D = "fdr", 
                      p.adjust.D = "global") 
# Plot significant associations
dev.new(title = "Fourth-corner analysis", noRStudioGD = TRUE)
plot(fourth.aravo.adj, alpha = 0.05, stat = "D2")

# Biplot combining RLQ and fourth-corner results
dev.new(
  title = "Combining RLQ and Fourth-corner results", 
  noRStudioGD = TRUE
)
plot(fourth.aravo.adj, 
  x.rlq = rlq.aravo, 
  alpha = 0.05, 
  stat = "D2", 
  type = "biplot"
)


###################################################################

# RLQ and fourth-corner analyses (Doubs data) ######

summary(fishtraits)
rownames(fishtraits)
names(spe)
names(fishtraits)
tra <- fishtraits[ , 6:15]
tra

# Preliminary analyses: CA, Hill-Smith and PCA
afcL.doubs <- dudi.coa(spe, scannf = FALSE)
acpR.doubs <- dudi.hillsmith(env3,
                             row.w = afcL.doubs$lw,
                             scannf = FALSE)
acpQ.doubs <- dudi.pca(tra, 
                       row.w = afcL.doubs$cw,
                       scannf = FALSE)

# RLQ analysis
rlq.doubs <- rlq(
  dudiR = acpR.doubs, 
  dudiL = afcL.doubs, 
  dudiQ = acpQ.doubs,
  scannf = FALSE)
dev.new(
  title = "RLQ",
  width = 12,
  height = 7,
  noRStudioGD = TRUE
)
plot(rlq.doubs)
# Traits by environment crossed table
rlq.doubs$tab


# Since the plots are crowded, one can plot them one by one 
# in large graphical windows:

dev.new(title = "RLQ - site (L) scores", noRStudioGD = TRUE)
s.label(rlq.doubs$lR, 
        plabels.boxes.draw = FALSE, 
        ppoints.alpha = 0,
        psub.text = "a",
        psub.cex = 2, 
        psub.position = "topleft"
)
dev.new(title = "RLQ - species abundances", noRStudioGD = TRUE)
s.label(rlq.doubs$lQ, 
        plabels.boxes.draw = FALSE, 
        ppoints.alpha = 0,
        psub.text = "b",
        psub.cex = 2, 
        psub.position = "topleft"
)
dev.new(title = "RLQ - environmental variables", 
        noRStudioGD = TRUE)
s.arrow(rlq.doubs$l1,
        psub.text = "c",
        psub.cex = 2, 
        psub.position = "topleft"
)
dev.new(title = "RLQ - species traits", noRStudioGD = TRUE)
s.arrow(rlq.doubs$c1,
        psub.text = "d",
        psub.cex = 2, 
        psub.position = "topleft"
)

# Global test
randtest(rlq.doubs, nrepet = 999, modeltype = 6)


# Fourth-corner analysis (takes time with 49999 permutations!) ######
?fourthcorner

fourth.doubs2 <- fourthcorner(
  tabR = env3,
  tabL = spe,
  tabQ = tra,
  modeltype = 2,
  p.adjust.method.G = "fdr",
  p.adjust.method.D = "fdr",
  nrepet = 49999
)
fourth.doubs2
summary(fourth.doubs2)

fourth.doubs <- fourthcorner(
  tabR = env2, 
  tabL = spe, 
  tabQ = tra,
  modeltype = 6,
  p.adjust.method.G = "none", 
  p.adjust.method.D = "none", 
  nrepet = 49999)
# Correction for multiple testing, here using FDR
fourth.doubs.adj <- p.adjust.4thcorner(
  fourth.doubs,
  p.adjust.method.G = "fdr", 
  p.adjust.method.D = "fdr", 
  p.adjust.D = "global") 

fourth.doubs.adj
summary(fourth.doubs.adj)

# Plot
dev.new(title = "Fourth-corner analysis", noRStudioGD = TRUE)
plot(fourth.doubs.adj, alpha = 0.05, stat = "D2")
plot(fourth.doubs2, stat = "D2")
plot(fourth.doubs2, stat = "G")

# Biplot combining RLQ and fourth-corner results
dev.new(
  title = "Combining RLQ and Fourth-corner results", 
  noRStudioGD = TRUE
)
plot(fourth.doubs.adj, 
     x.rlq = rlq.doubs, 
     alpha = 0.05, 
     stat = "D2", 
     type = "biplot"
)

plot(fourth.doubs2,
     x.rlq = rlq.doubs,
     alpha = 0.05,
     stat = "D2",
     type = "biplot"
)

   ```
