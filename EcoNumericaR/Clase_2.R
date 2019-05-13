### CLASE 2: MEDIDAS DE ASOCIACIÓN

library(ade4)
library(adespatial)
library(vegan)
library(gclus)
library(cluster)
library(FD)

source("funciones/coldiss.R")
source("funciones/panelutils.R")

load("Doubs.Rdata")
# Quitar un sitio vacio
spe <- spe[-8, ]
env <- env[-8, ]
spa <- spa[-8, ]


### Modo Q: matrices de disimilitud

## Modo Q medidas de disimilitud y distancia para datos (semi) cuantitativos

# Matriz de disimilitud de la diferencia de porcentajes (Bray-Curtis) on datos crudos
spe.db <- vegdist(spe)	# method = "bray" (default)
head(spe.db)
# En abundancias log-transformadas
spe.dbln <- vegdist(log1p(spe))
head(spe.dbln)
# Matriz de distancia de cuerda
spe.dc <- dist.ldc(spe, "chord")
   # En vegan:
   spe.norm <- decostand(spe, "nor")
   spe.dc <- dist(spe.norm)
head(spe.dc)
# Matriz de distancia Hellinger
spe.dh <- dist.ldc(spe) # Hellinger is the default distance
   # En vegan:
   spe.hel <- decostand(spe, "hel")
   spe.dh <- dist(spe.hel)
head(spe.dh)
# Matriz de distancia log-cuerda
spe.logchord <- dist.ldc(spe, "log.chord")
   # En vegan:
   spe.ln <- log1p(spe)
   spe.ln.norm <- decostand(spe.ln, "nor")
   spe.logchord <- dist(spe.ln.norm)
head(spe.logchord)

## _____ diapos ___________

## Modo Q medidas de disimilaridad para datos binarios

# Matriz de disimilitud Jaccard
spe.dj <- vegdist(spe, "jac", binary = TRUE)
head(spe.dj)
head(sqrt(spe.dj))
# function dist()
spe.dj2 <- dist(spe, "binary")
head(spe.dj2)
# function dist.binary()
spe.dj3 <- dist.binary(spe, method = 1)
head(spe.dj3)
# Matriz de disimilitud Sorensen function dist.ldc()
spe.ds <- dist.ldc(spe, "sorensen")
# function vegdist()
spe.ds2 <- vegdist(spe, method = "bray", binary = TRUE)
# function dist.binary()
spe.ds3 <- dist.binary(spe, method = 5)
head(spe.ds)
head(spe.ds2)
head(sqrt(spe.ds2))
head(spe.ds3)
# Matriz de disimilitud Ochiai
spe.och <- dist.ldc(spe, "ochiai")   # or
spe.och <- dist.binary(spe, method = 7)
head(spe.och)


## Muestra gráfica de las matrices de asociación

# Gráficos de color (heat maps, diagramas trellis)

## Compara las matrices de disimilitud y distancia obtenidas de los datos. Cuatro
## colores son usados con intervalos de igual longitud.

# Percentage difference (aka Bray-Curtis)
coldiss(spe.db, byrank = FALSE, diag = TRUE)

# Same but on log-transformed data
coldiss(spe.dbln, byrank = FALSE, diag = TRUE)

# Chord distance matrix
coldiss(spe.dc, byrank = FALSE, diag = TRUE)

# Hellinger distance matrix
coldiss(spe.dh, byrank = FALSE, diag = TRUE)

# log-chord distance matrix
coldiss(spe.logchord, byrank = FALSE, diag = TRUE)

# Jaccard distance matrix
coldiss(spe.dj, byrank = FALSE, diag = TRUE)

# Simple matching dissimilarity
# (called the Sokal and Michener index in ade4)
# El uso de un coeficiente de asociación simétrico (dobles ceros con importancia)
spe.s1 <- dist.binary(spe, method = 2)
coldiss(spe.s1 ^ 2, byrank = FALSE, diag = TRUE)


#_____________Diapos____________________________

# Quitar la variable "dfs"
env2 <- env[, -1]

# Matriz de distancia euclidiana en la base de datos estandarizada env2
env.de <- dist(scale(env2))
coldiss(env.de, nc = 16, diag = TRUE)

# Matriz de distancia Hellinger 
# Usar nc= 16 (número de colores)
coldiss(spe.dh, nc = 16, diag = TRUE)

# Matriz de distancia euclideana en coordenadas espaciales (2D)
spa.de <- dist(spa)
coldiss(spa.de, nc = 16, diag = TRUE)

# Matriz de distancia euclideana en la distancia al origen (1D)
dfs.df <- as.data.frame(env$dfs, row.names = rownames(env))
riv.de <- dist(dfs.df)
coldiss(riv.de, nc = 16, diag = TRUE)


## Ejemplos con datos artificiales

# Calcular cinco variables binarias con 30 objetos cada una
# Cada variable tiene números predefinidos de 0 y 1
# Variable 1: 10 x1 y 20 x 0; el orden es aleatorio
var1 <- sample(c(rep(1, 10), rep(0, 20)))
# Variable 2: 15 x 0 and 15 x 1, un bloque cada una
var2 <- c(rep(0, 15), rep(1, 15))
# Variable 3: alternación de 3 x 1 y 3 x 0 hasta 30 objetos
var3 <- rep(c(1, 1, 1, 0, 0, 0), 5)
# Variable 4: alternación de 5 x 1 y 10 x 0 hasta 30 objetos
var4 <- rep(c(rep(1, 5), rep(0, 10)), 2)
# Variable 5: 16 objetos con distribución aleatoria de 7 x 1 y 9 x 0, seguido por 4 x
# 0 y 10 x 1
var5.1 <- sample(c(rep(1, 7), rep(0, 9)))
var5.2 <- c(rep(0, 4), rep(1, 10))
var5 <- c(var5.1, var5.2)

# Juntar todas las bases de datos
(dat <- data.frame(var1, var2, var3, var4, var5))
dim(dat)

# Calculo de una matriz de coeficientes de coincidencia simple (Sokal and Michener 
# index in ade4)
dat.s1 <- dist.binary(dat, method = 2)
coldiss(dat.s1, diag = TRUE)

#_______________Diapos____________________

# Datos ficticios para el índice Gower (S15)
# Desviación normal aleatoria con media cero y unidad de desviación estandard
var.g1 <- rnorm(30, 0, 1)
# Desviación uniforme estandard de 0 a 5
var.g2 <- runif(30, 0, 5)
# Factor con 3 niveles (10 objetos cada una)
var.g3 <- gl(3, 10, labels = c("A", "B", "C"))
# Factor con 2 niveles, ortogonal a var.g3
var.g4 <- gl(2, 5, 30, labels = c("D", "E"))
(dat2 <- data.frame(var.g1, var.g2, var.g3, var.g4))
summary(dat2)

# Cálculo de una matriz de disimilitud Gower usando la función daisy()
dat2.S15 <- daisy(dat2, "gower")
range(dat2.S15)
coldiss(dat2.S15, diag = TRUE)

# Matriz de datos con los factores ortogonales solamente
dat2partial.S15 <- daisy(dat2[, 3:4], "gower")
coldiss(dat2partial.S15, diag = TRUE)
head(as.matrix(dat2partial.S15))

# ¿Cuales son los valores de disimilitud en la matriz dat2partial.S15?
levels(factor(dat2partial.S15))

# Cálculo de una matriz de disimilitud Gower usando la función gowdis() del paquete 
# FD
?gowdis
dat2.S15.2 <- gowdis(dat2)
range(dat2.S15.2)
coldiss(dat2.S15.2, diag = TRUE)

# Matriz de datos con los factores ortogonales solamente
dat2partial.S15.2 <- gowdis(dat2[, 3:4])
coldiss(dat2partial.S15.2, diag = TRUE)
head(as.matrix(dat2partial.S15.2))

levels(factor(dat2partial.S15.2))

#__________________Diapos________________--

#### Modo R matrices de disimilitud ####

spe.t <- t(spe)

# Pre-transformación chi cuadrada seguido por distancia Euclidiana
spe.t.chi <- decostand(spe.t, "chi.square")
spe.t.D16 <- dist(spe.t.chi)
coldiss(spe.t.D16, diag = TRUE)

# Índice Jaccar en presencia-ausencia de peces
spe.t.S7 <- vegdist(spe.t, "jaccard", binary = TRUE)
coldiss(spe.t.S7, diag = TRUE)

## Matrices de correlación modo R

# Correlación linear Pearson entre variables ambientales
env.pearson <- cor(env)	# default method = "pearson"
round(env.pearson, 2)

# Reordenar las variables antes de graficarlas
env.o <- order.single(env.pearson)

# pairs() es una función para graficar una matriz de gráficos bivariados.
# panelutils.R es un grupo de funciones que añaden características útiles a
# pairs():
#   upper.panel = panel.cor: para imprimir coeficientes de correlación en el
#     panel superior con nivel de significancia;
#   diag.panel = panel.hist: para graficar histogramas de las variables en la
#     diagonal.
# Especifica el método para elegir los coeficientes de correlación por default
# method = "pearson", otras opciones son "spearman" y "kendall".
# Para graficar en tonos grises en lugar de colores usa no.col = TRUE y 
# lower.panel = panel.smoothb.

pairs(
  env[, env.o],
  lower.panel = panel.smooth,
  upper.panel = panel.cor,
  diag.panel = panel.hist,
  main = "Pearson Correlation Matrix"
)

# 
# Correlación de rangos Kendall tau entre variables ambientales sin colores
env.ken <- cor(env, method = "kendall")
env.o <- order.single(env.ken)
pairs(
  env[, env.o],
  lower.panel = panel.smoothb,
  upper.panel = panel.cor,
  no.col = TRUE,
  method = "kendall",
  diag.panel = panel.hist,
  main = "Kendall Correlation Matrix"
)
