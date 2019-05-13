### ANÁLISIS EXPLORATORIO DE DATOS

# Cargar paquetes, funciones y datos ===============================
library(vegan)
library(RgoogleMaps)
library(googleVis)
library(labdsv)
library(leaflet)

# Funciones adicionales
source("funciones/panelutils.R")

# Cargar los datos
load("Doubs.RData")  

# El archivo Doubs.RData contiene los siguientes objetos:
#   spe: datos de especies (comunidad) (abundancia de peces)
#   env: datos ambientales
#   spa: datos espaciales coordenadas cartesianas
#   fishtraits: rasgos funcionales de especies de peces
#   latlong: datos espaciales - latitud y longitud


# Exploración de los datos usando las funciones básicas de R =============

spe                       # Muestra los datos completos
                          # ¡No se recomienda para muchos datos
spe[1:5, 1:10]            # Muestra solo 5 lineas y 10 columnas
head(spe)                 # Muestra las primeras 6 filas
tail(spe)                 # Muestra las últimas 6 filas
nrow(spe)                 # Número de filas (sitios)
ncol(spe)                 # Número de columnas (especies)
dim(spe)                  # Dimensiones de la base de datos (filas y columnas)
colnames(spe)             # Nombres de las columnas (especies)
rownames(spe)             # Nombres de las filas (sitios)
summary(spe)              # Estadísticas descriptivas para columnas


# Distribución general de abundancias (dominancia) ============

# Valores mínimos y máximos de abundancia en el set de datos
range(spe)
# Valores mínimos y máximos para cada especie
apply(spe, 2, range)
# Cuenta los casos para clase de abundancia
(ab <- table(unlist(spe)))

# Gráfica de barras de la distribución, todas las especies juntas
barplot(ab, 
  las = 1,
  xlab = "Clase de abundancia",
  ylab = "Frecuencia",
  col = gray(5 : 0 / 5)
)
# Número de ausencias
sum(spe == 0)
# Proporción de ceros en la base de datos
sum(spe == 0) / (nrow(spe) * ncol(spe))


# Mapa de los sitios ===============================

# Coordenadas geográficas x y y de la base de datos spa
plot(spa, 
  asp = 1, 
  type = "n", 
  main = "Site Locations", 
  xlab = "x coordinate (km)", 
  ylab = "y coordinate (km)"
)
# Agregar una linea azul conectando los sitios a lo largo del rio Dub
lines(spa, col = "light blue")
# Agregar los nombres de los sitios
text(spa, row.names(spa), cex = 1, col = "red")
text(68, 20, "Upstream", cex = 1.2, col = "red")
text(15, 35, "Downstream", cex = 1.2, col = "red")

# Sitios proyectados en Google Maps

# 1. Usando googlevis

nom <- latlong$Site
latlong2 <- paste(latlong$LatitudeN, latlong$LongitudeE, sep = ":")
df <- data.frame(latlong2, nom, stringsAsFactors = FALSE)

mymap1 <- gvisMap(df, 
  locationvar = "latlong2", 
  tipvar = "nom", 
  options = list(showTip = TRUE)
)

plot(mymap1)


# 2. Usando Leaflet

m <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(lng=latlong$LongitudeE, lat=latlong$LatitudeN)
m


# Mapas de algunas especies de peces =======================================

par(mfrow = c(2,2))
# Mapa de cuatro especies
plot(spa, 
  asp = 1, 
  cex.axis = 0.8, 
  col = "brown", 
  cex = spe$Satr, 
  main = "Brown trout", 
  xlab = "x coordinate (km)", 
  ylab = "y coordinate (km)"
)
lines(spa, col = "light blue")
plot(spa, 
  asp = 1, 
  cex.axis = 0.8, 
  col = "brown", 
  cex = spe$Thth, 
  main = "Grayling", 
  xlab = "x coordinate (km)", 
  ylab = "y coordinate (km)"
)
lines(spa, col = "light blue")
plot(spa, 
  asp = 1, 
  cex.axis = 0.8, 
  col = "brown", 
  cex = spe$Baba, 
  main = "Barbel", 
  xlab = "x coordinate (km)", 
  ylab = "y coordinate (km)"
)
lines(spa, col = "light blue")
plot(spa, 
  asp = 1, 
  cex.axis = 0.8, 
  col = "brown", 
  cex = spe$Abbr, 
  main = "Common bream", 
  xlab = "x coordinate (km)", 
  ylab = "y coordinate (km)"
)
lines(spa, col = "light blue")


# Compara especies: número de presencias ==========================

# Calcula el número de sitios donde cada especie está presente
# Suma por columnas, el segundo argumento de apply(), MARGIN, está en 2
spe.pres <- apply(spe > 0, 2, sum)
# Ordena los resultados en orden creciente
sort(spe.pres)
# Calcula el porcentaje de frecuencias
spe.relf <- 100 * spe.pres/nrow(spe)
# Redonde el resultados ordenados a 1 dígito
round(sort(spe.relf), 1)

par(mfrow = c(1,2))
hist(spe.pres, 
  main = "Species Occurrences", 
  right = FALSE, 
  las = 1, 
  xlab = "Number of occurrences", 
  ylab = "Number of species", 
  breaks = seq(0, 30, by = 5),
  col = "bisque"
)
hist(spe.relf, 
  main = "Species Relative Frequencies", 
  right = FALSE, 
  las = 1,
  xlab = "Frequency of occurrences (%)", 
  ylab = "Number of species",
  breaks = seq(0, 100, by = 10),
  col = "bisque"
)


# Compara sitios: riqueza de especies =================================

# Calcula el número de especies de cada sitio
# Para sumar las filas, el segundo argumento de apply(), MARGIN, está en 1
sit.pres <- apply(spe > 0, 1, sum)
# Ordena los resultados en orden ascendente
sort(sit.pres)

par(mfrow = c(1, 2))
# Grafica la riqueza de especies vs la posición de los sitios a lo largo del rio
plot(sit.pres,type = "s",
  las = 1, 
  col = "gray",
  main = "Species Richness vs. \n Upstream-Downstream Gradient",
  xlab = "Site numbers", 
  ylab = "Species richness"
)
text(sit.pres, row.names(spe), cex = .8, col = "red")

# Usa las coordenadas geográficas para graficar un mapa de burbujas
plot(spa, 
  asp = 1, 
  main = "Map of Species Richness", 
  pch = 21, 
  col = "white", 
  bg = "brown", 
  cex = 5 * sit.pres / max(sit.pres), 
  xlab = "x coordinate (km)", 
  ylab = "y coordinate (km)"
)
lines(spa, col = "light blue")


# Transformación y estandarización de los datos de especies ==========

?decostand

spe[1:5, 2:4]
# Transformarlas a presencia ausencia (1-0)
spe.pa <- decostand(spe, method = "pa")
spe.pa[1:5, 2:4]


## Estandarización por columnas (especies)

# Escalar las abundancias dividiendolas por el valor máximo de cada epecie
# Nota: MARGIN = 2 (columna, valor default) para argumento "max"
spe.scal <- decostand(spe, "max")
spe.scal[1:5, 2:4]
# Muestra el máximo de cada columna transformada
apply(spe.scal, 2, max)
# Escalar las abundancias dividiendolas por el total de especies (abundancia relativa
# por especies)
# Nota: Cambiamos el argumento MARGIN = 1 de "total"
spe.relsp <- decostand(spe, "total", MARGIN = 2)
spe.relsp[1:5, 2:4]
# Muestra la suma por columna
colSums(spe.relsp)

## Estandarización por filas (sitios)

# Escala las abundancias dividiendolas por el total de sitios (perfiles de abundancia
# relativa por sitio)
spe.rel <- decostand(spe, "total") # default MARGIN = 1
spe.rel[1:5, 2:4]
# Mostrar la suma de las filas para ver si el escalamiento funciona
rowSums(spe.rel)

# Dar una longitud (norma) de 1 a cada fila
# Esta es llamada la transformación de cuerda
spe.norm <- decostand(spe, "normalize") # default MARGIN = 1
spe.norm[1 : 5, 2 : 4]
# Verificar la norma de las filas transformadas
vec.norm <- function(x) sqrt(sum(x ^ 2))
apply(spe.norm, 1, vec.norm)

# Calcular la raiz cuadrada de las abundancias relativas por sitios
# Hellinger transformation
spe.hel <- decostand(spe, "hellinger")
spe.hel[1:5, 2:4]
apply(spe.hel, 1, vec.norm)

# Chord y Hellinger son las más recomendadas ver: Legendre and Gallagher (2001)

## Doble estandarización por columnas y filas

# Chi-square transformation
spe.chi <- decostand(spe, "chi.square")
spe.chi[1:5, 2:4]
# Que pasa con el sitio 8 donde no se encontraron especies
spe.chi[7:9, ]
# Se producen valores de 0 para 0/0 en lugar de NA 

# Wisconsin standardization
# Las abundancias son primero arregladas por la especie máxima y luego por el total 
# del sitio 
spe.wis <- wisconsin(spe)
spe.wis[1:5, 2:4]


## Boxplots de las abundancias transformadas de una especie comun (the stone loach, 
# species #4)
par(mfrow = c(2,2))
boxplot(spe$Babl,
  sqrt(spe$Babl), 
  log1p(spe$Babl),
  las = 1, 
  main = "Simple transformations",
  names = c("raw data", "sqrt", "log"), 
  col = "bisque"
)
boxplot(spe.scal$Babl, 
  spe.relsp$Babl,
  las = 1, 
  main = "Standardizations by species",
  names = c("max", "total"), 
  col = "lightgreen"
)
boxplot(spe.hel$Babl, 
  spe.rel$Babl, 
  spe.norm$Babl,
  las = 1, 
  main = "Standardizations by sites",
  names = c("Hellinger", "total", "norm"), 
  col = "lightblue"
)
boxplot(spe.chi$Babl, 
  spe.wis$Babl,
  las = 1, 
  main = "Double standardizations",
  names = c("Chi-square", "Wisconsin"), 
  col = "orange"
)

## Graficar las abundancias crudas y transformadas a lo largo del rio
par(mfrow = c(2, 2))
plot(env$dfs, 
  spe$Satr, 
  type = "l", 
  col = 4, 
  main = "Raw data",
  xlab = "Distance from the source [km]", 
  ylab = "Raw abundance code"
)
lines(env$dfs, spe$Thth, col = 3)
lines(env$dfs, spe$Baba, col = "orange")
lines(env$dfs, spe$Abbr, col = 2)
lines(env$dfs, spe$Babl, col = 1, lty = "dotted")

plot(env$dfs, 
  spe.scal$Satr, 
  type = "l", 
  col = 4, 
  main = "Species abundances ranged 
  by maximum",
  xlab = "Distance from the source [km]", 
  ylab = "Ranged abundance"
)
lines(env$dfs, spe.scal$Thth, col = 3)
lines(env$dfs, spe.scal$Baba, col = "orange")
lines(env$dfs, spe.scal$Abbr, col = 2)
lines(env$dfs, spe.scal$Babl, col = 1, lty = "dotted")

plot(env$dfs, 
  spe.hel$Satr, 
  type = "l", 
  col = 4, 
  main =  "Hellinger-transformed abundances", 
  xlab = "Distance from the source [km]", 
  ylab = "Standardized abundance"
)
lines(env$dfs, spe.hel$Thth, col = 3)
lines(env$dfs, spe.hel$Baba, col = "orange")
lines(env$dfs, spe.hel$Abbr, col = 2)
lines(env$dfs, spe.hel$Babl, col = 1, lty = "dotted")

plot(env$dfs, 
  spe.chi$Satr, 
  type = "l", 
  col = 4, 
  main = "Chi-square-transformed abundances", 
  xlab = "Distance from the source [km]", 
  ylab = "Standardized abundance"
)
lines(env$dfs, spe.chi$Thth, col = 3)
lines(env$dfs, spe.chi$Baba, col = "orange")
lines(env$dfs, spe.chi$Abbr, col = 2)
lines(env$dfs, spe.chi$Babl, col = 1, lty = "dotted")


## Conversion de las abundancias de peces usando una escala arbitraria
current <- c(0, 1, 2, 3, 4, 5)
converted <- c(0, 1 ,5, 10, 20, 50)
spe.conv <- vegtrans(spe, current, converted)


# Mapas de burbujas de algunas variables ambientales =====================
par(mfrow = c(2, 2))
plot(spa, 
  asp = 1, 
  cex.axis = 0.8, 
  main = "Elevation", 
  pch = 21, 
  col = "white", 
  bg = "red", 
  cex = 5 * env$ele / max(env$ele), 
  xlab = "x", 
  ylab = "y"
)
lines(spa, col = "light blue")
plot(spa, 
  asp = 1, 
  cex.axis = 0.8, 
  main = "Discharge", 
  pch = 21, 
  col = "white", 
  bg = "blue",
  cex = 5 * env$dis / max(env$dis),
  xlab = "x", 
  ylab = "y"
)
lines(spa, col = "light blue")
plot(spa, 
  asp = 1, 
  cex.axis = 0.8, 
  main = "Oxygen", 
  pch = 21, 
  col = "white", 
  bg = "green3",
  cex = 5 * env$oxy / max(env$oxy),
  xlab =  "x", 
  ylab = "y"
)
lines(spa, col = "light blue")
plot(spa, 
  asp = 1,
  cex.axis = 0.8, 
  main = "Nitrate", 
  pch = 21,
  col = "white", 
  bg = "brown",
  cex = 5 * env$nit / max(env$nit),
  xlab = "x", 
  ylab = "y"
)
lines(spa, col = "light blue")


# Line plots ======================================================
par(mfrow = c(2, 2))
plot(env$dfs, env$ele, 
  type = "l", 
  xlab = "Distance from the source (km)", 
  ylab = "Elevation (m)", 
  col = "red", main = "Elevation"
)
plot(env$dfs, env$dis, 
  type = "l", 
  xlab = "Distance from the source (km)", 
  ylab = "Discharge (m3/s)", 
  col = "blue", 
  main = "Discharge"
)
plot(env$dfs, env$oxy, 
  type = "l", 
  xlab = "Distance from the source (km)", 
  ylab = "Oxygen (mg/L)", 
  col = "green3", 
  main = "Oxygen"
)
plot(env$dfs, env$nit, 
  type = "l", 
  xlab = "Distance from the source (km)", 
  ylab = "Nitrate (mg/L)", 
  col = "brown", 
  main = "Nitrate"
)


# Gráficos de dispersión de todas las variables ambientales ==========

# Gráficos bivariados con histogramas en las diagonal y ajuste suavizado
pairs(env, 
  panel = panel.smooth, 
  diag.panel = panel.hist,
  main = "Bivariate Plots with Histograms and Smooth Curves"
)


# Transformación simple de una variable ambiental ==============

range(env$slo)
# Transformación log de la variable pendiente (y = ln(x))
# Compara los histogramas y los boxplots de los valores originales y transformados
par(mfrow = c(2, 2))
hist(env$slo, 
  col = "bisque", 
  right = FALSE
)
hist(log(env$slo), 
  col = "light green", 
  right = FALSE, 
  main = "Histogram of ln(env$slo)"
)
boxplot(env$slo, 
  col = "bisque", 
  main = "Boxplot of env$slo", 
  ylab = "env$slo"
)
boxplot(log(env$slo), 
  col = "light green", 
  main = "Boxplot of ln(env$slo)",
  ylab = "log(env$slo)"
)


# Estandarización de todas las variables ambientales ==================

# Centrar y escalar = estandarizar las variables (valor z)
env.z <- decostand(env, "standardize")
apply(env.z, 2, mean)	# means = 0
apply(env.z, 2, sd)	# standard deviations = 1

# Misma estandarización usando las función () que regresa una matriz
env.z <- as.data.frame(scale(env))
