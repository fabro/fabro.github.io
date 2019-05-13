# Este script instala o proporciona pautas para instalar todos los paquetes necesarios
# para ejecutar los códigos provisto del curso.

# Los pasos 1 a 2 deben ejecutarse solo una vez al instalar o actualizar R.
# El paso 3 no es obligatorio.


# 1. Actualizar paquetes instalados
#    -------------------------
update.packages(checkBuilt = TRUE, ask = FALSE)


# 2. Instalar paquetes desde el sitio principal de CRAN
#    ----------------------------------------

install.packages(
  c(
    "ade4",
    "adegraphics",
    "adespatial",
    "agricolae",
    "ape",
    "cluster",
    "cocorresp",
    "colorspace",
    "dendextend",
    "ellipse",
    "FactoMineR",
    "FD",
    "gclus",
    "ggplot2",
    "googleVis",
    "igraph",
    "indicspecies",
    "labdsv",
    "leaflet",
    "MASS",
    "missMDA",
    "mvpart",
    "MVPARTwrap",
    "picante",
    "pvclust",
    "RColorBrewer",
    "rgexf",
    "RgoogleMaps",
    "rioja",
    "rrcov",
    "SoDA",
    "spdep",
    "taxize",
    "vegan",
    "vegan3d",
    "vegclust",
    "vegetarian"
  ),
  dependencies = TRUE,
  type = "both"
)


# Instale mvpart y MVPARTwrap que ya no están disponibles en CRAN:
# En las máquinas con Windows, primero se deben instalar Rtools (3.4 y superiores). 
# Ir a: https://cran.r-project.org/bin/windows/Rtools/
# Después:
install.packages("devtools")
library(devtools)
install_github("cran/mvpart", force = TRUE)
install_github("cran/MVPARTwrap", force = TRUE)

#### A partir de aquí es opcional ####
# 3. OPCIONAL (para usuarios avanzados): instale todos los paquetes R de 
# Environmetrics, una vista de tareas CRAN para el análisis de datos ecológicos y 
# ambientales
# Ver: http://cran.r-project.org/web/views/Environmetrics.html
#    --------------------------------------------------------------------

install.packages("ctv")
library(ctv)
update.views("Environmetrics")

# Otras vistas de tareas CRAN potencialmente útiles...
update.views("Cluster")
update.views("Multivariate")
update.views("Spatial")
update.views("MachineLearning")
