######################
## MODELO ECTOTHERM ##
######################

#### Fun��es de transfer�ncia de calor y temperatura ####

Q_rad <- function(A, a, D, S) {A * a * (D + S)}  
# A: superficie exposta, a: absor��o de radia��o curta, D: radia��o direta, S: radia��o difusa

Q_IR <- function(A, epsilon, boltzmann, Tb) {-A * epsilon * boltzmann * Tb^4}
# epsilon: emissividade IR, bolzmann: constante Stefan-Botzmann; Tb: Temperatura corporal

Q_conv <- function(A, hc, Ta, Tb) {-A * hc * (Tb - Ta)}
# hc: coeficiente de convec��o, Ta: temperatura do ar

Q_net <- function(Q_absorcion_sol, Q_emision_IR, Q_conveccion) {Q_absorcion_sol + Q_emision_IR + Q_conveccion}

Tb_differential <- function(Q_net, M, C) {1/(C*M) * Q_net}
# Q_net: calor total, M: massa corporal, C: capacidade de calor


#### Par�metros do modelo ####
# Morfologia animal

M = 40                  # Massa corporal (g)
A = 1e-4 * 10 * M^(2/3) # Superficie (m^2)
l = A^(1/2)             # comprimento do corpo (m)
C = 0.3                 # Capacidade de calor (W �C^-1)

# Temperatura ambiente

Ta = 20  # Temperatura do ar (�C)
v = 2    # Velocidade do vento (m s^-1)

# Radia��o solar

a = 0.9       # Aabsorb�ncia radia��o onda solar curta
D = 400       # Radia��o direta (W m^-2)
s = D * 0.2   # Radia��o difusa (20%)

# Emiss�o de IR

epsilon = 0.9       # Emissividade IR
boltzmann = 5.67e-8 # Constante Stefan-Botzmann (W m^-2 �C^-4)

# Coeficiente de convec��o (c�lculo emp�rico)

hc = 6.77 * v^0.6 * l^-0.4 # Coeficiente de convec��o (rela��o v�lida para r�pteis, Spotila et al. 1992)

#### C�lculo de Tb ####

Tb = 20 # N�s come�amos com uma temperatura inicial do 20�C
q_rad <- Q_rad(A, a, D, s)  # Calculamos a energia absorvida da radia��o solar
q_IR <- Q_IR(A, epsilon, boltzmann, Tb) # Calculamos a energia t�rmica emitida, dada a corrente Tb
q_conv <- Q_conv(A, hc, Ta, Tb) # Calculamos a transfer�ncia (perda ou ganho) de calor por convec��o

q_net <- Q_net(q_rad, q_IR, q_conv) # Calculamos o calor total

# A temperatura ir� mudar (aumentar ou diminuir) em
dTb <- Tb_differential(q_net, M, C)
dTb

# Portanto, a temperatura no pr�ximo instante ser�
Tb <- Tb + dTb

#################################
#### C�lculo iterativo de Tb ####
#################################

t_total <- 240  # Tempo total de exposi��o (segundos)
Tb <- numeric(t_total) # Vector (de comprimento t_total) para armazenar os resultados
Tb[1] <- 20 # N�s armazenamos a temperatura inicial na primeira entrada do vetor

for(t in 2:t_total){ # Para cada entrada de vetor (come�ando com a segunda entrada)
  q_rad <- Q_rad(A, a, D, s)  # Repetimos os c�lculos anteriores
  q_IR <- Q_IR(A, epsilon, boltzmann, Tb[t-1]) 
  q_conv <- Q_conv(A, hc, Ta, Tb[t-1]) 
  q_net <- Q_net(q_rad, q_IR, q_conv) 
  
  dTb <- Tb_differential(q_net, M, C)
  Tb[t] <- Tb[t-1] + dTb # N�s armazenamos o novo valor de temperatura no vetor
}

# Representa��o gr�fica
plot(Tb, ylim = c(0,50), type = "l", ylab = "Temperatura corporal (�C)", xlab  ="Tempo", lwd = 2, col = "black")

 
#### Qual � a temperatura de equil�brio (Temperatura operativa) #### 
Te_estim <- function(hc, Ta, epsilon, boltzmann, a, D, s) {(hc*Ta+epsilon*boltzmann*Ta^4+a*(D+s)) / (hc + epsilon*boltzmann*Ta^3)}

Te <- Te_estim(hc, Ta, epsilon, boltzmann, a, D, s)
abline(h = Te, col = "red")


##############################################
## UPSCALING Te INTO GEOGRAPHICAL GRADIENTS ##
##############################################

##################### DADOS AMBIENTAIS ############################

#install.packages("raster"); install.packages("rgdal")
require(raster)
require(rgdal)

# Vamos carregar um mapa de refer�ncia
# ...primeiro!! -> movemos o file "Environmental data" para este diret�rio:
getwd()

shapeTotal <- raster("Environmental data/bio1.bil")
world <- aggregate(shapeTotal, fact=4/res(shapeTotal)[1])
regions0 <- shapefile("Environmental data/Shapes/newRealms.shp")
regions <- aggregate(rbind(regions0))

world <- mask(world, regions)

xy <- xyFromCell(world, 1:ncell(world))
cells <- cellFromPolygon(world, regions)[[1]]
xy.values <- xy[cells,]

# Selecionamos um m�s do ano 
mes = 4
load(file = paste0("Environmental data/data_sun_", mes, ".RData"))

## Temos um mapa com 935 c�lulas, cada uma delas cont�m um banco de dados:
# TIME: hora do dia; RAD: radia��o solar (direta + difusa), AIRT: Temp do ar,
# SOILT: Temp do solo, RH: umidade rel, v: velocidade do vento

plot(world); lines(regions)

# Visualizamos uma das c�lulas
points(x=xy.values[760,1], y=xy.values[760,2], pch = 20, col = "red")
head(data_sun[[760]])

##################### DADOS DAS ESP�CIES EM STUDO ############################

#### Par�metros do modelo para Anolis Auratus

# Carregamos a distribui��o de A. auratus
load(file="Environmental data/distAauratus.RData")
lines(distrib, col = "blue")

M = 5                   # Massa corporal (g)
A = 1e-4 * 10 * M^(2/3) # Superficie (m^2)
l = A^(1/2)             # comprimento do corpo (m)
C = 0.3                 # Capacidade de calor (W �C^-1)

a = 0.4                 # Aabsorb�ncia (cor da pele)
epsilon = 0.9           # Emissividade IR
boltzmann = 5.67e-8     # Constante Stefan-Botzmann (W m^-2 �C^-4)

hc_estim = function(v, l) 6.77 * v^0.6 * l^-0.4 # Coeficiente de convec��o

# Fun��o temperatura operativa
Te_estim <- function(hc, Ta, epsilon, boltzmann, a, D, s) {(hc*Ta+epsilon*boltzmann*Ta^4+a*(D+s)) / (hc + epsilon*boltzmann*Ta^3)}

##################### Proje��o para a �rea de distribui��o ############################

points_auratus <- match(cellFromPolygon(world,distrib)[[1]],cells) # c�lulas com presen�a de A. auratus

Te_Aauratus <- array(NA, dim=c(24,length(points_auratus))) # matriz para salvar os resultados
for(i in 1:length(points_auratus)){ # para cada uma dessas c�lulas
  data_environment <- data_sun[[points_auratus[i]]] # extrair os dados ambientais
  hc <- hc_estim(v=data_environment$V, l) # estimar o coeficiente de convec��o
  Te_Aauratus[,i] <- Te_estim(hc, data_environment$AIRT, epsilon, boltzmann, a, data_environment$RAD, s=0) # estimar temperatura operativa
}

head(Te_Aauratus) # linhas: hora do dia; colunas: c�lula do �rea da distribu��o
meanTe_Aauratus <- rowMeans(Te_Aauratus) # Temperatura m�dia ao longo do dia em toda a �rea de distribui��o

#meanTe_Aauratus_sol <- meanTe_Aauratus

# Temp operativa
plot(meanTe_Aauratus ~ data_environment$TIME, type = "o", col = "black", pch = 20, ylab = "Temperatura, �C", xlab = "Tempo, min")
# Temp do ar
points(data_environment$AIRT ~ data_environment$TIME, type = "o")


##################### proje��o global ############################

mapa_Te <- world
j=0
for(i in cells){ # para cada c�lula no mapa
  j=j+1
  data_environment <- data_sun[[j]]  # extrair os dados ambientais
  hc <- hc_estim(v=data_environment$V, l) # estimar o coeficiente de convec��o
  Te <- Te_estim(hc, data_environment$AIRT, epsilon, boltzmann, a, data_environment$RAD, s=0) 
  
  Te_max <- max(Te, na.rm=T) # selecionamos a vari�vel de interesse
  mapa_Te[i] <- Te_max
}

#install.packages("RColorBrewer")
require(RColorBrewer)
ramp <- seq(-50, 100, length.out = 20)
plot(mapa_Te, legend = F, breaks=ramp, col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(20))
lines(regions)



############### How to make a GIF in R ##################
require(RColorBrewer)
require(purrr)

for(i in 1:12){
  if (i < 10) {name = paste('000',i,'plot.png',sep='')}
  if (i < 100 && i >= 10) {name = paste('00',i,'plot.png', sep='')}
  if (i < 1000 && i >= 100) {name = paste('0', i,'plot.png', sep='')}
  if (i >= 1000) {name = paste(i,'plot.png', sep='')}
  png(file=name, width=400, height=400)
  
  mes = i
  load(file = paste0("Environmental data/data_sun_", mes, ".RData"))
  
  mapa_Te <- world
  j=0
  for(cell in cells){ # para cada c�lula no mapa
    j=j+1
    data_environment <- data_sun[[j]]  
    hc <- hc_estim(v=data_environment$V, l) 
    Te <- Te_estim(hc, data_environment$AIRT, epsilon, boltzmann, a, data_environment$RAD, s=0) 
    
    Te_max <- max(Te, na.rm=T) 
    mapa_Te[cell] <- Te_max
  }
  ramp <- seq(-50, 100, length.out = 20)
  plot(mapa_Te, legend = F, breaks=ramp, col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(20))
  lines(regions)
  dev.off()
}

# pegue todos os arquivos * .png no file Environmental data

list.files(path = "./Environmental data", pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=2) %>% # animates, can opt for number of loops
  image_write("map_seasonal.gif")



