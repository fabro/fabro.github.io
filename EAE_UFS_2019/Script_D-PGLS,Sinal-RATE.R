#Instalamos os pacotes e abrimos eles
install.packages("geomorph")
install.packages("ape") 
library("geomorph")
library("ape")

#Abrimos nossa árvore
arvore<-read.tree("Arvore_caranx.tree")

#Abrimos nossas coordenadas em 2D. Quando nos pergunte se queremos considerar os valores negativos como dados faltantes colocamos "n".
ladnmarks<-readland.tps("shapemean.tps", specID = c("ID"))

#Realizamos a superposição de procrustes
shape_GPA<- gpagen(ladnmarks)

#Abrimos co-variáveis (características das especies)
traits<-read.csv("nicho.csv", header=T, sep=";")

#Juntamos a árvore, as coordenadas superpostas e as características. Esta função ajuda a não dar erros posteriores.
gdf <- geomorph.data.frame(shape_GPA, phy = arvore,traits)

#Fazemos um D-PGLS para ver se a forma (coords) é afetada pela velocidade da agua e a produtividade

caranx.pgls <- procD.pgls(coords ~ velocidade+produtividade, phy = phy, data = gdf, iter = 999)
summary(caranx.pgls)

#Testamos a sinal filogenética da forma a partir do K-mult
sinal_shape_kmult<- physignal(gdf$coords,phy=arvore,iter=999)
sinal_shape_kmult

#Agora vamos a estimar as taxas de diversificação
#Instalamos e abrimos o pacote
install.packages("RRphylo")
library("RRphylo")

#Fazemos uma análise de Componentes Principais a partir das coordenadas para caracterizar a forma das especies
PCA<-plotTangentSpace(shape_GPA$coords, label=T)
Forma_PC<-PCA$pc.scores
Dados_diver<-as.matrix(Forma_PC)

#Realizamos a Phylogenetic Ridge Regression
Taxas<-RRphylo(tree=arvore,y=Dados_diver)

#ploteamos os shift encontrados ao longo da árvore. Temos que indicar a pasta onde queremos que salve a figura
search.shift(Taxas,auto.recognize="yes",test.single= "no",
status.type= "clade",foldername="C:/Users/saraiva/Documents/UFS/Workshop_Macroecologia/Pratica/Morfometria/")
