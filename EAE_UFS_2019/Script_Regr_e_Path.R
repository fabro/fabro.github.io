#Selecionamos nosso diretorio de trabalho, neste caso a pasta "Path_Analisis_phy"

#Nesta primeira parte vamos ver o efeito de controlar estatisticamente uma variável
#Abrimos nossa tabela de dados
Dados<-read.csv("Dados_analises.csv", header=T, sep=";")

#Observamos o cabeçalho dos dados
head(Dados)
#Gasto público em educação obtido do Banco Macro
#Taxa de homicidios de acordo com o  UNODC 
#https://www.indexmundi.com/g/r.aspx?v=69&l=pt

#Rodamos as regressões lineais simples e múltiples
modelo1<-lm(log10(Taxa_Homi)~Pobreza,Dados)
summary(modelo1)

modelo2<-lm(log10(Taxa_Homi)~log10(Inves_Educa+1),Dados)
summary(modelo2)

modelo3<-lm(log10(Taxa_Homi)~log10(Inves_Educa+1)+Pobreza,Dados)
summary(modelo3)


#Path Analysis
#instalamos o pacote phylopath
install.packages("phylopath")
#abrimos o pacote
library(phylopath)

#neste exemplo vamos utilizar dados hipotéticos de 100 especies de rinocerontes
#o objeto rhino são nossos dados com os traits, LS=Tamanho da ninhada, BM= Massa corporal, DD=Potencial de dispersão, NL=Tamanho da nariz e RS=Área Geográfica
rhino

#o objeto rhino_tree é nossa árvore filogenética, plotamos ela
plot(rhino_tree)

#Criamos os possiveis modelos competitivos
candidates <- list(
A = DAG(LS ~ BM, NL ~ BM, DD ~ BM, RS ~ BM + NL+ DD),
B = DAG(LS ~ BM, NL ~ BM,DD~BM, RS~DD),
C = DAG(LS ~ BM, NL ~ BM, DD ~ BM+NL, RS ~ BM + NL+ DD)
)

#Plotamos as hipóteses causais
plot_model_set(candidates)

#Avaliamos os modelos causais, se o valor de p é <0.05 o modelo não é confiável. Depois selecionamos pelo CICc (similar a AIC)
p <- phylo_path(candidates, rhino, rhino_tree)
summary(p)

#Do melhor modelo fazemos o path analysis
d_fitted <- est_DAG(candidates$C, rhino, rhino_tree, 'lambda',boot=100)

#Fazemos o plot do modelo causal
plot(d_fitted)

#Observamos os coeficientes estandarizados e seus intervalos de confiança do 95%
coef_plot(d_fitted, error_bar="ci")



