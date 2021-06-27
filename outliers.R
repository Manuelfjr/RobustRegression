### Trabalho: Red Wine #################################################################
# Aluno: Manuel Ferreira Junior
# Disciplina: Regressão II
#setwd('/home/manuel/Área de Trabalho/git/UFPB/REG2/Prova03/RobustRegression')
library(robustbase)
library(MASS)
library(sm)
library(quantreg)

# Carregando os dados
url = 'https://raw.githubusercontent.com/Manuelfjr/RobustRegression/main/data/winequality-red.csv'
data = read.csv(url)
#View(data)
alpha = 0.05
# Analisando os dados
str(data)
summary(data)

# Residual vs density -----------------------------------------------------
x = data$residual.sugar
y = data$density

plot(y~x)

# Ajustando modelo OLS
model.rd.o = lm(y~ x)
n=nrow(data)
p = model.rd.o$rank
s = summary(model.rd.o)
ard = ls.diag(model.rd.o)
devres = ard$std.res

# Pontos de Alavanca e Influentes
# Visualizacao Grafica

#Pontos aberrantes
plot(devres)
abline(h=-2, col = 2); abline(h=2, col=2)
abe = which(abs(devres)>2);abe #idenifica no grafico as observacoes aberrantes

#Identificacao de pontos de alavanca
plot(ard$hat)
abline(h=2*(p/n), col = 2)
ala = which(ard$hat > 2*(p/n));ala

#Identificacao de pontos de influencia
plot(ard$cook)
abline(h=qchisq(alpha,p)/p, col = 2)
inf = which(ard$cook > qchisq(alpha,p)/p);inf
#which(ard$cook > qchisq(0.05,p)/p) %in% which(ard$hat > 2*(p/n))
x.rd.no = x[-inf];y.rd.no = y[-inf] 
plot(x.rd.no,y.rd.no)

# Modelo sem outlier
model.rd.no = lm(y.rd.no~x.rd.no)
summary(model.rd.o);summary(model.rd.no)
abs(coef(model.rd.no) - coef(model.rd.o))/abs(coef(model.rd.o))


# Free vs Total -----------------------------------------------------------
x = data$free.sulfur.dioxide
y = data$total.sulfur.dioxide
plot(y~x)
# Ajustando modelo OLS
model.ft.o = lm(y~ x)
n=nrow(data)
p = model.ft.o$rank
s = summary(model.ft.o)
ard = ls.diag(model.ft.o)
devres = ard$std.res

# Pontos de Alavanca e Influentes
# Visualizacao Grafica

#Pontos aberrantes
plot(devres)
abline(h=-2, col = 2); abline(h=2, col=2)
abe = which(abs(devres)>2);abe #idenifica no grafico as observacoes aberrantes

#Identificacao de pontos de alavanca
plot(ard$hat)
abline(h=2*(p/n), col = 2)
ala = which(ard$hat > 2*(p/n));ala

#Identificacao de pontos de influencia
plot(ard$cook)
abline(h=qchisq(alpha,p)/p, col = 2)
inf = which(ard$cook > qchisq(alpha,p)/p);inf
#which(ard$cook > qchisq(0.05,p)/p) %in% which(ard$hat > 2*(p/n))
x.ft.no = x[-inf];y.ft.no = y[-inf] 

# Modelo sem outlier
model.ft.no = lm(y.ft.no~x.ft.no)
summary(model.ft.o);summary(model.ft.no)
plot(x.ft.no,y.ft.no)
abline(model.ft.o)
abline(model.ft.no,col='red')
abs(abs(coef(model.ft.no) - coef(model.ft.o))/abs(coef(model.ft.o)))

