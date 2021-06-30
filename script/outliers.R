### Trabalho: Red Wine #################################################################
# Aluno: Manuel Ferreira Junior
# Disciplina: Regressão II
#setwd('/RobustRegression')
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

#plot(y~x)

# Ajustando modelo OLS
model.rd.o = lm(y~ x)
n=nrow(data)
p = model.rd.o$rank
s = summary(model.rd.o)
ard = ls.diag(model.rd.o)
std.res = ard$std.res

# Pontos de Alavanca e Influentes
# Visualizacao Grafica

#Pontos aberrantes
plot(std.res)
abline(h=-2, col = 2); abline(h=2, col=2)
abe = which(abs(std.res)>2);abe #idenifica no grafico as observacoes aberrantes
x.rd.no = x[-abe];y.rd.no = y[-abe]

# Free vs Total -----------------------------------------------------------
x = data$free.sulfur.dioxide
y = data$total.sulfur.dioxide
# Ajustando modelo OLS
model.ft.o = lm(y ~ x)
n=nrow(data)
p = model.ft.o$rank
s = summary(model.ft.o)
ard = ls.diag(model.ft.o)
std.res = ard$std.res

# Pontos de Alavanca e Influentes
# Visualizacao Grafica

#Pontos aberrantes
plot(std.res)
abline(h=-2, col = 2); abline(h=2, col=2)
abe = which(abs(std.res)>2);abe #idenifica no grafico as observacoes aberrantes
x.ft.no = x[-abe];y.ft.no = y[-abe]