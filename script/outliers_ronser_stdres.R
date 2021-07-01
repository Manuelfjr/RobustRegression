### Trabalho: Red Wine #################################################################
# Aluno: Manuel Ferreira Junior e Marcos Antonio Bezerra
# Disciplina: Regressão II
#setwd('/RobustRegression')

# Pacotes
#install.packages("EnvStats")
library(EnvStats)
url = 'https://raw.githubusercontent.com/Manuelfjr/RobustRegression/main/data/winequality-red.csv'
data = read.csv(url)

# Residual vs density -----------------------------------------------------
x = data$residual.sugar
y = data$density

model.rd.o = lm(y~ x)
ard = ls.diag(model.rd.o)
std.res = ard$std.res

dim(data)
plot(x,y)
rd.rtest = rosnerTest(abs(std.res), k = 1400)
rd.rtest$n.outliers # Quantidade de outliers (TRUE)
non.out = rd.rtest$all.stats[rd.rtest$all.stats$Outlier,]$Obs.Num
x.rd.no = x[-non.out]; y.rd.no = y[-non.out]
par(mfrow=c(1,2))
plot(x, y, main='Residual vs Density',
     xlab = 'Residual Sugar',
     ylab='Density')
plot(x.rd.no, y.rd.no, main='Residual vs Density (Sem outlier)',
     xlab = 'Residual Sugar',
     ylab='')

# Free vs total -----------------------------------------------------------
x = data$free.sulfur.dioxide
y = data$total.sulfur.dioxide

model.rd.o = lm(y~ x)
ard = ls.diag(model.rd.o)
std.res = ard$std.res

dim(data)
plot(x,y)
rd.rtest = rosnerTest(abs(std.res), k = 1400)
ft.rtest$n.outliers # Quantidade de outliers (TRUE)
non.out = ft.rtest$all.stats[ft.rtest$all.stats$Outlier,]$Obs.Num
x.ft.no = x[-non.out]; y.ft.no = y[-non.out]
par(mfrow=c(1,2))
plot(x, y, main='Residual vs Density',
     xlab = 'Residual Sugar',
     ylab='Density')
plot(x.ft.no, y.ft.no, main='Residual vs Density (Sem outlier)',
     xlab = 'Residual Sugar',
     ylab='')

which(max(y) == y)

