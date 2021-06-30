### Trabalho: Red Wine #################################################################
# Aluno: Manuel Ferreira Junior e Marcos Antonio Bezerra
# Disciplina: Regressão II
#setwd('/home/manuel/Área de Trabalho/git/UFPB/REG2/Prova03/RobustRegression')
library(robustbase)
library(MASS)
library(sm)
library(quantreg)

## Obs:
# Rodar uma OLS e verificar residuos
# padronizado -3,3 
# Reajustar o OLS e outros modelos
# sem outliers
View(data)
# ----------------------------------------------
## functions required - ETKRR method ---------------

sigma2est1 = function(y1, y2, frac = .5) {
  n = length(y1)
  m = floor(n*frac)
  idx1 = sample(1:n, m, replace = T)
  idx2 = sample(1:n, m, replace = T)
  tmp = (y1[idx1] - y2[idx2])^2
  mean(quantile(tmp[tmp != 0], probs = c(.9, .1)))
}

sigma2est2 = function(y1, y2) {
  n = length(y1)
  tmp = c()
  for (i in 1:n) {
    for (j in 1:n) {
      tmp = c(tmp, (y1[i] - y2[j])^2)
    }
  }
  tmp = tmp[tmp != 0]
  median(tmp)
}

###### Gaussian kernel

gauss.kern = function(a, b, s)
{
  as.vector(exp(-(1/s)*(a-b)^2)) #Gaussian
}

#### ETKRR METHOD - According to the hyper-parameter

kernel.reg1 = function(x, y, tol = 1e-10, maxit = 100)
{
  x = as.matrix(x)
  n = nrow(x)
  # Initialization
  x = cbind(1, x)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  
  s2 = sigma2est1(y, yhat)
  
  K = gauss.kern(y, yhat, s2)
  S = sum(2-2*K)
  it = 1
  
  # Model Step
  repeat {
    it = it+1
    betahat = solve(t(x)%*%diag(K)%*%x)%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, s2)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coef = as.vector(betahat), fitted = as.vector(yhat), criterion = S, weigth = K))
}

kernel.reg2 = function(x, y, tol = 1e-10, maxit = 100)
{
  x = as.matrix(x)
  n = nrow(x)
  # Initialization
  x = cbind(1, x)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  
  s2 = sigma2est2(y, yhat)
  
  K = gauss.kern(y, yhat, s2)
  S = sum(2-2*K)
  it = 1
  
  # Model Step
  repeat {
    it = it+1
    betahat = solve(t(x)%*%diag(K)%*%x)%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, s2)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coef = as.vector(betahat), fitted = as.vector(yhat), criterion = S, weigth = K))
}

kernel.reg3 = function(x, y, tol = 1e-10, maxit = 100)
{
  x = as.matrix(x)
  n = nrow(x)
  # Initialization
  x = cbind(1, x)
  betahat = solve(t(x)%*%x)%*%t(x)%*%y
  yhat = x%*%betahat
  
  s2 = sum((y-yhat)^2)/(n-ncol(x))
  
  K = gauss.kern(y, yhat, s2)
  S = sum(2-2*K)
  it = 1
  
  # Model Step
  repeat {
    it = it+1
    betahat = solve(t(x)%*%diag(K)%*%x)%*%t(x)%*%diag(K)%*%y
    yhat = x%*%betahat
    K = gauss.kern(y, yhat, s2)
    S = c(S, sum(2-2*K))
    if (abs(S[it]-S[(it-1)]) <= tol || it >= maxit) break
  }
  (result = list(coef = as.vector(betahat), fitted = as.vector(yhat), criterion = S, weigth = K))
}

# Carregando dados --------------

# Carregando os dados
url = 'https://raw.githubusercontent.com/Manuelfjr/RobustRegression/main/data/winequality-red.csv'
data = read.csv(url)
#View(data)

# Analisando os dados
str(data)
summary(data)

############ OLS ----------------------------------------------
# Residual vs density -----------------------------------------------------
x = data$residual.sugar
y = data$density
#x = x.rd.no; y = y.rd.no
# Plotando gráficos
plot(x,y, xlab='Residual Sugar', ylab='Density')
plot(x.rd.no,y.rd.no, xlab='Residual Sugar', ylab='Density')
## Modelagem ---------------------------------------------------------------
# OLS
mod.ols = lm(y~x)
mod.ols.rd.no = lm(y.rd.no~x.rd.no)
fit.ols = fitted(mod.ols)
fit.ols.rd.no = fitted(mod.ols.rd.no)

# WLS
xmod = cbind(1,x)
sigma2 = (summary(lm(y~x))$sigma)^2
peso = 1/diag(sigma2*xmod%*%(solve(t(xmod)%*%xmod))%*%t(xmod))

## minimos quadrados ponderados
mod.wls = lm(y ~ x, weights = peso)
fit.wls = fitted (mod.wls)

## WLS - sem outlier
xmod = cbind(1,x.rd.no)
sigma2 = (summary(lm(y.rd.no~x.rd.no))$sigma)^2
peso = 1/diag(sigma2*xmod%*%(solve(t(xmod)%*%xmod))%*%t(xmod))

### minimos quadrados ponderados
mod.wls.rd.no = lm(y.rd.no ~ x.rd.no, weights = peso)
fit.wls.rd.no = fitted (mod.wls.rd.no)


# M-Estimator
mod.M = rlm(y~x)
fit.M = fitted(mod.M)

mod.M.rd.no = rlm(y.rd.no ~ x.rd.no)
fit.M.rd.no = fitted(mod.M.rd.no)

# MM-Estimator
mod.rob = lmrob(y ~ x)
fit.rob = fitted(mod.rob)

mod.rob.rd.no = lmrob(y.rd.no ~ x.rd.no)
fit.rob.rd.no = fitted(mod.rob.rd.no)

# Regressao L1 (mediana)
mod.l1 = rq(y~ x, tau = 0.5)
fit.l1 = fitted(mod.l1)

mod.l1.rd.no = rq(y.rd.no ~ x.rd.no, tau = 0.5)
fit.l1.rd.no = fitted(mod.l1.rd.no)

# ETKRR
mod.etkrrs1 = kernel.reg1(x, y)
#mod.etkrrs2 = kernel.reg2(x, y)
mod.etkrrs3 = kernel.reg3(x, y)

mod.etkrrs1.rd.no = kernel.reg1(x.rd.no, y.rd.no)
#mod.etkrrs2.rd.no = kernel.reg2(x.rd.no, y.rd.no)
mod.etkrrs3.rd.no = kernel.reg3(x.rd.no, y.rd.no)

dim(data)

##  Vizualiando e comparando os modelos (Com outlier)------------------------------------
colors.l = c("black",'red','green','deepskyblue3','darkorange', 'deeppink')
models = c("OLS", "WLS", "M",'MM', 'L1', 'ETKRRS1')
#par(mfrow=c(1,2))
plot(x,y, main = "Red Wine Data (Com outlier)",
     xlab='Residual Sugar',
     ylab='Density')#, ylim=c(0, 30))
lines(x, fit.ols, col = colors.l[1])
lines(x, fit.wls, col = colors.l[2])
lines(x, fit.M, col = colors.l[3])
lines(x, fit.rob, col = colors.l[4])
lines(x, fit.l1, col = colors.l[5])
lines(x, mod.etkrrs1$fitted, col = colors.l[6])
grid()
legend(13, 0.9935, legend=models,
              col=colors.l, lty=c(1,1,1), cex=0.5, text.font=2, bg='transparent')


##  Vizualiando e comparando os modelos (Sem outlier)------------------------------------
colors.l = c("black",'red','green','deepskyblue3','darkorange', 'deeppink')
models = c("OLS", "WLS", "M",'MM', 'L1', 'ETKRRS1')
plot(x.rd.no,y.rd.no, main = "Red Wine Data (Sem outlier)",
     xlab='Residual Sugar',
     ylab='Density')#, ylim=c(0, 30))
lines(x.rd.no, fit.ols.rd.no, col = colors.l[1])
lines(x.rd.no, fit.wls.rd.no, col = colors.l[2])
lines(x.rd.no, fit.M.rd.no, col = colors.l[3])
lines(x.rd.no, fit.rob.rd.no, col = colors.l[4])
lines(x.rd.no, fit.l1.rd.no, col = colors.l[5])
lines(x.rd.no, mod.etkrrs1.rd.no$fitted, col = colors.l[6])
grid()
legend(13, 0.9956, legend=models,
       col=colors.l, lty=c(1,1,1), cex=0.5, text.font=2, bg='transparent')



## DIFBETAS ----------------
diffbetas <- function(b,bo){
  abs((b - bo)/b)*100
}

# OLS
diffbetas.ols = diffbetas(coef(mod.ols), coef(mod.ols.rd.no))

# WLS
diffbetas.wls = diffbetas(coef(mod.wls), coef(mod.wls.rd.no))

# M-Estimator
diffbetas.M = diffbetas(coef(mod.M), coef(mod.M.rd.no))

# MM-Estimator
diffbetas.rob = diffbetas(coef(mod.rob), coef(mod.rob.rd.no))

# Regressao L1 (mediana)
diffbetas.l1 = diffbetas(coef(mod.l1),coef(mod.l1.rd.no))

# ETKRR1
ettkr1.betas.o = c('(Intercept)' = mod.etkrrs1$coef[1],'x' = mod.etkrrs1$coef[2])
ettkr1.betas.no = c('(Intercept)' = mod.etkrrs1.rd.no$coef[1],'x' = mod.etkrrs1.rd.no$coef[2])
diffbetas.ettkr1 = diffbetas(mod.etkrrs1$coef,mod.etkrrs1.rd.no$coef)

# Construindo tabela
diffsbetas=  data.frame(diffbetas.ols,
                        diffbetas.wls,
                        diffbetas.M,
                        diffbetas.rob,
                        diffbetas.l1,
                        diffbetas.ettkr1)
colnames(diffsbetas) = models
diffsbetas


# ---------------------------------------------
#  Free vs Total -----------------------------------------------------
x = data$free.sulfur.dioxide
y = data$total.sulfur.dioxide
#x = x.rd.no; y = y.rd.no
# Plotando gráficos
plot(x,y, xlab='Free Sulfur Dioxide', ylab='Total Sulfur Dioxide')
plot(x.ft.no,y.ft.no,xlab='Free Sulfur Dioxide', ylab='Total Sulfur Dioxide')
## Modelagem ---------------------------------------------------------------
# OLS
mod.ols = lm(y~x)
mod.ols.ft.no = lm(y.ft.no~x.ft.no)

fit.ols = fitted(mod.ols)
fit.ols.ft.no = fitted(mod.ols.ft.no)

# WLS
xmod = cbind(1,x)
sigma2 = (summary(lm(y~x))$sigma)^2
peso = 1/diag(sigma2*xmod%*%(solve(t(xmod)%*%xmod))%*%t(xmod))

## minimos quadrados ponderados
mod.wls = lm(y ~ x, weights = peso)
fit.wls = fitted (mod.wls)

## WLS - sem outlier
xmod = cbind(1,x.ft.no)
sigma2 = (summary(lm(y.ft.no~x.ft.no))$sigma)^2
peso = 1/diag(sigma2*xmod%*%(solve(t(xmod)%*%xmod))%*%t(xmod))

### minimos quadrados ponderados
mod.wls.ft.no = lm(y.ft.no ~ x.ft.no, weights = peso)
fit.wls.ft.no = fitted (mod.wls.ft.no)


# M-Estimator
mod.M = rlm(y~x)
fit.M = fitted(mod.M)

mod.M.ft.no = rlm(y.ft.no ~ x.ft.no)
fit.M.ft.no = fitted(mod.M.ft.no)

# MM-Estimator
mod.rob = lmrob(y ~ x)
fit.rob = fitted(mod.rob)

mod.rob.ft.no = lmrob(y.ft.no ~ x.ft.no)
fit.rob.ft.no = fitted(mod.rob.ft.no)

# Regressao L1 (mediana)
mod.l1 = rq(y~ x, tau = 0.5)
fit.l1 = fitted(mod.l1)

mod.l1.ft.no = rq(y.ft.no ~ x.ft.no, tau = 0.5)
fit.l1.ft.no = fitted(mod.l1.ft.no)

# ETKRR
mod.etkrrs1 = kernel.reg1(x, y)
#mod.etkrrs2 = kernel.reg2(x, y)
mod.etkrrs3 = kernel.reg3(x, y)

mod.etkrrs1.ft.no = kernel.reg1(x.ft.no, y.ft.no)
#mod.etkrrs2.ft.no = kernel.reg2(x.ft.no, y.ft.no)
mod.etkrrs3.ft.no = kernel.reg3(x.ft.no, y.ft.no)

#plot(x,y)
#lines(x, mod.etkrrs1$fitted, col = 2)
#lines(x, mod.etkrrs2$fitted, col = 3)
#lines(x, mod.etkrrs3$fitted, col = 4) ## vamos escolher o s3, por menor custo e facil intepretacao
dim(data)

##  Vizualiando e comparando os modelos (Com outlier)------------------------------------
colors.l = c("black",'red','green','deepskyblue3','darkorange', 'deeppink')
models = c("OLS", "WLS", "M",'MM', 'L1', 'ETKRRS1')
#par(mfrow=c(1,2))
plot(x,y, main = "Red Wine Data (Com outlier)",
     xlab='Free Sulfur Dioxide',
     ylab='Total Sulfur Dioxide')#, ylim=c(0, 30))
lines(x, fit.ols, col = colors.l[1])
lines(x, fit.wls, col = colors.l[2])
lines(x, fit.M, col = colors.l[3])
lines(x, fit.rob, col = colors.l[4])
lines(x, fit.l1, col = colors.l[5])
lines(x, mod.etkrrs1$fitted, col = colors.l[6])
grid()
legend(0, 280, legend=models,
       col=colors.l, lty=c(1,1,1), cex=0.5, text.font=2, bg='transparent')


##  Vizualiando e comparando os modelos (Sem outlier)------------------------------------
colors.l = c("black",'red','green','deepskyblue3','darkorange', 'deeppink')
models = c("OLS", "WLS", "M",'MM', 'L1', 'ETKRRS1')
plot(x.ft.no,y.ft.no, main = "Red Wine Data (Sem outlier)",
     xlab='Free Sulfur Dioxide',ylab='Total Sulfur Dioxide')#,
     #ylab='Total Sulfur Dioxide')#, ylim=c(0, 30))
lines(x.ft.no, fit.ols.ft.no, col = colors.l[1])
lines(x.ft.no, fit.wls.ft.no, col = colors.l[2])
lines(x.ft.no, fit.M.ft.no, col = colors.l[3])
lines(x.ft.no, fit.rob.ft.no, col = colors.l[4])
lines(x.ft.no, fit.l1.ft.no, col = colors.l[5])
lines(x.ft.no, mod.etkrrs1.ft.no$fitted, col = colors.l[6])
grid()
legend(53,47, legend=models,
       col=colors.l, lty=c(1,1,1), cex=0.5, text.font=2, bg='transparent')
## DIFBETAS ----------------
diffbetas <- function(b,bo){
  abs((b - bo)/b)*100
}

# OLS
diffbetas.ols = diffbetas(coef(mod.ols), coef(mod.ols.ft.no))

# WLS
diffbetas.wls = diffbetas(coef(mod.wls), coef(mod.wls.ft.no))

# M-Estimator
diffbetas.M = diffbetas(coef(mod.M), coef(mod.M.ft.no))

# MM-Estimator
diffbetas.rob = diffbetas(coef(mod.rob), coef(mod.rob.ft.no))

# Regressao L1 (mediana)
diffbetas.l1 = diffbetas(coef(mod.l1),coef(mod.l1.ft.no))

# ETKRR1
ettkr1.betas.o = c('(Intercept)' = mod.etkrrs1$coef[1],'x' = mod.etkrrs1$coef[2])
ettkr1.betas.no = c('(Intercept)' = mod.etkrrs1.ft.no$coef[1],'x' = mod.etkrrs1.ft.no$coef[2])
diffbetas.ettkr1 = diffbetas(mod.etkrrs1$coef,mod.etkrrs1.ft.no$coef)

# Construindo tabela
diffsbetas=  data.frame(diffbetas.ols,
                           diffbetas.wls,
                           diffbetas.M,
                           diffbetas.rob,
                           diffbetas.l1,
                           diffbetas.ettkr1)
colnames(diffsbetas) = models
diffsbetas


############# ROSNER ---------------------------------------------
# Residual vs density -----------------------------------------------------
x = data$residual.sugar
y = data$density
#x = x.rd.no; y = y.rd.no
# Plotando gráficos
plot(x,y, xlab='Residual Sugar', ylab='Density')
plot(x.rd.no,y.rd.no, xlab='Residual Sugar', ylab='Density')
## Modelagem ---------------------------------------------------------------
# OLS
mod.ols = lm(y~x)
mod.ols.rd.no = lm(y.rd.no~x.rd.no)
fit.ols = fitted(mod.ols)
fit.ols.rd.no = fitted(mod.ols.rd.no)

# WLS
xmod = cbind(1,x)
sigma2 = (summary(lm(y~x))$sigma)^2
peso = 1/diag(sigma2*xmod%*%(solve(t(xmod)%*%xmod))%*%t(xmod))

## minimos quadrados ponderados
mod.wls = lm(y ~ x, weights = peso)
fit.wls = fitted (mod.wls)

## WLS - sem outlier
xmod = cbind(1,x.rd.no)
sigma2 = (summary(lm(y.rd.no~x.rd.no))$sigma)^2
peso = 1/diag(sigma2*xmod%*%(solve(t(xmod)%*%xmod))%*%t(xmod))

### minimos quadrados ponderados
mod.wls.rd.no = lm(y.rd.no ~ x.rd.no, weights = peso)
fit.wls.rd.no = fitted (mod.wls.rd.no)


# M-Estimator
mod.M = rlm(y~x)
fit.M = fitted(mod.M)

mod.M.rd.no = rlm(y.rd.no ~ x.rd.no)
fit.M.rd.no = fitted(mod.M.rd.no)

# MM-Estimator
mod.rob = lmrob(y ~ x)
fit.rob = fitted(mod.rob)

mod.rob.rd.no = lmrob(y.rd.no ~ x.rd.no)
fit.rob.rd.no = fitted(mod.rob.rd.no)

# Regressao L1 (mediana)
mod.l1 = rq(y~ x, tau = 0.5)
fit.l1 = fitted(mod.l1)

mod.l1.rd.no = rq(y.rd.no ~ x.rd.no, tau = 0.5)
fit.l1.rd.no = fitted(mod.l1.rd.no)

# ETKRR
mod.etkrrs1 = kernel.reg1(x, y)
#mod.etkrrs2 = kernel.reg2(x, y)
mod.etkrrs3 = kernel.reg3(x, y)

mod.etkrrs1.rd.no = kernel.reg1(x.rd.no, y.rd.no)
#mod.etkrrs2.rd.no = kernel.reg2(x.rd.no, y.rd.no)
mod.etkrrs3.rd.no = kernel.reg3(x.rd.no, y.rd.no)

#plot(x,y)
#lines(x, mod.etkrrs1$fitted, col = 2)
#lines(x, mod.etkrrs2$fitted, col = 3)
#lines(x, mod.etkrrs3$fitted, col = 4) ## vamos escolher o s3, por menor custo e facil intepretacao
dim(data)

##  Vizualiando e comparando os modelos (Com outlier)------------------------------------
colors.l = c("black",'red','green','deepskyblue3','darkorange', 'deeppink')
models = c("OLS", "WLS", "M",'MM', 'L1', 'ETKRRS1')
#par(mfrow=c(1,2))
plot(x,y, main = "Red Wine Data (Com outlier)",
     xlab='Residual Sugar',
     ylab='Density')#, ylim=c(0, 30))
lines(x, fit.ols, col = colors.l[1])
lines(x, fit.wls, col = colors.l[2])
lines(x, fit.M, col = colors.l[3])
lines(x, fit.rob, col = colors.l[4])
lines(x, fit.l1, col = colors.l[5])
lines(x, mod.etkrrs1$fitted, col = colors.l[6])
grid()
legend(13, 0.9935, legend=models,
       col=colors.l, lty=c(1,1,1), cex=0.5, text.font=2, bg='transparent')


##  Vizualiando e comparando os modelos (Sem outlier)------------------------------------
colors.l = c("black",'red','green','deepskyblue3','darkorange', 'deeppink')
models = c("OLS", "WLS", "M",'MM', 'L1', 'ETKRRS1')
plot(x.rd.no,y.rd.no, main = "Red Wine Data (Sem outlier)",
     xlab='Residual Sugar',
     ylab='Density')#, ylim=c(0, 30))
lines(x.rd.no, fit.ols.rd.no, col = colors.l[1])
lines(x.rd.no, fit.wls.rd.no, col = colors.l[2])
lines(x.rd.no, fit.M.rd.no, col = colors.l[3])
lines(x.rd.no, fit.rob.rd.no, col = colors.l[4])
lines(x.rd.no, fit.l1.rd.no, col = colors.l[5])
lines(x.rd.no, mod.etkrrs1.rd.no$fitted, col = colors.l[6])
grid()
legend(4, 0.994, legend=models,
       col=colors.l, lty=c(1,1,1), cex=0.5, text.font=2, bg='transparent')



## DIFBETAS ----------------
diffbetas <- function(b,bo){
  abs((b - bo)/b)*100
}

# OLS
diffbetas.ols = diffbetas(coef(mod.ols), coef(mod.ols.rd.no))

# WLS
diffbetas.wls = diffbetas(coef(mod.wls), coef(mod.wls.rd.no))

# M-Estimator
diffbetas.M = diffbetas(coef(mod.M), coef(mod.M.rd.no))

# MM-Estimator
diffbetas.rob = diffbetas(coef(mod.rob), coef(mod.rob.rd.no))

# Regressao L1 (mediana)
diffbetas.l1 = diffbetas(coef(mod.l1),coef(mod.l1.rd.no))

# ETKRR1
ettkr1.betas.o = c('(Intercept)' = mod.etkrrs1$coef[1],'x' = mod.etkrrs1$coef[2])
ettkr1.betas.no = c('(Intercept)' = mod.etkrrs1.rd.no$coef[1],'x' = mod.etkrrs1.rd.no$coef[2])
diffbetas.ettkr1 = diffbetas(mod.etkrrs1$coef,mod.etkrrs1.rd.no$coef)

# Construindo tabela
diffsbetas=  data.frame(diffbetas.ols,
                        diffbetas.wls,
                        diffbetas.M,
                        diffbetas.rob,
                        diffbetas.l1,
                        diffbetas.ettkr1)
colnames(diffsbetas) = models
diffsbetas


# ---------------------------------------------
#  Free vs Total -----------------------------------------------------
x = data$free.sulfur.dioxide
y = data$total.sulfur.dioxide
#x = x.rd.no; y = y.rd.no
# Plotando gráficos
plot(x,y, xlab='Free Sulfur Dioxide', ylab='Total Sulfur Dioxide')
plot(x.ft.no,y.ft.no,xlab='Free Sulfur Dioxide', ylab='Total Sulfur Dioxide')
## Modelagem ---------------------------------------------------------------
# OLS
mod.ols = lm(y~x)
mod.ols.ft.no = lm(y.ft.no~x.ft.no)

fit.ols = fitted(mod.ols)
fit.ols.ft.no = fitted(mod.ols.ft.no)

# WLS
xmod = cbind(1,x)
sigma2 = (summary(lm(y~x))$sigma)^2
peso = 1/diag(sigma2*xmod%*%(solve(t(xmod)%*%xmod))%*%t(xmod))

## minimos quadrados ponderados
mod.wls = lm(y ~ x, weights = peso)
fit.wls = fitted (mod.wls)

## WLS - sem outlier
xmod = cbind(1,x.ft.no)
sigma2 = (summary(lm(y.ft.no~x.ft.no))$sigma)^2
peso = 1/diag(sigma2*xmod%*%(solve(t(xmod)%*%xmod))%*%t(xmod))

### minimos quadrados ponderados
mod.wls.ft.no = lm(y.ft.no ~ x.ft.no, weights = peso)
fit.wls.ft.no = fitted (mod.wls.ft.no)


# M-Estimator
mod.M = rlm(y~x)
fit.M = fitted(mod.M)

mod.M.ft.no = rlm(y.ft.no ~ x.ft.no)
fit.M.ft.no = fitted(mod.M.ft.no)

# MM-Estimator
mod.rob = lmrob(y ~ x)
fit.rob = fitted(mod.rob)

mod.rob.ft.no = lmrob(y.ft.no ~ x.ft.no)
fit.rob.ft.no = fitted(mod.rob.ft.no)

# Regressao L1 (mediana)
mod.l1 = rq(y~ x, tau = 0.5)
fit.l1 = fitted(mod.l1)

mod.l1.ft.no = rq(y.ft.no ~ x.ft.no, tau = 0.5)
fit.l1.ft.no = fitted(mod.l1.ft.no)

# ETKRR
mod.etkrrs1 = kernel.reg1(x, y)
#mod.etkrrs2 = kernel.reg2(x, y)
mod.etkrrs3 = kernel.reg3(x, y)

mod.etkrrs1.ft.no = kernel.reg1(x.ft.no, y.ft.no)
#mod.etkrrs2.ft.no = kernel.reg2(x.ft.no, y.ft.no)
mod.etkrrs3.ft.no = kernel.reg3(x.ft.no, y.ft.no)

#plot(x,y)
#lines(x, mod.etkrrs1$fitted, col = 2)
#lines(x, mod.etkrrs2$fitted, col = 3)
#lines(x, mod.etkrrs3$fitted, col = 4) ## vamos escolher o s3, por menor custo e facil intepretacao
dim(data)

##  Vizualiando e comparando os modelos (Com outlier)------------------------------------
colors.l = c("black",'red','green','deepskyblue3','darkorange', 'deeppink')
models = c("OLS", "WLS", "M",'MM', 'L1', 'ETKRRS1')
#par(mfrow=c(1,2))
plot(x,y, main = "Red Wine Data (Com outlier)",
     xlab='Free Sulfur Dioxide',
     ylab='Total Sulfur Dioxide')#, ylim=c(0, 30))
lines(x, fit.ols, col = colors.l[1])
lines(x, fit.wls, col = colors.l[2])
lines(x, fit.M, col = colors.l[3])
lines(x, fit.rob, col = colors.l[4])
lines(x, fit.l1, col = colors.l[5])
lines(x, mod.etkrrs1$fitted, col = colors.l[6])
grid()
legend(0, 280, legend=models,
       col=colors.l, lty=c(1,1,1), cex=0.5, text.font=2, bg='transparent')


##  Vizualiando e comparando os modelos (Sem outlier)------------------------------------
colors.l = c("black",'red','green','deepskyblue3','darkorange', 'deeppink')
models = c("OLS", "WLS", "M",'MM', 'L1', 'ETKRRS1')
plot(x.ft.no,y.ft.no, main = "Red Wine Data (Sem outlier)",
     xlab='Free Sulfur Dioxide',ylab='Total Sulfur Dioxide')#,
#ylab='Total Sulfur Dioxide')#, ylim=c(0, 30))
lines(x.ft.no, fit.ols.ft.no, col = colors.l[1])
lines(x.ft.no, fit.wls.ft.no, col = colors.l[2])
lines(x.ft.no, fit.M.ft.no, col = colors.l[3])
lines(x.ft.no, fit.rob.ft.no, col = colors.l[4])
lines(x.ft.no, fit.l1.ft.no, col = colors.l[5])
lines(x.ft.no, mod.etkrrs1.ft.no$fitted, col = colors.l[6])
grid()
legend(53,47, legend=models,
       col=colors.l, lty=c(1,1,1), cex=0.5, text.font=2, bg='transparent')
## DIFBETAS ----------------
diffbetas <- function(b,bo){
  abs((b - bo)/b)*100
}

# OLS
diffbetas.ols = diffbetas(coef(mod.ols), coef(mod.ols.ft.no))

# WLS
diffbetas.wls = diffbetas(coef(mod.wls), coef(mod.wls.ft.no))

# M-Estimator
diffbetas.M = diffbetas(coef(mod.M), coef(mod.M.ft.no))

# MM-Estimator
diffbetas.rob = diffbetas(coef(mod.rob), coef(mod.rob.ft.no))

# Regressao L1 (mediana)
diffbetas.l1 = diffbetas(coef(mod.l1),coef(mod.l1.ft.no))

# ETKRR1
ettkr1.betas.o = c('(Intercept)' = mod.etkrrs1$coef[1],'x' = mod.etkrrs1$coef[2])
ettkr1.betas.no = c('(Intercept)' = mod.etkrrs1.ft.no$coef[1],'x' = mod.etkrrs1.ft.no$coef[2])
diffbetas.ettkr1 = diffbetas(mod.etkrrs1$coef,mod.etkrrs1.ft.no$coef)

# Construindo tabela
diffsbetas=  data.frame(diffbetas.ols,
                        diffbetas.wls,
                        diffbetas.M,
                        diffbetas.rob,
                        diffbetas.l1,
                        diffbetas.ettkr1)
colnames(diffsbetas) = models
diffsbetas
