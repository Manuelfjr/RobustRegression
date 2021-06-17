### Trabalho: Red Wine #################################################################
# Aluno: Manuel Ferreira Junior
# Disciplina: Regressão II
#setwd('/home/manuel/Área de Trabalho/git/UFPB/REG2/Prova03/RobustRegression')

# Carregando os dados
url = file.choose()
data = read.csv(url)
View(data)

# Analisando os dados
str(data)
summary(data)

x = data$free.sulfur.dioxide
y = data$total.sulfur.dioxide
#x = data$residual.sugar
#y = data$density

# Plotando gráficos
plot(x,y, xlab='Free Sulfur Dioxide', ylab='Total Sulfur Dioxide')


# Modelagem ---------------------------------------------------------------
## functions required - ETKRR method

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

# OLS
mod.ols = lm(y~x)
fit.ols = fitted(mod.ols)

# WLS
xmod = cbind(1,x)
sigma2 = (summary(lm(y~x))$sigma)^2
peso = 1/diag(sigma2*xmod%*%(solve(t(xmod)%*%xmod))%*%t(xmod))

# minimos quadrados ponderados
mod.wls = lm(y ~ x, weights = peso)
fit.wls = fitted (mod.wls)

# M-Estimator
mod.M = rlm(y~x)
fit.M = fitted(mod.M)

# MM-Estimator
mod.rob = lmrob(y~x)
fit.rob = fitted(mod.rob)

# Regressao L1 (mediana)
mod.l1 = rq(y~x, tau = 0.5)
fit.l1 = fitted(mod.l1)

# ETKRR
mod.etkrrs1 = kernel.reg1(x, y)
mod.etkrrs2 = kernel.reg2(x, y)
mod.etkrrs3 = kernel.reg3(x, y)

plot(x,y)
lines(x, mod.etkrrs1$fitted, col = 2)
lines(x, mod.etkrrs2$fitted, col = 3)
lines(x, mod.etkrrs3$fitted, col = 4) ## vamos escolher o s3, por menor custo e facil intepretacao

#  Vizualiando e comparando os modelos ------------------------------------
colors.l = c("black",'red','green','deepskyblue3','darkorange', 'deeppink')
models = c("OLS", "WLS", "M",'ROB', 'L1', 'ETKRRS1')
plot(x,y, main = "Red Wine Data",
     xlab='Free Sulfur Dioxide',
     ylab='Total Sulfur Dioxide')#, ylim=c(0, 30))
lines(x, fit.ols, col = colors.l[1])
lines(x, fit.wls, col = colors.l[2])
lines(x, fit.M, col = colors.l[3])
lines(x, fit.rob, col = colors.l[4])
lines(x, fit.l1, col = colors.l[5])
lines(x, mod.etkrrs1$fitted, col = colors.l[6])
grid()
legend(0, 290, legend=models,
       col=colors.l, lty=c(1,1,1), cex=0.54, text.font=2, bg='transparent')
#legend(12.6, 0.996, legend=models,
#              col=colors.l, lty=c(1,1,1), cex=0.5, text.font=2, bg='transparent')
