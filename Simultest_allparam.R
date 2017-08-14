### Simulation for Pooling All Parameters
### Copy Right by Hao Zhou
### Please Contact hzhou@stat.wisc.edu for Questions and Bugs

### 2 Sets Example
library(MASS)

### Define Parameters
k = 2;p = 3;n = 2^6;noise = c(3,0.5)

### Generate True Signals
beta1 = runif(p)*4;beta2 = beta1+0.1*rep(1,p)

### Define True Covariance Matrices
Sigma1 = 0.5*diag(rep(1,p))+0.5*rep(1,p)%*%t(rep(1,p))
Sigma2 = 0.5*diag(rep(1,p))+0.5*rep(1,p)%*%t(rep(1,p))

### Generate Features (Predictors)
X1 = mvrnorm(n,rep(0,p),Sigma1)
X2 = mvrnorm(n,rep(0,p),Sigma2)

### Generate Noise and Responses
error1 = rnorm(n,0,noise[1])
y1 = X1%*%beta1 + error1
error2 = rnorm(n,0,noise[2])
y2 = X2%*%beta2 + error2

### Get High-Level Information from Each Set
Fitols <- function(X,y){
  fit = lm(y~X)
  coeff = fit$coefficients[2:(ncol(X)+1)]
  noiest = sqrt(sum(fit$res^2/fit$df.res))
  n = nrow(X)
  sacov = cov(X)*(nrow(X)-1)/nrow(X)
  Inform = list(beta_est=coeff,noise_est=noiest,samplecov=sacov,n=n)
  return(Inform)
}

Set1_Inform = Fitols(X1,y1) 
Set2_Inform = Fitols(X2,y2)
ksets_Inform = list(Set1=Set1_Inform,Set2=Set2_Inform)

### Conduct the Hypothesis Test
pvalue = Hypotest_allparam(k,ksets_Inform)

### Fit k Sets Model and Check Performance
tau_est = Set1_Inform$noise_est/Set2_Inform$noise_est
y=rbind(y1,tau_est*y2)
X=rbind(X1,tau_est*X2)
fitk = lm(y~X)
betak_est = fitk$coefficients[2:(ncol(X)+1)]
beta1_est = Set1_Inform$beta_est
performance = c(n,pvalue,sum((beta1_est-beta1)^2),sum((betak_est-beta1)^2))
names(performance) = c('sample size','p-value','Set1 One-Time Squared-Error','k-Sets One-Time Squared-Error')
performance
