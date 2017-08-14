### Simulation for Pooling Sub Parameters with Confounding Z
### Copy Right by Hao Zhou
### Please Contact hzhou@stat.wisc.edu for Questions and Bugs

### 2 Sets Example
library(MASS)

### Define Parameters
k = 2;p = 3;q = 5;n = 2^6;noise = c(3,0.5)

### Generate True Signals
beta1 = runif(p)*4;beta2 = beta1+0.1*rep(1,p)
gamma1 = c(rep(1,2),rep(2,q-2));gamma2 = c(rep(2,2),rep(1,q-2))

### Define True Covariance Matrices
Sigma1_X = 0.5*diag(rep(1,p))+0.5*rep(1,p)%*%t(rep(1,p));
Sigma2_X = 0.5*diag(rep(1,p))+0.5*rep(1,p)%*%t(rep(1,p));
Sigma1_Z = 0.8*diag(rep(1,q))+0.2*rep(1,q)%*%t(rep(1,q));
Sigma2_Z = 0.8*diag(rep(1,q))+0.2*rep(1,q)%*%t(rep(1,q));
Sigma1_XZ = 0.2*rep(1,p)%*%t(rep(1,q));
Sigma2_XZ = 0.2*rep(1,p)%*%t(rep(1,q));
Sigma1 = rbind(cbind(Sigma1_X,Sigma1_XZ),cbind(t(Sigma1_XZ),Sigma1_Z));
Sigma2 = rbind(cbind(Sigma2_X,Sigma2_XZ),cbind(t(Sigma2_XZ),Sigma2_Z));


### Generate Features (Predictors)
Predictor1 = mvrnorm(n,rep(0,p+q),Sigma1)
X1 = Predictor1[,1:p]
Z1 = Predictor1[,(p+1):(p+q)]
Predictor2 = mvrnorm(n,rep(0,p+q),Sigma2)
X2 = Predictor2[,1:p]
Z2 = Predictor2[,(p+1):(p+q)]

### Generate Noise and Responses
error1 = rnorm(n,0,noise[1])
y1 = X1%*%beta1 + Z1%*%gamma1 + error1
error2 = rnorm(n,0,noise[2])
y2 = X2%*%beta2 + Z2%*%gamma2 + error2

### Get High-Level Information from Each Set
Fitols_confounding <- function(X,Z,y){
  fit = lm(y~X+Z)
  coeff = fit$coefficients[2:(1+ncol(X)+ncol(Z))]
  noisest = sqrt(sum(fit$res^2/fit$df.res))
  n = nrow(X)
  p = ncol(X)
  q = ncol(Z)
  samplecov_raw = cov(cbind(X,Z))*(nrow(X)-1)/nrow(X)
  samplecov_x = samplecov_raw[1:p,1:p]
  samplecov_z = samplecov_raw[(p+1):(p+q),(p+1):(p+q)]
  samplecov_xz = samplecov_raw[1:p,(p+1):(p+q)]
  samcov_cond = samplecov_x-samplecov_xz%*%solve(samplecov_z,t(samplecov_xz))  
  Inform = list(beta_est=coeff[1:ncol(X)],gamma_est=coeff[(ncol(X)+1):(ncol(X)+ncol(Z))],
                noise_est=noisest,samplecov_cond=samcov_cond,n=n)
  return(Inform)
}

Set1_Inform = Fitols_confounding(X1,Z1,y1) 
Set2_Inform = Fitols_confounding(X2,Z2,y2)
ksets_Inform = list(Set1=Set1_Inform,Set2=Set2_Inform)

### Conduct the Hypothesis Test
pvalue = Hypotest_subparam(k,ksets_Inform)

### Fit k Sets Model and Check Performance
tau_est = Set1_Inform$noise_est/Set2_Inform$noise_est
y=rbind(y1,tau_est*y2)
X=rbind(X1,tau_est*X2)
bdiag <- Matrix::bdiag
Z=as.matrix(bdiag(Z1,tau_est*Z2))
fitk = lm(y~X+Z)
betak_est = fitk$coefficients[2:(ncol(X)+1)]
beta1_est = Set1_Inform$beta_est

y1new=y1-X1%*%betak_est
fitkgamma1 = lm(y1new~Z1)
gammak1_est = fitkgamma1$coefficients[2:(1+ncol(Z1))]
gamma1_est = Set1_Inform$gamma_est
performance = c(n,pvalue,sum((beta1_est-beta1)^2),sum((betak_est-beta1)^2),sum((gamma1_est-gamma1)^2),sum((gammak1_est-gamma1)^2))
names(performance) = c('sample size','p-value','Set1 One-Time Squared-Error','k-Sets One-Time Squared-Error',
                       'Set1 Confounding gamma One-Time Squared-Error',
                       'k-sets Confounding gamma One-Time Squared-Error')
performance
