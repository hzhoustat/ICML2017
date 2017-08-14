### Hypothesis Test for Pooling All Parameters
### Copy Right by Hao Zhou
### Please Contact hzhou@stat.wisc.edu for Questions and Bugs

### Define the Base Set to be Set1
### ksets_Inform = list(Set1,Set2,...,Setk)
### Set1 = list(beta_est= ,noise_est= ,samplecov= ,n= )
### Set2 = list(beta_est= ,noise_est= ,samplecov= ,n= )
### ...
### Setk = list(beta_est= ,noise_est= ,samplecov= ,n= )

Hypotest_allparam <- function(k,ksets_Inform){
  basemat = matrix(rep(1,(k-1)^2),ncol=(k-1))
  Set1_Inform = ksets_Inform[[1]]
  n1 = Set1_Inform$n;noisest1 = Set1_Inform$noise_est
  samplecov1 = Set1_Inform$samplecov
  Set1mat = noisest1^2/n1*solve(samplecov1)
  G = kronecker(basemat,Set1mat)
  p = nrow(Set1mat)
  for (i in 2:k){
    Seti_Inform = ksets_Inform[[i]]
    ni = Seti_Inform$n;noisesti =Seti_Inform$noise_est
    samplecovi = Seti_Inform$samplecov
    Setimat = noisesti^2/ni*solve(samplecovi)
    stpos = (i-2)*p+1
    edpos = (i-1)*p
    G[stpos:edpos,stpos:edpos] = G[stpos:edpos,stpos:edpos]+Setimat
  }
  basevec = rep(1,(k-1))
  beta1_est = ksets_Inform[[1]]$beta_est
  dbeta_est = kronecker(basevec,beta1_est)
  for (i in 2:k){
    betai_est = ksets_Inform[[i]]$beta_est
    stpos =(i-2)*p+1
    edpos = (i-1)*p
    dbeta_est[stpos:edpos] = dbeta_est[stpos:edpos]-betai_est
  }
  test_stat = t(dbeta_est)%*%solve(G,dbeta_est)
  pvalue=1-pchisq(test_stat,df=(k-1)*p,ncp=1);
  return(pvalue)
}