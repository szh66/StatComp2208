## -----------------------------------------------------------------------------
obj_lasso<-function(beta,x,y,lambda){  
  return(sum((y-x%*%beta)^2)/(2*n)+lambda*sum(abs(beta)))
}

lasso<-function(x,y,lambda){  
  p<-dim(x)[2]
  init<-rep(0,p)
  return(optim(init,obj_lasso,method = "BFGS" ,x=x,y=y,lambda=lambda)$par)
}

obj_scalelasso<-function(beta_sigma,x,y,lambda){ 
  n<-dim(x)[1]
  p<-dim(x)[2]
  return(sum((y-x%*%beta_sigma[1:p])^2)/(2*n*beta_sigma[p+1])+
           lambda*sum(abs(beta_sigma[1:p]))+beta_sigma[p+1]/2)
}
scale_lasso<-function(x,y,lambda){ 
  p<-dim(x)[2]
  init_beta<-rep(0,p)
  init_sigma<-1
  return(optim(c(init_beta,init_sigma),obj_scalelasso,method = "BFGS"
               ,x=x,y=y,lambda=lambda)$par)
} 
beta_sigma_init<-function(x,y){ 
  p<-dim(x)[2]
  n<-dim(x)[1]
  lambda_univ<-sqrt(log(p)*2/n)
  coef<-scale_lasso(x,y,lambda = lambda_univ)
  return(coef)
}

## ----echo=FALSE---------------------------------------------------------------
library(MASS)
roi<-1/5;alpha<-2
n<-30;p<-500
repl<-1
lambda_univ<-sqrt(log(p)*2/n)
beta<-numeric(p)
for(i in 1:p){
  beta[i]<-3*lambda_univ/(i^alpha)
}
beta[250]<-beta[300]<-beta[350]<-beta[400]<-beta[450]<-beta[500]<-3*lambda_univ
sigma<-matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    sigma[i,j]<-roi^abs(j-i)
  }
}
X_hat<-mvrnorm(repl*n,rep(0,p),sigma)
for(i in 1:repl){
  for(j in 1:p){
    xij<-X_hat[(n*(i-1)+1):(n*i),j]
    X_hat[(n*(i-1)+1):(n*i),j]<-sqrt(n)*xij/sqrt(sum(xij^2))
  } 
}
X<-X_hat
e<-rnorm(repl*n)
Y<-X%*%beta+e



## -----------------------------------------------------------------------------

beta_hat<-lasso(X,Y,lambda_univ)
mean_all=mean(beta_hat-beta)
round(c(mean_bias=mean_all),4)

## -----------------------------------------------------------------------------

beta_hat<-scale_lasso(X,Y,lambda_univ)
mean_all=mean(beta_hat[1:p]-beta)
round(c(mean_bias=mean_all),4)

## -----------------------------------------------------------------------------
z_generate<-function(x){ 
  eta_star<-sqrt(2*log(p))
  k0<-1/4
  n<-dim(x)[1]
  p<-dim(x)[2]
  z<-matrix(0,n,p)
  for(j in 1:p){
    f<-glmnet::glmnet(x[,-j],x[,j],nlambda = 100)  
    all_lmabda<-f$lambda
    all_beta<-f$beta
    k<-length(all_lmabda)
    lambda_star<-all_lmabda[k] 
    eta_j<-numeric(k)
    tao_j<-numeric(k)
    zj<-matrix(0,n,k)
    star<-k
    for(s in 1:k){
      beta_j<-all_beta[,s]
      zj[,s]<-x[,j]-x[,-j]%*%beta_j
      m<-max(t(zj[,s])%*%x[,-j])
      eta_j[s]<-m/sqrt(sum(zj[,s]^2))
      tao_j[s]<-sqrt(sum(zj[,s]^2))/abs(t(zj[,s])%*%x[,j])
    }
    if(eta_j[k]>=eta_star){
      star<-k
      
    }
    else{
      if(prod(!(eta_j>=eta_star))){
        star<-k
      }
      else{
        tao_star<-(1+k0)*min(tao_j[eta_j>=eta_star]) 
        star<-(1:k)[eta_j==min(eta_j[tao_j<=tao_star])]
      }
    }
    z[,j]<-zj[,star]
  }
  return(z)
}

## -----------------------------------------------------------------------------

Z<-z_generate(X)
print(Z[,1])

## -----------------------------------------------------------------------------
LDPE<-function(x,y){
  n<-dim(x)[1]
  p<-dim(x)[2]
  beta_sigma_initial<-beta_sigma_init(x,y)
  beta_init<-beta_sigma_initial[1:p]
  sigma_init<-beta_sigma_initial[p+1]
  z<-z_generate(x)
  
  beta_hat<-numeric(p)  
  tao<-numeric(p) 
  
  for(j in 1:p){
    beta_hat[j]<-beta_init[j]+t(z[,j])%*%(y-x%*%beta_init)/(t(z[,j])%*%x[,j])
    tao[j]<-sqrt(sum(z[,j]^2))/abs(t(z[,j])%*%x[,j])
  }
  return(list(beta_hat=beta_hat,sigma_hat=tao*sigma_init))
}

## -----------------------------------------------------------------------------
repl<-1
result1<-list()
result1[[1]]<-LDPE(X,Y)
mean_all=mean(result1[[1]]$beta_hat-beta)
round(c(mean_bias=mean_all),4)
a<-0.05
for (i  in 1:p) {
  intervals1<-matrix(0,2,repl)
  In<-numeric(repl)
  
  for(j in 1:repl){
    
    intervals1[1,j]<-(result1[[j]]$beta_hat)[i]+10*qnorm(1-a/2)*(result1[[j]]
                                                              $sigma_hat)[i]
    intervals1[2,j]<-(result1[[j]]$beta_hat)[i]-10*qnorm(1-a/2)*(result1[[j]]
                                                              $sigma_hat)[i]
    if(beta[i]<=intervals1[1,j] & beta[i]>=intervals1[2,j]){
      In[j]<-1
    }
  }
}
print(intervals1)

## -----------------------------------------------------------------------------
RLDPE_z_generate<-function(x,m){ 
  
  eta_star<-sqrt(2*log(p))
  k0<-1/4
  n<-dim(x)[1]
  p<-dim(x)[2]
  z<-matrix(0,n,p)
  for(j in 1:p){
    idx<-numeric(m)
    product<-as.array(t(x[,j])%*%x[,-j])
    First4<-sort(product,decreasing = T)[1:m]
    for(l in 1:m){
      idx[l]<-(1:(p-1))[product==First4[l]]
    }
    
    xx<-(x[,-j])[,idx]
    P<-xx%*%solve(t(xx)%*%xx)%*%t(xx)
    x_j<-x[,j]-P%*%x[,j]
    x_nj<-x[,-j]-P%*%x[,-j]
    f<-glmnet::glmnet(x_nj,x_j,nlambda = 30)  
    all_lmabda<-f$lambda
    all_beta<-f$beta
    k<-length(all_lmabda)
    lambda_star<-all_lmabda[k] 
    eta_j<-numeric(k)
    tao_j<-numeric(k)
    zj<-matrix(0,n,k)
    star<-k
    for(s in 1:k){
      beta_j<-all_beta[,s]
      zj[,s]<-x_j-x_nj%*%beta_j
      ma<-max(t(zj[,s])%*%x[,-j])
      eta_j[s]<-ma/sqrt(sum(zj[,s]^2))
      tao_j[s]<-sqrt(sum(zj[,s]^2))/abs(t(zj[,s])%*%x[,j])
    }
    if(eta_j[k]>=eta_star)star<-k
    else{
      if(prod(!(eta_j>=eta_star))){
        star<-k
      }
      else{
        tao_star<-(1+k0)*min(tao_j[eta_j>=eta_star]) 
        star<-(1:k)[eta_j==min(eta_j[tao_j<=tao_star])]
      }
      
      
    }
    z[,j]<-zj[,star]
  }
  
  return(z)
}

## -----------------------------------------------------------------------------
m<-4
Z<-RLDPE_z_generate(X,m)
print(Z[,1])

## -----------------------------------------------------------------------------
RLDPE<-function(x,y,m){
  n<-dim(x)[1]
  p<-dim(x)[2]
  beta_sigma_initial<-beta_sigma_init(x,y)
  beta_init<-beta_sigma_initial[1:p]
  sigma_init<-beta_sigma_initial[p+1]
  
  z<-RLDPE_z_generate(x,m)
  
  beta_hat<-numeric(p)  
  tao<-numeric(p) 
  
  for(j in 1:p){
    beta_hat[j]<-beta_init[j]+t(z[,j])%*%(y-x%*%beta_init)/(t(z[,j])%*%x[,j])
    tao[j]<-sqrt(sum(z[,j]^2))/abs(t(z[,j])%*%x[,j])
  }
  return(list(beta_hat=beta_hat,sigma_hat=tao*sigma_init))
}

## -----------------------------------------------------------------------------
repl<-1
result1_R<-list()
result1_R[[1]]<-RLDPE(X,Y,4)
mean_all=mean(result1_R[[1]]$beta_hat-beta)
round(c(mean_bias=mean_all),4)
intervals1_R<-matrix(0,2,p)
rate1_R<-numeric(p)
width1_R<-matrix(0,repl,p)  #RLDPE的区间长度
for (i  in 1:p) { #对每个变量计算：
  intervals1_R<-matrix(0,2,repl)
  In<-numeric(repl)
  
  for(j in 1:repl){
    intervals1_R[1,j]<-(result1_R[[j]]$beta_hat)[i]+qnorm(1-a/2)*(result1_R[[j]]
                                                                  $sigma_hat)[i]
    intervals1_R[2,j]<-(result1_R[[j]]$beta_hat)[i]-qnorm(1-a/2)*(result1_R[[j]]
                                                                  $sigma_hat)[i]
  }
}


