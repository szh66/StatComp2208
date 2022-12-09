#' @import MASS
#' @import boot
#' @import DAAG
#' @import glmnet


#' @title design matrix
#' @name X
#' @description design matrix 
#' @examples
#' \dontrun{
#' data(X)
#' cor(t(X))[1,2]
#' }
NULL

#' @title response vector
#' @name Y
#' @description response vector
#' @examples
#' \dontrun{
#' data(Y)
#' print(Y)
#' }
NULL

#' @title real coefficient
#' @name beta
#' @description real coefficient
#' @examples
#' \dontrun{
#' data(beta)
#' 
#' print(beta)
#' }
NULL

#' @title punishment coefficient 
#' @name lambda_univ
#' @description punishment coefficient of Lasso and scale lasso
#' @examples
#' \dontrun{
#' data(lambda_univ)
#' print(lambda_univ)
#' }
NULL



obj_lasso<-function(beta,x,y,lambda){
  n<-dim(x)[1]
  p<-dim(x)[2]
  return(sum((y-x%*%beta)^2)/(2*n)+lambda*sum(abs(beta)))
}


#' @title Raw lasso method using optim 
#' @description given punishment coefficient lambda and initial preconceived beta, this function return estimator of beta
#' @param x the design matrix
#' @param y response vector
#' @param lambda the punishment coefficient lambda
#' @return the estimator of beta 
#' @importFrom stats optim 
#' @examples
#' \dontrun{
#'library(MASS)
#'roi<-1/5;alpha<-2
#'n<-30;p<-500
#'repl<-1
#'lambda_univ<-sqrt(log(p)*2/n)
#'beta<-numeric(p)
#'for(i in 1:p){
#'  beta[i]<-3*lambda_univ/(i^alpha)
#'}
#'beta[250]<-beta[300]<-beta[350]<-beta[400]<-beta[450]<-beta[500]<-3*lambda_univ
#'sigma<-matrix(0,p,p)
#'for(i in 1:p){
#'  for(j in 1:p){
#'    sigma[i,j]<-roi^abs(j-i)
#'  }
#'}
#'X_hat<-mvrnorm(repl*n,rep(0,p),sigma)
#'for(i in 1:repl){
#'  for(j in 1:p){
#'    xij<-X_hat[(n*(i-1)+1):(n*i),j]
#'    X_hat[(n*(i-1)+1):(n*i),j]<-sqrt(n)*xij/sqrt(sum(xij^2))
#'  } 
#'}
#'X<-X_hat
#'e<-rnorm(repl*n)
#'Y<-X%*%beta+e
#'beta_hat<-lasso(X,Y,lambda_univ)
#'}
#' @export
lasso<-function(x,y,lambda){  
  p<-dim(x)[2]
  init<-rep(0,p)
  return(optim(init,obj_lasso,method = "BFGS" ,x=x,y=y,lambda=lambda)$par)
}


obj_scalelasso<-function(beta_sigma,x,y,lambda){ 
  n<-dim(x)[1]
  p<-dim(x)[2]
  return(sum((y-x%*%beta_sigma[1:p])^2)/(2*n*beta_sigma[p+1])+lambda*sum(abs(beta_sigma[1:p]))+beta_sigma[p+1]/2)
}

#' @title scale lasso method using optim 
#' @description given punishment coefficient lambda and initial preconceived beta, this function return estimator of beta
#' @param x the design matrix
#' @param y response vector
#' @param lambda the punishment coefficient lambda
#' @return the estimator of beta and sigma 
#' @importFrom stats optim
#' @examples
#' \dontrun{
#'library(MASS)
#'roi<-1/5;alpha<-2
#'n<-30;p<-500
#'repl<-1
#'lambda_univ<-sqrt(log(p)*2/n)
#'beta<-numeric(p)
#'for(i in 1:p){
#'  beta[i]<-3*lambda_univ/(i^alpha)
#'}
#'beta[250]<-beta[300]<-beta[350]<-beta[400]<-beta[450]<-beta[500]<-3*lambda_univ
#'sigma<-matrix(0,p,p)
#'for(i in 1:p){
#'  for(j in 1:p){
#'    sigma[i,j]<-roi^abs(j-i)
#'  }
#'}
#'X_hat<-mvrnorm(repl*n,rep(0,p),sigma)
#'for(i in 1:repl){
#'  for(j in 1:p){
#'    xij<-X_hat[(n*(i-1)+1):(n*i),j]
#'    X_hat[(n*(i-1)+1):(n*i),j]<-sqrt(n)*xij/sqrt(sum(xij^2))
#'  } 
#'}
#'X<-X_hat
#'e<-rnorm(repl*n)
#'Y<-X%*%beta+e
#'beta_hat<-scale_lasso(X,Y,lambda_univ)
#'}
#' @export
scale_lasso<-function(x,y,lambda){ 
  p<-dim(x)[2]
  init_beta<-rep(0,p)
  init_sigma<-1
  return(optim(c(init_beta,init_sigma),obj_scalelasso,method = "BFGS" ,x=x,y=y,lambda=lambda)$par)
} 


beta_sigma_init<-function(x,y){ 
  p<-dim(x)[2]
  n<-dim(x)[1]
  lambda_univ<-sqrt(log(p)*2/n)
  coef<-scale_lasso(x,y,lambda = lambda_univ)
  return(coef)
}


#' @title Generate score vector z_j for each column x_j of design X
#' @description Generate score vector z_j for each column x_j of design X
#' @param x the design matrix
#' @return score vector z 
#' @importFrom  glmnet glmnet
#' @examples
#' \dontrun{
#'library(MASS)
#'roi<-1/5;alpha<-2
#'n<-30;p<-500
#'repl<-1
#'lambda_univ<-sqrt(log(p)*2/n)
#'sigma<-matrix(0,p,p)
#'for(i in 1:p){
#'  for(j in 1:p){
#'    sigma[i,j]<-roi^abs(j-i)
#'  }
#'}
#'X_hat<-mvrnorm(repl*n,rep(0,p),sigma)
#'for(i in 1:repl){
#'  for(j in 1:p){
#'    xij<-X_hat[(n*(i-1)+1):(n*i),j]
#'    X_hat[(n*(i-1)+1):(n*i),j]<-sqrt(n)*xij/sqrt(sum(xij^2))
#'  } 
#'}
#'X<-X_hat
#'Z<-z_generate(X)
#'}
#' @export
z_generate<-function(x){ 
  eta_star<-sqrt(2*log(p))
  k0<-1/4
  n<-dim(x)[1]
  p<-dim(x)[2]
  z<-matrix(0,n,p)
  for(j in 1:p){
    init <- -beta[-j]/beta[j]
    f<-glmnet(x[,-j],x[,j],nlambda = 100)  
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

#' @title Low dimensional projection estimator 
#' @description use LDPE method to give a smaller bias estimator for beta_j
#' @param x the design matrix
#' @param y response vector
#' @return the estimators of beta and sd(sigma) of beta 
#' @examples
#' \dontrun{
#'library(MASS)
#'roi<-1/5;alpha<-2
#'n<-30;p<-500
#'repl<-1
#'lambda_univ<-sqrt(log(p)*2/n)
#'beta<-numeric(p)
#'for(i in 1:p){
#'  beta[i]<-3*lambda_univ/(i^alpha)
#'}
#'beta[250]<-beta[300]<-beta[350]<-beta[400]<-beta[450]<-beta[500]<-3*lambda_univ
#'sigma<-matrix(0,p,p)
#'for(i in 1:p){
#'  for(j in 1:p){
#'    sigma[i,j]<-roi^abs(j-i)
#'  }
#'}
#'X_hat<-mvrnorm(repl*n,rep(0,p),sigma)
#'for(i in 1:repl){
#'  for(j in 1:p){
#'    xij<-X_hat[(n*(i-1)+1):(n*i),j]
#'    X_hat[(n*(i-1)+1):(n*i),j]<-sqrt(n)*xij/sqrt(sum(xij^2))
#'  } 
#'}
#'X<-X_hat
#'e<-rnorm(repl*n)
#'Y<-X%*%beta+e
#'result<-LDPE(X,Y)
#'}
#' @export
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


#' @title Generate score vector z_j for each column x_j of design X using RLPDE method
#' @description Generate score vector z_j for each column x_j of design X using RLPDE method
#' @param x the design matrix
#' @param m the selected number of maximum correlated beta_i of beta_j
#' @return score vector z 
#' @importFrom  glmnet glmnet
#' @examples
#' \dontrun{
#'library(MASS)
#'roi<-1/5;alpha<-2
#'n<-30;p<-500
#'repl<-1
#'lambda_univ<-sqrt(log(p)*2/n)
#'sigma<-matrix(0,p,p)
#'for(i in 1:p){
#'  for(j in 1:p){
#'    sigma[i,j]<-roi^abs(j-i)
#'  }
#'}
#'X_hat<-mvrnorm(repl*n,rep(0,p),sigma)
#'for(i in 1:repl){
#'  for(j in 1:p){
#'    xij<-X_hat[(n*(i-1)+1):(n*i),j]
#'    X_hat[(n*(i-1)+1):(n*i),j]<-sqrt(n)*xij/sqrt(sum(xij^2))
#'  } 
#'}
#'X<-X_hat
#'m<-4
#'Z<-RLDPE_z_generate(X,m)
#'}
#' @export
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
    f<-glmnet(x_nj,x_j,nlambda = 30)  
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

#' @title restricted low dimensional projection estimator 
#' @description use RLDPE method to give a smaller bias estimator for beta_j
#' @param x the design matrix
#' @param y response vector
#' @param m the selected number of maximum correlated beta_i of beta_j
#' @return the estimators of beta and sd(sigma) of beta \code{n}
#' @examples
#' \dontrun{
#'library(MASS)
#'roi<-1/5;alpha<-2
#'n<-30;p<-500
#'repl<-1
#'lambda_univ<-sqrt(log(p)*2/n)
#'beta<-numeric(p)
#'for(i in 1:p){
#'  beta[i]<-3*lambda_univ/(i^alpha)
#'}
#'beta[250]<-beta[300]<-beta[350]<-beta[400]<-beta[450]<-beta[500]<-3*lambda_univ
#'sigma<-matrix(0,p,p)
#'for(i in 1:p){
#'  for(j in 1:p){
#'    sigma[i,j]<-roi^abs(j-i)
#'  }
#'}
#'X_hat<-mvrnorm(repl*n,rep(0,p),sigma)
#'for(i in 1:repl){
#'  for(j in 1:p){
#'    xij<-X_hat[(n*(i-1)+1):(n*i),j]
#'    X_hat[(n*(i-1)+1):(n*i),j]<-sqrt(n)*xij/sqrt(sum(xij^2))
#'  } 
#'}
#'X<-X_hat
#'e<-rnorm(repl*n)
#'Y<-X%*%beta+e
#'m<-4
#'result<-RLDPE(X,Y,m)
#'}
#' @export
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


