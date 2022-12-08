## -----------------------------------------------------------------------------
library(MASS)
data("Cars93")
plot(Cars93$Manufacturer)
plot(Cars93$EngineSize)
hist(Cars93$Price,col='pink',xlab='price/1k ',main='Cars93 price')

## -----------------------------------------------------------------------------
knitr::kable(Cars93[1:12,1:6])
rmarkdown::paged_table(Cars93[1:18,1:6])

## ----echo=FALSE---------------------------------------------------------------
a<-2;b<-2;n=100
f<-function(x){
  (a/b)*(b/x)^(a+1)
}
inverse_F<-function(u){
  b*(1-u)^(-1/a)
}
u<-runif(n)
x<-inverse_F(u)
hist(x,breaks=60,freq = FALSE,xlim=c(0,20),ylim=c(0,1.2),ylab=NULL)
par(new=TRUE)
curve(f,2,30,xlim = c(0,20),ylim=c(0,1.2),ylab=NULL,col="blue")


## ----echo=FALSE---------------------------------------------------------------
a<-3;b<-2
f<-function(x){
  1/beta(a,b)*x^(a-1)*(1-x)^(b-1)
}
c=f((a-1)/(a+b-2))

n <- 1000;k<-0;sample <- numeric(n)
while (k < n) {
  u <- runif(1)
  x <- runif(1) 
  if (dbeta(x,3,2)/c > u) {
    k <- k + 1 
    sample[k] <- x
  }
}  
hist(sample,xlim=c(0,1),ylim=c(0,2),breaks = 50,freq = FALSE,ylab = NULL,xlab = NULL)
par(new=TRUE)
curve(f,ylim=c(0,2),col='red',ylab=NULL)
abline(h=c,col='pink')



## ----echo=FALSE---------------------------------------------------------------
n <- 1000; r <- 4; beta <- 3
lambda <- rgamma(n, r, beta)
x <- rexp(n,lambda)
x<-data.frame(x)
rmarkdown::paged_table(x)

## ----echo=FALSE---------------------------------------------------------------
n <- 1000; r <- 4; beta <- 2
lambda <- rgamma(n, r, beta)
x <- rexp(n,lambda)
hist(x,breaks = 50,xlim=c(0,15),ylim=c(0,1.3),freq=FALSE,ylab = NULL)

f<-function(y){
  (r/beta)*(beta/(beta+y))^(r+1)
}
par(new=TRUE)
curve(f,xlim = c(0,15),ylim=c(0,1.3),col='red')

## ----echo=FALSE---------------------------------------------------------------
n<-c(10^4,2*10^4,4*10^4,6*10^4,8*10^4)
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}
}
test<-lapply(1:length(n), function(i)sample(1:n[i]))
#rmarkdown::paged_table(data.frame(test[[1]]))
sorted_test<-lapply(1:length(test),function(i)quick_sort(test[[i]]))

## ----echo=FALSE---------------------------------------------------------------
rmarkdown::paged_table(data.frame('Original_list'=test[[1]],'Sorted_list'=sorted_test[[1]]))

## ----echo=FALSE---------------------------------------------------------------
rmarkdown::paged_table(data.frame('Original_list'=test[[2]],'Sorted_list'=sorted_test[[2]]))

## ----echo=FALSE---------------------------------------------------------------
rmarkdown::paged_table(data.frame('Original_list'=test[[3]],'Sorted_list'=sorted_test[[3]]))

## ----echo=FALSE---------------------------------------------------------------
rmarkdown::paged_table(data.frame('Original_list'=test[[4]],'Sorted_list'=sorted_test[[4]]))

## ----echo=FALSE---------------------------------------------------------------
rmarkdown::paged_table(data.frame('Original_list'=test[[5]],'Sorted_list'=sorted_test[[5]]))

## ----echo=FALSE---------------------------------------------------------------
a<-numeric(length(n))
sim_num<-100 # the number of simulations for each n 
for(i in 1:length(n)){
  t_j<-sapply(1:sim_num,function(j){# For each n, we did 100 times simulation.
    test<-sample(1:n[i])
    return(system.time(quick_sort(test))[[1]])
  })# All the time costs of 100 simulations formed vector t_j
  a[i]<-mean(t_j) 
}
for(i in 1:length(n)){
  cat('\nThe computation time averaged over 100 simulations of n =',n[i],'is',a[i])
}

## ----echo=FALSE---------------------------------------------------------------
t<-sapply(n,function(i)i*log(i))
my_lm<-lm(t~a)
plot(a,t,main='a~t regression')
abline(my_lm,lwd=2,col='red')

## ----echo=FALSE---------------------------------------------------------------
MC.Phi <- function(R = 1, antithetic = FALSE) {
  u <- runif(R)
  if (antithetic){
    v <- 1 - u 
    u <- c(u, v)
  } 
  g <- exp(u)
  cdf <- mean(g) 
  cdf
}
m <- 1000
MC1 <- MC2 <- numeric(m) #MC2[i]=(exp(u)+exp(1-u))/2; MC1[i]=exp(u)
for (i in 1:m) {
  MC1[i] <- MC.Phi(R = 1, antithetic = FALSE)
  MC2[i] <- MC.Phi(R = 1, antithetic = TRUE)
}
p<-1-2*var(MC2)/var(MC1) 
print(p)

## ----echo=FALSE---------------------------------------------------------------
g<-function(x){
  x^2*exp(-x^2/2)/sqrt(2*pi)
}
curve(g,1,10)

## -----------------------------------------------------------------------------
set.seed(12)
n<-10;m<-10000
mu<-0;sigma<-1
low<-high<-numeric(m) # left end point and right end point
for(i in 1:m){
  x<-rnorm(n,mu,sigma)
  x_bar<-mean(x)
  S<-sd(x)
  low[i]<-x_bar-S*qt(0.975,n-1)/sqrt(n)
  high[i]<-x_bar-S*qt(0.025,n-1)/sqrt(n)
}
mean((0>low)*(0<high)) #empirical estimate of the confidence level


## -----------------------------------------------------------------------------
set.seed(123)
n<-c(10,100,1000)
m<-10000
sigma1 <- 1
sigma2 <- 1.5
# generate samples under H1 to estimate power
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
for (i in 1:3){
  power <- mean(replicate(m, expr={
    x <- rnorm(n[i], 0, sigma1)
    y <- rnorm(n[i], 0, sigma2)
    count5test(x, y)
    }))
  print(power)
}





## -----------------------------------------------------------------------------
alpha<-0.055
set.seed(123)
F_hat<-numeric(m)
for(i in 1:3){
  for(j in 1:m){
  x <- rnorm(n[i], 0, sigma1)
  y <- rnorm(n[i], 0, sigma2)
  F_hat[j] = var(x)/var(y)
  }
  low<-qf(alpha/2,n[i]-1,n[i]-1);high<-qf(1-alpha/2,n[i]-1,n[i]-1) # left end point and right end point of accept interval 
  print(1-mean((F_hat<high)*(F_hat>low))) # empirical estimate of power
}


## -----------------------------------------------------------------------------
library(boot)
x<-c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
stat<-function(x,i){1/mean(x[i])}
lambda_mle<-12/sum(x)
obj<-boot(data=x,statistic = stat,R=99)
round(c(original=obj$t0,bias=mean(obj$t)-obj$t0,se=sd(obj$t)),6)

## -----------------------------------------------------------------------------
ci<- boot.ci(obj,type=c("norm","basic","perc",'bca'))
cat('normal:\n' ,ci$norm[2:3],'\n basic:\n',ci$basic[4:5], '\n percentile:\n',ci$percent[4:5],'\n BCa:\n',ci$bca[4:5] )

## -----------------------------------------------------------------------------
set.seed(1234)
mu<-0
n<-10;m<-100
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)
stat<-function(x,i)mean(x[i])
for(i in 1:m){
  x<-rnorm(n,mu)
  obj <- boot(data=x,statistic=stat, R = 99)
  ci <- boot.ci(obj,type=c("norm","basic","perc"))
  ci.norm[i,]<-ci$norm[2:3]
  ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5]
}
cat('proportion that CI miss on the left:  norm=',mean(ci.norm[,1]>=mu),'basic =',mean(ci.basic[,1]>=mu ),'perc =',mean(ci.perc[,1]>=mu))
cat('proportion that CI miss on the right:  norm =',mean(ci.norm[,2]<=mu ),'basic =',mean(ci.basic[,2]<=mu ),'perc =',mean(ci.perc[,2]<=mu))


## -----------------------------------------------------------------------------
library(bootstrap)
n=dim(scor)[1]
theta.jack<-numeric(n)
lambda<-eigen(cor(scor))$val
theta.hat<-lambda[1]/sum(lambda)
for(i in 1:n){
  jack_df<-scor[-i,]
  Sigma<-matrix(nrow = 5,ncol=5)
  for(j in 1:5){
    for(k in j:5){
      Sigma[j,k]=Sigma[k,j]=cor(jack_df[,k],jack_df[,j])
    }
  }
  lambda_hat<-eigen(Sigma)$val
  theta.jack[i]<-lambda_hat[1]/sum(lambda_hat)
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))
round(c(bias.jack=bias.jack,
se.jack=se.jack),6)


## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1)/2)
i=0
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:(n-1)) {
  for(j in (k+1):n){
    i=i+1
    y <- magnetic[-c(k,j)]
    x <- chemical[-c(k,j)]
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(k,j)]
    e1[i]<- sum((magnetic[c(k,j)] - yhat1)^2)
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(k,j)] +
    J2$coef[3] * chemical[c(k,j)]^2
    e2[i]<- sum((magnetic[c(k,j)] - yhat2)^2)
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(k,j)]
    yhat3 <- exp(logyhat3)
    e3[i]<- sum((magnetic[c(k,j)] - yhat3)^2)
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(k,j)])
    yhat4 <- exp(logyhat4)
    e4[i]<- sum((magnetic[c(k,j)] - yhat4)^2)
  }
}
round(c(mean_e1=mean(e1), mean_e2=mean(e2), mean_e3=mean(e3), mean_e4=mean(e4)),3)

## -----------------------------------------------------------------------------
n=20
set.seed(123)
x<-rnorm(n,1,2)
y<-runif(n,2,10)
R <- 999;z <- c(x, y);K <- 1:(2*n)
reps <- numeric(R);t0 <- cor(x, y,method='spearman')
for (i in 1:R) {
  k <- sample(K, size = n, replace = FALSE)
  x1 <- z[k]; y1 <- z[-k] #complement of x1
  reps[i] <- cor(x1, y1,method='spearman')
}
p <- mean(abs(c(t0, reps)) >= abs(t0))
round(c(permutation_p=p,cor.test_p=cor.test(x,y)$p.value),3)

## -----------------------------------------------------------------------------
rw.Metropolis <- function(sigma, x0, N){ 
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (exp(-abs(y) )/exp(-abs(x[i-1])) ))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
    }
return(list(x=x, k=k))
}

Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic 
  return(r.hat)
  
}

set.seed(123)
sigma <- c(1,2,3,4) #different sigma
k <- 4 #number of chains to generate for each sigma
n <- 6000 #length of chains
b <- 500 #burn-in length
#choose overdispersed initial values for each sigma
x0 <- c(-10, -5, 5, 10)
#generate the chains
S <- list() #generate a list,each element corresponds to different sigma and contains 4 chains 
R <- numeric(length(sigma)) # R
accept_rate <- matrix(0,length(sigma),k)
for(j in 1:length(sigma)){
  
  X <- matrix(0, nrow=k, ncol=n)
  for (i in 1:k){
    results <- rw.Metropolis(sigma[j], x0[i], n)
    X[i, ] <- results[[1]] 
    accept_rate[j,i] <- (n+b-results[[2]])/(n+b)
  }
  #compute diagnostic statistics
  psi <- t(apply(X, 1, cumsum))
  for (i in 1:nrow(psi)){
    psi[i,] <- psi[i,] / (1:ncol(psi))
  }
  S[[j]] <- psi #save diagnostic stat from each sigma as a matrix
  R[j] <- Gelman.Rubin(psi)
}


#plot the sequence of R-hat statistics
for(i in 1:length(sigma)){
  psi<-S[[i]]
  rhat <- rep(0, n)
  for (j in 1:n) rhat[j] <- Gelman.Rubin(psi[,1:j])
  plot(rhat[b:n], type="l", xlab="", ylab="R",main=bquote(sigma ==.(sigma[i])))
  abline(h=1.2, lty=2)
}


## -----------------------------------------------------------------------------
#data prepare
set.seed(12345)
n<-20
x<-runif(n,1,10) 
alpha<-c(0,0,1)
beta<-c(0,1,0)
lambda<-rep(1,3)
paramas<-cbind(alpha,beta,lambda)
a_m<-a_y<-1

#linear model 
#simulate M&Y
M<-a_m+x%*%t(alpha)+rnorm(n*3,0,0.1)
Y<-a_y+M%*%diag(beta)+x%*%t(lambda)+rnorm(n*3,0,0.1)

## -----------------------------------------------------------------------------
#First: for H0: a=0
###for each of the three paramates vectors (alpha,beta)
###use permutation method to compute p-value of this test.
p<-numeric(3)
R <- 999;K <- 1:(2*n)
for(i in 1:3){
  alpha_hat<-numeric(R)
  my_fit<-lm(M[,i]~x)
  my_fit<-summary(my_fit)
  alpha_0<-my_fit$coefficients[2,1]
  z <- c(M[,i], x)
  for (j in 1:R) {
    k <- sample(K, size = n, replace = FALSE)
    Mj<- z[k]; xj<- z[-k] 
    my_fit<-lm(Mj~xj)
    my_fit<-summary(my_fit)
    alpha_hat[j]<-my_fit$coefficients[2,1]
  }
  p[i]<-mean(abs(c(alpha_0, alpha_hat)) >= abs(alpha_0))
}

## -----------------------------------------------------------------------------
round(p,3)

## -----------------------------------------------------------------------------
#First: for H0: a=0
###
###

for(i in 1:3){
  beta_hat<-numeric(R)
  my_fit<-lm(Y[,i]~M[,i]+x)
  my_fit<-summary(my_fit)
  beta_0<-my_fit$coefficients[2,1]
  z <- c(M[,i], Y[,i])
  for (j in 1:R) {
    k <- sample(K, size = n, replace = FALSE)
    Yj<- z[k]; Mj<- z[-k] 
    my_fit<-lm(Yj~Mj+x)
    my_fit<-summary(my_fit)
    beta_hat[j]<-my_fit$coefficients[2,1]
  }
  p[i]<-mean(abs(c(beta_0, beta_hat)) >= abs(beta_0))
}
round(p,3)

## -----------------------------------------------------------------------------
#First: for H0: a=0
###for each of the three paramates vectors (alpha,beta,lambda)
###use permutation method to compute p-value of this test.

for(i in 1:3){
  alpha_hat<-beta_hat<-numeric(R)
  my_fit<-lm(Y[,i]~M[,i]+x)
  my_fit<-summary(my_fit)
  beta_0<-my_fit$coefficients[2,1]
  alpha_hat<-numeric(R)
  my_fit<-lm(M[,i]~x)
  my_fit<-summary(my_fit)
  alpha_0<-my_fit$coefficients[2,1]
  z1 <- c(M[,i], Y[,i])
  z2 <- c(M[,i], x)
  for (j in 1:R) {
    k <- sample(K, size = n, replace = FALSE)
    Yj<- z1[k]; Mj<- z1[-k]
    Mj2<-z2[k]; xj<- z2[-k]
    my_fit1<-lm(Yj~Mj);my_fit2<-lm(Mj2~xj)
    my_fit1<-summary(my_fit1);my_fit2<-summary(my_fit2)
    beta_hat[j]<-my_fit1$coefficients[2,1]
    alpha_hat[j]<-my_fit2$coefficients[2,1]
  }
  c<-(abs(c(beta_0, beta_hat)) >= abs(beta_0))*(abs(c(alpha_0, alpha_hat)) >= abs(alpha_0))
  p[i]<-mean(c)
}
round(p,3)

## -----------------------------------------------------------------------------
solve_alpha<-function(N,b1,b2,b3,f0){ 
  x1 <- rpois(N,1); x2 <- rexp(N);x3<-rbinom(N,1,0.5)
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
  }
  solution <- uniroot(g,c(-20,20))
  return(solution$root)
}

## -----------------------------------------------------------------------------
N <- 1e6; b1 <- 0; b2 <- 1; b3 <- -1
f0<-c(0.1,0.01,0.001,0.0001)
alpha<-numeric(4)
for (i in 1:4){
  alpha[i]<-solve_alpha(N,b1,b2,b3,f0[i])
}
plot(alpha,f0)

## -----------------------------------------------------------------------------

u<-c(11,8,27,13,16,0,23,10,24,2)
v<-c(12,9,28,14,17,1,24,11,25,3)

## -----------------------------------------------------------------------------
# log likehood function
l<-function(lambda){
  log(prod(exp(-lambda*u)-exp(-lambda*v)))
}
f<-function(x){
  x^2-9
}
solve<-optimize(l,lower=0,upper=1,maximum = TRUE)
lambda_mle<-solve$maximum
round(lambda_mle,3)

## -----------------------------------------------------------------------------
it <- 7000  
n<-length(u)
lambda_em<-numeric(it)  
lambda_em[1]<-1  
for(i in 2:it){  
  l<-(u*exp(-lambda_em[i-1]*u)-v*exp(-lambda_em[i-1]*v))/(exp(-lambda_em[i-1]*u)-exp(-lambda_em[i-1]*v))
  s<-sum(l)
  lambda_em[i]<-n/(n/lambda_em[i-1]+s)  
}  
round(lambda_em[it],3)

## -----------------------------------------------------------------------------
l<-list(1,2,3)
a=as.vector(l)
print(unlist(l))
typeof(l[1])

## -----------------------------------------------------------------------------
print(1=='1')
print(-1<FALSE)
print('one'<2)

## -----------------------------------------------------------------------------
v<-c(1,2,3)
dim(v)

## -----------------------------------------------------------------------------
df <- data.frame(x = 1:3, y = c("a", "b", "c"))
names(attributes(df))

## -----------------------------------------------------------------------------
df <- data.frame(x = 1:3, y = c("a", "b", "c"))
as.matrix(df)

## -----------------------------------------------------------------------------
df<-data.frame()
dim(df)

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}
X<-data.frame(x=6:10,y=1:5,z=c('a','b','c','d','e'))
print(X)
apply(X[,sapply(X,is.numeric)], 2, scale01)

## -----------------------------------------------------------------------------
X<-data.frame(matrix(runif(40,1,10),ncol = 5))
X
vapply(X, sd, numeric(1))

## -----------------------------------------------------------------------------
X<-data.frame(matrix(runif(24,1,10),ncol = 3),z=strsplit('abcdefgh',split = ''))
colnames(X)<-c(1,2,3,4)
vapply(X[,vapply(X, is.numeric, logical(1))], sd, numeric(1))  

## -----------------------------------------------------------------------------
library(Rcpp)
sourceCpp('../src/StatCompC.cpp')
N=5000;mu1=0;mu2=0;sigma1=1;sigma2=1;rho=0.9
xy_cpp<-cppgibbs(N,mu1,mu2,sigma1,sigma2,rho)
cor(xy_cpp[,1],xy_cpp[,2])

## -----------------------------------------------------------------------------
N=5000;mu1=0;mu2=0;sigma1=1;sigma2=1;rho=0.9
gibbsR<-function(N=5000,mu1=0,mu2=0,sigma1=1,sigma2=1,rho=0.9){
  xy<-matrix(0,N,2)
  xy[1,]<-c(mu1,mu2)
  s1 <- sqrt(1-rho^2)*sigma1
  s2 <- sqrt(1-rho^2)*sigma2
  xy[1,]<-c(mu1,mu2)
  for (i in 2:N) {
    x2 <- xy[i-1, 2]
    m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
    xy[i, 1] <- rnorm(1, m1, s1)
    x1 <- xy[i, 1]
    m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
    xy[i, 2] <- rnorm(1, m2, s2)
  }
  return(xy)
}
xy<-gibbsR()
cor(xy[,1],xy[,2])
qqplot(xy_cpp,xy)

