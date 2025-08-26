
######################################
###    Simulation experiment 1     ###
######################################

rm(list=ls())

start.time <- Sys.time()
set.seed(123)

source('AuxFunctions.R')
source('MQ2.R')

library(MASS)

###########################
### Simulation scenario ###
###########################

B=200 # number of iterations
D=40 # number of areas
R=500 # number of bootstrap resamples
grid.k=seq(0,10,by=.001) # grid for c

Nd=rep(100,D)
N=sum(Nd)

nd=rep(5,D)
n=sum(nd)

regioncode<-rep(1:D,Nd)

# Continuous Part
xdj2<-rlnorm(N,1,0.5)

# Sample and non-sample subsets
s<-NULL
for (d in 1:D){ s<-sort(c(s, sample((1:N)[regioncode==d], nd[d]))) }  

x.s<-xdj2[s]
regioncode.s<-regioncode[s]

x.r<-xdj2[-s]
regioncode.r<-regioncode[-s]

# Quantiles on the grid
tau=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))

# Outliers 
mean.ui<-9
mean.e<-20
ki<-as.integer(.1 * D) 

# Prediction results, bias and MSE estimation
k.d <- mqo.d <- matrix(nrow=B, ncol=D)
res.s <- res.theta.d <- A1.50 <- A2.50 <- A3.50 <- A1.d <- A2.d <- A3.d <- list()

for (b in 1:B){
  set.seed(b)
  print(b)
  
  # Outliers 
  ud<-rnorm(D,0,sqrt(3))
  # % contamination in ui
  out.ui<-1
  uui<-ud
  u1 <- rnorm(ki, mean.ui, sqrt(20))
  uui[(D-ki+1):D]<-u1
  out.ui <- rep(out.ui, D)
  ud <- ifelse(out.ui > 0, uui, ud)
  
  # % contamination in e
  out.e<-0.03
  n1<-rbinom(N,1,out.e)
  edj <- (1-n1)*rnorm(N,0,sqrt(6))+n1*rnorm(N, mean.e, sqrt(150))
  
  ydj <-100+5*xdj2+rep(ud,Nd)+edj
  
  # Sample observations
  y.s<-ydj[s]
  data.s <- data.frame(regioncode.s, y.s, x.s)
  names(data.s) <- c('regioncode.s', 'y.s', 'x.s')
  
  sum.s <- aggregate(data.s$y.s,by=list(data.s$regioncode.s),sum)$x
  
  ############################
  ###   2. Area-level MQ   ###
  ############################
  
  mod<-QRLM(x=cbind(1,x.s), y=y.s, q=tau, maxit=30, k = 1.345)
  qo <- matrix(c(gridfitinter(y.s,mod$fitted.values,mod$q.values)),nrow=n,ncol=1)
  qmat <- data.frame(qo, data.s$regioncode.s)
  mqo <- aggregate(qmat[,1],by=list(qmat[,2]),mean)[,2]
  
  mod.SAE<-QRLM(x=cbind(1,x.s), y=y.s, q=mqo, maxit=30, k = 1.345)
  
  for(d in 1:D){ res.s[[d]] <- y.s-cbind(1,x.s)%*%mod.SAE$coef[,d] }
  sd <- sapply(res.s, fun.MAD)
  
  MQ.results <- MQpred.mse(mqo, data.s, x.r, sd, mod.SAE, Nd, nd, regioncode.s,
  				 regioncode.r, seq(0,10,by=0.01))
  
  mqo.d[b, ] <- mqo
  k.d[b,] <- MQ.results$k.d
  
  #######################################
  ###  Test to detect atypical areas  ###
  #######################################
 
  d.check <- 37

  ###   Bootstrap distribution of c_d(0.5)  
  theta.d <- 0.50
  
  mod.Naive<-QRLM(x=cbind(1,x.s), y=y.s, q=theta.d, maxit=30, k = 1.345)
  weig <- diag(as.numeric(mod.Naive$q.weights))
  H <- cbind(1,data.s$x.s)%*%t(weig%*%cbind(1,data.s$x.s)%*%
       solve(t(cbind(1,data.s$x.s))%*%weig%*%cbind(1,data.s$x.s)))
  
  for(d in 1:D){  res.theta.d[[d]] <- y.s[regioncode.s==d]-cbind(1,x.s)[regioncode.s==d, ]%*%mod.Naive$coef }
  sigma.theta.d <- fun.MAD(unlist(res.theta.d))
  
  A1.50[[b]] <- f.c.d.a1(res=res.theta.d[[d.check]], sigma.d=sigma.theta.d, nd[d.check], R, grid.k)
  
  A2.50[[b]] <- f.c.d.a2(mu.d.vec=unlist(res.theta.d)+y.s, sigma.d=sigma.theta.d, 
               theta.d=theta.d, n=n, nd=nd[d.check], d.check, R, grid.k, H=H)
  
  A3.50[[b]] <- f.c.d.a3(theta.d=theta.d, nd[d.check], grid=seq(-50,50,by=.005), R, grid.k)

  ###   Bootstrap distribution of c_d(theta_d) 
  
  A1.d[[b]] <- f.c.d.a1(res=res.s[[d.check]], sigma.d=sd[d.check], nd[d.check], R, grid.k)
  
  A2.d[[b]] <- f.c.d.a2(mu.d.vec=unlist(res.s[[d.check]])+y.s, sigma.d=sd[d.check], 
               theta.d=mqo[d.check], n=n, nd=nd[d.check], d.check, R, grid.k, H=MQ.results$H[[d.check]])
  
  A3.d[[b]] <- f.c.d.a3(theta.d=mqo[d.check], nd[d.check], grid=seq(-50,50,by=.005), R, grid.k)
}  


(z1 <- mean(k.d[,d.check] > sapply(A1.50, quantile, 0.95)))
(z2 <- mean(k.d[,d.check] > sapply(A2.50, quantile, 0.95)))
(z3 <- mean(k.d[,d.check] > sapply(A3.50, quantile, 0.95)))

mean(sapply(A1.50, quantile, 0.95))
mean(sapply(A2.50, quantile, 0.95))
mean(sapply(A3.50, quantile, 0.95))

# save.image(paste0("ResultsTest/Results_for_[e,u]_d=", d.check, "B=", B, ".Rdata"))

# Figure 1
par(mar=c(4, 4.5, 4, 4), xpd=F)
plot(density(sapply(A1.d,unlist)), type='l', lwd=3, ylab='Density', xlim=c(0,4),
     main='[0,u]', ylim=c(0,2.6),
     xlab=expression(hat(c)[phi][d]), cex.lab=2, cex.axis=2, cex.main=2)	
lines(density(sapply(A2.d,unlist)), lwd=3, col='red')
lines(density(sapply(A3.d,unlist)), lwd=3, col='blue')

lines(density(sapply(A1.d,unlist), bw=0.1), lwd=3, lty=3)
lines(density(sapply(A2.d,unlist),bw=0.2), lwd=3, lty=3, col='red')
lines(density(sapply(A3.d,unlist), bw=0.1), lwd=3, lty=3, col='blue')

legend("topright", legend = c("Alg. 1.1", "Alg. 1.2", "Alg. 1.3"), 
       col=c('black', 'red', 'blue'), bty = "n", cex = 2,  lwd=3, horiz=F) 


# Figure 2
par(mar=c(4, 5.4, 4, 4), xpd=F)
plot(k.d[,d.check], ylim=c(0,4.3), pch=20, xlab='', ylab=expression(hat(c)[phi][d]),
     cex.lab=2, cex.axis=2, cex.main=2, xaxt='n', main='[0,u]    d=37')
lines(sapply(A1.d, quantile, 0.95), lwd=1.5, col='black')
lines(sapply(A2.50, quantile, 0.95), lwd=2, col='red')
lines(sapply(A3.50, quantile, 0.95), lwd=2, col='blue')
