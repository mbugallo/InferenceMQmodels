
### MEDIAN ABSOLUTE DEVIATION

fun.MAD <- function(x){ median(abs( x-median(x)) )/0.6745 }

### HUBER FUNCTIONS

hub.psi <- function(x, k){ ifelse(abs(x) <= k, x, sign(x) * k) }
der.hub.psi <- function(x, k){ ifelse(abs(x) <= k, 1, 0) }

### MQ MODELS (MQ and BMQ prediction)

MQpred.mse <- function(mqo, data.s, x.r, sd, mod.SAE, Nd, nd, 
                       regioncode.s, regioncode.r, grid){
  
  n <- dim(data.s)[1]
  p <- dim(cbind(1,x.s))[2]
  D <- length(unique(data.s$regioncode.s))
  
  MQ<-BMQ<-pred.s<-pred.r<-res.s<-mse1<-mse2<-bias1<-bias2<-fun.dif<-phi<-H<-list()
  k.d<-var.k.d<-rep(NA, D)
  
  for (d in 1:D){
    xd.r.sum<-colSums(cbind(1,x.r)[regioncode.r==d,])
    
    weig.d<-diag(mod.SAE$q.weights[,d])
    H[[d]] <- weig.d%*%cbind(1,data.s$x.s)%*%solve(t(cbind(1,data.s$x.s))%*%weig.d%*%cbind(1,data.s$x.s))
    uu.d <- H[[d]]%*%xd.r.sum
    ww.d<-(uu.d+(regioncode.s==d))
    
    H[[d]] <- cbind(1,data.s$x.s)%*%t(H[[d]])
    
    pred.s[[d]]<-cbind(1,data.s$x.s)[regioncode.s==d,]%*%mod.SAE$coef[,d]
    pred.r[[d]]<-cbind(1,x.r)[regioncode.r==d,]%*%mod.SAE$coef[,d]
    
    # MQ and BMQ predictors
    MQ[[d]]<-c(ww.d)%*%y.s/Nd[d]
    
    suma.bias=0
    for (d.k in 1:D){
      index <- (regioncode.s==d.k)
      suma.bias=suma.bias+sum(ww.d[index]*cbind(1,data.s$x.s)[index,]%*%mod.SAE$coef[,d.k])
    }
  
    bias1[[d]]<-(suma.bias-(sum(pred.s[[d]])+sum(pred.r[[d]])))/Nd[d]
    
    res.s[[d]]<-y.s[regioncode.s==d]- pred.s[[d]]
    Nn.rob<-(Nd[d]-nd[d])*sd[d]/(Nd[d]*nd[d])
  
    for (i in 1:length(grid)){
      phi[[i]]<-Nn.rob*sum(hub.psi(res.s[[d]]/sd[d], k=grid[i]))
      fun.dif[[i]]<-Nn.rob^2*sum(hub.psi(res.s[[d]]/sd[d], k=grid[i])^2) + (bias1[[d]] + phi[[i]])^2
    }
    
    k.d[d]<-grid[which.min(unlist(fun.dif))]
    i<-(which(grid==k.d[d]))
    
    BMQ[[d]] <-MQ[[d]]+phi[[i]]
    bias2[[d]]<-bias1[[d]]+phi[[i]]
  }
  
  return(list(MQ=unlist(MQ), BMQ=unlist(BMQ), bias1=unlist(bias1), 
              bias2=unlist(bias2), k.d=k.d, H=H))
}


###  A Working Likelihood for M-quantiles: The Generalised Asymmetric 
### Least Informative (GALI) distribution

rho <- function(x, k){ 2*ifelse(abs(x) <= k, x^2/2, abs(x)*k-k^2/2) }
rho.q <- function(u, q){ ifelse(u < 0, (q-1)*u, q*u) * rho(u, k=1.345) }

### Option 1: Algorithm 1.1	

f.c.d.a1 <- function(res, sigma.d, nd, R, grid.k){
  
  c.d.boot <- rep(NA, R)	
  
  # Bootstrap process
  for (r in 1:R){
    boot.sample <- sample(res/sigma.d, nd, replace=TRUE)
    
    # Optimization of k
    c.prob<-rep(NA,length(grid.k))
    for (i in 1:length(grid.k)){
      c.prob[i] <- sum(hub.psi(boot.sample,k=grid.k[i])^2)+ 
        (sum(hub.psi(boot.sample,k=grid.k[i]))-sum(boot.sample))^2
    }
    c.d.boot[r]<-grid.k[which.min(c.prob)]      
  } 	
  return(c.d.boot)	
}

### Option 2: Algorithm 1.2
  
f.c.d.a2 <- function(mu.d.vec, sigma.d, theta.d, n, nd, d, R, grid.k, H){
	
  GALI.dist.emp <- list()
  for(l in 1:n){
    mu.d <- mu.d.vec[l]
    grid <- seq(mu.d-100, mu.d+100,by=.005)	
    tol <- max(50*sigma.d, 50)
    
    # Preliminary functions
    integrando <- function(u, mu.d, sigma.d, theta.d){ 1/sigma.d*exp(-rho.q((u-mu.d)/sigma.d, theta.d)) }
    B <- integrate(integrando, mu.d, sigma.d=sigma.d, theta.d=theta.d, lower = mu.d-tol, upper = mu.d+tol)[[1]] 
    integrando <- function(u, mu.d, sigma.d, theta.d, B){ exp(-rho.q((u-mu.d)/sigma.d, theta.d))/(sigma.d*B) }
    
    # GALI.distribucion
    GALI.distribucion <- function(x, mu.d, sigma.d, theta.d, B){ 
      integrate(integrando, mu.d=mu.d, sigma.d=sigma.d, theta.d=theta.d,  B=B, lower = mu.d-tol, upper = x)[[1]] }
    
    GALI.dist.emp[[l]] <- sapply(grid,  GALI.distribucion, mu.d=mu.d, sigma.d=sigma.d, theta.d=theta.d, B=B)	
    
  }
  c.d.boot <- rep(NA, R)	
  
  # Bootstrap process
  for (r in 1:R){
    boot.sample <- ((diag(1,n,n)-H)%*%sapply(GALI.dist.emp, 
                   function(x){ grid[which.min(abs(runif(1)-x))] }))
    
    # boot.sample <- boot.sample[(nd*(d-1)+1):(nd*d)]/sigma.d
    boot.sample <- boot.sample[(nd*(d-1)+1):(nd*d)]/fun.MAD(boot.sample)
    
    # Optimization of k
    c.prob<-rep(NA,length(grid.k))
    for (i in 1:length(grid.k)){
      c.prob[i] <- sum(hub.psi(boot.sample,k=grid.k[i])^2)+ 
        (sum(hub.psi(boot.sample,k=grid.k[i]))-sum(boot.sample))^2
    }
    c.d.boot[r]<-grid.k[which.min(c.prob)]      
  } 	
  return(c.d.boot)	
}

### Option 3: Algorithm 1.3
f.c.d.a3 <- function(theta.d, n, grid, R, grid.k){
  
  # Preliminary functions
  integrando <- function(u, theta.d){ exp(-rho.q(u, theta.d)) }
  B <- integrate(integrando, theta.d=theta.d, lower = -Inf, upper = Inf)[[1]] 
  integrando <- function(u, theta.d, B){ exp(-rho.q(u, theta.d))/B }
  
  # GALI01.distribucion
  GALI01.distribucion <- function(x, theta.d, B){ integrate(integrando, theta.d=theta.d, 
                                                            B=B, lower = -Inf, upper = x)[[1]] }
  
  GALI01.dist.emp <- sapply(grid,  GALI01.distribucion, theta.d=theta.d, B=B)	
  
  c.d.boot <- rep(NA, R)	
  
  # Bootstrap process
  for (r in 1:R){
    boot.sample <- sapply(runif(n), function(x){ grid[which.min(abs(x-GALI01.dist.emp))] })
    
    # Optimization of k
    c.prob<-rep(NA,length(grid.k))
    for (i in 1:length(grid.k)){
      c.prob[i] <- sum(hub.psi(boot.sample,k=grid.k[i])^2)+ 
        (sum(hub.psi(boot.sample,k=grid.k[i]))-sum(boot.sample))^2
    }
    c.d.boot[r]<-grid.k[which.min(c.prob)]      
  } 	
  return(c.d.boot)	
}

