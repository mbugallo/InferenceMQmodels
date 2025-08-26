

###  A Working Likelihood for M-quantiles: The Generalised Asymmetric 
### Least Informative (GALI) distribution

source('AuxFunctions.R')

# Initial parameters

theta.d <- 0.5
n <- 10
grid <- seq(-10,10,by=.001)
R <- 2000
grid.k <- seq(0,10,by=.001)

f.mean <- function(theta.d,  grid){
	
	# Preliminary functions
	integrando <- function(u, theta.d){ exp(-rho.q(u, theta.d)) }
	B <- integrate(integrando, theta.d=theta.d, lower = -Inf, upper = Inf)[[1]] 
	integrando <- function(u, theta.d, B){ exp(-rho.q(u, theta.d))/B }
	
	# GALI01.distribucion
	GALI01.distribucion <- function(x, theta.d, B){ integrate(integrando, theta.d=theta.d, 
						B=B, lower = -Inf, upper = x)[[1]] }
						
	# Cumulative distribution function
	# plot(grid, sapply(grid, GALI01.distribucion, theta.d=theta.d, B=B), ylim=c(0,1))
	# abline(h=0); abline(h=1)

	# Density function
	# plot(grid, sapply(grid, integrando, theta.d=theta.d, B=B))
	# abline(v=0)
						
	media <- mean(grid*sapply(grid, integrando, theta.d=theta.d, B=B))
	return(media)
}

media.res <- sapply(seq(0.01, 0.99, by=0.01),  f.mean, grid=grid)

par(mar=c(4, 5.5, 4, 4), xpd=F)
plot(seq(0.01, 0.99, by=0.01), media.res, xlab=expression(hat(theta)[d]), lwd=3,
			ylab=expression(paste('E[',e[psi][d][j], ']')), cex.lab=2, 
			cex.axis=2, cex.main=2, type='l')
abline(h=0, col='red', lwd=2, lty=2)
points(0.5, 0, pch=20, cex=2.5)

###

result.n.5 <- lapply(c(0.5,0.75,0.9,0.95), f.c.d.a3, n=5, grid=grid, R=R, grid.k=grid.k)
result.n.10 <- lapply(c(0.5,0.75,0.9,0.95), f.c.d.a3, n=10, grid=grid, R=R, grid.k=grid.k)
result.n.15 <- lapply(c(0.5,0.75,0.9,0.95), f.c.d.a3, n=15, grid=grid, R=R, grid.k=grid.k)
result.n.25 <- lapply(c(0.5,0.75,0.9,0.95), f.c.d.a3, n=25, grid=grid, R=R, grid.k=grid.k)
result.n.50 <- lapply(c(0.5,0.75,0.9,0.95), f.c.d.a3, n=50, grid=grid, R=R, grid.k=grid.k)
result.n.100 <- lapply(c(0.5,0.75,0.9,0.95), f.c.d.a3, n=100, grid=grid, R=R, grid.k=grid.k)
result.n.1000 <- lapply(c(0.5,0.75,0.9,0.95), f.c.d.a3, n=1000, grid=grid, R=R, grid.k=grid.k)

round(sapply(result.n.5, median),3)

df <- cbind(result.n.5[[1]], result.n.10[[1]], result.n.15[[1]], result.n.25[[1]], 
		result.n.50[[1]], result.n.100[[1]], result.n.1000[[1]])
colnames(df) <- c('5', '10', '15', '25', '50', '100', '1000')		

par(mar=c(4, 5.2, 4, 0), xpd=F)
boxplot(df,ylim=c(0,5.4), xlab=expression(n[d]), ylab=expression(hat(c)[phi][d]),
		main=expression(paste(hat(theta)[d],'=0.50')), cex.lab=2, cex.axis=2, cex.main=2)

# save.image('ResultsTest/Results_for_n.Rdata')

####

c.d.boot <- lapply(seq(0.01, 0.99, by=0.01),  f.c.d.a3, n=n, grid=grid, R=R, grid.k=grid.k)

# save.image('ResultsTest/Results_for_theta_d.Rdata')

par(mar=c(4, 5.5, 4, 4), xpd=F)
plot(seq(0.01, 0.99, by=0.01), sapply(c.d.boot, mean), xlab=expression(hat(theta)[d]), lwd=3,
			ylab=expression(hat(c)[phi][d]), cex.lab=2, cex.axis=2, cex.main=2, type='l', ylim=c(0,5))
abline(h=min(sapply(c.d.boot,mean)), col='red', lwd=2, lty=2)	
points(0.5, min(sapply(c.d.boot,mean)), pch=20, cex=2.5)	

index <- which(seq(0.01, 0.99, by=0.01)==0.50)
par(mar=c(4, 4.5, 4, 4), xpd=F)
plot(density(c.d.boot[[index]]), type='l', lwd=3, ylab='Density', xlim=c(0,4),
	main=expression(paste(hat(theta)[d],'=0.50 and ', n[d], '=10')), ylim=c(0,2.1),
	xlab=expression(hat(c)[phi][d]), cex.lab=2, cex.axis=2, cex.main=2)	

abline(v=median(c.d.boot[[index]]), col='blue', lwd=3, lty=5)
abline(v=quantile(c.d.boot[[index]],0.95), col='red', lwd=4, lty=3)

q95 <- quantile(c.d.boot[[index]],0.95); q95 # 1.18905

# legend("topright", legend = c(expression(q[0.5]), expression(q[0.95])),
#  	col=c('blue', 'red'), lwd=c(3,4), lty=c(5,3), cex = 2,  horiz=F)


####
theta.d<-0.50
mu.d <- 0
grid <- seq(mu.d-100, mu.d+100,by=.005)	

GALI.dist.emp <- list(); l <- 0
for(sigma.d in c(0.5, 1, 5, 10, 15)){
	l <- l+1
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

par(mar=c(4, 4.5, 4, 4), xpd=F)
plot(grid, GALI.dist.emp[[1]], xlim=c(-100,50), type='l', lwd=2, cex.lab=2, cex.axis=2, ylim=c(0,1),
		cex.main=2, xlab=' ', ylab='Distribution', main='q=0.50')
lines(sort(grid),sort(GALI.dist.emp[[2]]), lwd=2)
lines(sort(grid),sort(GALI.dist.emp[[3]]), lwd=2)
lines(sort(grid),sort(GALI.dist.emp[[4]]), lwd=2)
lines(sort(grid),sort(GALI.dist.emp[[5]]), lwd=2)

# legend("topleft", legend = c(expression(paste(sigma[q], '=0.5')), expression(paste(sigma[q], '=1')), 
			# expression(paste(sigma[q], '=5')), expression(paste(sigma[q], '=10')), expression(paste(sigma[q], '=15'))), 
      		# bty = "n", cex = 2,  horiz=F)



	