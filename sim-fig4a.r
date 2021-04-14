################################################################################
# LIBRARIES
library(MASS)

EXPORT_PATH = "PATH/TO/EXPORT/FOLDER/"

################################################################################
# PART III ESSAY - UNMEASURED CONFOUNDING IN HIGH-DIMENSIONAL DATA
# - SIMULATIONS
# - FILE N.1
################################################################################

################################################################################
# SIMULATION

# constants
n = 300
theta = 8

# list of correlations
cor_gamma <- list()
cor_b <- list()

# number of iterations per case
h = 1000
span = seq(50, 1500, by=50)

for (p in span){
  print(paste("Beginning: p=", as.character(p), sep=""))
  
  for (i in 1:h){
    print(paste("- Iteration: ", as.character(i), "/", as.character(h), sep=""))
    
    # sample gamma from a N(O, 1) distribution
    gamma = matrix(rnorm(p, mean=0, sd=1), p, 1)
    gamma = (1/sqrt(sum(gamma^2))) * gamma
    gamma = sqrt(theta) * gamma
    
    # compute sigma
    sigma = (gamma %*% t(gamma)) + diag(p)
  
    # sample delta from a N(O, 1) distribution
    delta = rnorm(1, mean=0, sd=1)
  
    # compute b
    b = delta * (solve(sigma) %*% gamma)
  
    # sample X from a N(0, sigma)
    X = matrix(mvrnorm(n, mu=rep(0, p), Sigma=sigma), n, p)
  
    # compute correlation between:
    # - gamma and the first right singular vector of X
    # - beta and the first right singular vector of X
    cor_gamma[[as.character(p)]][[as.character(i)]] <- c(cor_gamma[[as.character(p)]][[as.character(i)]], cor(gamma, svd(X)$v[, 1]))
    cor_b[[as.character(p)]][[as.character(i)]] <- c(cor_b[[as.character(p)]][[as.character(i)]], cor(b, svd(X)$v[, 1]))
  
  }
  
}

# plot
png(file=paste(EXPORT_PATH, "fig4a.png", sep=""), width=500, height=450)
plot(x=rep(50, h), y=abs(unlist(cor_gamma$"50")), type="p", lwd=2, xlab="Number of covariates (p)", ylab="Correlation in absolute value", col="lightgray", xlim=c(0, 1500), ylim=c(0, 1))
counter = 1
for (p in span[-1]){
  counter = counter + 1
  lines(x=rep(p, h), y=abs(unlist(cor_gamma[[counter]])), col="lightgray", type="p", lwd=2)
}
lines(x=span, y=rowMeans(t(matrix(abs(unlist(cor_gamma)), h, length(span)))), col="royalblue1", type="o", lwd=2)
grid (NULL,NULL, lty = 1,lwd=0.5, col = "lightgray") 
dev.off()
