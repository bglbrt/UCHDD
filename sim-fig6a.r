################################################################################
# LIBRARIES
library(MASS)
library(glmnet)
library(corpcor)

EXPORT_PATH = "PATH/TO/EXPORT/FOLDER/"

################################################################################
# PART III ESSAY - UNMEASURED CONFOUNDING IN HIGH-DIMENSIONAL DATA
# - SIMULATIONS
# - FILE N.5
################################################################################

################################################################################
# DATA

# generates data from the linear confounding model with:
# - n the number of samples 
# - p the number of covariates
# - r the number of confounders
# - s the number of non-zero coefficients of beta
# - b the value of the non-zero coefficients of beta
# - sigma_e the scalar for for the covariance matrix of the rows of E
# - kappa the scalar for the covariance matrix of the rows of Gamma
generate_data <- function(n=300, p=400, r=6, s=5, b=1, sigma_e=1){
  
  # finding min(n, p)
  m = min(n, p)
  
  # generate the coefficient vectors beta and delta
  # - where beta is sparse with the first s entries non null
  # - and delta is sampled from a N(d, 1) distribution
  beta = c(rep(b, s), rep(0, p-s))
  delta = rnorm(r, mean=0, sd=1)
  
  # sample the rows of gamma from a N(O, kappa*Ip) distribution
  gamma = matrix(rnorm(r*p, mean=0, sd=1), r, p)
  
  # sample nu from a N(0, 1) distribution
  nu = rnorm(n, mean=0, sd=1)
  
  # sample E from a N(0, sigma_e*Ip) distribution
  E = matrix(mvrnorm(n, mu=rep(0, p), Sigma = (sigma_e**2)*diag(p)), n, p)
    
  # sample Z from a N(0, Ir) distribution
  Z = mvrnorm(n, mu=rep(0, r), Sigma = diag(r))
  
  # compute X and Y as in the linear confounding model
  X = Z %*% gamma + E
  Y = X %*% beta + Z %*% delta + nu
  
  # output dataset as a list
  data <- list()
  # - with all constants
  data$n <- n
  data$p <- p
  data$r <- r
  data$s <- s
  data$m <- m
  data$b <- b
  data$sigma_e <- sigma_e
  # - all parameters
  data$beta <- beta
  data$delta <- delta
  data$gamma <- gamma
  # - all errors
  data$nu <- nu
  data$E <- E
  # - and X, Z, Y
  data$X <- X
  data$Z <- Z
  data$Y <- Y
    
  return(data)
}

# example of data generation
data <- generate_data()

################################################################################
# SPECTRAL TRANSFORMATIONS

# outputs the transformed singular values for computing the Lava transform
d_lava <- function(di, l, n) {
  return(sqrt(n*l*(di**2)/(n*l + (di**2))))
}

############################################
# SPECTRAL TRANSFORMS (RETURNS THE F MATRIX)
# METHODS AVAILABLE :
# - LASSO
# - Oracle PCA
# - Trim
# - Puffer
# - Lava

# outputs F the spectral transform given:
# - X the original design matrix
# - the method of computation to choose from:
#   - LASSO (transform is equal to UU^T)
#   - Oracle PCA (PCA adjustment where the number of confounders r is known)
#   - Trim (Trim transform)
#   - Puffer (Puffer transform)
#   - Lava (Lava transform)
# - r the number of confounders
# - l the tuning parameter for the Lava transform
# - svd the singular value decomposition of X (svd(X)) if already computed
spectral_transform <- function(X, method, r, l=NULL, svd=NULL) {
  
  # finds n, p and m
  n <- nrow(X)
  p <- ncol(X)
  m <- min(n, p)
  
  # computes the SVD of X if not already given
  if (is.null(svd)) {
    X.svd <- svd(X)
  } else {
    X.svd <- svd
  }
  
  # decomposes the SVD of X
  d <- X.svd$d
  D <- diag(X.svd$d)
  U <- X.svd$u
  V <- X.svd$v
  
  # computes F the spectral transform depending on the chosen method
  if (method=="LASSO") {
    # LASSO
    F_spectral_transform <- U %*% diag(rep(1, m)) %*% t(U)
    
  } else if (method=="Oracle PCA") {
    # Oracle PCA
    d_tilde <- c(rep(0, r), tail(d, m-r))
    D_tilde = diag(d_tilde / d)
    F_spectral_transform <- U %*% D_tilde %*% t(U)
    
  } else if (method == "Trim") {
    # Trim
    tau = d[floor(m/2)]
    d_tilde <- c(rep(tau, floor(m/2)), tail(d, m-floor(m/2)))
    D_tilde = diag(d_tilde / d)
    F_spectral_transform <- U %*% D_tilde %*% t(U)
    
  } else if (method == "Puffer") {
    # Puffer
    dmin <- min(d)
    d_tilde = rep(dmin, m)
    D_tilde <- diag(d_tilde / d)
    F_spectral_transform <- U %*% D_tilde %*% t(U)
    
  } else if (method == "Lava") {
    # Lava
    if (is.null(l)){
      l = (1/n)*(d[floor(m/2)])**2
    }
    
    d_tilde <- unlist(lapply(d, FUN=d_lava, l=l, n=n))
    D_tilde <- diag(d_tilde / d)
    F_spectral_transform <- U %*% D_tilde %*% t(U)
    
  } else {
    print("ERROR: Please provide a valid method!")
    
  }
  
  return(F_spectral_transform)
}

# example of transformation matrix generation
F <- spectral_transform(data$X, "Trim", 5)

################################################################################
# SIMULATION

# generates data from the linear confounding model with:
# - n the number of samples 
# - p the number of covariates
# - r the number of confounders
# - s the number of non-zero coefficients of beta
# - b the value of the non-zero coefficients of beta
# - sigma_e the scalar for for the covariance matrix of the rows of E
# - kappa the scalar for the covariance matrix of the rows of Gamma
generate_data_mod1 <- function(n=300, p=400, r=6, s=5, b=1, sigma_e=1, theta=1){
  
  # finding min(n, p)
  m = min(n, p)
  
  # generate the coefficient vectors beta and delta
  # - where beta is sparse with the first s entries non null
  # - and delta is sampled from a N(d, 1) distribution
  beta = c(rep(b, s), rep(0, p-s))
  delta = rnorm(r, mean=0, sd=1)
  
  # sample the rows of gamma from a N(O, kappa*Ip) distribution
  gamma = matrix(rnorm(r*p, mean=0, sd=1), r, p)
  gamma = (1/sqrt(sum(gamma^2))) * gamma
  gamma = sqrt(theta) * gamma
  
  # sample nu from a N(0, 1) distribution
  nu = rnorm(n, mean=0, sd=1)
  
  # sample E from a N(0, sigma_e*Ip) distribution
  E = matrix(mvrnorm(n, mu=rep(0, p), Sigma = (sigma_e**2)*diag(p)), n, p)
  
  # sample Z from a N(0, Ir) distribution
  Z = mvrnorm(n, mu=rep(0, r), Sigma = diag(r))
  
  # compute X and Y as in the linear confounding model
  X = Z %*% gamma + E
  Y = X %*% beta + Z %*% delta + nu
  
  # output dataset as a list
  data <- list()
  # - with all constants
  data$n <- n
  data$p <- p
  data$r <- r
  data$s <- s
  data$m <- m
  data$b <- b
  data$sigma_e <- sigma_e
  data$kappa <- kappa
  # - all parameters
  data$beta <- beta
  data$delta <- delta
  data$gamma <- gamma
  # - all errors
  data$nu <- nu
  data$E <- E
  # - and X, Z, Y
  data$X <- X
  data$Z <- Z
  data$Y <- Y
  
  return(data)
}

# constants
p = 600
n = 300
r = 1
s = 5
b = 1
sigma_e = 1

# list of errors
lasso_errors <- c()
pca_errors <- c()
trim_errors <- c()

# number of iterations per case
h = 1000
span = 10^(seq(-2,3,1/10))

for (theta in span){
  print(paste("Beginning: theta=", as.character(theta), sep=""))
  
  step <- list()
  step["LASSO"] <- 0.0
  step["Oracle PCA"] <- 0.0
  step["Trim"] <- 0.0
  
  for (i in 1:h){
    print(paste("Iteration: ", as.character(i), "/", as.character(h), sep=""))
    
    #data generation
    data <- generate_data_mod1(n=n, p=p, r=r, s=s, b=b, sigma_e=sigma_e, theta=theta)
    data$svd <- fast.svd(data$X)
    
    for (method in c("LASSO", "Oracle PCA", "Trim")){
      
      # transformation matrix computation
      F <- spectral_transform(X=data$X, method=method, r=data$r, svd=data$svd)
      
      # X_tilde and Y_tilde computation
      X_tilde <- F %*% data$X
      Y_tilde <- F %*% data$Y
      
      # lasso regression of X_tilde on Y_tilde
      data$reg <- glmnet(X_tilde, Y_tilde, alpha=1, lambda=sqrt(log(p)/n))
      
      # append the l1 errors
      l1_error <- sum(abs(data$beta - data$reg$beta))
      step[[method]] <- step[[method]] + l1_error
      
    }
    
  }
  
  lasso_errors <- c(lasso_errors, step[["LASSO"]]/h)
  pca_errors <- c(pca_errors, step[["Oracle PCA"]]/h)
  trim_errors <- c(trim_errors, step[["Trim"]]/h)
  
}

# plot
png(file=paste(EXPORT_PATH, "fig6a.png", sep=""), width=500, height=450)
plot(x=span, y=lasso_errors, type="l", lty=5, lwd=2, xaxt="n", log="x", xlab=expression(theta), ylab="L1 estimation error", col="gray", xlim=c(10^-2, 10^3), ylim=c(0, 3))
axis(1, at = 10^(seq(-2,3,1)), las=0)
grid (NULL,NULL, lty = 1,lwd=0.5, col = "lightgray") 
lines(x=span, y=pca_errors, col="seagreen3", lwd=2)
lines(x=span, y=trim_errors, col="royalblue2", lwd=2)
abline(v = sqrt(2), col="tomato3", lwd=2, lty=1)
legend("topright", legend=c("Lasso", "Oracle PCA", "Trim"), col=c("gray", "seagreen3", "royalblue2"), lty=c(2, 1, 1), cex=0.8)
dev.off()