################################################################################
# LIBRARIES
library(MASS)
library(glmnet)

EXPORT_PATH = "PATH/TO/EXPORT/FOLDER/"
DATA_PATH = "PATH/TO/DATA/FOLDER/"

################################################################################
# PART III ESSAY - UNMEASURED CONFOUNDING IN HIGH-DIMENSIONAL DATA
# - CASE STUDY
# - FILE N.7
################################################################################

################################################################################
# DATA

# read documentation CSV file
var <- read.csv(paste(DATA_PATH, "raw_data/urb_codlab.csv", sep=""))

# read data CSV file
df <- read.csv(paste(DATA_PATH, "data.csv", sep=""))[,-1]

# set (for each city)
# - Y to be:
#   - total nights spent in tourist accommodation establishments per 1000 inhabitants
# - X to be:
#   - a dataframe of all other covariates
#   - except the ones directly linked to the number of nights spent in tourist accommodation establishments
Y_dataframe <- 1000 * (df$CR2001V / df$DE1001V)
X_dataframe <- df[ , -which(names(df) %in% c("NAME", "CR2001V", "CR2009V", "CR2010I", "CR2011V", "CR2011I"))]

# filter out constant columns in the design matrix
X_dataframe <- Filter(function(x) sd(x, na.rm = TRUE) != 0, X_dataframe)

# use matrix representation for Y and X
Y <- as.matrix(Y_dataframe)
X <- as.matrix(X_dataframe)

X <- scale(X, scale=TRUE)
Y <- scale(Y, scale=TRUE)

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

################################################################################
# NAIVE ANALYSIS
# LASSO REGRESSION WITHOUT ADJUSTING FOR CONFOUNDING

# Lasso regression of X on Y
naive_regression <- list()
naive_regression$lasso <- glmnet(X, Y, alpha=1, lambda.min.ratio=1e-4, nlambda=1000)

# plot
plot(naive_regression$lasso, xvar="dev", label=TRUE)
plot(naive_regression$lasso, xvar="lambda", label=TRUE)

# Lasso regression of X on Y with cross validation
naive_regression$lambdamat <- c()
for (i in 1:10){
  print(paste("Iteration: ", as.character(i), "/", as.character(10), sep=""))
  
  naive_regression$cvlasso <- cv.glmnet(X, Y, alpha=1, lambda=naive_regression$lasso$lambda, nfolds=10)
  naive_regression$lambdamat <- c(naive_regression$lambdamat, naive_regression$cvlasso$lambda.min)
}
naive_regression$lambda <- mean(naive_regression$lambdamat)

# results for optimal penalty term
naive_regression$lambda

# Lasso regression of X on Y with optimal penalty term
naive_regression$finallasso <- glmnet(X, Y, alpha=1, lambda=naive_regression$lambda)

# number of non-zero coefficients
sum(naive_regression$finallasso$beta != 0)

# coefficient estimation
naive_regression$finallasso$beta

# list remaining covariates
var$LABEL[var$CODE %in% c(rownames(which(naive_regression$finallasso$beta != 0, arr.ind=TRUE)))]

################################################################################
# ADJUSTING FOR CONFOUNDING
# LASSO REGRESSION WITH ADJUSTMENT USING TRIM SPECTRAL TRANSFORMATION

# visualisation of the singular values of X
png(file=paste(EXPORT_PATH, "fig7a.png", sep=""), width=500, height=450)
plot(svd(X)$d, type="h", lty=1, lwd=1.5, xlab="Ordered singular values", ylab="", col="royalblue2", ylim=c(0, 120))
lines(svd(X)$d, type="p", pch=16, lty=1, lwd=1, cex = 0.6, col="royalblue2")
dev.off()

# transformation matrix computation
F <- spectral_transform(X=X, method="Trim", r=3)

# X_tilde and Y_tilde computation
X_tilde <- F %*% X
Y_tilde <- F %*% Y

png(file=paste(EXPORT_PATH, "fig7b.png", sep=""), width=500, height=450)
plot(svd(X_tilde)$d, type="h", lty=1, lwd=1.5, xlab="Ordered singular values", ylab="", col="royalblue2", ylim=c(0, 10))
lines(svd(X_tilde)$d, type="p", pch=16, lty=1, lwd=1, cex = 0.6, col="royalblue2")
dev.off()

# Lasso regression of X_tilde on Y_tilde
regression <- list()
regression$lasso <- glmnet(X_tilde, Y_tilde, alpha=1, lambda.min.ratio=1e-6, nlambda=1000)

# plot
plot(regression$lasso, xvar="dev", label=TRUE)
plot(regression$lasso, xvar="lambda", label=TRUE)

# Lasso regression of X on Y with cross validation
regression$lambdamat <- c()
for (i in 1:10){
  print(paste("Iteration: ", as.character(i), "/", as.character(10), sep=""))
  
  regression$cvlasso <- cv.glmnet(X_tilde, Y_tilde, alpha=1, lambda=regression$lasso$lambda, nfolds=10)
  regression$lambdamat <- c(regression$lambdamat, regression$cvlasso$lambda.min)
}
regression$lambda <- mean(regression$lambdamat)

# results for optimal penalty term
regression$lambda

# Lasso regression of X on Y with optimal penalty term
regression$finallasso <- glmnet(X_tilde, Y_tilde, alpha=1, lambda=regression$lambda)

# number of non-zero coefficients
sum(regression$finallasso$beta != 0)

# coefficient estimation
regression$finallasso$beta

# list remaining covariates
var$LABEL[var$CODE %in% c(rownames(which(regression$finallasso$beta != 0, arr.ind=TRUE)))]