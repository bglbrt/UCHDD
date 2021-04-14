################################################################################
# LIBRARIES
library(MASS)
library(cate)

EXPORT_PATH = "PATH/TO/EXPORT/FOLDER/"
DATA_PATH = "PATH/TO/DATA/FOLDER/"

################################################################################
# PART III ESSAY - UNMEASURED CONFOUNDING IN HIGH-DIMENSIONAL DATA
# - SIMULATIONS
# - FILE N.8
################################################################################

################################################################################
# DATA

# read documentation CSV file
var <- read.csv(paste(DATA_PATH, "raw_data/urb_codlab.csv", sep=""))

# read data CSV file
df <- read.csv(paste(DATA_PATH, "data.csv", sep=""))[,-1]

# set (for each city)
# - X to be:
#   - total nights spent in tourist accommodation establishments per 1000 inhabitants
# - Y to be:
#   - a dataframe of all other covariates
#   - except the ones directly linked to the number of nights spent in tourist accommodation establishments
X_dataframe <- data.frame(X = df$CR2011I)#1000 * (df$CR2001V / df$DE1001V))
Y_dataframe <- df[ , -which(names(df) %in% c("NAME", "CR2001V", "CR2009V", "CR2010I", "CR2011V", "CR2011I", "df$DE1001V"))]

# filter out constant columns in the outcomes matrix
Y_dataframe <- Filter(function(x) sd(x, na.rm = TRUE) != 0, Y_dataframe)

# use matrix representation for Y and X
X <- as.matrix(X_dataframe)
Y <- as.matrix(Y_dataframe)

X_dataframe <- as.data.frame(scale(X_dataframe, scale=TRUE))
Y_dataframe <- as.data.frame(scale(Y_dataframe, scale=TRUE))

X <- scale(X, scale=TRUE)
Y <- scale(Y, scale=TRUE)

################################################################################
# NAIVE ANALYSIS

results <- list()

# naive linear regression of X on each column of Y
results$naive <- cate(~ X | 1, X.data=X_dataframe, Y, r=0, adj.method="naive", calibrate = FALSE)

# remaining covariates when controlling the FWER at 5% error rate
var$LABEL[match(colnames(Y)[which(p.adjust(results$naive$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)], var$CODE)]

# corresponding coefficients
results$naive$beta[which(p.adjust(results$naive$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)]

# corresponding p-values
p.adjust(results$naive$beta.p.value, "bonferroni")[which(p.adjust(results$naive$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)]

# output table
output_naive_1 <- colnames(Y)[which(p.adjust(results$naive$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)]
output_naive_2 <- var$LABEL[match(colnames(Y)[which(p.adjust(results$naive$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)], var$CODE)]
output_naive_3 <- results$naive$beta[which(p.adjust(results$naive$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)]
output_naive_4 <- p.adjust(results$naive$beta.p.value, "bonferroni")[which(p.adjust(results$naive$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)]
output_naive <- data.frame("CODE" = output_naive_1, "Outcome variable" = output_naive_2, "Coefficient" = output_naive_3, "p-value" = output_naive_4)
write.csv(output_naive, file=paste(EXPORT_PATH, "output_naive.csv", sep=""))

################################################################################
# NUMBER OF CONFOUNDERS TO ADJUST FOR

# estimate number of factors to adjust for using BCV
results$r <- est.confounder.num(~ X | 1, X_dataframe, Y, method = "bcv", bcv.plot = FALSE, rmax = 30, nRepeat = 100)$r

# estimate number of factors
results$r

################################################################################
# ADJUSTING FOR CONFOUNDING USING CATE

# robust regression with adjusting for confounding (from CATE package)
results$cate <- cate(~ X | 1, X.data=X_dataframe, Y, r=results$r, fa.method="ml", adj.method="rr", calibrate=FALSE)

# test for confounding
results$cate$alpha.p.value

# remaining covariates when controlling the FWER at 5% error rate
var$LABEL[match(colnames(Y)[which(p.adjust(results$cate$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)], var$CODE)]

# corresponding coefficients
results$cate$beta[which(p.adjust(results$cate$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)]

# corresponding p-values
p.adjust(results$cate$beta.p.value, "bonferroni")[which(p.adjust(results$cate$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)]

# investigate the coefficients in Gamma
var$LABEL[match(names(sort(abs(results$cate$Gamma[,1]))), var$CODE)][1:10]

# output table
output_cate_1 <- colnames(Y)[which(p.adjust(results$cate$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)]
output_cate_2 <- var$LABEL[match(colnames(Y)[which(p.adjust(results$cate$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)], var$CODE)]
output_cate_3 <- results$cate$beta[which(p.adjust(results$cate$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)]
output_cate_4 <- p.adjust(results$cate$beta.p.value, "bonferroni")[which(p.adjust(results$cate$beta.p.value, "bonferroni") < 0.05, arr.ind=TRUE)]
output_cate <- data.frame("CODE" = output_cate_1, "Outcome variable" = output_cate_2, "Coefficient" = output_cate_3, "p-value" = output_cate_4)
write.csv(output_cate, file=paste(EXPORT_PATH, "output_naive.csv", sep=""))
