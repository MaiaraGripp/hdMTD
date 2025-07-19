## Ensure needed packages are downloaded
pkgs <- c("hdMTD")
install.packages(pkgs[!pkgs %in% installed.packages()])

require(hdMTD)

####################################
## Section 5: Using hdMTD
## 5.1 Data generation
####################################

## Generating Model

set.seed(11)
Lambda <- c(1, 15, 30)
A <- c(0, 1)
lam0 <- 0.01
lamj <- c(0.39, 0.3, 0.3)
p0 <- c(0.5, 0.5)
MTD <- MTDmodel(Lambda = Lambda, A = A, lam0 = lam0, lamj = lamj, p0 = p0)
MTD

## Sampling from invariant distribution

X <- perfectSample(MTD, N = 1000)

####################################
## 5.2 Estimation
####################################

## Estimation methods of Lambda

hdMTD_FS(X, d = 40, l = 4)
# Output: 30 15  1 27
# Equivalent: hdMTD(X, d = 40, method = "FS", l=4)

# WARNING: hdMTD_BIC with d = 40 bellow ~30min (on i7-1255U, 10 cores)
hdMTD_BIC(X, d = 40, minl = 4, maxl = 4)
# 1 15 17 30
# Equivalent: hdMTD(X, d = 40, method = "BIC", minl = 4, maxl = 4)

hdMTD_BIC(X, d = 40,
          S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40),
          minl = 4, maxl = 4)
# 1 15 17 30
# Equivalent: hdMTD(X,d=40,method="BIC",S=c(1,5,10,15,17,20,27,30,35,40),minl=4,maxl=4)

hdMTD_BIC(X, d = 40,
          S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40),
          minl = 1, maxl = 4)
# 30

hdMTD_BIC(X, d = 40,
          S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40),
          minl = 1, maxl = 4,
          byl = TRUE, BICvalue = TRUE)
#       30        15,30      1,15,30   1,15,17,30 smallest: 30 
# 644.4959     648.0111     649.4950     650.2869     644.4959 

hdMTD_BIC(X, d = 40,
          S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40),
          minl = 1, maxl = 4,
          byl = TRUE, BICvalue = TRUE,
          xi = 0.4)
#       30      15,30    1,15,30   1,15,17,30  smallest: 1,15,17,30
# 641.7328   643.1757   642.5873     641.3069              641.3069  

hdMTD_BIC(X, d = 40,
          S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40),
          minl = 1, maxl = 4,
          byl = TRUE, BICvalue = TRUE,
          single_matrix = TRUE, indep_part = FALSE)
#       30      15,30    1,15,30   1,15,17,30  smallest: 1,15,17,30 
# 637.5881   634.1956   628.7718     622.6559              622.6559


# WARNING: hdMTD_CUT bellow ~2.5min (on i7-1255U, 10 cores)
hdMTD_CUT(X, d = 40, S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40))
# 40 35 30 27 20 17 15 10  5  1
# Equivalent: hdMTD(X,d=40,S=c(1,5,10,15,17,20,27,30,35,40),method="CUT")

# WARNING: hdMTD_CUT bellow ~2.5min (on i7-1255U, 10 cores)
hdMTD_CUT(X, d = 40, S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40), alpha = 0.13)
# 35 27  5  1

hdMTD_CUT(X, d = 40, S = c(1, 5, 17, 27, 30, 35), alpha = 0.13)
# 30 17  5  1

hdMTD_FSC(X, d = 40, l = 4, alpha = 0.1)
# 30 24
# Equivalent: hdMTD(X,d=40,method="FSC",l=4,alpha = 0.1)

hdMTD_FS(X[1:500], d = 40, l = 4)
# 11 30  7 24

## Estimation of transition probabilities

probs(X, S = c(1, 15, 30))
probs(X, S = c(1, 15, 30), matrixform = TRUE)


## Estimating the oscillations

oscillation(MTD)
#        -1       -15       -30
# 0.1233648 0.1017517 0.1846987

oscillation(X, S = c(1, 15, 30))
#        -1       -15       -30
# 0.1076339 0.1166363 0.1675360


## Estimating MTD parameters through the EM algorithm

# Initial parameters for EM method
init <- list(
  'lambdas'= c(0.01, 0.33, 0.33, 0.33),
  'p0' = c(0.5, 0.5),
  'pj' = rep(list(matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2, nrow = 2)), 3)
)

MTDest(X, S = c(1, 15, 30), init = init, iter = TRUE)

MTDest(X, S = c(1, 15, 30), M = NULL, nIter = 9, init = init, oscillations = TRUE)


## Estimating P (transition probabilities)

estParam <- MTDest(X, S = c(1, 15, 30), init = init)

estMTD <- MTDmodel(Lambda, A, lam0 = estParam$lambdas[1],
         lamj = estParam$lambdas[-1], p0 = estParam$p0,
         pj = estParam$pj)

estMTD$P


