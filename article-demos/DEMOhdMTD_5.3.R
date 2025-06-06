## Ensure needed packages are downloaded
pkgs <- c("hdMTD", "dplyr", "ggplot2")
install.packages(pkgs[!pkgs %in% installed.packages()])

require(hdMTD)
require(dplyr)
require(ggplot2)

####################################
## Section 5: Using hdMTD
## 5.3 Testing hdMTD
####################################

set.seed(123)
# True model specification
Lambda <- c(1, 5)
A <- c(0, 1)
lam0 <- 0.01
p0 <- c(0.5, 0.5)
MTD <- MTDmodel(Lambda, A, lam0, p0 = p0) # Generates an MTD model

# Simulation parameters:
n <- 100      # Number of replications
N <- 10000    # Full sample size
m <- c(1000, 1500, 2000, 2500, 3000, 5000, 10000) # Subsample sizes
d <- 100      # Max order for FS and Oracle
dNaive <- 5   # Max order for Naive
pairList <- t(combn(100, 2)) # All possible pairs with digits from 1 to 100
npairs <- nrow(pairList)

# Initialize results storage
FSP <- matrix(0, ncol = length(m),nrow = n)
FS <- matrix(0, ncol = length(m),nrow = n)
NaiveP <- matrix(0, ncol = length(m), nrow = n)
Naive <- matrix(0, ncol = length(m), nrow = n)
OracleP <- matrix(0, ncol = length(m),nrow = n)
Oracle <- matrix(0, ncol = length(m), nrow = n)
SFS <- matrix(0, ncol = length(m)*2, nrow = n)
ZOracle <- matrix(0, ncol = length(m)*2, nrow = n)

# WARNING: this next loop takes 5 to 6 DAYS to complete
for (i in seq_len(n)) {
    X <- perfectSample(MTD, N = N) # Generate full sample

    for (k in seq_along(m)) {
        Y <- X[seq_len(m[k])] # Take sub sample
        ct <- countsTab(Y, d = d) # Counts of size d+1 sequences from Y

    # FS:
        S <- hdMTD_FS(Y, d = d, l = 2)
        SFS[i, ( (k * 2 - 1):(k * 2) )] <- S # Stores S
        p_FS <- freqTab(S = S, A = A, countsTab = ct)$qax_Sj[1]
        FS[i, k] <- abs(p_FS - MTD$P[1, 1])
        FSP[i, k] <- abs(p_FS - MTD$P[1, 1])/min(MTD$P[1, 1], MTD$P[1, 2])

    # NAIVE:
        ct_dNaive  <- countsTab(Y, dNaive) # Counts of size dNaive+1 sequences from Y
        p_Naive <- freqTab(S = seq_len(dNaive), A = A, countsTab = ct_dNaive)$qax_Sj[1]
        Naive[i, k] <- abs(p_Naive - MTD$P[1, 1])
        NaiveP[i, k] <- abs(p_Naive - MTD$P[1, 1])/min(MTD$P[1, 1], MTD$P[1, 2])

    # ORACLE:
        p_pairs <- numeric(npairs)
        for (s in seq_len(npairs)) {
          p_pairs[s] <- freqTab(S = pairList[s, ], A = A, countsTab = ct)$qax_Sj[1]
        }
        minpos <- which.min(abs(p_pairs - MTD$P[1,1]))
        ZOracle[i, (( k * 2 - 1):(k * 2))] <- pairList[minpos, ] # Stores lags set with minimal error (Oracle set)
        p_Oracle <- p_pairs[minpos]
        Oracle[i, k] <- abs(p_Oracle - MTD$P[1, 1])
        OracleP[i, k] <- abs(p_Oracle - MTD$P[1, 1] )/min(MTD$P[1, 1],MTD$P[1, 2])
    }
}


# Article Table 1 (Mean error of estimators (d = 5))
round(apply(FS, 2, mean), 5)
round(apply(Oracle, 2, mean), 5)
round(apply(Naive, 2, mean), 5)
round(apply(FSP, 2, mean), 5)
round(apply(OracleP, 2, mean), 5)
round(apply(NaiveP, 2, mean), 5)

# How many times SFS != ZOracle for m = 5000.
sum(apply(SFS[, c(11, 12)], 1, sort) !=  apply(ZOracle[, c(11, 12)], 1, sort))

# How many times SFS != ZOracle for m = 10000.
sum(apply(SFS[, c(13, 14)], 1, sort) != apply(ZOracle[, c(13, 14)], 1, sort))




## Plots (Generating Article Figure 1)
# Data arrangement
# FS
tab <- FS
FStab <- rbind(apply(tab, 2, summary),'sd'=apply(tab, 2, sd))
FStab <- rbind(FStab,'sdLo'=FStab[4,]-FStab[7,],'sdUp'=FStab[4,]+FStab[7,])
Fmean <- FStab[4,]
FsdUp <- FStab[9,]
FsdLo <- FStab[8,]
Fq1 <- FStab[2,]
Fq2 <- FStab[3,]
Fq3 <- FStab[5,]

# NAIVE
tab <- Naive
Naivetab <- rbind(apply(tab, 2, summary),'sd'=apply(tab, 2, sd))
Naivetab <- rbind(Naivetab,'sdLo'=Naivetab[4,]-Naivetab[7,],'sdUp'=Naivetab[4,]+Naivetab[7,])
Nmean <- Naivetab[4,]
NsdUp <- Naivetab[9,]
NsdLo <- Naivetab[8,]
Nq1 <- Naivetab[2,]
Nq2 <- Naivetab[3,]
Nq3 <- Naivetab[5,]

# ORACLE
tab <- Oracle
Oracletab <- rbind(apply(tab, 2, summary),'sd'=apply(tab, 2, sd))
Oracletab <- rbind(Oracletab,'sdLo'=Oracletab[4,]-Oracletab[7,],'sdUp'=Oracletab[4,]+Oracletab[7,])
Omean <- Oracletab[4,]
OsdUp <- Oracletab[9,]
OsdLo <- Oracletab[8,]
Oq1 <- Oracletab[2,]
Oq2 <- Oracletab[3,]
Oq3 <- Oracletab[5,]

# Plotting
# pdf("trio_plot.pdf", width = 9, height = 6)
par(mfrow=c(1,2))
plot(m/100, Fmean, type = "l", col = "#377EB8",
     xlab = "m (x100)", ylab = "Mean error", ylim = c(0, 0.06),lwd=3,
     frame.plot = FALSE, xaxt="n", xlim = c(10,100))
lines(m/100, Omean, type = "l", col = "#E41A1C",lwd=3)
lines(m/100, Nmean, type = "l", col = "#4DAF4A",lwd=3)
points(m/100, Fmean, col = "#377EB8", pch=19,cex=0.7)
points(m/100, Omean, col = "#E41A1C", pch=19,cex=0.7)
points(m/100, Nmean, col = "#4DAF4A", pch=19,cex=0.7)
lines(m/100, FsdUp, type = "l", col = "#377EB8", lty = 2)
lines(m/100, FsdLo, type = "l", col = "#377EB8", lty = 2)
lines(m/100, OsdUp, type = "l", col = "#E41A1C", lty = 2)
lines(m/100, OsdLo, type = "l", col = "#E41A1C", lty = 2)
lines(m/100, NsdUp, type = "l", col = "#4DAF4A", lty = 2)
lines(m/100, NsdLo, type = "l", col = "#4DAF4A", lty = 2)
axis(side = 1, at = m/100, labels = m/100)
legend("topright",
       legend = c(expression(bar(Delta) ~ "FS"),
                  expression(bar(Delta) ~ "FS" %+-% "sd"),
                  expression(bar(Delta) ~ "Oracle"),
                  expression(bar(Delta) ~ "Oracle" %+-% "sd"),
                  expression(bar(Delta) ~ "Naive"),
                  expression(bar(Delta) ~ "Naive" %+-% "sd")),
       col = c("#377EB8","#377EB8","#E41A1C","#E41A1C","#4DAF4A","#4DAF4A"), lty = c(1,3,1,3,1,3), bty = "n")

plot(m/100, Fq2, type = "l", col = "#377EB8",
     xlab = "m (x100)", ylab = "Quartiles of mean error", ylim = c(0, 0.06),lwd=3,
     frame.plot = FALSE, xaxt="n", xlim = c(10,100), lty=1)
lines(m/100, Oq2, type = "l", col = "#E41A1C",lwd=3, lty=1)
lines(m/100, Nq2, type = "l", col = "#4DAF4A",lwd=3, lty=1)
points(m/100, Fq2, col = "#377EB8", pch=19,cex=0.7)
points(m/100, Oq2, col = "#E41A1C", pch=19,cex=0.7)
points(m/100, Nq2, col = "#4DAF4A", pch=19,cex=0.7)
lines(m/100, Fq1, type = "l", col = "#377EB8", lty = 2)
lines(m/100, Fq3, type = "l", col = "#377EB8", lty = 2)
lines(m/100, Oq1, type = "l", col = "#E41A1C", lty = 2)
lines(m/100, Oq3, type = "l", col = "#E41A1C", lty = 2)
lines(m/100, Nq1, type = "l", col = "#4DAF4A", lty = 2)
lines(m/100, Nq3, type = "l", col = "#4DAF4A", lty = 2)
axis(side = 1, at = m/100, labels = m/100)
legend("topright",
       legend = c(expression("Med  " ~ bar(Delta) ~ "FS"),
                  expression("q1,q3" ~ bar(Delta) ~ "FS"),
                  expression("Med  " ~ bar(Delta) ~ "Oracle"),
                  expression("q1,q3" ~ bar(Delta) ~ "Oracle"),
                  expression("Med  " ~ bar(Delta) ~ "Naive"),
                  expression("q1,q3" ~ bar(Delta) ~ "Naive")),
       col = c("#377EB8","#377EB8","#E41A1C","#E41A1C","#4DAF4A","#4DAF4A"), lty = c(1,3,1,3,1,3), bty = "n")
# dev.off()

