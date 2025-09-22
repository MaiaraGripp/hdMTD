####################################################
### R code to reproduce the submitted manuscript
###
### "hdMTD:
###  An R Package for High-Dimensional Mixture Transition
###  Distribution Models"
###
### Authors: Maiara Gripp, Guilherme Ost, Giulio Iacobelli, Daniel Y. Takahashi
### Date: July 2025
####################################################

# Required packages:
# install.packages(c("hdMTD", "dplyr", "ggplot2", "lubridate", "purrr", "tidyr"))

## Load packages
library("hdMTD")
library("dplyr")
library("ggplot2")
library("lubridate")
library("purrr")
library("tidyr")

#'
#' ## Section 5: Using hdMTD
#' ### 5.1 Data generation
#'
#' 1. Generate MTD model:
#'
#' Parameters: $\Lambda = \{-30,-15,-1\}$, $\mathcal{A} = \{0,1\}$,
#' $\lambda_0= \{0.01\}$, $\lambda_{-1} = 0.39$, $\lambda_{-15} = \lambda_{-30} = 0.3$,
#' $p_0(0)=p_0(1)=0.5$, and transition matrices $p_j$, $j\in\Lambda$, sampled uniformly.
set.seed(11)
Lambda <- c(1, 15, 30)
A <- c(0, 1)
lam0 <- 0.01
lamj <- c(0.39, 0.3, 0.3)
p0 <- c(0.5, 0.5)
MTD <- MTDmodel(Lambda = Lambda, A = A, lam0 = lam0, lamj = lamj, p0 = p0)
MTD
#'
#' 2. Sample from the invariant distribution
#'
X <- perfectSample(MTD, N = 1000)
#'
#' ### 5.2 Estimation
#'
#' 3. Estimate relevant lags using FS method
#'
hdMTD_FS(X, d = 40, l = 4)
#' Equivalent call: hdMTD(X, d = 40, method = "FS", l=4)
#'
#' 4. Estimate relevant lags using BIC method
#'
#' **WARNING**: next line takes ~30min (on i7-1255U, 10 cores). Uncomment to compute. <br>
#' hdMTD_BIC(X, d = 40, minl = 4, maxl = 4)
#'
#' Equivalent call: hdMTD(X, d = 40, method = "BIC", minl = 4, maxl = 4)
#'
#' Custom subset S
hdMTD_BIC(X, d = 40,
          S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40),
          minl = 4, maxl = 4)
#'
#' Varying number of lags to be selected
hdMTD_BIC(X, d = 40,
          S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40),
          minl = 1, maxl = 4)
#'
#' With BIC values by number of lags
hdMTD_BIC(X, d = 40,
          S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40),
          minl = 1, maxl = 4,
          byl = TRUE, BICvalue = TRUE)
#'
#' Setting $\xi=0.4$
hdMTD_BIC(X, d = 40,
          S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40),
          minl = 1, maxl = 4,
          byl = TRUE, BICvalue = TRUE,
          xi = 0.4)
#'
#' All matrices $p_j$ are equal and $\lambda_0=0$
hdMTD_BIC(X, d = 40,
          S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40),
          minl = 1, maxl = 4,
          byl = TRUE, BICvalue = TRUE,
          single_matrix = TRUE, indep_part = FALSE)
#'
#' 5. Estimate relevant lags using CUT method
#'
#' **WARNING**: next line takes ~2.5min (on i7-1255U, 10 cores). Uncomment to compute. <br>
#' hdMTD_CUT(X, d = 40, S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40))
#'
#' Equivalent call: hdMTD(X,d=40,S=c(1,5,10,15,17,20,27,30,35,40),method="CUT")
#'
#' Setting $\alpha = 0.13$ <br>
#' **WARNING**: next line takes ~2.5min (on i7-1255U, 10 cores). Uncomment to compute. <br>
#' hdMTD_CUT(X, d = 40, S = c(1, 5, 10, 15, 17, 20, 27, 30, 35, 40), alpha = 0.13)
#'
#' Custom subset S
hdMTD_CUT(X, d = 40, S = c(1, 5, 17, 27, 30, 35), alpha = 0.13)
#'
#' 6. Estimate relevant lags using FSC method
#'
hdMTD_FSC(X, d = 40, l = 4, alpha = 0.1)
#' Equivalent call: hdMTD(X, d = 40, method = "FSC", l = 4, alpha = 0.1)
#'
#' FS method with halved sample
hdMTD_FS(X[1:500], d = 40, l = 4)
#'
#' 7. Estimating transition probabilities
#'
empirical_probs(X, S = c(1, 15, 30))
empirical_probs(X, S = c(1, 15, 30), matrixform = TRUE)
#'
#' 8. Oscillations
#'
#' Computing from MTD
oscillation(MTD)
#' Estimating from sample
oscillation(X, S = c(1, 15, 30))
#'
#' 9. Estimating MTD parameters through the EM algorithm
#'
#' Initial parameters for EM method
init <- list(
  'lambdas'= c(0.01, 0.33, 0.33, 0.33),
  'p0' = c(0.5, 0.5),
  'pj' = rep(list(matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2, nrow = 2)), 3)
)
#' Run EM
MTDest(X, S = c(1, 15, 30), init = init, iter = TRUE)
#' Stops after $9$ iterations
MTDest(X, S = c(1, 15, 30), M = NULL, nIter = 9, init = init, oscillations = TRUE)
#'
#' Estimating P (transition probabilities) with EM
#'
#' Estimate MTD parameters with EM
estParam <- MTDest(X, S = c(1, 15, 30), init = init)
#' Set MTD model with estimations
estMTD <- MTDmodel(Lambda, A, lam0 = estParam$lambdas[1],
                   lamj = estParam$lambdas[-1], p0 = estParam$p0,
                   pj = estParam$pj)
#' Return estimated transition matrix
estMTD$P
#'
#' ### 5.3 Testing hdMTD
#'
#' 1. MTD model specification:
#'
set.seed(123)
Lambda <- c(1, 5)
A <- c(0, 1)
lam0 <- 0.01
p0 <- c(0.5, 0.5)
MTD <- MTDmodel(Lambda, A, lam0, p0 = p0) # Generates an MTD model
#'
#' 2. Simulation settings and results
#'
#' Simulation parameters:
n <- 100      # Number of replications
N <- 10000    # Full sample size
m <- c(1000, 1500, 2000, 2500, 3000, 5000, 10000) # Subsample sizes
d <- 100      # Max order for FS and Oracle
dNaive <- 5   # Max order for Naive
pairList <- t(combn(100, 2)) # All possible pairs with digits from 1 to 100
npairs <- nrow(pairList)
#'
#' **WARNING**: The following simulation takes approximately 5–6 DAYS to complete.
#' For reproducibility and ease of access, we use a precomputed `.rds` file (simulated_data.rds).
#' To regenerate the results manually, set `recompute <- TRUE` bellow:
recompute <- FALSE
if (recompute) {
  # Initialize results storage
  FSP <- matrix(0, ncol = length(m),nrow = n)
  FS <- matrix(0, ncol = length(m),nrow = n)
  NaiveP <- matrix(0, ncol = length(m), nrow = n)
  Naive <- matrix(0, ncol = length(m), nrow = n)
  OracleP <- matrix(0, ncol = length(m),nrow = n)
  Oracle <- matrix(0, ncol = length(m), nrow = n)
  SFS <- matrix(0, ncol = length(m)*2, nrow = n)
  ZOracle <- matrix(0, ncol = length(m)*2, nrow = n)
  # WARNING: the next loop takes 5 to 6 DAYS to complete
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
} else {
  #' Load precomputed results from `simulated_data.rds`
  simulated_data <- readRDS("simulated_data.rds")
  FS <- simulated_data$FS
  FSP <- simulated_data$FSP
  Oracle <- simulated_data$Oracle
  OracleP <- simulated_data$OracleP
  Naive <- simulated_data$Naive
  NaiveP <- simulated_data$NaiveP
  SFS <- simulated_data$SFS
  ZOracle <- simulated_data$ZOracle
}
#'
#' 3. Generate Table 1
means_table <- rbind(
  round(apply(FS, 2, mean), 5),
  round(apply(Oracle, 2, mean), 5),
  round(apply(Naive, 2, mean), 5),
  round(apply(FSP, 2, mean), 5),
  round(apply(OracleP, 2, mean), 5),
  round(apply(NaiveP, 2, mean), 5)
)
rownames(means_table) <- c(
  "Delta_FS(m)",
  "Delta_Oracle(m)",
  "Delta_Naive,5(m)",
  "std Delta_FS(m)",
  "std Delta_Oracle(m)",
  "std Delta_Naive,5(m)"
)
colnames(means_table) <- m
#'
#' Table 1: Mean error of estimators
#' ```{r mean-error-table, results='asis', message=FALSE, warning=FALSE}
#' knitr::kable(means_table, caption = "Mean error of estimators", format = "markdown", booktabs = TRUE, escape = FALSE)
#' ```
#'
#' 4. Compute how often the FS output differs from Oracle by subsample size
m_index <- seq(1, length(m)*2, by = 2)
names(m_index) <- as.character(m)
SFS_vs_ZOracle_diff <- sapply(m_index, function(idx) {
  SFS_set <- apply(SFS[, c(idx, idx + 1)], 1, function(x) paste(sort(x), collapse = "-"))
  ZOracle_set <- apply(ZOracle[, c(idx, idx + 1)], 1, function(x) paste(sort(x), collapse = "-"))
  sum(SFS_set != ZOracle_set)
})
SFS_vs_ZOracle_diff
#'
#' 5. Generate Figure 1:
#'
#' Data arrangement
#'
#' FS
tab <- FS
FStab <- rbind(apply(tab, 2, summary),'sd'=apply(tab, 2, sd))
FStab <- rbind(FStab,'sdLo'=FStab[4,]-FStab[7,],'sdUp'=FStab[4,]+FStab[7,])
Fmean <- FStab[4,]
FsdUp <- FStab[9,]
FsdLo <- FStab[8,]
Fq1 <- FStab[2,]
Fq2 <- FStab[3,]
Fq3 <- FStab[5,]
#'
#' NAIVE
tab <- Naive
Naivetab <- rbind(apply(tab, 2, summary),'sd'=apply(tab, 2, sd))
Naivetab <- rbind(Naivetab,'sdLo'=Naivetab[4,]-Naivetab[7,],'sdUp'=Naivetab[4,]+Naivetab[7,])
Nmean <- Naivetab[4,]
NsdUp <- Naivetab[9,]
NsdLo <- Naivetab[8,]
Nq1 <- Naivetab[2,]
Nq2 <- Naivetab[3,]
Nq3 <- Naivetab[5,]
#'
#' ORACLE
tab <- Oracle
Oracletab <- rbind(apply(tab, 2, summary),'sd'=apply(tab, 2, sd))
Oracletab <- rbind(Oracletab,'sdLo'=Oracletab[4,]-Oracletab[7,],'sdUp'=Oracletab[4,]+Oracletab[7,])
Omean <- Oracletab[4,]
OsdUp <- Oracletab[9,]
OsdLo <- Oracletab[8,]
Oq1 <- Oracletab[2,]
Oq2 <- Oracletab[3,]
Oq3 <- Oracletab[5,]
#'
#' ### Plot Figure 1: Estimators mean error across $N_{rep}=100$ replications.
#' ```{r fig1-plot, fig.width=10, fig.height=5, message=FALSE, warning=FALSE}
#' par(mfrow=c(1,2))
#'
#' # --- Left panel: Mean error with standard deviation bands ---
#'  plot(m/100, Fmean, type = "l", col = "#377EB8",
#'       xlab = "m (x100)", ylab = "Mean error", ylim = c(0, 0.06),lwd=3,
#'       frame.plot = FALSE, xaxt="n", xlim = c(10,100))
#'  lines(m/100, Omean, type = "l", col = "#E41A1C",lwd=3)
#'  lines(m/100, Nmean, type = "l", col = "#4DAF4A",lwd=3)
#'  points(m/100, Fmean, col = "#377EB8", pch=19,cex=0.7)
#'  points(m/100, Omean, col = "#E41A1C", pch=19,cex=0.7)
#'  points(m/100, Nmean, col = "#4DAF4A", pch=19,cex=0.7)
#'  lines(m/100, FsdUp, type = "l", col = "#377EB8", lty = 2)
#'  lines(m/100, FsdLo, type = "l", col = "#377EB8", lty = 2)
#'  lines(m/100, OsdUp, type = "l", col = "#E41A1C", lty = 2)
#'  lines(m/100, OsdLo, type = "l", col = "#E41A1C", lty = 2)
#'  lines(m/100, NsdUp, type = "l", col = "#4DAF4A", lty = 2)
#'  lines(m/100, NsdLo, type = "l", col = "#4DAF4A", lty = 2)
#'  axis(side = 1, at = m/100, labels = m/100)
#'  legend("topright",
#'         legend = c(expression(bar(Delta) ~ "FS"),
#'                    expression(bar(Delta) ~ "FS" %+-% "sd"),
#'                    expression(bar(Delta) ~ "Oracle"),
#'                    expression(bar(Delta) ~ "Oracle" %+-% "sd"),
#'                    expression(bar(Delta) ~ "Naive"),
#'                    expression(bar(Delta) ~ "Naive" %+-% "sd")),
#'         col = c("#377EB8","#377EB8","#E41A1C","#E41A1C","#4DAF4A","#4DAF4A"), lty = c(1,3,1,3,1,3), bty = "n")
#'  # --- Right panel: Median and quartiles ---
#'  plot(m/100, Fq2, type = "l", col = "#377EB8",
#'       xlab = "m (x100)", ylab = "Quartiles of mean error", ylim = c(0, 0.06),lwd=3,
#'       frame.plot = FALSE, xaxt="n", xlim = c(10,100), lty=1)
#'  lines(m/100, Oq2, type = "l", col = "#E41A1C",lwd=3, lty=1)
#'  lines(m/100, Nq2, type = "l", col = "#4DAF4A",lwd=3, lty=1)
#'  points(m/100, Fq2, col = "#377EB8", pch=19,cex=0.7)
#'  points(m/100, Oq2, col = "#E41A1C", pch=19,cex=0.7)
#'  points(m/100, Nq2, col = "#4DAF4A", pch=19,cex=0.7)
#'  lines(m/100, Fq1, type = "l", col = "#377EB8", lty = 2)
#'  lines(m/100, Fq3, type = "l", col = "#377EB8", lty = 2)
#'  lines(m/100, Oq1, type = "l", col = "#E41A1C", lty = 2)
#'  lines(m/100, Oq3, type = "l", col = "#E41A1C", lty = 2)
#'  lines(m/100, Nq1, type = "l", col = "#4DAF4A", lty = 2)
#'  lines(m/100, Nq3, type = "l", col = "#4DAF4A", lty = 2)
#'  axis(side = 1, at = m/100, labels = m/100)
#'  legend("topright",
#'         legend = c(expression("Med  " ~ bar(Delta) ~ "FS"),
#'                    expression("q1,q3" ~ bar(Delta) ~ "FS"),
#'                    expression("Med  " ~ bar(Delta) ~ "Oracle"),
#'                    expression("q1,q3" ~ bar(Delta) ~ "Oracle"),
#'                    expression("Med  " ~ bar(Delta) ~ "Naive"),
#'                    expression("q1,q3" ~ bar(Delta) ~ "Naive")),
#'         col = c("#377EB8","#377EB8","#E41A1C","#E41A1C","#4DAF4A","#4DAF4A"), lty = c(1,3,1,3,1,3), bty = "n")
#' ```
#'
#' ### 5.4 Analysis of Real-World Data
#'
data("tempdata")
#'
#' 1. Treat NA data:
#'
#' Removing days before "2010-08-05"
tempdata <- hdMTD::tempdata %>% filter(DATE >= "2010-08-05")
#' Identify remaining $155$ NA positions
posNA <- which(is.na(tempdata$MAXTEMP))
#' Fill short sequences of NAs ( $≤6$ NAs) using nearest neighbors mean
for (i in posNA) {
  if(!is.na(tempdata$MAXTEMP[i - 1]) && !all(is.na(tempdata$MAXTEMP[(i + 1):(i + 6)]))) {
    aux <- which(!is.na(tempdata$MAXTEMP[(i + 1):(i + 6)]))[1]
    tempdata$MAXTEMP[i] <- mean(tempdata$MAXTEMP[i - 1],tempdata$MAXTEMP[i + aux])
  }
}
#' Identify remaining $88$ NA positions
posNA <- which(is.na(tempdata$MAXTEMP))
#' Fill remaining NA with the mean of the previous hour temperature, next hour
#' temperature, and same hour of previous day temperature.
for (i in posNA) {
  tempdata$MAXTEMP[i] <- mean(c(tempdata$MAXTEMP[i - 1], tempdata$MAXTEMP[i + 1],
                                tempdata$MAXTEMP[i - 24]),na.rm = TRUE)
}
#'
#' 2. Compute mean daily maximum temperatures:
#'
temp <- tempdata %>%
  group_by(DATE) %>%
  summarize(MAXTEMP = mean(MAXTEMP), .groups = 'drop')
head(temp, 4)
#'
#' ### Plot Figure 2: Time series with quarterly mean of daily maximum temperatures
#' ```{r temp-plot, fig.width=9, fig.height=6, message=FALSE, warning=FALSE}
#' TRIM_DATA <- temp %>%
#'   mutate(
#'     Y_TRIMESTER = paste0(year(DATE), "-T", quarter(DATE))
#'   ) %>%
#'   group_by(Y_TRIMESTER) %>%
#'   summarise(
#'     MEAN_TEMP = mean(MAXTEMP),
#'     DATA_REF = min(DATE)
#'   ) %>%
#'   ungroup() %>%
#'   arrange(DATA_REF)
#' TRIM_DATA <- TRIM_DATA[-c(1, nrow(TRIM_DATA)),]
#'
#' ggplot(TRIM_DATA, aes(x = DATA_REF, y = MEAN_TEMP)) +
#'   geom_line(color = "steelblue", linewidth = 0.5) +
#'   geom_point(color = "steelblue", size = 1.2) +
#'   scale_x_date(
#'     date_breaks = "1 year",
#'     date_labels = "%Y",
#'     minor_breaks = NULL
#'   ) +
#'   theme_minimal() +
#'   theme(
#'     axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
#'     axis.text.y = element_text(size = 12),
#'     panel.grid.major = element_line(color = "gray90"),
#'     axis.title.x = element_text(size = 15),
#'     axis.title.y = element_text(size = 15),
#'     plot.title = element_text(hjust = 0.5, size = 17, face = "bold")
#'   ) +
#'   labs(
#'     title = "Quarterly mean of daily maximum temperatures across the years",
#'     x = "Year",
#'     y = "Mean Temperature (°C)"
#'   )
#' ```
#'
#' 3. Create categories of temperature:
#'
xn <- max(temp$MAXTEMP)
x1 <- min(temp$MAXTEMP)
maxAmp <- xn - x1
temp$MAXTEMP1  <- ifelse(temp$MAXTEMP < x1 + maxAmp/2, 1, 2)
head(temp, 4)
prop.table(table(temp$MAXTEMP1)) # frequency of thermal regimes
#'
#' 4. Run FS for temp:
#'
Temp12 <- rev(temp$MAXTEMP1)
#' hdMTD functions assume the sample is sorted from the latest observations to
#' oldest.
#'
#' **WARNING**: next line takes ~7min (on i7-1255U, 10 cores). Uncomment to compute. <br>
#' hdMTD_FS(Temp12, d = 400, l = 3)
#'
#' Reduce maximum order to improve estimation <br>
#' Note: This code line is mentioned in the article but without a CodeChunk <br>
#' **WARNING**: next line takes ~6min (on i7-1255U, 10 cores). Uncomment to compute. <br>
#' hdMTD_FS(Temp12, d = 364, l = 3)
#'
#' 5. Split sample in Train and Test data:
#'
ndays <- nrow(temp %>%
                filter(DATE >= "2023-09-01")) # 366 days in the latest year of the sample
Temp12_Train <- Temp12[-seq_len(ndays)] # Training data with 4775 days
Temp12_Test <- Temp12[seq_len(ndays)] # Test data
#'
#' 6. Rerun FS for Train data:
#'
#' **WARNING**: next line takes ~5min (on i7-1255U, 10 cores). Uncomment to compute. <br>
#' hdMTD_FS(Temp12_Train, d = 364, l = 3)
#'
#' 7. Trim out irrelevant lags:
#'
#' With CUT method
hdMTD_CUT(Temp12_Train, d = 364, S = c(1, 364, 6))
#' With BIC method
hdMTD_BIC(Temp12_Train, d = 364, S = c(1, 364, 6), minl = 1, maxl = 3,
          byl = TRUE, BICvalue = TRUE )
#'
#' 8. Lag selection with FSC method:
#'
#' **WARNING**: next line takes ~3min (on i7-1255U, 10 cores). Uncomment to compute. <br>
#' hdMTD_FSC(Temp12_Train, d = 364, l = 3)
#'
#' 9. Estimated transition matrix for FS method output:
#'
P_FS <- empirical_probs(Temp12_Train, S = c(1, 6, 364), matrixform = T)
P_FS
#'
#' ### Classic method for choosing relevant lag set:
#'
#' 10. Compute models:
#'
ct <- countsTab(Temp12_Train, d = 6) # Table with size 6 sequence counts
head(ct,4)
#'
#' MC1
ft <- freqTab(S = 1, A = c(1, 2), countsTab = ct)
LL <- sum(log(ft$qax_Sj) * ft$Nxa_Sj)
freeParam <- 2 * 1
BICMC1 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
BICMC1
#' Comparable BIC if the model is a Markov chain of order $1$: $1869.162$
#'
#' MC2
ft <- freqTab(S = c(1, 2), A = c(1, 2), countsTab = ct)
LL <- sum(log(ft$qax_Sj) * ft$Nxa_Sj)
freeParam <- 2^2 * 1
BICMC2 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
BICMC2
#' Comparable BIC if the model is a Markov chain of order $2$: $1850.598$
#'
#' MC3. Used as example in the article
ft <- freqTab(S = c(1, 2, 3), A = c(1, 2), countsTab = ct)
head(ft, 4)
LL <- sum(log(ft$qax_Sj) * ft$Nxa_Sj)
freeParam <- 2^3 * 1
BICMC3 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
BICMC3
#' Comparable BIC if the model is a Markov chain of order $3$: $1854.029$
#'
#' MC4
ft <- freqTab(S = c(1, 2, 3, 4), A = c(1, 2), countsTab = ct)
LL <- sum(log(ft$qax_Sj) * ft$Nxa_Sj)
freeParam <- 2^4 * 1
BICMC4 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
BICMC4
#' Comparable BIC if the model is a Markov chain of order $4$: $1877.888$
#'
#' MC5
#'
ft <- freqTab(S = c(1, 2, 3, 4, 5), A = c(1, 2), countsTab = ct)
pos <- which(ft$Nxa_Sj > 0)
LL <- sum(log(ft$qax_Sj[pos]) * ft$Nxa_Sj[pos])
freeParam <- 2^5 * 1
BICMC5 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
BICMC5
#' Comparable BIC if the model is a Markov chain of order $5$: $1925.962$
#'
#' MC6
ft <- freqTab(S = c(1, 2, 3, 4, 5, 6), A = c(1, 2), countsTab = ct)
pos <- which(ft$Nxa_Sj > 0)
LL <- sum(log(ft$qax_Sj[pos]) * ft$Nxa_Sj[pos])
freeParam <- 2^6 * 1
BICMC6 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
BICMC6
#' Comparable BIC if the model is a Markov chain of order $6$: $2031.679$
#'
#' 11. Comparing models:
#'
BIC_vals <-  c(BICMC1, BICMC2, BICMC3, BICMC4, BICMC5, BICMC6)
model_names <- paste0("MC", 1:6)
minBIC_idx <- which.min(BIC_vals)
BIC_fmt <- sprintf("%.3f", BIC_vals)
BIC_fmt[minBIC_idx] <- paste0("**", BIC_fmt[minBIC_idx], "**")
minBIC_idx
#' The classic method chooses order $2$ ($S=\{-2,-1\}$).
#'
#' Table 2: BIC values
bic_matrix <- data.frame(matrix(BIC_fmt, nrow = 1))
colnames(bic_matrix) <- model_names
rownames(bic_matrix) <- "BIC"
#' ```{r bic-table, results='asis', message=FALSE, warning=FALSE}
#' knitr::kable(bic_matrix, caption = "BIC values computed for classical Markov chain models of different orders.")
#' ```
#'
#' Estimated matrix for $S=\{-2,-1\}$
P_MC2 <- empirical_probs(Temp12_Train, S = c(1, 2), matrixform = TRUE)
P_MC2
#' Independent model distribution:
P_Ind <- prop.table(table(Temp12_Train))
P_Ind
#'
#' ### Comparing methods
#'
#' 12. Computing values for Table 3.
#'
Days1 <- which(Temp12_Test == 1)
lenDays1 <- length(Days1)
lenDays1 / ndays # frequency of low temperature days
Temp12_Test <- c(Temp12_Test, Temp12_Train[seq_len(364)])
#'
set.seed(1)
hitInd <- numeric(1000)
hitMC2 <- numeric(1000)
hitFS <- numeric(1000)
T1Ind <- numeric(1000)
T1MC2 <- numeric(1000)
T1FS <- numeric(1000)
F1Ind <- numeric(1000)
F1MC2 <- numeric(1000)
F1FS <- numeric(1000)
pasts2 <- rownames(P_MC2)
pastsFS <- rownames(P_FS)
#'
for (j in seq_len(1000)){
  u <- runif(ndays)
  predInd <- numeric(ndays)
  predMC2 <- numeric(ndays)
  predFS <- numeric(ndays)

  for (i in ndays:1) {
    predInd[i] <- ifelse(u[i] <= P_Ind[1], 1, 2)
    pastRow <- which(pasts2 == paste0(rev(Temp12_Test[c(i+1, i+2)]), collapse = ""))
    predMC2[i] <- ifelse(u[i] <= P_MC2[pastRow, 1], 1, 2)
    pastRow <- which(pastsFS == paste0(rev(Temp12_Test[c(i+1, i+6, i+364)]), collapse = ""))
    predFS[i] <- ifelse(u[i] <= P_FS[pastRow, 1], 1, 2)
  }

  hitInd[j] <- sum(predInd == Temp12_Test[seq_len(ndays)])
  hitMC2[j] <- sum(predMC2 == Temp12_Test[seq_len(ndays)])
  hitFS[j] <- sum(predFS == Temp12_Test[seq_len(ndays)])
  T1Ind[j] <- sum(predInd[Days1] == 1)
  T1MC2[j] <- sum(predMC2[Days1] == 1)
  T1FS[j] <- sum(predFS[Days1] == 1)
  F1Ind[j] <- sum(predInd[-Days1] == 1)
  F1MC2[j] <- sum(predMC2[-Days1] == 1)
  F1FS[j] <- sum(predFS[-Days1] == 1)
}

#' ### Accuracy
AccInd <- mean(hitInd)/ndays
AccMC2 <- mean(hitMC2)/ndays
AccFS <- mean(hitFS)/ndays
AccInd; AccMC2; AccFS

#' ### Precision
PrecInd <- mean(T1Ind/(T1Ind + F1Ind))
PrecMC2 <- mean(T1MC2/(T1MC2 + F1MC2))
PrecFS <- mean(T1FS/(T1FS + F1FS))
PrecInd; PrecMC2; PrecFS

#' ### Sensitivity (Recall)
SensInd <- mean(T1Ind)/lenDays1
SensMC2 <- mean(T1MC2)/lenDays1
SensFS <- mean(T1FS)/lenDays1
SensInd; SensMC2; SensFS

#' ### Specificity
SpecInd <- 1 - mean(F1Ind)/(ndays - lenDays1)
SpecMC2 <- 1 - mean(F1MC2)/(ndays - lenDays1)
SpecFS <- 1 - mean(F1FS)/(ndays - lenDays1)
SpecInd; SpecMC2; SpecFS

#' ### F1-Score
F1ScoreInd <- 2 * (PrecInd * SensInd) / (PrecInd + SensInd)
F1ScoreMC2 <- 2 * (PrecMC2 * SensMC2) / (PrecMC2 + SensMC2)
F1ScoreFS <- 2 * (PrecFS * SensFS) / (PrecFS + SensFS)
F1ScoreInd; F1ScoreMC2; F1ScoreFS

#' Table 3: Model performance metrics
metric <- c("Accuracy", "Precision", "Sensitivity (Recall)", "Specificity", "F1-Score")
formula <- c("(TP+TN)/(TP+TN+FP+FN)",
             "TP/(TP+FP)",
             "TP/(TP+FN)",
             "TN/(TN+FP)",
             "2(PPV*Recall)/(PPV+Recall)")
performance_table <- data.frame(
  Metric = metric,
  Formula = formula,
  indc = round(c(AccInd, PrecInd, SensInd, SpecInd, F1ScoreInd) * 100, 2),
  mc2c = round(c(AccMC2, PrecMC2, SensMC2, SpecMC2, F1ScoreMC2) * 100, 2),
  fsc = round(c(AccFS, PrecFS, SensFS, SpecFS, F1ScoreFS) * 100, 2),
  check.names = FALSE
)
names(performance_table) <- c("Metric", "Formula", "Ind (\\%)", "MC2 (\\%)", "FS (\\%)")
#' ```{r performance-table, results='asis', message=FALSE, warning=FALSE}
#' knitr::kable(performance_table, align = "l", caption = "Model performance metrics.")
#' ```
#'
#' ### Plot Figure 3: Exploratory analysis of accuracies
#'
#' ```{r accuracy-plot, fig.width=9, fig.height=6, message=FALSE, warning=FALSE}
#' accuracy_data <- data.frame(
#'   MC2 = hitMC2 / ndays,
#'   FS = hitFS / ndays
#' ) %>%
#'   pivot_longer(
#'     everything(),
#'     names_to = "Model",
#'     values_to = "Accuracy"
#'   )
#'
#' ggplot(accuracy_data, aes(x = Model, y = Accuracy, fill = Model)) +
#'   geom_boxplot() +
#'   labs(
#'     title = "Accuracy distribution (1000 replications)",
#'     x = "Model",
#'     y = "Accuracy"
#'   ) +
#'   theme_minimal() +
#'   scale_fill_brewer(palette = "Paired") +
#'   theme(
#'     plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
#'     axis.title = element_text(size = 16),
#'     axis.text = element_text(size = 14),
#'     legend.title = element_text(size = 16),
#'     legend.text = element_text(size = 14)
#'   )
#' ```
#'
#' ### Empirical $\nu$ Analysis
#'
#' 13. FS sequential selection based on $\hat{\nu}_{n,j,S}$ values:
#'
run_sequential_lag_selection <- function(Temp12_Train, d = 364) {
  # Initialization
  A <-  sort(unique(Temp12_Train))
  lenA <- length(A)
  lenX <- length(Temp12_Train)
  A_pairs <- matrix(A, ncol = 2) # All unique state pairs
  ct <- countsTab(X = Temp12_Train, d = 364) # Sequence counts table

  # Initialize storage
  results <- list(
    nuj1 = numeric(d),
    nuj2 = numeric(d-1),
    nuj3 = numeric(d-2),
    selected_lags = numeric(3)
  )

  # Helper function for empirical distribution calculation
  PI <- function(S, groupTab, x_S, lenX, d) {
    if (length(S) > 0) {
      filtr_S <- paste0("x", S)
      groupTab <- groupTab %>%
        dplyr::mutate(match = purrr::pmap_lgl(dplyr::pick(dplyr::all_of(filtr_S)),
                                              ~ all(c(...) == x_S))) %>%
        dplyr::filter(match) %>%
        dplyr::select(-match)
    }
    PI <- matrix(groupTab$Nx_Sj/(lenX - d),ncol = 1)
    PI
  }

  # Sequential Selection Process
  cat("=== Starting Sequential Lag Selection ===\n")

  # Step 1: Initial selection (S = ∅)
  cat("\n[Step 1] Selecting first lag (S = ∅)...\n")
  S <- NULL
  Sc <- sort(setdiff(seq_len(d), S), decreasing = TRUE) # Complement of S in 1:d

  for (z in seq_along(Sc)) { # Runs across all available lags
    j <- Sc[z]
    # Frequency tables
    b_Sja <- freqTab(S = S, j = j, A = A, countsTab = ct)
    b_Sj <- b_Sja %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(paste0("x", j)))) %>%
      dplyr::summarise(Nx_Sj = sum(Nxa_Sj), .groups = "drop")

    # Compute νj
    PIs <- PI(S = S, groupTab = b_Sj, x_S = S, lenX = lenX, d = 364)
    dTVs <- dTV_sample(S = S, j = j, lenA = lenA, base = b_Sja,
                       A_pairs = A_pairs, x_S = S)
    results$nuj1[z] <- prod(PIs) * dTVs
  }
  results$selected_lags[1] <- Sc[which.max(results$nuj1)]
  cat(sprintf("Selected: j = %d (ν = %.4f)\n",
              results$selected_lags[1], max(results$nuj1)))

  # Step 2: Second selection (S = {1})
  cat(sprintf("\n[Step 2] Selecting second lag (S = {%d})...\n",
              results$selected_lags[1]))
  S <- results$selected_lags[1]
  Sc <- sort(setdiff(seq_len(d), S), decreasing = TRUE)

  for (z in seq_along(Sc)) {
    j <- Sc[z]
    Sj <- sort(c(S, j), decreasing = TRUE)

    # Frequency tables
    b_Sja <- freqTab(S = S, j = j, A = A, countsTab = ct)
    b_Sj <- b_Sja %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(paste0("x", Sj)))) %>%
      dplyr::summarise(Nx_Sj = sum(Nxa_Sj), .groups = "drop")
    b_S <- b_Sja %>%
      dplyr::dplyr::group_by(dplyr::across(dplyr::all_of(paste0("x", S)))) %>%
      dplyr::summarise(Nx_Sj = sum(Nxa_Sj), .groups = "drop")

    subx <- b_S[, -ncol(b_S)]

    # Compute ν_j
    for (t in which(b_S$Nx_Sj > 0)) {
      PIs <- PI(S = S, groupTab = b_Sj, x_S = subx[t, ],
                lenX = lenX, d = d)
      dTVs <- dTV_sample(S = S, j = j, lenA = lenA, base = b_Sja,
                         A_pairs = A_pairs, x_S = subx[t, ])
      PI_xS <- as.numeric(b_S[t, ncol(b_S)]/(lenX - d))
      results$nuj2[z] <- results$nuj2[z] + prod(PIs) * dTVs/PI_xS
    }
  }
  results$selected_lags[2] <- Sc[which.max(results$nuj2)]
  cat(sprintf("Selected: j = %d (ν = %.4f)\n",
              results$selected_lags[2], max(results$nuj2)))

  # Step 3: Third selection (S = {1,364})
  cat(sprintf("\n[Step 3] Selecting third lag (S = {%d,%d})...\n",
              results$selected_lags[1], results$selected_lags[2]))
  S <- c(S, results$selected_lags[2])
  Sc <- sort(setdiff(seq_len(d), S), decreasing = TRUE)

  for (z in seq_along(Sc)) {
    j <- Sc[z]
    Sj <- sort(c(S, j), decreasing = TRUE)
    dec_S <- rev(S) # S in decreasing order

    # Frequency tables
    b_Sja <- freqTab(S = dec_S, j = j, A = A, countsTab = ct)
    b_Sj <- b_Sja %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(paste0("x", Sj)))) %>%
      dplyr::summarise(Nx_Sj = sum(Nxa_Sj), .groups = "drop")
    b_S <- b_Sja %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(paste0("x", dec_S)))) %>%
      dplyr::summarise(Nx_Sj = sum(Nxa_Sj), .groups = "drop")

    subx <- b_S[, -ncol(b_S)]

    # Compute ν_j
    for (t in which(b_S$Nx_Sj > 0)) {
      PIs <- PI(S = dec_S, groupTab = b_Sj, x_S = subx[t, ],
                lenX = lenX, d = d)
      dTVs <- dTV_sample(S = dec_S, j = j, lenA = lenA, base = b_Sja,
                         A_pairs = A_pairs, x_S = subx[t, ])
      PI_xS <- as.numeric(b_S[t, ncol(b_S)]/(lenX - d))
      results$nuj3[z] <- results$nuj3[z] + prod(PIs) * dTVs/PI_xS
    }
  }
  results$selected_lags[3] <- Sc[which.max(results$nuj3)]
  cat(sprintf("Selected: j = %d (ν = %.4f)\n",
              results$selected_lags[3], max(results$nuj3)))

  # Final Results
  cat("\n=== Final Selection Results ===\n")
  print(data.frame(Step = 1:3, Selected_Lag = results$selected_lags,
                   nu = c(max(results$nuj1), max(results$nuj2), max(results$nuj3))))

  return(results)
}
#'
#' **WARNING**: The next code line takes ~6 min (on i7-1255U, 10 cores).
#' For reproducibility and ease of access, we use a precomputed `.rds` file
#' (results_sequential_selection.rds). To regenerate the results manually, set
#' `recompute <- TRUE` bellow: <br>
recompute <- FALSE
if (recompute) {
  results <- run_sequential_lag_selection(Temp12_Train)
} else {
  results <- readRDS("results_sequential_selection.rds") # using precomputed results
}
#'
#' ### Plot Figure 4: FS sequential step analysis through $\hat{\nu}_{n,j,S}$.
#' ```{r nu-plot, fig.width=11, fig.height=6, message=FALSE, warning=FALSE}
#' par(mfrow = c(1, 3), mar = c(5, 6, 4, 2), oma = c(0, 0, 4, 0))
#' palette <- c("#E41A1C", "#377EB8", "#4DAF4A")
#' with(results, {
#'   # Graph 1
#'   Sc <-  364:1
#'   plot(1:364, rev(nuj1), type = "p", pch = 19, cex = 0.8, col = "gray70",
#'        ylab = "", xlab = "Lag (-j)", cex.lab = 1.8, cex.axis = 1.3,
#'        ylim = c(0,0.13), main = "", panel.first = grid())
#'   title(main = expression(paste("S = ", Ø)), cex.main = 1.5, font.main = 1)
#'   title(ylab = expression(widehat(nu)[n*","*j*","*S]), line = 4, cex.lab = 1.8)
#'   points(Sc[which.max(nuj1)], nuj1[which.max(nuj1)], pch = 21, bg = palette[1], cex = 1.5, lwd = 1)
#'   text(Sc[which.max(nuj1)], nuj1[which.max(nuj1)], labels = paste0(Sc[which.max(nuj1)]),
#'        pos = 3, col = palette[1], font = 2, cex = 1.4)
#'   # Graph 2
#'   Sc <- 364:2
#'   plot(2:364, rev(nuj2), type = "p", pch = 19, cex = 0.8, col = "gray70",
#'        ylab = "", xlab = "Lag (-j)", cex.lab = 1.8, cex.axis = 1.3,
#'        ylim = c(0,0.025),  main = "", panel.first = grid())
#'   title(main = "With S = {-1}", cex.main = 1.5, font.main = 1)
#'   title(ylab = expression(widehat(nu)[n*","*j*","*S]), line = 4, cex.lab = 1.8)
#'   points(Sc[which.max(nuj2)], nuj2[which.max(nuj2)], pch = 21, bg = palette[2], cex = 1.5, lwd = 1)
#'   text(Sc[which.max(nuj2)]-5, nuj2[which.max(nuj2)], labels = paste0(Sc[which.max(nuj2)]),
#'        pos = 3, col = palette[2], font = 2, cex = 1.4)
#'   # Graph 3
#'   Sc <- 363:2
#'   plot(2:363, rev(nuj3), type = "p", pch = 19, cex = 0.8, col = "gray70",
#'        ylab = "", xlab = "Lag (-j)", cex.lab = 1.8, cex.axis = 1.3,
#'        ylim = c(0,0.025), main = "", panel.first = grid())
#'   title(main = "With S = {-364, -1}", cex.main = 1.5, font.main = 1)
#'   title(ylab = expression(widehat(nu)[n*","*j*","*S]), line = 4, cex.lab = 1.8)
#'   points(Sc[which.max(nuj3)], nuj3[which.max(nuj3)], pch = 21, bg = palette[3], cex = 1.5, lwd = 1)
#'   text(Sc[which.max(nuj3)], nuj3[which.max(nuj3)], labels = paste0(Sc[which.max(nuj3)]),
#'        pos = 3, col = palette[3], font = 2, cex = 1.4)
#'
#'   mtext(expression(paste("Sequential lag selection based on ", widehat(nu)[n*","*j*","*S])),
#'         outer = TRUE, cex = 1.6, font = 2, line = 1.4)
#' })
#' ```
#'
#'
#' ## Session info
sessionInfo()
