
pkgs <- c("hdMTD", "dplyr", "lubridate", "purrr", "ggplot2", "tidyr")
install.packages(pkgs[!pkgs %in% installed.packages()])

require(hdMTD)
require(dplyr)
require(ggplot2)
require(lubridate)
require(purrr)
require(tidyr)

####################################
## Section 5: Using hdMTD
## 5.4 Analysis of Real-World Data.
####################################

#################################
## Temperatures in Brazil data ##
#################################

## Treat NA data:
  tempdata <- hdMTD::tempdata %>% filter(DATE >= "2010-08-05") # Remove days before "2010-08-05"
  # Identify remaining 155 NA positions
  posNA <- which(is.na(tempdata$MAXTEMP))
  # Fill short sequences of NAs ( ≤6 NAs) using nearest neighbors mean
  for (i in posNA) {
      if(!is.na(tempdata$MAXTEMP[i-1]) && !all(is.na(tempdata$MAXTEMP[(i+1):(i+6)]))) {
          aux <- which(!is.na(tempdata$MAXTEMP[ (i+1):(i+6) ]))[1]
          tempdata$MAXTEMP[i] <- mean(tempdata$MAXTEMP[i-1],tempdata$MAXTEMP[i+aux])
      }
  }
  # Identify remaining 88 NA positions (sequences >=7 NAs)
  posNA <- which(is.na(tempdata$MAXTEMP))
  # Fill remaining NA with the mean of the previous hour temperature, next hour 
  # temperature, and same hour of previous day temperature (any is this mean NA
  # will be ignored).
  for (i in posNA) {
      tempdata$MAXTEMP[i] <- mean(c(tempdata$MAXTEMP[i - 1], tempdata$MAXTEMP[i + 1],
          tempdata$MAXTEMP[i - 24]),na.rm = TRUE)
  }
  
  
## Compute mean daily maximum temperatures:
  temp <- tempdata %>%
      group_by(DATE) %>%
      summarize(MAXTEMP = mean(MAXTEMP), .groups = 'drop')
  head(temp, 4)


## Create categories of temperature:
  # Transform the continuous range of temperature data into categorical data by
  # dividing the range into equal intervals.
  xn <- max(temp$MAXTEMP)
  x1 <- min(temp$MAXTEMP)
  maxAmp <- xn - x1
  temp$MAXTEMP1  <- ifelse(temp$MAXTEMP < x1 + maxAmp/2, 1, 2)
  head(temp, 4)
  prop.table(table(temp$MAXTEMP1))

## Run FS for Temp12:
  Temp12 <- rev(temp$MAXTEMP1)
  # hdMTD functions assume the sample is sorted from the latest observations to
  # oldest.
  
  # WARNING: hdMTD_FS bellow ~7min (on i7-1255U, 10 cores)
  hdMTD_FS(Temp12, d = 400, l = 3)
  # 1 364 6

  # Reduce maximum order to improve estimation
  # WARNING: hdMTD_FS bellow ~6min (on i7-1255U, 10 cores)
  hdMTD_FS(Temp12, d = 364, l = 3)
  # 1 364 6

## Split sample in Train and Test data:
  ndays <- nrow(temp %>%
      filter(DATE >= "2023-09-01")) # 366 days in the latest year of the sample
  Temp12_Train <- Temp12[-seq_len(ndays)] # Training data with 4775 days
  Temp12_Test <- Temp12[seq_len(ndays)] # Test data

## Rerun FS for Train data:
  # WARNING: hdMTD_FS bellow ~6min (on i7-1255U, 10 cores)
  hdMTD_FS(Temp12_Train, d = 364, l = 3)
  # 1 364 6
  
  # Trim out irrelevant lags:
  hdMTD_CUT(Temp12_Train, d = 364, S = c(1, 364, 6))
  # 364, 6, 1
  
  hdMTD_BIC(Temp12_Train, d = 364, S = c(1, 364, 6), minl = 1, maxl = 3,
            byl = TRUE, BICvalue = TRUE )
  #        1             1,364           1,6,364 smallest: 1,6,364 
  # 1720.801          1690.543          1674.080          1674.080
  
  # WARNING: hdMTD_FSC bellow ~3min (on i7-1255U, 10 cores)
  hdMTD_FSC(Temp12_Train, d = 364, l = 3)
  # 364  6   1
  
  P_FS <- probs(Temp12_Train, S = c(1, 6, 364), matrixform = T)
  P_FS
  #              1         2
  # 111 0.86626140 0.1337386
  # 112 0.24736842 0.7526316
  # 121 0.77157360 0.2284264
  # 122 0.13318777 0.8668122
  # 211 0.78846154 0.2115385
  # 212 0.10972569 0.8902743
  # 221 0.57506361 0.4249364
  # 222 0.07283555 0.9271645

## Classic method for choosing relevant lag set:
  ct <- countsTab(Temp12_Train, d = 6) # Table with size 6 sequence counts
  head(ct,4)

  # MC(1)
  ft <- freqTab(S = 1, A = c(1, 2), countsTab = ct)
  LL <- sum(log(ft$qax_Sj)*ft$Nxa_Sj)
  freeParam <- 2 * 1
  BICMC1 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
  # Comparable BIC if the model is a Markov chain of order 1: 1869.162
  
  # MC(1,2)
  ft <- freqTab(S = c(1, 2), A = c(1, 2), countsTab = ct)
  LL <- sum(log(ft$qax_Sj)*ft$Nxa_Sj)
  freeParam <- 2^2 * 1
  BICMC2 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
  # Comparable BIC if the model is a Markov chain of order 2: 1850.598
  
  # MC(1,2,3)
  ft <- freqTab(S = c(1, 2, 3), A = c(1, 2), countsTab = ct)
  LL <- sum(log(ft$qax_Sj)*ft$Nxa_Sj)
  freeParam <- 2^3 * 1
  BICMC3 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
  # Comparable BIC if the model is a Markov chain of order 3: 1854.029
  
  # MC(1,2,3,4)
  ft <- freqTab(S = c(1, 2, 3, 4), A = c(1, 2), countsTab = ct)
  LL <- sum(log(ft$qax_Sj)*ft$Nxa_Sj)
  freeParam <- 2^4 * 1
  BICMC4 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
  # Comparable BIC if the model is a Markov chain of order 4: 1877.888
  
  # Helper function
  prodinf <- function(x, y){
    # This function simply allows the product: -inf * 0 = 0, which will be useful
    # for the next likelihood calculus. With it, log(ft$qax_Sj) * ft$Nxa_Sj = 0
    # when ft$Nxa_Sj = 0.
    prinf <- numeric(length(x))
    for (i in seq_len(length(x))) {
      if (is.infinite(x[i]) && y[i] == 0) {
        prinf[i] <- 0
      } else {
        prinf[i] <- x[i] * y[i]
      }
    }
    prinf
  }
  
  # MC(1,2,3,4,5)
  ft <- freqTab(S = c(1, 2, 3, 4, 5), A = c(1, 2), countsTab = ct)
  LL <- sum(prodinf(log(ft$qax_Sj), ft$Nxa_Sj))
  freeParam <- 2^5 * 1
  BICMC5 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
  # Comparable BIC if the model is a Markov chain of order 5: 1925.962
  
  # MC(1,2,3,4,5,6)
  ft <- freqTab(S = c(1, 2, 3, 4, 5, 6), A = c(1, 2), countsTab = ct)
  LL <- sum(prodinf(log(ft$qax_Sj), ft$Nxa_Sj))
  freeParam <- 2^6 * 1
  BICMC6 <- -LL + 0.5 * log(length(Temp12_Train)) * freeParam
  # Comparable BIC if the model is a Markov chain of order 6: 2031.679
  
  which.min(c(BICMC1, BICMC2, BICMC3, BICMC4, BICMC5, BICMC6))
  # The classic method chooses order 2 (S=c(1, 2)).

  P_MC2 <- probs(Temp12_Train, S = c(1, 2), matrixform = TRUE)
  P_MC2
  #             1         2
  # 11 0.77813505 0.2218650
  # 12 0.16326531 0.8367347
  # 21 0.60233918 0.3976608
  # 22 0.09064976 0.9093502

# Independent model:
  P_Ind <- prop.table(table(Temp12_Train))
  P_Ind
  # Temp12_Train
  #         1         2
  # 0.2672251 0.7327749

## Comparing methods:
  # Store information for computing performance metrics
  Days1 <- which(Temp12_Test == 1) 
  lenDays1 <- length(Days1) # 64 days
  lenDays1/ndays # 0.1748634 = frequency o 1 days in test sample
    
  Temp12_Test <- c(Temp12_Test, Temp12_Train[seq_len(364)]) 
  # Add latest 364 pasts (of Training sample) to test data to facilitate
  # computations
  
  # Monte Carlo Simulation
  set.seed(1)
  
  # Initiation
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

  
  for (j in seq_len(1000)){ # 1000 replications
    u <- runif(ndays) # 366 uniforms for each replication
  
    predInd <- numeric(ndays) # Stores predictions of Ind method
    predMC2 <- numeric(ndays) # Stores predictions of MC2 method
    predFS <- numeric(ndays) # Stores predictions of FS method
    
    for (i in ndays:1) {
      # Model Ind:
        predInd[i] <- ifelse(u[i] <= P_Ind[1], 1, 2)
      # Model MC2:
        pastRow <- which(pasts2 == paste0(rev(Temp12_Test[c(i+1, i+2)]),
            collapse = ""))
        # Why reverse? Temp12_Test most recent data is to the left, but the in the
        # matrices the sequences most recent symbol is to the right.
        predMC2[i] <- ifelse(u[i] <= P_MC2[pastRow, 1], 1, 2)
      # Model FS:
        pastRow <- which(pastsFS == paste0(rev(Temp12_Test[c(i+1, i+6, i+364)]),
            collapse = ""))
        predFS[i] <- ifelse(u[i] <= P_FS[pastRow, 1], 1, 2)
    }

    # Store data for computing metrics:
    # TP + TN (True Positives + True Negatives = total hits)
    hitInd[j] <- sum(predInd == Temp12_Test[seq_len(ndays)])
    hitMC2[j] <- sum(predMC2 == Temp12_Test[seq_len(ndays)])
    hitFS[j] <- sum(predFS == Temp12_Test[seq_len(ndays)])
  
    # TP
    T1Ind[j] <- sum(predInd[Days1] == 1)
    T1MC2[j] <- sum(predMC2[Days1] == 1)
    T1FS[j] <- sum(predFS[Days1] == 1)
  
    # FP (False Positives)
    F1Ind[j] <- sum(predInd[-Days1] == 1)
    F1MC2[j] <- sum(predMC2[-Days1] == 1)
    F1FS[j] <- sum(predFS[-Days1] == 1)
  }

  # Computing metrics:
  # Mean of Accuracies
  AccInd <- mean(hitInd)/ndays # 0.6509262
  AccMC2 <- mean(hitMC2)/ndays # 0.8236585
  AccFS <- mean(hitFS)/ndays   # 0.8349153
  # ndays = TN + TP + FN + FP

  # Mean of Precisions
  PrecInd <- mean(T1Ind/(T1Ind + F1Ind)) # 0.1745824
  PrecMC2 <- mean(T1MC2/(T1MC2 + F1MC2)) # 0.4981026
  PrecFS <- mean(T1FS/(T1FS + F1FS)) # 0.5249133
  
  # Mean of Sensitivities//Recall
  SensInd <- mean(T1Ind)/lenDays1 # 0.2671875
  SensMC2 <- mean(T1MC2)/lenDays1 # 0.5870313
  SensFS <- mean(T1FS)/lenDays1 # 0.6294531
  # lenDays1 = TP + FN
  
  # Mean of Specificity
  SpecInd <- 1 - mean(F1Ind)/(ndays - lenDays1) # 0.7322483
  SpecMC2 <- 1 - mean(F1MC2)/(ndays - lenDays1) # 0.8738046
  SpecFS <- 1 - mean(F1FS)/(ndays - lenDays1) # 0.878457
  
  # F1 - Score
  F1ScoreInd <- 2 * (PrecInd * SensInd) / (PrecInd + SensInd) # 0.2111789
  F1ScoreMC2 <- 2 * (PrecMC2 * SensMC2) / (PrecMC2 + SensMC2) # 0.538923
  F1ScoreFS <- 2 * (PrecFS * SensFS) / (PrecFS + SensFS) # 0.5724496

  
## FS sequential selection based on ν_{n,j,S} values:
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
          dplyr::group_by_at(paste0("x", j)) %>%
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
          dplyr::group_by_at(paste0("x", Sj)) %>%
          dplyr::summarise(Nx_Sj = sum(Nxa_Sj), .groups = "drop")
        b_S <- b_Sja %>%
          dplyr::group_by_at(paste0("x", S)) %>%
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
          dplyr::group_by_at(paste0("x", Sj)) %>%
          dplyr::summarise(Nx_Sj = sum(Nxa_Sj), .groups = "drop")
        b_S <- b_Sja %>%
          dplyr::group_by_at(paste0("x", dec_S)) %>%
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
  
  # WARNING: run_sequential_lag_selection bellow ~6min (on i7-1255U, 10 cores)
  results <- run_sequential_lag_selection(Temp12_Train)

  
## Plots:
  # Time series with quarterly mean of daily maximum temperatures
  # (Generating Article Figure 2)
  # Data arrangement
  TRIM_DATA <- temp %>%
    mutate(
      Y_TRIMESTER = paste0(year(DATE), "-T", quarter(DATE))
    ) %>%
    group_by(Y_TRIMESTER) %>%
    summarise(
      MEAN_TEMP = mean(MAXTEMP),
      DATA_REF = min(DATE)  # For order in graph
    ) %>%
    ungroup() %>%
    arrange(DATA_REF)
  TRIM_DATA <- TRIM_DATA[-c(1, nrow(TRIM_DATA)),]
  
# Plotting
  #pdf("Timetemp_plot.pdf", width = 9, height = 6)
  ggplot(TRIM_DATA, aes(x = DATA_REF, y = MEAN_TEMP)) +
    geom_line(color = "steelblue", linewidth = 0.5) +
    geom_point(color = "steelblue", size = 1.2) +
    scale_x_date(
      date_breaks = "1 year",
      date_labels = "%Y",
      minor_breaks = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15), 
      plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),
    ) +
    labs(
      title = "Quarterly mean of daily maximum temperatures across the years",
      x = "Year",
      y = "Mean Temperature (°C)"
    )
  #dev.off()
  
  # Empirical Nu analysis
  # (Generating Article Figure 3)
  #pdf("nu_plot.pdf", width = 11, height = 6)
  par(mfrow = c(1, 3), mar = c(5, 6, 4, 2), oma = c(0, 0, 4, 0))
  palette <- c("#E41A1C", "#377EB8", "#4DAF4A")
  with(results, {
    # Graph 1
    Sc <-  364:1
    plot(1:364, rev(nuj1), type = "p", pch = 19, cex = 0.8, col = "gray70",
         ylab = "", xlab = "Lag (-j)", cex.lab = 1.8, cex.axis = 1.3,
         ylim = c(0,0.13), main = "", panel.first = grid())
    title(main = expression(paste("S = ", Ø)), cex.main = 1.5, font.main = 1)
    title(ylab = expression(widehat(nu)[n*","*j*","*S]), line = 4, cex.lab = 1.8) 
    points(Sc[which.max(nuj1)], nuj1[which.max(nuj1)], pch = 21, bg = palette[1], cex = 1.5, lwd = 1)
    text(Sc[which.max(nuj1)], nuj1[which.max(nuj1)], labels = paste0(Sc[which.max(nuj1)]),
         pos = 3, col = palette[1], font = 2, cex = 1.4)
    # Graph 2
    Sc <- 364:2
    plot(2:364, rev(nuj2), type = "p", pch = 19, cex = 0.8, col = "gray70",
         ylab = "", xlab = "Lag (-j)", cex.lab = 1.8, cex.axis = 1.3,
         ylim = c(0,0.025),  main = "", panel.first = grid())
    title(main = "With S = {-1}", cex.main = 1.5, font.main = 1)
    title(ylab = expression(widehat(nu)[n*","*j*","*S]), line = 4, cex.lab = 1.8)
    points(Sc[which.max(nuj2)], nuj2[which.max(nuj2)], pch = 21, bg = palette[2], cex = 1.5, lwd = 1)
    text(Sc[which.max(nuj2)]-5, nuj2[which.max(nuj2)], labels = paste0(Sc[which.max(nuj2)]),
         pos = 3, col = palette[2], font = 2, cex = 1.4)
    # Graph 3
    Sc <- 363:2
    plot(2:363, rev(nuj3), type = "p", pch = 19, cex = 0.8, col = "gray70",
         ylab = "", xlab = "Lag (-j)", cex.lab = 1.8, cex.axis = 1.3,
         ylim = c(0,0.025), main = "", panel.first = grid())
    title(main = "With S = {-364, -1}", cex.main = 1.5, font.main = 1)
    title(ylab = expression(widehat(nu)[n*","*j*","*S]), line = 4, cex.lab = 1.8)
    points(Sc[which.max(nuj3)], nuj3[which.max(nuj3)], pch = 21, bg = palette[3], cex = 1.5, lwd = 1)
    text(Sc[which.max(nuj3)], nuj3[which.max(nuj3)], labels = paste0(Sc[which.max(nuj3)]),
         pos = 3, col = palette[3], font = 2, cex = 1.4)
    
    mtext(expression(paste("Sequential lag selection based on ", widehat(nu)[n*","*j*","*S])),
          outer = TRUE, cex = 1.6, font = 2, line = 1.4)
  })
  #dev.off()
  