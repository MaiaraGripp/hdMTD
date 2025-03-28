
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hdMTD

## Overview

hdMTD is an R programming language package for the estimation of
parameters in Mixture Transition Distribution (MTD) models. An MTD is a
Markov chain with transition probabilities represented as a convex
combination of conditional distributions. Given a sample from an MTD
chain, hdMTD can estimate the relevant past set using various methods,
such as Bayesian Information Criterion (BIC) and the Forward Stepwise
and CUT (FSC) algorithm, which is efficient in estimating the set of
relevant lags even in high-dimension, i.e., when dependence extends to
pasts so distant that they may be close to the sample size. The package
also computes the Maximum Likelihood Estimate (MLE) for transition
probabilities, estimates oscillations, determines MTD parameters through
the Expectation-Maximization (EM) algorithm, and can also simulate an
MTD sample from its invariant distribution using the perfect sample
algorithm.

## Installation

``` r
remotes::install_github("MaiaraGripp/hdMTD")
```

## **Usage**

Given a sample from an MTD chain, the `hdMTD()` function estimates the
relevant lag set $\Lambda$ using a specified method and suitable
parameters. The available methods include:

- The *FSC* algorithm developed by [Ost and Takahashi
  (2023)](http://jmlr.org/papers/v24/22-0266.html).  
- The first step of *FSC*, referred to as the *FS* method.  
- The second step of *FSC*, known as the *CUT* method.  
- A Bayesian Information Criterion ( *BIC* ) approach.

Additionally, the package provides the following functionalities:

### **Model Definition and Sampling**

- Create an `MTD` object with all parameters necessary to define an MTD
  model using `MTDmodel()`.  
- Generate a perfect sample from the invariant distribution of an MTD
  Markov chain using `perfectSample()`.

### **Estimation and Computation**

- Compute the oscillations of an `MTD` object or estimate them from a
  chain sample using `oscillations()`.  
- Obtain the Maximum Likelihood Estimates (MLE) of transition
  probabilities for a given set of lags using `probs()`.  
- Compute the absolute frequency of all observed sequences of length $d$
  in a sample using `countsTab()`. Aggregate these frequencies to
  consider only sequences indexed by subsets of $\{1, \dots, d\}$ using
  `freqTab()`.

### **Parameter Estimation**

- Estimate the parameters of an MTD model from a sample using
  `MTDest()`, which implements an Expectation-Maximization (EM)
  algorithm based on [Lèbre and Bourguinon
  (2008)](https://doi.org/10.1080/00949650701266666).

``` r

library(hdMTD)
set.seed(1234)

## 1. Simulating an MTD Markov Chain:

Lambda <- c(1, 5, 10) # Set of relevant lags Λ = {-10, -5, -1}
A <- c(0, 1) # State space

# Create an MTD model
MTD <- MTDmodel(Lambda = Lambda, A = A)
MTD # An MTD class object
#> $P
#>             0         1
#> 000 0.3796866 0.6203134
#> 001 0.4474859 0.5525141
#> 010 0.2728461 0.7271539
#> 011 0.3406455 0.6593545
#> 100 0.4389599 0.5610401
#> 101 0.5067593 0.4932407
#> 110 0.3321195 0.6678805
#> 111 0.3999188 0.6000812
#> 
#> $lambdas
#>      lam0     lam-1     lam-5    lam-10 
#> 0.2228608 0.2280200 0.3149060 0.2342131 
#> 
#> $pj
#> $pj$`p-1`
#>            0         1
#> 0 0.01405572 0.9859443
#> 1 0.31139528 0.6886047
#> 
#> $pj$`p-5`
#>           0         1
#> 0 0.7104103 0.2895897
#> 1 0.3711330 0.6288670
#> 
#> $pj$`p-10`
#>           0         1
#> 0 0.5052655 0.4947345
#> 1 0.7583400 0.2416600
#> 
#> 
#> $p0
#>     p0(0)     p0(1) 
#> 0.1544877 0.8455123 
#> 
#> $Lambda
#> [1]  1  5 10
#> 
#> $A
#> [1] 0 1

# Compute oscillations for the MTD model
oscillation(MTD)
#>         -1         -5        -10 
#> 0.06779938 0.10684048 0.05927337

# Generate a sample from the MTD
X <- perfectSample(MTD = MTD, N = 1000)
X[1:10] # Display the last 10 sampled states (X[1] is the latest sampled state).
#>  [1] 1 0 0 0 1 0 1 0 0 1

## 2. Inference on MTD Markov Chains

oscillation(X, S = c(1, 10)) # Estimated oscillations for lags -10 and -1
#>         -1        -10 
#> 0.04735179 0.06575918

# Estimate the set of relevant lags using the "Forward Stepwise" (FS) method
S1 <- hdMTD(X = X, d = 15, method = "FS", l = 3)
S1 # Estimated lags (size l=3)
#> [1]  5 10  1

# Alternative equivalent function:
# S1 <- hdMTD_FS(X = X, d = 15, l = 3)

# Refining the estimated lags using "BIC"
S2 <- hdMTD(X, d = max(S1), method = "BIC", S = S1, minl = 1, maxl = 3)
S2 # Estimated set of relevant lags with the "BIC" method using "FS" output.
#> [1] 5

# Alternative equivalent function:
# S2 <- hdMTD_BIC(X, d = max(S1), S = S1, minl = 1, maxl = 3)

# "BIC" estimation for sets of relevant lags
S3 <- hdMTD(X, d = 12, method = "BIC", minl = 1, maxl = 3, byl = TRUE, BICvalue = TRUE)
S3 # Sets of relevant lags (subsets of 1:12) with sizes from size 1 to 3 with lowest BIC.
#>           5        5,10      5,7,10 smallest: 5 
#>    668.7065    675.4400    682.0869    668.7065

# "CUT" estimation for sets of relevant lags given S1
S4 <- hdMTD(X, d = 20, method = "CUT", S = S1, alpha = 0.1)
S4 # Estimated set of relevant lags with the "CUT" method using "FS" output.
#> [1] 10  5

# # "FSC" method combining FS and CUT
S5 <- hdMTD(X, d = 20, method = "FSC", l=3, alpha = 0.01);
S5 # Estimated set of relevant lags with the "FSC".
#> [1] 5

# Validation: S5 should match the result of running "FS" on the first half of the sample
# and "CUT" on the second half
all.equal(S5, hdMTD_CUT(X[501:1000], d = 20,
                    S = hdMTD_FS(X[1:500], d = 20, l=3), alpha = 0.01))
#> [1] TRUE

## 3. Estimating Transition Probabilities

# Estimate the transition probability matrix given a sample and set of relevant lags
p <- probs(X, S = S4, matrixform = TRUE)
p # MLE given the CUT output as set of relevant lags
#>            0         1
#> 00 0.3779070 0.6220930
#> 01 0.3254717 0.6745283
#> 10 0.5070423 0.4929577
#> 11 0.3638677 0.6361323

## 4. Estimating MTD Parameters with the Expectation-Maximization (EM) Algorithm

# Initial parameter values for the EM algorithm
init <- list(
  'lambdas'= c(0.05,0.3,0.3,0.35),
  'p0' = c(0.5,0.5),
  'pj' = list(
    matrix(c(0.5,0.5,
             0.5,0.5),ncol=2,nrow = 2),
    matrix(c(0.5,0.5,
             0.5,0.5),ncol=2,nrow = 2),
    matrix(c(0.5,0.5,
             0.5,0.5),ncol=2,nrow = 2)
  )
)

# Estimate the MTD model parameters
estParam <- MTDest(X,S=c(1,5,10),init=init, iter = TRUE); estParam
#> $lambdas
#>      lam-0      lam-1      lam-5     lam-10 
#> 0.04889442 0.29599651 0.30684622 0.34826285 
#> 
#> $pj
#> $pj$`p_-1`
#>           0         1
#> 0 0.2986150 0.7013850
#> 1 0.4432366 0.5567634
#> 
#> $pj$`p_-5`
#>           0         1
#> 0 0.5928177 0.4071823
#> 1 0.2647329 0.7352671
#> 
#> $pj$`p_-10`
#>           0         1
#> 0 0.2658844 0.7341156
#> 1 0.4675134 0.5324866
#> 
#> 
#> $p0
#>    p_0(0)    p_0(1) 
#> 0.3845184 0.6154816 
#> 
#> $iterations
#> [1] 9
#> 
#> $distlogL
#>  [1] 28.863219344  2.133782102  1.082393935  0.549555520  0.279169153
#>  [6]  0.141869452  0.072120503  0.036675669  0.018657448  0.009494794

# Build the estimated MTD model using the inferred parameters
estMTD <- MTDmodel(Lambda, A,
                   lam0 = estParam$lambdas[1],
                   lamj = estParam$lambdas[-1],
                   p0 = estParam$p0,
                   pj = estParam$pj)
# Display the estimated transition probability matrix
estMTD$P
#>             0         1
#> 000 0.3816913 0.6183087
#> 001 0.4244988 0.5755012
#> 010 0.2810198 0.7189802
#> 011 0.3238273 0.6761727
#> 100 0.4519112 0.5480888
#> 101 0.4947187 0.5052813
#> 110 0.3512397 0.6487603
#> 111 0.3940471 0.6059529
```

## Data sets

This package includes **three real-world data sets** acquired from
external sources and **one simulated data set** (`testChains`), which
was generated using the `perfectSample` function.

- **`raindata`**: This dataset was obtained from
  [Kaggle](https://www.kaggle.com/datasets/jsphyg/weather-dataset-rattle-package),
  but its original source is the **Australian Bureau of Meteorology**  
  ([Data Source](http://www.bom.gov.au/climate/dwo/), [Climate
  Data](http://www.bom.gov.au/climate/data)).  
  **Copyright:** Commonwealth of Australia 2010, Bureau of Meteorology.

- **`sleepscoring`**: This dataset was collected at the **Haaglanden
  Medisch Centrum (HMC, The Netherlands) Sleep Center** and is publicly
  available on **PhysioNet**  
  ([DOI: 10.13026/t79q-fr32](https://doi.org/10.13026/t79q-fr32)).

- **`tempdata`**: This dataset was obtained from **INMET** (National
  Institute of Meteorology, Brazil) and is available at [INMET Data
  Portal](https://bdmep.inmet.gov.br/).

------------------------------------------------------------------------

## References

The methods implemented in this package are based on the following
works:

- **“FS”, “CUT”, and “FSC” methods:**  
  These functions are based on: **Ost G, Takahashi DY (2023).** *Sparse
  Markov Models for High-Dimensional Inference.*  
  *Journal of Machine Learning Research, 24(279), 1–54.*  
  [URL:
  http://jmlr.org/papers/v24/22-0266.html](http://jmlr.org/papers/v24/22-0266.html)

- **“BIC” method:**  
  This function is based on the well-known **Bayesian Information
  Criterion (BIC) method** for model selection:  
  **Imre Csiszár & Paul C. Shields (2000).** *The Consistency of the BIC
  Markov Order Estimator.*  
  *The Annals of Statistics, 28(6), 1601-1619.*  
  [DOI: 10.1214/aos/1015957472](https://doi.org/10.1214/aos/1015957472)

- **`MTDest` function (Expectation-Maximization method):**  
  The `MTDest` function applies the **Expectation-Maximization (EM)
  algorithm** for parameter estimation in **Mixture Transition
  Distribution (MTD) models**, based on the following paper:  
  **Lebre, Sophie & Bourguignon, Pierre-Yves (2008).** *An EM Algorithm
  for Estimation in the Mixture Transition Distribution Model.*  
  *Journal of Statistical Computation and Simulation, 78.*  
  [DOI:
  10.1080/00949650701266666](https://doi.org/10.1080/00949650701266666)
