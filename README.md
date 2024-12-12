
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hdMTD

## Overview

hdMTD is a package designed for conducting inference in Mixture
Transition Distribution (MTD) models, which can represent higher-order
Markov models as a convex mixture of single-step Markov chains. The
algorithms employed here are primarily built upon the work of [Ost and
Takahashi (2023)](http://jmlr.org/papers/v24/22-0266.html) in their
paper titled “Sparse Markov Models for High-Dimensional Inference.”

## Installation

``` r
remotes::install_github("MaiaraGripp/hdMTD")
```

## Usage

Given a sample from an MTD chain, the hdMTD() function estimates the
relevant lag set $\Lambda$ using a specified method and a suitable set
of parameters. The available methods include the *FSC* algorithm
developed by [Ost and Takahashi
(2023)](http://jmlr.org/papers/v24/22-0266.html), FSC’s first step,
referred to as *FS* method, its second step, referred to as *CUT*
method, and a Bayesian Information Criterion (*BIC*) method.
Furthermore, the package can perform the following tasks:

- Create an object of class MTD with all the parameters necessary to
  define a proper MTD model using the `MTDmodel()` function. Such
  objects can be used to sample an MTD Markov Chain from its invariant
  distribution with the `perfectSample()` function.
- Calculate the oscillations of an MTD object or estimate them from a
  chain sample using the `oscillations()` function.
- Compute the Maximum Likelihood Estimates (MLE) of the transition
  probabilities for a chain sample given a relevant lag set with the
  `probs()` function.
- Determine the absolute frequency of each sequence of size $d$ that
  appears in a sample using the `countsTab()` function.
- Group such frequencies to consider only sequences time-indexed by a
  subset of ${-1, \dots, -d}$ with the `freqTabSj()` function.
- Lastly, `MTDest()` function uses a Expectation Maximization (EM)
  algorithm based on [Lèbre and Bourguinon
  (2008)](https://doi.org/10.1080/00949650701266666), to estimate the
  parameter set of an MTD model given a sample and a set of initial
  parameters.

``` r

library(hdMTD)

set.seed(1234)
Lambda <- c(1,5,10) #relevant lag set
A <- c(0,1) #state space
MTD <- MTDmodel(Lambda = Lambda,A=A);MTD # MTD object
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

oscillation(MTD)
#>         -1         -5        -10 
#> 0.06779938 0.10684048 0.05927337

X <- perfectSample(MTD=MTD,N=1000);X[1:10] # MTD chain sampled from invariant dist.
#>  [1] 1 0 0 0 1 0 1 0 0 1

S1 <- hdMTD(X=X,d=15,method = "FS",l=3);S1 # estimate a relevant lag set of size 3 with the "Foward Stepwise" method.
#> [1]  5 10  1

S2 <- hdMTD(X,d=max(S1),method = "BIC",S=S1, minl=1, maxl=3);S2# estimate a relevant lag set with the "BIC" method using FS output.
#> [1] 5

S3 <- hdMTD(X,d=12,method = "BIC", minl=1, maxl=3, byl=TRUE, BICvalue=TRUE);S3
#>           5        5,10      5,7,10 smallest: 5 
#>    668.7065    675.4400    682.0869    668.7065
# estimate a relevant lag set with the "BIC" method.

S4 <- hdMTD(X,d=20,method = "CUT",S=S1);S4# estimate a relevant lag set with the "CUT" method using FS output.
#> [1] 10  5  1

p <- probs(X,S=S4, matrixform = TRUE);p #estimates the MLE transition probabilities given the CUT output.
#>             0         1
#> 000 0.4062500 0.5937500
#> 001 0.3611111 0.6388889
#> 010 0.3086420 0.6913580
#> 011 0.3358779 0.6641221
#> 100 0.4444444 0.5555556
#> 101 0.5528455 0.4471545
#> 110 0.3245033 0.6754967
#> 111 0.3884298 0.6115702

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

MTDest(X,S=c(1,5,10),init=init, iter = TRUE)
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
## Estimating P (transition probabilities)
estParam <- MTDest(X,S=c(1,15,30),init=init)
estMTD <- MTDmodel(Lambda,A,lam0=estParam$lambdas[1],
         lamj = estParam$lambdas[-1],p0=estParam$p0,
         pj=estParam$pj)
estMTD$P
#>             0         1
#> 000 0.3205668 0.6794332
#> 001 0.3759126 0.6240874
#> 010 0.4008207 0.5991793
#> 011 0.4561665 0.5438335
#> 100 0.2868808 0.7131192
#> 101 0.3422266 0.6577734
#> 110 0.3671347 0.6328653
#> 111 0.4224804 0.5775196
```

## Data sets

This package includes three real-data sets that where acquired
externaly.

The raindata set was acquired at
[Kaggle](https://www.kaggle.com/datasets/jsphyg/weather-dataset-rattle-package)
but it’s original source is Data source:\[Australian Bureau of
Meteorology\] (<http://www.bom.gov.au/climate/dwo/>) and
(<http://www.bom.gov.au/climate/data>). Copyright Commonwealth of
Australia 2010, Bureau of Meteorology.

The sleepscoring data was acquired at the Haaglanden Medisch Centrum
(HMC, The Netherlands) sleep center, and is available at Physionet.org.
URL <https://doi.org/10.13026/t79q-fr32>.

The tempdata data set was acquired at INMET, the National Institute of
Meteorology in Brazil <https://bdmep.inmet.gov.br/>.

## References

If you use our package in your work, please cite the following articles:

Ost G, Takahashi DY (2023). “Sparse Markov Models for High-dimensional
Inference.” Journal of Machine Learning Research, 24(279), 1–54. URL
<http://jmlr.org/papers/v24/22-0266.html>

Function MTDest was based on the ideas of the following article:

Lebre, Sophie & Bourguignon, Pierre-Yves. (2008). An EM algorithm for
estimation in the Mixture Transition Distribution model. Journal of
Statistical Computation and Simulation. 78. [Link to
Article](https://doi.org/10.1080/00949650701266666)
