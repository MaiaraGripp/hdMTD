
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hdMTD

## Overview

hdMTD is a package designed for conducting inference in Mixture
Transition Distribution (MTD) models, which can represent higher-order
Markov models as a convex mixture of single-step Markov chains. The
algorithms employed here are primarily built upon the work of [Ost and
Takahashi (2023)](https://arxiv.org/abs/2202.08007) in their paper
titled “Sparse Markov Models for High-Dimensional Inference.”

## Installation

``` r
remotes::install_github("MaiaraGripp/hdMTD")
```

## Usage

Given a sample from an MTD chain, the hdMTD() function estimates the
relevant lag set $\Lambda$ using a specified method and a suitable set
of parameters. The available methods include the *FSC* algorithm
developed by [Ost and Takahashi
(2023)](https://arxiv.org/abs/2202.08007), FSC’s first step, referred to
as *FS* method, its second step, referred to as *CUT* method, and a
Bayesian Information Criterion (*BIC*) method. Furthermore, the package
can perform the following tasks:

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
  (2008)](https://arxiv.org/abs/0803.0525), to estimate the parameter
  set of an MTD model given a sample and a set of initial parameters.

``` r

library(hdMTD)

set.seed(1234)
Lambda <- c(1,10) #relevant lag set
A <- c(0,1) #state space
MTD <- MTDmodel(Lambda = Lambda,A=A);MTD # MTD object
#> $P
#>            0         1
#> 00 0.4630341 0.5369659
#> 01 0.2487901 0.7512099
#> 10 0.5555227 0.4444773
#> 11 0.3412787 0.6587213
#> 
#> $lambdas
#>      lam0     lam-1    lam-10 
#> 0.2910220 0.2977591 0.4112189 
#> 
#> $p_j
#> $p_j$`p_-1`
#>            0         1
#> 0 0.73357674 0.2664233
#> 1 0.01405572 0.9859443
#> 
#> $p_j$`p_-10`
#>           0         1
#> 0 0.4854971 0.5145029
#> 1 0.7104103 0.2895897
#> 
#> 
#> $p0
#>    p_0(0)    p_0(1) 
#> 0.1544877 0.8455123 
#> 
#> $Lambda
#> [1]  1 10
#> 
#> $A
#> [1] 0 1

oscillation(MTD)
#> Calculating d<U+2096> = <U+03BB><U+2096>*max{b,c in A: dTV[p<U+2096>(.|b),p<U+2096>(.|c)]}, 
#> for each k in Lambda: 
#> 
#>         -1        -10 
#> 0.21424395 0.09248858

X <- perfectSample(MTD=MTD,N=1000);X[1:10] # MTD chain sampled from invariant dist.
#>  [1] 1 0 0 0 1 1 1 1 1 1

S1 <- hdMTD(X=X,d=15,method = "FS",l=3);S1 # estimate a relevant lag set of size 3 with the "Foward Stepwise" method.
#> [1]  1 10 13
S2 <- hdMTD(X,d=max(S1),method = "BIC",S=S1, minl=1, maxl=3, xi=0.3);S2# estimate a relevant lag set with the "BIC" method.
#> [1]  1 10

probs(X,S=S2) #estimates the MLE transition probabilities given S2.
#>   past_{ -10,-1 } a_0 p(a|past)
#> 1              00   0 0.5167785
#> 2              00   1 0.4832215
#> 3              01   0 0.2195122
#> 4              01   1 0.7804878
#> 5              10   0 0.5582329
#> 6              10   1 0.4417671
#> 7              11   0 0.3728324
#> 8              11   1 0.6271676
```

## Data sets

This package includes two real-data sets that where acquired externaly.
The raindata set was acquired at
[Kaggle](https://www.kaggle.com/datasets/jsphyg/weather-dataset-rattle-package)
but it’s original source is Data source:\[Australian Bureau of
Meteorology\] (<http://www.bom.gov.au/climate/dwo/>) and
(<http://www.bom.gov.au/climate/data>). Copyright Commonwealth of
Australia 2010, Bureau of Meteorology.

## References

If you use our package in your work, please cite the following article:

Function MTDest was based on the ideas of the following article:

Lebre, Sophie & Bourguignon, Pierre-Yves. (2008). An EM algorithm for
estimation in the Mixture Transition Distribution model. Journal of
Statistical Computation and Simulation. 78. [Link to
Article](https://doi.org/10.1080/00949650701266666)
