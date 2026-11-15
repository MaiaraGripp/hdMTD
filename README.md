# hdMTD

**Inference for High-Dimensional Mixture Transition Distribution Models**

`hdMTD` is an R package for inference in mixture transition distribution (MTD) models, a parsimonious class of high-order Markov chains. It provides methods for selecting relevant past lags, estimating transition probabilities and model parameters, and generating exact samples from the stationary distribution.

In particular, its forward stepwise (FS) and forward stepwise and cut (FSC) methods support lag selection when dependencies extend far into the past.

## Installation

Install the current CRAN release:

```r
install.packages("hdMTD")
```

Load the package:

```r
library(hdMTD)
```

## Main features

- **Relevant lag selection:** Bayesian information criterion (BIC), forward stepwise (FS), cut (CUT), and forward stepwise and cut (FSC) methods, implemented in `hdMTD_BIC()`, `hdMTD_FS()`, `hdMTD_CUT()`, and `hdMTD_FSC()`. The `hdMTD()` function offers a common interface.
- **Model specification:** Define an MTD model using `MTDmodel()`.
- **Exact simulation:** Generate stationary MTD samples using `perfectSample()`.
- **Transition probability estimation:** Compute empirical conditional probabilities using `empirical_probs()`.
- **Parameter estimation:** Estimate MTD parameters with the expectation-maximization algorithm using `MTDest()`.
- **Lag relevance:** Compute or estimate lag oscillations using `oscillation()`.

## Quick example

```r
library(hdMTD)

set.seed(11)

# Define an MTD model with relevant lags -1, -15, and -30.
model <- MTDmodel(
  Lambda = c(1, 15, 30),
  A = c(0, 1),
  lam0 = 0.01,
  lamj = c(0.39, 0.30, 0.30),
  p0 = c(0.5, 0.5)
)

# Generate an exact sample from its stationary distribution.
X <- perfectSample(model, N = 1000)

# Select four candidate lags among the 40 most recent lags.
hdMTD_FS(X, d = 40, l = 4)
```

The input `Lambda` uses **positive lag distances** (for example, `30` represents lag `-30`).

## Documentation and replication materials

- [Replication Materials for hdMTD (CRAN vignette)](https://cran.r-project.org/web/packages/hdMTD/vignettes/replication-materials-hdMTD.html) — reproducible code and supporting results for the examples and analyses in Section 5 of the article.
- [Package reference manual (CRAN)](https://cran.r-project.org/web/packages/hdMTD/hdMTD.pdf) — documentation of package functions and arguments.
- [CRAN package page](https://cran.r-project.org/package=hdMTD) — releases and installation information.
- [Article preprint (arXiv)](https://arxiv.org/abs/2509.01808) — *hdMTD: An R Package for High-Dimensional Mixture Transition Distribution Models*.

## License

GPL-3.
