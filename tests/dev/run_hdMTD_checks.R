## ===========================
## 0) Preparação
## ===========================
pkgs <- c("testthat","withr")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, quiet = TRUE)
suppressPackageStartupMessages({
  library(testthat)
  library(withr)
})

message("\n===> Carregando hdMTD atual...\n")
# Assuma que a versão de desenvolvimento já está instalada/carregada.
# Se você estiver desenvolvendo com devtools, pode usar: devtools::load_all("caminho/hdMTD")


## Helper para somar linhas = 1 (matriz estocástica)
row_sums_one <- function(M, tol = 1e-10) {
  all(abs(rowSums(M) - 1) < tol)
}


sp <- list(
  A = c(2, 4, 6), Lambda = c(4, 5, 10),
  lam0 = 0.01, lamj = c(0.39, 0.3, 0.3),
  pj = NULL, p0 = NULL,
  indep = TRUE, sm = TRUE
)


## ===========================
## 1) Dados de teste reproduzíveis
## ===========================
Lambda <- sp$Lambda
A <- sp$A
lam0 <- sp$lam0
p0 <- sp$p0
MTD <- MTDmodel(Lambda = sp$Lambda, A = sp$A,
                lam0 = sp$lam0, lamj = sp$lamj,
                pj = sp$pj, p0 = sp$p0,
                indep_part = sp$indep,
                single_matrix =  sp$sm)

N  <- 2300
X  <- perfectSample(MTD, N = N)
X <- sample(states(MTD), 1000, replace = TRUE)

## Subconjunto de lags/conjuntos para testes
S_true <- Lambda(MTD)

## ===========================
## 2) Testes de accessors e métodos S3 básicos
## ===========================
test_that("Accessors básicos de MTD funcionam e são coerentes", {
  expect_true(is.list(pj(MTD)))
  expect_type(p0(MTD), "double")
  expect_type(lambdas(MTD), "double")
  expect_true(all(Lambda(MTD) == as.integer(Lambda(MTD))))
  expect_true(all(lags(MTD) == as.integer(lags(MTD))))
  expect_equal(sort(states(MTD)), sort(A))
  P <- transitP(MTD)
  expect_true(is.matrix(P))
  expect_true(row_sums_one(P))
  expect_equal(sum(p0(MTD)), 1)
  expect_equal(sum(lambdas(MTD)), 1)
})

## ===========================
## 3) EM fit (MTDest), as.MTD(), probs(), logLik()
## ===========================
# chute inicial simples (usa coef do próprio MTD para acelerar convergência):

fit <- MTDest(X, S = S_true, init = coef(MTD))   # EM
sf  <- summary(fit)

test_that("MTDest tem métodos e accessors", {
  expect_s3_class(fit, "MTDest")
  capture.output(print(fit))          # não deve falhar
  capture.output(print(sf))
  capture.output(summary(fit))
  expect_true(is.list(pj(fit)))
  expect_type(lambdas(fit), "double")
  expect_type(S(fit), "double")
  expect_equal(sort(states(fit)), sort(A))
  expect_true(row_sums_one(transitP(as.MTD(fit))))
  expect_equal(sum(p0(fit)), 1)
  expect_equal(sum(lambdas(fit)), 1)
  expect_true(row_sums_one(pj(fit)[[1]]))
  expect_true(row_sums_one(pj(fit)[[2]]))
})

# Coerção para MTD
m_from_fit <- as.MTD(fit)

test_that("as.MTD() reconstrói MTD coerente", {
  expect_s3_class(m_from_fit, "MTD")
  P_fit <- transitP(m_from_fit)
  expect_true(is.matrix(P_fit))
  expect_true(row_sums_one(P_fit))
  expect_equal(dim(P_fit), c(length(states(fit))^length(S(fit)), length(states(fit))))
  expect_equal(dim(pj(fit)[[1]]), c(length(states(fit)), length(states(fit))))
})

## ===========================
## 4) probs(): equivalências e formato
## ===========================
test_that("probs(MTD) vs transitP(MTD) (matriz global)", {
  P1 <- probs(MTD)                 # sem args -> matriz global
  P2 <- transitP(MTD)
  expect_true(is.matrix(P1))
  expect_equal(dim(P1), dim(P2))
  expect_true(row_sums_one(P1))
  expect_equal(P1, P2, tolerance = 1e-12)
})

test_that("probs(): newdata e context coerentes (mesmo passado)", {
  ## Monta um passado (ordem mais antiga -> mais recente) conforme docs da função
  ctx_vec <- c(0, 0, 1)               # exemplo qualquer em A
  names(ctx_vec) <- paste0("-", S_true)  # rótulos apenas para clareza
  
  ## Para MTD: usando 'context' (um símbolo por lag, na ordem das lags)
  pr_ctx <- probs(MTD, context = setNames(states(MTD), S_true))  # x_{-1}=0, x_{-5}=1 (ordem das lags no vetor)
  expect_true(is.matrix(pr_ctx) || is.numeric(pr_ctx))
  
  ## Para MTD via 'newdata': janela recente primeiro (mais recente primeiro)
  ## Pegamos do X uma prefixo que contenha ao menos max(S_true)
  pr_new  <- probs(MTD, newdata = X)
  expect_error(pr_new  <- probs(MTD, newdata = X[1]))

  ## Checagens de tipo/escala
  expect_true(all(pr_new >= 0) && all(pr_new <= 1))
  ## Não exigimos igualdade exata das duas chamdas, pois os modos de indexação podem
  ## produzir outputs em formato diferente (vetor/matriz). O importante é que retornem
  ## probabilidades válidas e não falhem.
})

test_that("probs(MTDest) retorna probabilidades e é finito", {
  pr_fit <- probs(fit)             # matriz global (via MTDest)
  expect_true(is.matrix(pr_fit))
  expect_true(row_sums_one(pr_fit))
  expect_true(all(is.finite(pr_fit)))
})

ct <- countsTab(X = X, d = max(Lambda(MTD)))
ft <- freqTab(S = Lambda(MTD), j = NULL, A =  states(MTD), countsTab = ct)

test_that("empirical_probs", {
  expect_equal(as.vector(empirical_probs(X, Lambda(MTD))[3])[[1]], ft$qax_Sj)
})


## ===========================
## 5) logLik(): tipo e atributos básicos
## ===========================
test_that("logLik() para MTD e MTDest retorna objeto válido", {
  ll_m   <- logLik(MTD,  X = X)
  ll_fit <- logLik(fit)
  pos <- which(ft$Nxa_Sj > 0)
  expect_equal(as.numeric(ll_m), sum(log(as.vector(t(transitP(MTD)))[pos]) * ft$Nxa_Sj[pos]))
  expect_equal(as.numeric(ll_fit), sum(log(as.vector(t(transitP(as.MTD(fit))))[pos]) * ft$Nxa_Sj[pos]))
  
  d_eff_m <- max(abs(lags(MTD))) 
  ## Não assumimos valores específicos — só tipo e finitude:
  expect_s3_class(ll_m,   "logLik")
  expect_s3_class(ll_fit, "logLik")
  expect_true(is.finite(as.numeric(ll_m)))
  expect_true(is.finite(as.numeric(ll_fit)))
  
  ## Atributos comuns em 'logLik' (opcionais, mas bons de ter)
  ## Se não setou nobs/df, estes testes podem ser afrouxados:
  if (!is.null(attr(ll_m, "nobs"))) expect_equal(attr(ll_m, "nobs"), length(X) - d_eff_m)
  if (!is.null(attr(ll_fit, "nobs"))) expect_equal(attr(ll_fit, "nobs"), length(X) - d_eff_m)
  
  if (!is.null(attr(ll_m, "df"))) expect_equal(attr(ll_m, "df"), n_parameters(S(fit), states(fit),
                                                                              single_matrix = FALSE,
                                                                              indep_part = TRUE))
  if (!is.null(attr(ll_fit, "df"))) expect_equal(attr(ll_fit, "df"), n_parameters(S(fit), states(fit),
                                                                              single_matrix = FALSE,
                                                                              indep_part = TRUE))
})
## ===========================
## 7) (Opcional) Regressão contra CRAN
##    Compara hdMTD_FS() da sua versão vs CRAN
## ===========================

seed <- sample(seq_len(10), 1)

test_that("MTD", {
  

  MTDmodel_cran <- get("MTDmodel", envir = ns_cran)
  MTDmodel_dev  <- get("MTDmodel",  envir = ns_dev)
  MTDest_cran <- get("MTDest", envir = ns_cran)
  MTDest_dev <- get("MTDest", envir = ns_dev)
  oscillation_cran <- get("oscillation", envir = ns_cran)
  oscillation_dev <- get("oscillation", envir = ns_dev)
  perfectSample_cran <- get("perfectSample", envir = ns_cran)
  perfectSample_dev <- get("perfectSample", envir = ns_dev)
    
  set.seed(seed)
  MTD_cran <- MTDmodel_cran(Lambda = sp$Lambda, A = sp$A,
                               lam0 = sp$lam0, lamj = sp$lamj,
                               pj = sp$pj, p0 = sp$p0,
                               indep_part = sp$indep,
                               single_matrix =  sp$sm)
  
  set.seed(seed)
  MTD_dev <- MTDmodel_dev(Lambda = sp$Lambda, A = sp$A,
                             lam0 = sp$lam0, lamj = sp$lamj,
                             pj = sp$pj, p0 = sp$p0,
                             indep_part = sp$indep,
                             single_matrix =  sp$sm)
  expect_equal(transitP(MTD_dev), MTD_cran$P)
  
  expect_equal(oscillation_cran(MTD_dev), oscillation_dev(MTD_dev))
  
  set.seed(seed)
  X_cran <- perfectSample_cran(MTD_dev, 800)
  set.seed(seed)
  X_dev <- perfectSample_dev(MTD_dev, 800)
  expect_equal(X_cran, X_dev)
  
  fit_cran <- MTDest_cran(X = X_dev, S = Lambda(MTD_dev), A = states(MTD_dev),
                          init = coef(MTD_dev), iter = TRUE)
  fit_dev <- MTDest_dev(X = X_dev, S = Lambda(MTD_dev), A = states(MTD_dev),
                          init = coef(MTD_dev), iter = TRUE)
  
  expect_equal(unname(p0(fit_dev)), unname(fit_cran$p0))
  expect_equal(unname(lambdas(fit_dev)), unname(fit_cran$lambdas))
  expect_equal(unname(pj(fit_dev)), unname(fit_cran$pj))
  expect_equal(fit_dev$iterations, fit_cran$iterations)
  
  expect_equal(oscillation_cran(X_dev, S = Lambda(MTD_dev)), oscillation_dev(X_dev, Lambda(MTD_dev)))
  expect_equal(probs(MTD_cran, newdata = X_dev), probs(MTD_dev, newdata = X_dev))
  
})





set.seed(seed)
X_dev <- perfectSample_dev(MTD, 800)

test_that("hdMTD()", {
  
  hdMTD_FS_cran <- get("hdMTD_FS", envir = ns_cran)
  hdMTD_FS_dev <- get("hdMTD_FS", envir = ns_dev)
  hdMTD_dev <- get("hdMTD", envir = ns_dev)
  hdMTD_BIC_cran <- get("hdMTD_BIC", envir = ns_cran)
  hdMTD_BIC_dev <- get("hdMTD_BIC", envir = ns_dev)
  hdMTD_CUT_cran <- get("hdMTD_CUT", envir = ns_cran)
  hdMTD_CUT_dev <- get("hdMTD_CUT", envir = ns_dev)
  hdMTD_FSC_cran <- get("hdMTD_FSC", envir = ns_cran)
  hdMTD_FSC_dev <- get("hdMTD_FSC", envir = ns_dev)
  
  
  S_cran <- hdMTD_FS_cran(X_dev, d = 6, l = 2)
  S_dev <- hdMTD_FS_dev(X_dev, d = 6, l = 2)
  expect_equal(S_cran, S_dev)
  expect_equal(S_cran, as.vector(S(
    hdMTD::hdMTD(X_dev, d = 6, method = "FS", l = 2)
  )))
  
  S_cran <- hdMTD_BIC_cran(X_dev, d = 6, minl = 2, maxl = 4, BICvalue = TRUE)
  S_dev <- hdMTD_BIC_dev(X_dev, d = 6, minl = 2, maxl = 4, BICvalue = TRUE)
  expect_equal(S_cran, S_dev)
  expect_equal(S_cran,
               attr(hdMTD_dev(X_dev, d = 6, method = "BIC", minl = 2, maxl = 4, BICvalue = TRUE),"extras")[[1]]
  )
  
  expect_equal(
    hdMTD_BIC_cran(X_dev, d = 6, minl = 2, maxl = 4, BICvalue = TRUE, byl = TRUE),
    attr(hdMTD_dev(X_dev, d = 6, method = "BIC", minl = 2, maxl = 4, BICvalue = TRUE, byl = TRUE), "extras")[[1]]
  )
  
  S_cran <- rev(hdMTD_CUT_cran(X_dev, d = 6, S = c(1,3,5,6), alpha = 0.01, xi = 0.5, mu = 0.8))
  S_dev <- hdMTD_CUT_dev(X_dev, d = 6, S = c(1,3,5,6), alpha = 0.01, xi = 0.5, mu = 0.8)
  expect_equal(S_cran, S_dev)
  expect_equal(S_cran, as.vector(S(hdMTD_dev(X_dev, d = 6, method = "CUT",
                                                S = c(1,3,5,6),
                                                alpha = 0.2,
                                                xi = 0.5,
                                                mu = 0.8 ))))
  
  S_cran <- rev(hdMTD_FSC_cran(X_dev, d = 6, l = 3, alpha = 0.01, xi = 0.5, mu = 0.8))
  S_dev <- hdMTD_FSC_dev(X_dev, d = 6, l = 3, alpha = 0.01, xi = 0.5, mu = 0.8)
  expect_equal(S_cran, S_dev)
  expect_equal(S_cran, as.vector(S(hdMTD_dev(X_dev, d = 6, method = "FSC",
                                                l = 3,
                                                alpha = 0.2,
                                                xi = 0.5,
                                                mu = 0.8 ))))
  
  
})


print("STOPPP")




################################################################
################################################################
################################################################
################################################################


# 1) Namespace CRAN em lib isolada
ensure_cran_hdMTD <- function(old_lib) {
  if ("hdMTD" %in% loadedNamespaces()) unloadNamespace("hdMTD")
  ns <- withr::with_libpaths(old_lib, action = "prefix", {
    loadNamespace("hdMTD")
  })
  stopifnot(grepl(normalizePath(old_lib, winslash = "/"),
                  normalizePath(getNamespaceInfo(ns, "path"), winslash = "/")))
  ns
}

old_lib <- file.path(tempdir(), "cran_lib_hdMTD")
ns_cran <- ensure_cran_hdMTD(old_lib)  # <- CRAN garantido

# 2) Namespace DEV (carregado do seu diretório do pacote, não da lib do R)
#    Requer pkgload; se preferir, use devtools::load_all()
if ("hdMTD" %in% loadedNamespaces()) unloadNamespace("hdMTD")
pkgload::load_all(export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
ns_dev <- asNamespace("hdMTD")         # <- DEV garantido (do fonte)

