# tests/codezin.R
suppressPackageStartupMessages({
  library(hdMTD)
  library(dplyr)
})

# Função auxiliar para arredondar numerics dentro de listas/matrizes/data.frames
.round_all <- function(x, digits = 8) {
  if (is.numeric(x)) return(round(x, digits))
  if (is.matrix(x))  { storage.mode(x) <- "double"; return(round(x, digits)) }
  if (is.data.frame(x)) {
    y <- x
    for (i in seq_along(y)) if (is.numeric(y[[i]])) y[[i]] <- round(y[[i]], digits)
    return(y)
  }
  if (is.list(x)) return(lapply(x, .round_all, digits = digits))
  x
}

# Coleta resultados "rápidos" e determinísticos
codezin_collect <- function() {
  # --- Seção 5.1: gerar modelo e amostra (reduzida) ---
  set.seed(11)
  Lambda <- c(1, 15, 30)
  A <- c(0, 1)
  lam0 <- 0.01
  lamj <- c(0.39, 0.30, 0.30)
  p0 <- c(0.5, 0.5)
  MTD <- MTDmodel(Lambda = Lambda, A = A, lam0 = lam0, lamj = lamj, p0 = p0)

  set.seed(11)
  X <- perfectSample(MTD, N = 800)  # 800 em vez de 1000 p/ ganhar tempo

  # --- Seção 5.2: inferência (subset reduzido para ser rápido) ---
  set.seed(11)                      # FS desempata com sample(), fixe a seed!
  fs_lags <- hdMTD_FS(X, d = 30, l = 3)

  bic_fix <- hdMTD_BIC(
    X, d = 30,
    S = c(1, 5, 10, 15, 20, 25, 30),
    minl = 3, maxl = 3
  )
  bic_range <- hdMTD_BIC(
    X, d = 30,
    S = c(1, 5, 10, 15, 20, 25, 30),
    minl = 1, maxl = 3, byl = TRUE, BICvalue = TRUE
  )

  cut_set <- hdMTD_CUT(X, d = 30, S = c(1, 5, 15, 30), alpha = 0.1)

  set.seed(11)
  fsc_set <- hdMTD_FSC(X, d = 30, l = 3, alpha = 0.1)

  # --- empirical_probs e Oscillation (rápidos) ---
  P_df  <- empirical_probs(X, S = c(1, 15, 30), matrixform = FALSE)
  P_mat <- empirical_probs(X, S = c(1, 15, 30), matrixform = TRUE)

  osc_model  <- oscillation(MTD)
  osc_sample <- oscillation(X, S = c(1, 15, 30))

  # --- EM (curto) ---
  init <- list(
    lambdas = c(0.01, 0.33, 0.33, 0.33),
    p0 = c(0.5, 0.5),
    pj = rep(list(matrix(c(0.5, 0.5, 0.5, 0.5), ncol = 2, nrow = 2)), 3)
  )
  set.seed(11)
  fit_em <- MTDest(X, S = c(1, 15, 30), init = init, iter = TRUE, nIter = 20)

  # --- Monta lista, padroniza para comparações estáveis ---
  out <- list(
    fs_lags   = sort(fs_lags),
    bic_fix   = sort(as.numeric(bic_fix)),      # BIC=min(lags); se por nomes, parseie
    bic_range = bic_range,                      # aqui pode vir nomeado; manter “cru” para sinal de mudança
    cut_set   = sort(cut_set),
    fsc_set   = sort(fsc_set),
    P_df_head = head(P_df, 12),                 # só um “slice” para ser leve
    P_mat_blk = P_mat[1:min(8, nrow(P_mat)), , drop = FALSE],
    osc_model  = osc_model,
    osc_sample = osc_sample,
    em_lambdas = fit_em$lambdas,
    em_iters   = fit_em$iterations,
    em_logLik  = fit_em$loglik
  )

  .round_all(out, digits = 8)
}

# Salva baseline (rode 1x antes das mudanças)
codezin_save_baseline <- function(path = "tests/codezin_baseline.rds") {
  b <- codezin_collect()
  saveRDS(b, path)
  message("Baseline salva em: ", path)
  invisible(b)
}

# Compara com baseline
codezin_check <- function(path = "tests/dev/codezin_baseline.rds", tol = 1e-8) {
  if (!file.exists(path)) stop("Baseline não encontrada em ", path, ". Rode codezin_save_baseline() primeiro.")
  base <- readRDS(path)
  cur  <- codezin_collect()

  # compara chave a chave
  keys <- union(names(base), names(cur))
  diffs <- list()
  for (k in keys) {
    a <- base[[k]]; b <- cur[[k]]
    ok <- isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))
    if (!ok) diffs[[k]] <- list(base = a, current = b, diff = all.equal(a, b, tolerance = tol, check.attributes = FALSE))
  }

  if (length(diffs) == 0) {
    message("✅ codezin_check: tudo estável vs baseline.")
    return(invisible(TRUE))
  } else {
    message("⚠️ codezin_check: diferenças detectadas nas chaves: ", paste(names(diffs), collapse = ", "))
    return(invisible(diffs))
  }
}

# Uso típico:
# devtools::load_all()
# codezin_save_baseline()  # rode 1x antes de começar alterações
# codezin_check()          # rode após cada alteração
