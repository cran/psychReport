# Internal sphericity functions
# These compute Mauchly's test for sphericity and GG/HF epsilon corrections
# using base R linear algebra (no external dependencies required).

# Build a (k x k-1) contrast matrix.
# Rows 1..k-1 have a 1 on the diagonal, row k has -1 in every column.
.contrastMatrix <- function(k) {
  P <- matrix(0, nrow = k, ncol = k - 1)
  for (i in 1:(k - 1)) {
    P[i, i] <- 1
    P[k, i] <- -1
  }
  return(P)
}

# Build a subject x condition matrix of cell means from long-format data.
# Returns an N x k numeric matrix (N subjects, k conditions).
.subjectConditionMatrix <- function(data, dv, id, factors) {

  # Build condition grid (all factor-level combinations)
  # Use rev() so first factor varies slowest, matching kronecker() convention
  level_list <- lapply(factors, function(f) sort(unique(data[[f]])))
  cond_grid  <- expand.grid(rev(level_list))
  cond_grid  <- cond_grid[, rev(seq_len(ncol(cond_grid))), drop = FALSE]
  names(cond_grid) <- factors

  # Create a single interaction key for each row
  if (length(factors) == 1) {
    data_key <- as.character(data[[factors]])
    cond_key <- as.character(cond_grid[[factors]])
  } else {
    data_key <- do.call(paste, c(data[factors], sep = ":"))
    cond_key <- do.call(paste, c(cond_grid[factors], sep = ":"))
  }

  # Compute cell means using tapply (vectorized, no R loops)
  cell_means <- tapply(data[[dv]], list(data[[id]], data_key), mean)

  # Reorder columns to match the condition grid ordering
  # and rows to sorted subjects
  subjects <- sort(unique(data[[id]]))
  Y <- cell_means[as.character(subjects), cond_key, drop = FALSE]

  return(unname(Y))
}

# Compute the SSPE (sum of squares and products of errors) matrix.
.computeSSPE <- function(Y) {
  Y_centered <- sweep(Y, 2, colMeans(Y))
  crossprod(Y_centered)
}


# Build contrast matrix for a (possibly interaction) effect.
# For single factors, returns .contrastMatrix(k).
# For interactions, returns the Kronecker product of per-factor contrast matrices.
.buildContrastMatrixForEffect <- function(factors, data) {
  n_per_factor <- sapply(factors, function(f) length(unique(data[[f]])))

  if (length(factors) == 1) {
    return(.contrastMatrix(n_per_factor))
  }

  P <- .contrastMatrix(n_per_factor[1])
  for (i in 2:length(n_per_factor)) {
    P <- kronecker(P, .contrastMatrix(n_per_factor[i]))
  }
  return(P)
}


# Compute all sphericity components for a given set of factors.
# Returns list with u_matrix, n_contrasts, df, n_conditions, W, and p (Mauchly).
# Returns NULL if < 3 cells (sphericity is trivially satisfied for 2 levels).
.computeEffectSphericity <- function(data, dv, id, factors) {
  n_cells <- prod(sapply(factors, function(f) length(unique(data[[f]]))))
  if (n_cells < 3) return(NULL)

  Y    <- .subjectConditionMatrix(data, dv, id, factors)
  SSPE <- .computeSSPE(Y)
  P    <- .buildContrastMatrixForEffect(factors, data)

  transformed_sspe <- crossprod(P, SSPE %*% P)
  p_transpose_p    <- crossprod(P)
  u_matrix         <- solve(p_transpose_p, transformed_sspe)

  n_contrasts  <- ncol(transformed_sspe)
  n_conditions <- nrow(P)
  df           <- nrow(Y) - 1

  # Compute Mauchly's W and p directly (avoid redundant recomputation)
  W <- 1.0
  p <- 1.0
  if (n_contrasts >= 2) {
    trace_u <- sum(diag(u_matrix))
    logW    <- log(det(u_matrix)) - n_contrasts * log(trace_u / n_contrasts)
    W       <- exp(logW)

    # Bartlett correction (from R's stats:::mauchly.test.SSD)
    rho <- 1.0 - (2 * n_contrasts^2 + n_contrasts + 2) / (6 * n_contrasts * df)

    # Second-order correction factor
    w2 <- (n_contrasts + 2) * (n_contrasts - 1) * (n_contrasts - 2) *
      (2 * n_contrasts^3 + 6 * n_contrasts^2 + 3 * n_conditions + 2) /
      (288 * (df * n_contrasts * rho)^2)

    # Test statistic
    z <- -df * rho * logW
    f <- n_contrasts * (n_contrasts + 1) / 2 - 1
    if (f > 0) {
      Pr1 <- 1.0 - stats::pchisq(z, df = f)
      Pr2 <- 1.0 - stats::pchisq(z, df = f + 4)
      p   <- max(0.0, min(1.0, Pr1 + w2 * (Pr2 - Pr1)))
    }
  }

  # Compute GG epsilon
  trace_u       <- sum(diag(u_matrix))
  sum_squared   <- sum(u_matrix * t(u_matrix))
  gg_eps        <- trace_u^2 / (n_contrasts * sum_squared)

  # Compute HF epsilon
  n_subjects    <- df + 1
  hf_numerator  <- n_subjects * n_contrasts * gg_eps - 2
  hf_denom      <- n_contrasts * (n_subjects - 1 - n_contrasts * gg_eps)
  hf_eps        <- min(1.0, hf_numerator / hf_denom)

  return(list(
    W = W, p = p,
    gg_eps = gg_eps, hf_eps = hf_eps,
    n_contrasts = n_contrasts
  ))
}


# Top-level function: compute sphericity data for all effects in an aov object.
#
# Takes the raw data and the ANOVA table (from aovTidyTable), and returns
# a list with two components:
#   $"Mauchly's Test for Sphericity" — data.frame with Effect, W, p, p<.05
#   $"Sphericity Corrections"        — data.frame with Effect, GGe, p[GG], HFe, p[HF]
#
# Only effects with DFn > 1 are included (DFn == 1 means 2 levels = no sphericity issue).
.computeSphericity <- function(aovTable, data) {

  effects <- aovTable$ANOVA

  # find effects with DFn > 1
  sph_rows <- which(effects$DFn > 1)
  if (length(sph_rows) == 0) return(NULL)

  # Identify columns
  effect_names <- effects$Effect
  all_factors  <- unique(unlist(strsplit(effect_names, ":")))

  # psychReport: VP is the subject ID
  numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
  factor_cols  <- names(data)[vapply(data, function(x) is.factor(x) | is.character(x), logical(1))]

  # DV is the numeric column that is NOT a factor in the ANOVA
  dv_candidates <- setdiff(numeric_cols, all_factors)
  if (length(dv_candidates) == 0) {
    dv_candidates <- setdiff(numeric_cols, c(all_factors, factor_cols))
  }

  # ID is the factor column not in the ANOVA effects
  id_candidates <- setdiff(factor_cols, all_factors)

  if (length(dv_candidates) == 0 || length(id_candidates) == 0) return(NULL)

  dv_name <- dv_candidates[1]
  id_name <- id_candidates[1]

  # Pre-allocate result vectors
  n_effects <- length(sph_rows)
  m_effect <- character(n_effects)
  m_W      <- numeric(n_effects)
  m_p      <- numeric(n_effects)
  c_effect <- character(n_effects)
  c_GGe    <- numeric(n_effects)
  c_pGG    <- numeric(n_effects)
  c_HFe    <- numeric(n_effects)
  c_pHF    <- numeric(n_effects)

  for (k in seq_along(sph_rows)) {
    row_idx <- sph_rows[k]
    effect  <- effects$Effect[row_idx]
    factors <- strsplit(effect, ":")[[1]]

    # Single call computes everything: Mauchly + GG + HF
    comp <- .computeEffectSphericity(data, dv_name, id_name, factors)

    if (is.null(comp)) {
      m_W[k] <- 1.0; m_p[k] <- 1.0
      c_GGe[k] <- 1.0; c_HFe[k] <- 1.0
    } else {
      m_W[k] <- comp$W; m_p[k] <- comp$p
      c_GGe[k] <- comp$gg_eps; c_HFe[k] <- comp$hf_eps
    }

    m_effect[k] <- effect
    c_effect[k] <- effect

    # Corrected p-values
    DFn  <- effects$DFn[row_idx]
    DFd  <- effects$DFd[row_idx]
    Fval <- effects$F[row_idx]
    c_pGG[k] <- 1.0 - stats::pf(Fval, DFn * c_GGe[k], DFd * c_GGe[k])
    c_pHF[k] <- 1.0 - stats::pf(Fval, DFn * c_HFe[k], DFd * c_HFe[k])
  }

  mauchly_df <- data.frame(
    Effect  = m_effect,
    W       = m_W,
    p       = m_p,
    "p<.05" = ifelse(m_p < 0.05, "*", ""),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(mauchly_df) <- m_effect

  corrections_df <- data.frame(
    Effect  = c_effect,
    GGe     = c_GGe,
    "p[GG]" = c_pGG,
    HFe     = c_HFe,
    "p[HF]" = c_pHF,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(corrections_df) <- c_effect

  return(list(
    "Mauchly's Test for Sphericity" = mauchly_df,
    "Sphericity Corrections"        = corrections_df
  ))
}
