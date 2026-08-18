# ============================================================
# SURVEY-WEIGHTED BASE-MEASURE FPS SENSITIVITY ANALYSIS
# ============================================================

fps_weights_base <- function(A,
                             C,
                             base_weights,
                             ridge = 0,
                             maxit = 5000,
                             reltol = 1e-10) {
  A <- as.matrix(A)
  C <- as.matrix(C)
  base_weights <- as.numeric(base_weights)
  
  if (
    nrow(A) != nrow(C) ||
    length(base_weights) != nrow(A)
  ) {
    stop(
      "A, C and base_weights have incompatible dimensions."
    )
  }
  
  if (
    any(!is.finite(A)) ||
    any(!is.finite(C)) ||
    any(!is.finite(base_weights)) ||
    any(base_weights < 0) ||
    sum(base_weights) <= 0
  ) {
    stop(
      "Invalid A, C or base_weights."
    )
  }
  
  q <- base_weights / sum(base_weights)
  
  weighted_column_mean <- function(x) {
    sum(q * x)
  }
  
  weighted_column_sd <- function(x) {
    mu <- weighted_column_mean(x)
    sqrt(sum(q * (x - mu)^2))
  }
  
  standardize_q <- function(M) {
    M <- as.matrix(M)
    
    center <- apply(
      M,
      2,
      weighted_column_mean
    )
    
    scale_value <- apply(
      M,
      2,
      weighted_column_sd
    )
    
    keep <-
      is.finite(scale_value) &
      scale_value > 1e-10
    
    M <- M[, keep, drop = FALSE]
    
    M <- sweep(
      M,
      2,
      center[keep],
      "-"
    )
    
    M <- sweep(
      M,
      2,
      scale_value[keep],
      "/"
    )
    
    M
  }
  
  Az <- standardize_q(A)
  Cz <- standardize_q(C)
  
  colnames(Az) <- paste0(
    "A_",
    seq_len(ncol(Az))
  )
  
  colnames(Cz) <- paste0(
    "C_",
    seq_len(ncol(Cz))
  )
  
  AC <- do.call(
    cbind,
    lapply(
      seq_len(ncol(Cz)),
      function(j) {
        result <- sweep(
          Az,
          1,
          Cz[, j],
          "*"
        )
        
        colnames(result) <- paste0(
          colnames(Az),
          ":",
          colnames(Cz)[j]
        )
        
        result
      }
    )
  )
  
  G <- cbind(Az, Cz, AC)
  
  weighted_sd_G <- apply(
    G,
    2,
    weighted_column_sd
  )
  
  keep <-
    is.finite(weighted_sd_G) &
    weighted_sd_G > 1e-10
  
  G <- G[, keep, drop = FALSE]
  weighted_sd_G <- weighted_sd_G[keep]
  
  Gs <- sweep(
    G,
    2,
    weighted_sd_G,
    "/"
  )
  
  objective <- function(theta) {
    eta <- -drop(Gs %*% theta)
    eta_max <- max(eta)
    
    eta_max +
      log(
        sum(q * exp(eta - eta_max))
      ) +
      0.5 * ridge * sum(theta^2)
  }
  
  gradient <- function(theta) {
    eta <- -drop(Gs %*% theta)
    
    w <- q * exp(
      eta - max(eta)
    )
    
    w <- w / sum(w)
    
    -drop(crossprod(Gs, w)) +
      ridge * theta
  }
  
  optimization <- stats::optim(
    par = rep(0, ncol(Gs)),
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(
      maxit = maxit,
      reltol = reltol
    )
  )
  
  eta <- -drop(
    Gs %*% optimization$par
  )
  
  weights <- q * exp(
    eta - max(eta)
  )
  
  weights <- weights / sum(weights)
  
  weighted_moments <- drop(
    crossprod(G, weights)
  )
  
  result <- list(
    weights = weights,
    theta = optimization$par,
    convergence = optimization$convergence,
    message = optimization$message,
    objective = optimization$value,
    weighted_moments = weighted_moments,
    maximum_imbalance = max(
      abs(weighted_moments)
    ),
    effective_sample_size =
      1 / sum(weights^2),
    optimization = optimization,
    base_weights = q
  )
  
  class(result) <- "fps_weights"
  result
}
