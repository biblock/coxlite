# R version of src/cox_derivs.cpp (engine = "R"). Same inputs/outputs.
# Rows sorted by stratum then decreasing time, x centred.
cox_derivs_r <- function(time, status, x, strata, beta, efron) {
  n <- nrow(x)
  p <- ncol(x)
  eta <- drop(x %*% beta)
  risk <- exp(eta)

  # block = rows with the same stratum and time. time is decreasing so the
  # risk set for a block is everything up to its last row
  new_block <- c(TRUE, time[-1L] != time[-n] | strata[-1L] != strata[-n])
  block <- cumsum(new_block)
  block_last <- c(which(new_block)[-1L] - 1L, n)

  # risk set sums S0, S1, S2
  rx <- risk * x
  rxx <- row_outer(rx, x)
  sums <- group_cumsum(cbind(risk, rx, rxx), strata)[block_last, , drop = FALSE]
  s0 <- sums[, 1L]
  s1 <- sums[, 1L + seq_len(p), drop = FALSE]
  s2 <- sums[, 1L + p + seq_len(p * p), drop = FALSE]

  # same sums over only the events (D0, D1, D2)
  dead <- status == 1L
  dead_block <- block[dead]
  dsum <- rowsum(cbind(1, risk, rx, rxx)[dead, , drop = FALSE], dead_block,
                 reorder = FALSE)
  ublock <- unique(dead_block)
  nd <- dsum[, 1L]
  d0 <- dsum[, 2L]
  d1 <- dsum[, 2L + seq_len(p), drop = FALSE]
  d2 <- dsum[, 2L + p + seq_len(p * p), drop = FALSE]

  # one row per event. efron weight is l/d for the l-th of d ties, breslow is 0
  idx <- rep(seq_along(ublock), nd)
  frac <- if (efron) (sequence(nd) - 1) / nd[idx] else 0
  bi <- ublock[idx]
  den <- s0[bi] - frac * d0[idx]
  a1 <- (s1[bi, , drop = FALSE] - frac * d1[idx, , drop = FALSE]) / den
  a2 <- (s2[bi, , drop = FALSE] - frac * d2[idx, , drop = FALSE]) / den

  loglik <- sum(eta[dead]) - sum(log(den))
  score <- colSums(x[dead, , drop = FALSE]) - colSums(a1)
  info <- colSums(a2) - colSums(row_outer(a1, a1))

  return(list(loglik = loglik,
              score = unname(score),
              info = matrix(info, p, p)))
}

# row i is as.vector(a[i, ] %o% b[i, ])
row_outer <- function(a, b) {
  p <- ncol(a)
  return(a[, rep(seq_len(p), times = p), drop = FALSE] *
           b[, rep(seq_len(p), each = p), drop = FALSE])
}

# cumsum down each column, restarting at each group (groups must be contiguous)
group_cumsum <- function(m, g) {
  n <- nrow(m)
  cs <- matrix(apply(m, 2L, cumsum), nrow = n)
  starts <- c(TRUE, g[-1L] != g[-n])
  first <- cummax(ifelse(starts, seq_len(n), 0L))
  return(cs - rbind(0, cs)[first, , drop = FALSE])
}
