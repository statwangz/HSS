irwls <- function(y, L2, update_x, weights, intercept = 1,
                  M, N, N_bar, fix_intercept = T, seperator) {

  # calculate new weights
  weights <- weights / sum(weights)
  x <- L2 * N / N_bar

  if (fix_intercept) {
    y <- y - intercept
    for (i in 1:2) {
      wy <- y * weights
      wx <- as.matrix(x * weights)

      fit <- lm(wy ~ wx + 0)

      h2 <- coef(fit)[1] * M / N_bar
      h2 <- max(h2, 0)
      h2 <- min(h2, 1)

      # update weights
      weights <- get_weights(h2, intercept, L2, N, M, update_x)
      weights <- weights / sum(weights)
    }

    # preweight LD and chi:
    wLD <- as.matrix(x * weights)
    wCHI <- as.matrix(y * weights)
  } else {
    for (i in 1:2) {
      wy <- y * weights
      wx <- as.matrix(cbind(x * weights, weights))

      fit <- lm(wy ~ wx + 0)

      h2 <- coef(fit)[1] * M / N_bar
      h2 <- max(h2, 0)
      h2 <- min(h2, 1)
      intercept <- coef(fit)[2]

      # update weights
      weights <- get_weights(h2, intercept, L2, N, M, update_x)
      weights <- weights / sum(weights)
    }

    # preweight LD and chi:
    wLD <- as.matrix(cbind(x * weights, weights))
    wCHI <- as.matrix(y * weights)
  }


  # Perfrom analysis
  n_blocks <- length(seperator) - 1
  n_snps <- length(y)
  p <- ncol(wLD)
  xty_block_values <- matrix(data = NA, nrow = n_blocks, ncol = p)
  xtx_block_values <- matrix(data = NA, nrow = p * n_blocks, ncol = p)

  from <- seperator
  to <- c(seperator[2:n_blocks] - 1, n_snps)

  rep_from <- seq(from = 1, to = p * n_blocks, by = p)
  rep_to <- seq(from = p, to = p * n_blocks, by = p)

  colnames(xty_block_values) <- colnames(xtx_block_values) <- colnames(wLD)

  for (i in 1:n_blocks) {
    xty_block_values[i, ] <- t(t(wLD[from[i]:to[i], ]) %*% wCHI[from[i]:to[i], ])
    xtx_block_values[rep_from[i]:rep_to[i], ] <- as.matrix(t(wLD[from[i]:to[i], ]) %*%
      wLD[from[i]:to[i], ])
  }

  xty <- as.matrix(colSums(xty_block_values))
  xtx <- matrix(data = NA, nrow = p, ncol = p)
  colnames(xtx) <- colnames(wLD)

  for (i in 1:p) {
    xtx[i, ] <- t(colSums(as.matrix(xtx_block_values[seq(from = i, to = p * n_blocks, by = p), ])))
  }

  reg <- solve(xtx) %*% xty

  if (p == 1) {
    intercept_est <- intercept
    h2_est <- reg[1] / N_bar * M
  }

  if (p == 2) {
    intercept_est <- reg[2]
    h2_est <- reg[1] / N_bar * M
  }

  return(list(h2 = h2_est, intercept = intercept_est))
}
