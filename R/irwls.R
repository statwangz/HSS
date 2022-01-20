irwls <- function(y, L2, update_x, weights, intercept = 1,
                  M, N, N_bar, fix_intercept = T, seperator, jackknife = T) {

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


  # Perform analysis
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

  if (jackknife) {
    # perform jackknife
    delete_from <- seq(from = 1, to = p * n_blocks, by = p)
    delete_to <- seq(from = p, to = p * n_blocks, by = p)
    delete_values <- matrix(data = NA, nrow = n_blocks, ncol = p)
    colnames(delete_values) <- colnames(wLD)

    for (i in 1:n_blocks) {
      xty_delete <- xty - xty_block_values[i, ]
      xtx_delete <- xtx - xtx_block_values[delete_from[i]:delete_to[i], ]
      delete_values[i, ] <- solve(xtx_delete) %*% xty_delete
    }

    delete_values <- as.matrix(delete_values[, 1:p])

    pseudo_values <- matrix(data = NA, nrow = n_blocks, ncol = p)
    colnames(pseudo_values) <- colnames(wLD)

    for (i in 1:n_blocks) {
      pseudo_values[i, ] <- (n_blocks * reg) - ((n_blocks - 1) * delete_values[i, ])
    }

    jackknife_cov <- cov(pseudo_values) / n_blocks
    jackknife_se <- sqrt(diag(jackknife_cov))

    coef_cov <- jackknife_cov[1, 1] / (N_bar^2) * M^2
    h2_se <- sqrt(coef_cov)

    if (p == 1) {
      intercept_est <- intercept
      h2_est <- reg[1] / N_bar * M
      intercept_se <- NA
    }

    if (p == 2) {
      intercept_est <- reg[2]
      h2_est <- reg[1] / N_bar * M
      intercept_se <- jackknife_se[length(jackknife_se)]
    }

    return(list(
      h2 = h2_est, h2_se = h2_se,
      intercept = intercept_est,
      intercept_se = intercept_se,
      delete_values = delete_values,
      jk_est = reg[1]
    ))
  } else {
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
}
