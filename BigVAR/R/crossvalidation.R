crossValidation <- function (
  Y, ZFull, h, recursive, window.size, cvtype, T1, T2, group, beta, lambda, tol, p, m, k1, k, s, s1,
  MN, C, intercept, separate_lambdas, dual, activeset, starting_eigvals, gran2,
  groups, compgroups, VARXI, alpha, palpha, gamma, RVAR, loss, delta, verbose
) {
    if (!dual) {
        if (separate_lambdas) {
            if (!VARXI) {
                MSFE <- array(0, dim = c(T2 - T1 + 1, gran2, k))
            } else {
                MSFE <- array(0, dim = c(T2 - T1 + 1, gran2, k1))
            }
        } else {
            MSFE <- matrix(0, nrow = T2 - T1 + 1, ncol = gran2)
            lambda <- as.matrix(lambda)
        }
    } else {
        nalpha <- length(alpha)
        MSFE <- matrix(0, nrow = T2 - T1 + 1, ncol = gran2 * nalpha)
    }

  if (verbose) {
    pb <- txtProgressBar(min = T1, max = T2, style = 3)
    cat("Validation Stage:", group, "\n")
  }
  YT <- Y[1:T2, , drop = FALSE]
  betaWS <- beta
  for (v in (T1 - h + 1):T2) {
    if (v + h - 1 > T2) {
      break
    }
    if (cvtype == "Rolling") {
      if (h > 1 & !recursive) {
        if (window.size != 0) {
          ws1 <- max(c(v - window.size - h, 1))
          index_y <- (ws1 + h-1):(v - 1)
          index_z <- (ws1 ):(v - h)
          trainY <- ZFull$Y[index_y, , drop = FALSE]
          trainZ <- ZFull$Z[, index_z,drop=FALSE]
        } else {
          index_y <- (h):(v - 1)
          index_z <- 1:(v - h)
          trainY <- ZFull$Y[index_y, , drop = FALSE]
          trainZ <- ZFull$Z[, index_z, drop = F]
        }
      } else {
        if (window.size != 0) {
          ws1 <- max(c(v - window.size, 1))
          trainY <- ZFull$Y[(ws1):(v - 1), , drop = FALSE]
          trainZ <- ZFull$Z[, (ws1):(v - 1), drop = FALSE]
        } else {
          trainY <- ZFull$Y[(1):(v - 1), , drop = FALSE]
          trainZ <- ZFull$Z[, (1):(v - 1), drop = FALSE]
        }
      }
    } else {
      if (VARXI) {
        YT2 <- YT[-v, , drop = FALSE]
        Y1 <- YT2[, 1:k1, drop = FALSE]
        X <- YT2[, (ncol(YT2) - m + 1):ncol(YT2), drop = FALSE]
        trainZ <- VARXCons(Y1, X, k1, p, m, s, contemp = contemp)
        trainZ <- trainZ[2:nrow(trainZ), , drop = FALSE]
        trainY <- YT2[(max(c(p, s)) + 1):nrow(YT2), 1:k1, drop = FALSE]
      } else {
        YT2 <- YT[-v, , drop = FALSE]
        Z1 <- VARXCons(YT2, matrix(0, nrow = nrow(YT2)), k, p, 0, 0)
        trainZ <- Z1[2:nrow(Z1), ]
        trainY <- YT2[(p + 1):nrow(YT2), , drop = FALSE]
      }
    }
    temp <- .BigVAR.fit(group, betaWS, trainZ, trainY, lambda, tol, p, m, k1, k, s, s1, MN, C, intercept, separate_lambdas,
                        dual, activeset, starting_eigvals, groups, compgroups, VARXI, alpha, palpha, gamma)
    beta <- temp$beta
    betaWS <- temp$beta
    if (MN) {
      for (i in 1:dim(betaWS)[3]) {
        submat <- adrop(betaWS[1:dim(betaWS)[1], 1:dim(betaWS)[1], i, drop = F], drop = 3)
        diag(submat) <- diag(submat) - C
        betaWS[1:dim(betaWS)[1], 1:dim(betaWS)[1], i] <- submat
      }
    }
    activeset <- temp$activeset
    q1a <- temp$q1a
    eZ <- c(1, ZFull$Z[, v])
    msfe_index <- v - (T1 - h)

    temp_results <- refine_and_forecast(beta, as.matrix(eZ), trainZ, trainY, ZFull$Y[v + h - 1, , drop = F], lambda = lambda,
                                        h = h, recursive = recursive, MN = MN, RVAR = RVAR, refit_fraction = refit_fraction, separate_lambdas = separate_lambdas,
                                        inds = NULL, loss = loss, delta = delta, k = k, p = p, k1 = k1, s = s, oos = FALSE)
    if (separate_lambdas) {
      MSFE[msfe_index, , ] <- temp_results$MSFE
    } else {
      MSFE[msfe_index, ] <- temp_results$MSFE
    }
    if (verbose) {
      setTxtProgressBar(pb, v)
    }
  }

  return(MSFE)
}