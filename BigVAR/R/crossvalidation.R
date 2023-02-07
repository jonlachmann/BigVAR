crossValidation <- function (
  Y, ZFull, h, recursive, window.size, cvtype, T1, T2, group, beta, lambda, tol, p, m, k1, k, s, s1,
  MN, C, intercept, separate_lambdas, dual, activeset, starting_eigvals, gran2,
  groups, compgroups, VARXI, alpha, palpha, gamma, RVAR, loss, delta, verbose
) {
  msfe_dim <- c(T2 - T1 + 1, gran2)
  if (!dual) {
    if (separate_lambdas) {
        if (!VARXI) {
          msfe_dim <- c(msfe_dim, k)
        } else {
          msfe_dim <- c(msfe_dim, k1)
        }
    } else {
      lambda <- as.matrix(lambda)
    }
  } else {
    msfe_dim[2] <- msfe_dim[2] * length(alpha)
  }
  MSFE <- array(0, dim = msfe_dim)

  if (verbose) {
    pb <- txtProgressBar(min = T1, max = T2 - h + 1, style = 3)
    cat("Validation Stage:", group, "\n")
  }
  YT <- Y[1:T2, , drop = FALSE]

  for (v in T1:T2 - h + 1) {
    cvsample <- getCrossValidationSample(v, YT, ZFull, VARXI, cvtype, h, k, k1, p, m, s, recursive, window.size)

    temp <- .BigVAR.fit(group, beta, cvsample$trainZ, cvsample$trainY, lambda, tol, p, m, k1, k, s, s1, MN, C, intercept, separate_lambdas,
                        dual, activeset, starting_eigvals, groups, compgroups, VARXI, alpha, palpha, gamma)
    beta <- temp$beta

    activeset <- temp$activeset
    eZ <- c(1, ZFull$Z[, v])

    temp_results <- refine_and_forecast(beta, as.matrix(eZ), cvsample$trainZ, cvsample$trainY, ZFull$Y[v + h - 1, , drop = F], lambda = lambda,
                                        h = h, recursive = recursive, MN = MN, RVAR = RVAR, refit_fraction = refit_fraction, separate_lambdas = separate_lambdas,
                                        inds = NULL, loss = loss, delta = delta, k = k, p = p, k1 = k1, s = s, oos = FALSE)
    if (separate_lambdas) {
      MSFE[v - (T1 - h), , ] <- temp_results$MSFE
    } else {
      MSFE[v - (T1 - h), ] <- temp_results$MSFE
    }
    if (verbose) {
      setTxtProgressBar(pb, v)
    }
  }

  return(MSFE)
}

applyMinnesotaPrior <- function (betaWS, C) {
  for (i in 1:dim(betaWS)[3]) {
    submat <- adrop(betaWS[1:dim(betaWS)[1], 1:dim(betaWS)[1], i, drop = F], drop = 3)
    diag(submat) <- diag(submat) - C
    betaWS[1:dim(betaWS)[1], 1:dim(betaWS)[1], i] <- submat
  }
  return(betaWS)
}

getCrossValidationSample <- function (v, YT, ZFull, VARXI, cvtype, h, k, k1, p, m, s, recursive, window.size) {
  if (cvtype == "Rolling") {
    if (h > 1 & !recursive) {
      if (window.size != 0) {
        ws1 <- max(c(v - window.size - h, 1))
        index_y <- (ws1 + h - 1):(v - 1)
        index_z <- (ws1 ):(v - h)
      } else {
        index_y <- (h):(v - 1)
        index_z <- 1:(v - h)
      }
    } else {
      if (window.size != 0) {
        ws1 <- max(c(v - window.size, 1))
        index_y <- (ws1):(v - 1)
        index_z <- index_y
      } else {
        index_y <- (1):(v - 1)
        index_z <- index_y
      }
    }
    trainY <- ZFull$Y[index_y, , drop = FALSE]
    trainZ <- ZFull$Z[, index_z, drop = F]
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
  return(list(trainY = trainY, trainZ = trainZ))
}