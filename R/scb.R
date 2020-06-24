scb <-
function (object, drv = 0, level = 0.95, pred = 500, div = 1000, 
          pages = 0) 
{
  if (identical(as.numeric(object$sp), numeric(0))) 
    stop("Do not run scb without any smooth term")
  if (drv > 2) 
    stop("Only first and second derivatives are supported")
  if (class(object)[1] != "pam") 
    stop("Object must be of class pam")
  y_orig <- object$model[, 1]
  gls <- object$gls
  model <- object$model
  sigma2 <- object$sig2
  coef <- object$coefficients
  index_data <- object$index_data
  if (!is.null(object$full.sp)) 
    sp <- object$full.sp
  else sp <- object$sp
  smi <- object$smi
  df <- object$df.residual
  if ("(weights)" %in% colnames(model)) {
    weights <- object$weights
    model <- object$model[, -which(colnames(object$model) == 
                                     "(weights)")]
  }
  else weights = NULL
  if (!is.null(weights)) {
    if (!is.null(smi)) 
      smi <- .sparseDiagonal(n = length(weights), x = weights) %*% 
        smi
    else smi <- .sparseDiagonal(n = length(weights), x = weights)
  }
  num_factors <- sum((lapply(model, is.factor) == TRUE))
  if (num_factors > 0) {
    B <- matrix(0, nrow(model), 0)
    col_factors <- c()
    for (i in 2:(ncol(model))) if (is.factor(model[, i])) {
      M <- model.matrix(~-1 + model[, i])[, -1, drop = FALSE]
      colnames(M) <- paste(colnames(model)[i], levels(model[, 
                                                            i])[-1], sep = "")
      B <- cbind(B, M)
      col_factors <- c(col_factors, i)
    }
    model <- cbind(model[, 1, drop = FALSE], B, model[, -c(1, 
                                                           col_factors), drop = FALSE])
    rm(B, M, col_factors)
  }
  n <- dim(model)[1]
  nsmooth <- 1:length(sp)
  if (pages > 1) 
    pages = 1
  if (pages < 0) 
    pages = 0
  if (pages == 1) 
    par(mfrow = c(length(sp), 1))
  else par(mfrow = c(1, 1))
  SMAT <- KNOTS <- M <- list()
  for (i in nsmooth) {
    SMAT[[i]] <- object$smooth[[i]]$S[[1]]
    KNOTS[[i]] <- object$smooth[[i]]$knots
    M[[i]] <- object$smooth[[i]]$m
  }
  nlin <- ncol(model) - length(nsmooth) - 1
  x = as.matrix(model[, -c(1:(nlin + 1)), drop = FALSE])
  if (nlin != 0) {
    xlin = (model[, 2:(1 + nlin), drop = FALSE])
    xlin_dif <- pdata.frame(cbind(index_data, xlin), index = names(index_data))
    max <- length(colnames(xlin_dif))
    formula <- pFormula(as.formula(paste(c(paste(colnames(xlin_dif)[1], 
                                                 " ~ -1"), colnames(xlin_dif)[-c(1:2)]), collapse = "+")))
    xlin_dif <- mod.matrix(pFormula(formula), data = xlin_dif, 
                           effect = "individual", model = "fd")
  }
  else {
    xlin = NULL
    xlin_dif = NULL
  }
  splinemodel = data.frame(y = y_orig, x = x)
  names(splinemodel)[2:ncol(splinemodel)] = colnames(model)[-c(1:(nlin + 
                                                                    1))]
  rm(xlin)
  dimen <- nlin
  rm(object, model)
  crit = Stdev.fit = ucb = lcb = fitted = k0 = seqx = list()
  for (j in nsmooth) {
    Cnj = Znj = numeric()
    relength <- numeric()
    for (i in nsmooth) {
      data.grid = data.frame(x = x[, i])
      names(data.grid) <- names(splinemodel)[1 + i]
      xx = as.matrix(x)[, i]
      fitX = splineDesign(knots = KNOTS[[i]], x = xx, ord = M[[i]][1] + 
                            2, derivs = rep(0, length(xx)))
      C <- apply(fitX, 2, sum)
      C <- matrix(C, 1)
      qrc <- qr(t(C))
      Z <- qr.Q(qrc, complete = TRUE)[, (nrow(C) + 1):ncol(C)]
      Xx <- fitX %*% Z
      Xx <- pdata.frame(cbind(index_data, Xx), index = names(index_data))[, 
                                                                          -1]
      max <- length(colnames(Xx))
      colnames(Xx) <- c(colnames(Xx)[1], paste("V", colnames(Xx)[2:(max)], 
                                               sep = ""))
      formula <- pFormula(as.formula(c(paste(colnames(Xx)[1], 
                                             " ~ -1 +"), c(paste(colnames(Xx)[2:(max - 1)], 
                                                                 "+ ", sep = " "), colnames(Xx)[max]))))
      Xx <- mod.matrix(pFormula(formula), data = Xx, effect = "individual", 
                       model = "fd")
      Straf <- SMAT[[i]]
      eig <- eigen(Straf)
      U <- eig$vectors
      D <- diag(eig$values)
      fixed_coef <- M[[i]][2]
      if (fixed_coef >= 2) {
        dis <- (dim(D)[1] - fixed_coef + 1):dim(D)[1]
        D_p <- D[-dis, -dis]
        U_r <- U[, -dis]
        U_f <- U[, dis]
        X_f <- Xx %*% U_f
      }
      else {
        D_p <- D
        U_r <- U
        U_f <- NULL
        X_f <- matrix(0, dim(Xx)[1], 0)
      }
      X_r <- Xx %*% U_r
      Z <- X_r %*% sqrt(solve(D_p))
      relength[i] <- dim(Z)[2]
      if (i == j) {
        if (i == 1) 
          firstZ = 1
        else firstZ = ncol(Znj) + 1
        CZj <- list()
        CZj$Z <- Z
        CZj$knots <- KNOTS[[j]]
        CZj$C <- X_f
        dimen <- dimen + dim(Xx)[2]
        lastZ = firstZ + ncol(CZj$Z) - 1
        rm(X_r, X_f, Z, Xx)
      }
      else {
        CZ.temp <- list()
        CZ.temp$Z <- Z
        CZ.temp$knots <- KNOTS[[j]]
        CZ.temp$C <- X_f
        if (M[[i]][1] == 0) 
          CZ.temp$C = matrix(0, nrow(CZ.temp$Z), 0)
        Cnj = cbind(Cnj, CZ.temp$C)
        Znj = cbind(Znj, CZ.temp$Z)
        rm(Z, X_r, U_r, X_f, U_f, Xx)
      }
    }
    Zj = CZj$Z
    Cj = as.matrix(CZj$C[, drop = F])
    if (M[[j]][1] == 0) 
      Cj = matrix(0, nrow(Zj), 0)
    Xj = cbind(Cj, Zj)
    Cnj = cbind(xlin_dif, Cnj)
    Xnj = cbind(Cnj, Znj)
    Pen <- numeric()
    for (i in nsmooth) {
      Pen <- c(Pen, rep(sp[i], relength[i]))
    }
    Lambdaj = diag(c(rep(0, ncol(Cj)), rep(sp[j], length(firstZ:lastZ))))
    suppressWarnings(rm(Znj, Zj, Cj, CZ.temp))
    if (!gls) {
      if (ncol(x) > 1 | nlin != 0) {
        SnjXj = Xnj %*% (tcrossprod(solve(crossprod(Xnj) + 
                                            diag(c(rep(0, ncol(Cnj)), Pen[-((firstZ:lastZ))]), 
                                                 ncol = ncol(Xnj))), Xnj) %*% Xj)
        WjXj = Xj - SnjXj
      }
      else WjXj = Xj
    }
    else {
      if (ncol(x) > 1 | nlin != 0) {
        smiXnj <- smi %*% Xnj
        Sigma <- t(smi) %*% smi
        SnjXj = crossprod(t(Xnj), (tcrossprod(solve(crossprod(smiXnj) + 
                                                      diag(c(rep(0, ncol(Cnj)), Pen[-((firstZ:lastZ))]), 
                                                           ncol = ncol(Xnj))), crossprod(t(Sigma), Xnj)) %*% 
                                     Xj))
      }
      else {
        Sigma <- t(smi) %*% smi
        SnjXj = matrix(0, dim(Xj)[1], dim(Xj)[2])
      }
      WjXj = crossprod(t(Sigma), (Xj - SnjXj))
    }
    cov.coef = solve(crossprod(Xj, WjXj) + Lambdaj)
    if (gls) 
      WjXj = crossprod(t(smi), (Xj - SnjXj))
    suppressWarnings(rm(SnjXj))
    cc.ev <- eigen(cov.coef)
    cov.coef12 <- cc.ev$vectors %*% diag(sqrt(cc.ev$values)) %*% 
      t(cc.ev$vectors)
    rm(cc.ev, Xj, Xnj)
    integ <- function(xx, diffe) {
      data.grid <- data.frame(x = xx)
      names(data.grid) <- names(splinemodel)[1 + j]
      xhilf = as.matrix(x)[, j]
      hilf = splineDesign(knots = KNOTS[[j]], x = xhilf, 
                          ord = M[[j]][1] + 2, derivs = rep(0, length(xhilf)))
      C <- apply(hilf, 2, sum)
      C <- matrix(C, 1)
      qrc <- qr(t(C))
      Z <- qr.Q(qrc, complete = TRUE)[, (nrow(C) + 1):ncol(C)]
      if (drv == 0) 
        Xp <- splineDesign(knots = KNOTS[[j]], ord = M[[j]][1] + 
                             2, x = xx, derivs = rep(0, length(xx)))
      else Xp <- splineDesign(knots = KNOTS[[j]], ord = M[[j]][1] + 
                                2, x = xx, derivs = rep(drv, length(xx)))
      Xp <- (Xp %*% Z)
      Straf <- SMAT[[j]]
      eig <- eigen(Straf)
      U <- eig$vectors
      D <- diag(eig$values)
      fixed_coef <- M[[j]][2]
      if (fixed_coef >= 2) {
        dis <- (dim(D)[1] - fixed_coef + 1):dim(D)[1]
        D_p <- D[-dis, -dis]
        U_r <- U[, -dis]
        U_f <- U[, dis]
        X_f <- Xp %*% U_f
      }
      else {
        D_p <- D
        U_r <- U
        U_f <- NULL
        X_f <- matrix(0, dim(Xp)[1], 0)
      }
      X_r <- Xp %*% U_r
      X <- cbind(X_f, X_r)
      Z <- X_r %*% sqrt(solve(D_p))
      C <- X_f
      Cj.grid <- C[-diffe, , drop = F]
      Zj.grid = Z[-diffe, ]
      if (M[[j]][1] == 0) 
        Cj.grid = matrix(0, nrow(Zj.grid), 0)
      Xj.grid = cbind(Cj.grid, Zj.grid)
      SX <- tcrossprod(cov.coef12, (Xj.grid))
      SX.norm = sqrt(apply(SX^2, 2, sum))
      SX/SX.norm
    }
    sx = seq(min(x[, j]), max(x[, j]), length = div)
    k0[[j]] = sum(sqrt(apply((integ(sx, div) - integ(sx, 
                                                     1))^2, 2, sum)))
    crit[[j]] <- .C("scritval", k0 = as.numeric(c(k0[[j]], 
                                                  1)), d = as.integer(1), cov = as.numeric(level), 
                    m = as.integer(2), rdf = as.numeric(df), x = numeric(1), 
                    k = as.integer(1), PACKAGE = "pamfe")$x
    seqx[[j]] <- seq(min(x[, j]), max(x[, j]), length = pred)
    xx = as.matrix(x)[, j]
    fitX = splineDesign(knots = KNOTS[[j]], x = xx, ord = M[[j]][1] + 
                          2, derivs = rep(0, length(xx)))
    C <- apply(fitX, 2, sum)
    C <- matrix(C, 1)
    qrc <- qr(t(C))
    Z <- qr.Q(qrc, complete = TRUE)[, (nrow(C) + 1):ncol(C)]
    Xp <- splineDesign(knots = KNOTS[[j]], x = seqx[[j]], 
                       ord = M[[j]][1] + 2, derivs = rep(drv, length(seqx[[j]])))
    Xp <- (Xp %*% Z)
    Straf <- SMAT[[j]]
    eig <- eigen(Straf)
    U <- eig$vectors
    D <- diag(eig$values)
    fixed_coef <- M[[j]][2]
    if (fixed_coef >= 2) {
      dis <- (dim(D)[1] - fixed_coef + 1):dim(D)[1]
      D_p <- D[-dis, -dis]
      U_r <- U[, -dis]
      U_f <- U[, dis]
      X_f <- Xp %*% U_f
    }
    else {
      D_p <- D
      U_r <- U
      U_f <- NULL
      X_f <- matrix(0, dim(Xp)[1], 0)
    }
    X_r <- Xp %*% U_r
    X <- cbind(X_f, X_r)
    Z <- X_r %*% sqrt(solve(D_p))
    if (M[[j]][1] == 0) 
      X_f = matrix(0, nrow(Z), 0)
    Xj.grid = cbind(X_f, Z)
    Dimj <- dim(U_r)[1]
    rm(Z, X_f, X_r, U_r, U_f, D, U, Straf, C, qrc, D_p)
    fitted[[j]] <- Xp %*% coef[(dimen - Dimj + 1):dimen]
    Stdev.fit[[j]] <- sqrt(rowSums((Xj.grid %*% tcrossprod(cov.coef, 
                                                           (WjXj)))^2) * (sigma2))
    ucb[[j]] <- fitted[[j]] + crit[[j]] * Stdev.fit[[j]]
    lcb[[j]] <- fitted[[j]] - crit[[j]] * Stdev.fit[[j]]
    ylab <- paste("sfe", paste(rep("'", drv), collapse = "", 
                               "", sep = ""), "(", names(splinemodel)[j + 1], ")", 
                  sep = "")
    xlab <- names(splinemodel)[j + 1]
    plot(seqx[[j]], fitted[[j]], type = "l", xlab = xlab, 
         ylab = ylab, ylim = c(min(lcb[[j]]) - 0.1 * (max(lcb[[j]]) - 
                                                        min(lcb[[j]])), max(ucb[[j]]) + 0.1 * (max(ucb[[j]]) - 
                                                                                                 min(ucb[[j]]))))
    polygon(c(seqx[[j]], rev(seqx[[j]])), c(lcb[[j]], rev(ucb[[j]])), 
            col = grey(0.85), border = NA)
    points(seqx[[j]], lcb[[j]], type = "l", lty = "dotted", 
           col = grey(0.55))
    points(seqx[[j]], ucb[[j]], type = "l", lty = "dotted", 
           col = grey(0.55))
    lines(seqx[[j]], fitted[[j]], type = "l")
    if (drv > 0) 
      points(seqx[[j]], rep(0, length(seqx[[j]])), col = 1, 
             lwd = 1, type = "l")
  }
  scb <- list(crit = crit, seqx = seqx, Stdev = Stdev.fit, 
              sigma2 = sigma2, drv = drv, fitted = fitted, lcb = lcb, 
              ucb = ucb)
  scb
}
