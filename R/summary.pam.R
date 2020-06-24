summary.pam <-
function(object, freq = FALSE, ...) {
        p.type = 0
        dispersion = NULL
        pinv <- function(V, M, rank.tol = 1e-06) {
            D <- eigen(V, symmetric = TRUE)
            M1 <- length(D$values[D$values > rank.tol * D$values[1]])
            if (M > M1)
                M <- M1
            if (M + 1 <= length(D$values))
                D$values[(M + 1):length(D$values)] <- 1
            D$values <- 1/D$values
            if (M + 1 <= length(D$values))
                D$values[(M + 1):length(D$values)] <- 0
            res <- D$vectors %*% (D$values * t(D$vectors))
            attr(res, "rank") <- M
            res
        }
        if (is.null(object$R)) {
            warning("p-values for any terms that can be penalized to zero will be unreliable:
                    refit model to fix this.")
            useR <- FALSE
        } else useR <- TRUE
        p.table <- pTerms.table <- s.table <- NULL
        if (freq)
            covmat <- object$Ve else covmat <- object$Vp
        name <- names(object$edf)
        dimnames(covmat) <- list(name, name)
        covmat.unscaled <- covmat/object$sig2
        est.disp <- object$scale.estimated
        if (!is.null(dispersion)) {
            covmat <- dispersion * covmat.unscaled
            object$Ve <- object$Ve * dispersion/object$sig2
            object$Vp <- object$Vp * dispersion/object$sig2
            est.disp <- FALSE
        } else dispersion <- object$sig2
        se <- diag(covmat)^0.5
        residual.df <- length(object$y) - sum(object$edf)
        if (sum(object$nsdf) > 0) {
            if (length(object$nsdf) > 1) {
                pstart <- attr(object$nsdf, "pstart")
                ind <- rep(0, 0)
                for (i in 1:length(object$nsdf)) if (object$nsdf[i] > 0)
                    ind <- c(ind, pstart[i]:(pstart[i] + object$nsdf[i] - 1))
            } else {
                pstart <- 1
                ind <- 1:object$nsdf
            }
            p.coeff <- object$coefficients[ind]
            p.se <- se[ind]
            p.t <- p.coeff/p.se
            if (!est.disp) {
                p.pv <- 2 * pnorm(abs(p.t), lower.tail = FALSE)
                p.table <- cbind(p.coeff, p.se, p.t, p.pv)
                dimnames(p.table) <- list(names(p.coeff), c("Estimate", "Std. Error", "z value",
                                                            "Pr(>|z|)"))
            } else {
                p.pv <- 2 * pt(abs(p.t), df = residual.df, lower.tail = FALSE)
                p.table <- cbind(p.coeff, p.se, p.t, p.pv)
                dimnames(p.table) <- list(names(p.coeff), c("Estimate", "Std. Error", "t value",
                                                            "Pr(>|t|)"))
            }
        } else {
            p.coeff <- p.t <- p.pv <- array(0, 0)
        }
        pterms <- if (is.list(object$pterms))
            object$pterms else list(object$pterms)
        if (!is.list(object$assign))
            object$assign <- list(object$assign)
        npt <- length(unlist(lapply(pterms, attr, "term.labels")))
        if (npt > 0)
            pTerms.df <- pTerms.chi.sq <- pTerms.pv <- array(0, npt)
        term.labels <- rep("", 0)
        k <- 0
        for (j in 1:length(pterms)) {
            tlj <- attr(pterms[[j]], "term.labels")
            nt <- length(tlj)
            if (j > 1 && nt > 0)
                tlj <- paste(tlj, j - 1, sep = ".")
            term.labels <- c(term.labels, tlj)
            if (nt > 0) {
                np <- length(object$assign[[j]])
                ind <- pstart[j] - 1 + 1:np
                Vb <- covmat[ind, ind, drop = FALSE]
                bp <- array(object$coefficients[ind], np)
                for (i in 1:nt) {
                    k <- k + 1
                    ind <- object$assign[[j]] == i
                    b <- bp[ind]
                    V <- Vb[ind, ind]
                    if (length(b) == 1) {
                        V <- 1/V
                        pTerms.df[k] <- nb <- 1
                        pTerms.chi.sq[k] <- V * b * b
                    } else {
                        V <- pinv(V, length(b), rank.tol = .Machine$double.eps^0.5)
                        pTerms.df[k] <- nb <- attr(V, "rank")
                        pTerms.chi.sq[k] <- t(b) %*% V %*% b
                    }
                    if (!est.disp)
                        pTerms.pv[k] <- pchisq(pTerms.chi.sq[k], df = nb, lower.tail = FALSE) else 
                            pTerms.pv[k] <- pf(pTerms.chi.sq[k]/nb, df1 = nb, df2 = residual.df,
                                               lower.tail = FALSE)
                }
            }
        }
        if (npt) {
            attr(pTerms.pv, "names") <- term.labels
            if (!est.disp) {
                pTerms.table <- cbind(pTerms.df, pTerms.chi.sq, pTerms.pv)
                dimnames(pTerms.table) <- list(term.labels, c("df", "Chi.sq", "p-value"))
            } else {
                pTerms.table <- cbind(pTerms.df, pTerms.chi.sq/pTerms.df, pTerms.pv)
                dimnames(pTerms.table) <- list(term.labels, c("df", "F", "p-value"))
            }
        } else {
            pTerms.df <- pTerms.chi.sq <- pTerms.pv <- array(0, 0)
        }
        m <- length(object$smooth)
        df <- edf1 <- edf <- s.pv <- chi.sq <- array(0, m)
        if (m > 0) {
            if (p.type < 5) {
                if (useR)
                    X <- object$R else {
                        sub.samp <- max(1000, 2 * length(object$coefficients))
                        if (nrow(object$model) > sub.samp) {
                            seed <- try(get(".Random.seed", envir = .GlobalEnv), silent = TRUE)
                            if (inherits(seed, "try-error")) {
                                runif(1)
                                seed <- get(".Random.seed", envir = .GlobalEnv)
                            }
                            kind <- RNGkind(NULL)
                            RNGkind("default", "default")
                            set.seed(11)
                            ind <- sample(1:nrow(object$model), sub.samp, replace = FALSE)
                            X <- predict(object, object$model[ind, ], type = "lpmatrix")
                            RNGkind(kind[1], kind[2])
                            assign(".Random.seed", seed, envir = .GlobalEnv)
                        } else {
                            X <- model.matrix(object)
                        }
                        X <- X[!is.na(rowSums(X)), ]
                    }
            }
            for (i in 1:m) {
                start <- object$smooth[[i]]$first.para
                stop <- object$smooth[[i]]$last.para
                V <- object$Vp[start:stop, start:stop, drop = FALSE]
                p <- object$coefficients[start:stop]
                edf1[i] <- edf[i] <- sum(object$edf[start:stop])
                if (!is.null(object$edf1))
                    edf1[i] <- sum(object$edf1[start:stop])
                Xt <- X[, start:stop, drop = FALSE]
                if (object$smooth[[i]]$null.space.dim == 0 && !is.null(object$R)) {
                    res <- do.call("reTest", list(object, i), envir = environment(summary.gam))
                } else {
                    df[i] <- min(ncol(Xt), edf1[i])
                    if (est.disp)
                        rdf <- residual.df else rdf <- -1
                        res <- do.call("testStat", list(p, Xt, V, df[i], type = p.type, res.df = rdf),
                                       envir = environment(summary.gam))
                }
                df[i] <- res$rank
                chi.sq[i] <- res$stat
                s.pv[i] <- res$pval
                names(chi.sq)[i] <- gsub("s", "sfe", object$smooth[[i]]$label)
            }
            if (!est.disp) {
                s.table <- cbind(edf, df, chi.sq, s.pv)
                dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "Chi.sq", "p-value"))
            } else {
                s.table <- cbind(edf, df, chi.sq/df, s.pv)
                dimnames(s.table) <- list(names(chi.sq), c("edf", "Ref.df", "F", "p-value"))
            }
        }
        w <- as.numeric(object$prior.weights)
        mean.y <- sum(w * object$y)/sum(w)
        w <- sqrt(w)
        nobs <- object$n
        r.sq <- if (inherits(object$family, "general.family") || !is.null(object$family$no.r.sq))
            NULL else 1 - var(w * (as.numeric(object$y) - object$fitted.values)) * (nobs - 1)/
            (var(w * (as.numeric(object$y) - mean.y)) * residual.df)
        if (object$method %in% c("REML", "ML"))
            object$method <- paste("-", object$method, sep = "")
        ret <- list(p.coeff = p.coeff, se = se, p.t = p.t, p.pv = p.pv,
                    residual.df = residual.df, m = m, chi.sq = chi.sq, s.pv = s.pv,
                    scale = dispersion, r.sq = r.sq, family = object$family,
                    formula = object$formula, n = nobs, edf = edf, pTerms.pv = pTerms.pv,
                    pTerms.chi.sq = pTerms.chi.sq, pTerms.df = pTerms.df,
                    cov.unscaled = covmat.unscaled, cov.scaled = covmat, p.table = p.table,
                    pTerms.table = pTerms.table, s.table = s.table, method = object$method,
                    sp.criterion = object$gcv.ubre, np = length(object$coefficients))
        class(ret) <- "summary.gam"
        ret
    }
