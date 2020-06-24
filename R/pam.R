pam <-
function (formula, data = list(), weights = NULL, method = "REML", 
              knots = NULL, optimizer = c("outer", "newton"), control = list(), 
              sp = NULL, gls = TRUE, corMatrix = list(), ...) 
    {
        control <- do.call("gam.control", control)
        if (class(data)[1] != "pdata.frame") 
            stop("Data must be of class pdata.frame from the plm package")
        if (!method %in% c("REML", "ML", "P-ML", "P-REML")) 
            stop("Unknown smoothing parameter estimation method")
        family <- gaussian()
        scale <- 0
        gamma <- 1
        H <- NULL
        na.action <- na.omit
        drop.unused.levels <- TRUE
        in.out <- NULL
        G <- NULL
        char <- as.character(formula)
        formula_orig <- as.formula(paste(char[2], char[1], char[3]), 
                                   env = .GlobalEnv)
        formula <- as.formula(paste(char[2], char[1], gsub("sfe", 
                                                           "s", char[3])), env = .GlobalEnv)
        gp <- interpret.gam(formula)
        lss <- length(gp$smooth.spec)
        if (lss > 0) {
            for (i in 1:lss) class(gp$smooth.spec[[i]]) <- "ps.smooth.spec"
        }
        cl <- match.call()
        mf <- match.call(expand.dots = FALSE)
        mf$formula <- gp$fake.formula
        mf$family <- mf$gls <-mf$corMatrix <- mf$control <- mf$scale <- mf$knots <- mf$sp <- mf$min.sp <- mf$H <- mf$gamma <- mf$method <- mf$fit <- mf$paraPen <- mf$G <- mf$optimizer <- mf$in.out <- mf$... <- NULL
        mf$drop.unused.levels <- drop.unused.levels
        mf[[1]] <- as.name("model.frame")
        pmf <- mf
        mf <- eval(mf, parent.frame())
        if (nrow(mf) < 2) 
            stop("Not enough (non-NA) data to do anything meaningful")
        terms <- attr(mf, "terms")
        vars <- all.vars(gp$fake.formula[-2])
        index_data <- index(data)
        if (any(is.na(data[, vars]))) {
            ind <- which(rowSums(is.na(data[, vars, drop = F])) != 
                             0)
            data <- data[-ind, , drop = F]
            index_data <- index_data[-ind, , drop = F]
        }
        col_fac <- which((lapply(data[, vars], is.factor) == TRUE))
        if (length(col_fac) != 0) {
            for (i in 1:length(col_fac)) A <- model.matrix(~-1 + 
                                                               data[, vars[col_fac[i]]])
            con_var <- which(as.numeric(apply(aggregate(A, by = list(index_data[, 
                                                                                1]), var)[, -1, drop = FALSE], 2, function(x) sum(as.numeric(x)))) < 
                                 1e-12)
            if (!identical((con_var), integer(0))) 
                stop(paste("Exclude variable", colnames(data[, vars[col_fac], 
                                                             drop = FALSE])[i], "from you analysis as it does not have any within variation."))
            vars <- vars[-col_fac]
        }
        con_var <- which(as.numeric(apply(aggregate(data[, vars], 
                                                    by = list(index_data[, 1]), var)[, -1, drop = FALSE], 
                                          2, function(x) sum(as.numeric(x)))) < 1e-12)
        if (length(con_var) != 0) 
            stop(paste("Exclude variable", colnames(data[, vars, 
                                                         drop = FALSE])[con_var[1]], "from you analysis as it does not have any within variation."))
        if (is.list(formula)) {
            environment(formula) <- environment(formula[[1]])
            pterms <- list()
            tlab <- rep("", 0)
            for (i in 1:length(formula)) {
                pmf$formula <- gp[[i]]$pf
                pterms[[i]] <- attr(eval(pmf, parent.frame()), "terms")
                tlabi <- attr(pterms[[i]], "term.labels")
                if (i > 1 && length(tlabi) > 0) 
                    tlabi <- paste(tlabi, i - 1, sep = ".")
                tlab <- c(tlab, tlabi)
            }
            attr(pterms, "term.labels") <- tlab
        }
        else {
            pmf$formula <- gp$pf
            pmf <- eval(pmf, parent.frame())
            pterms <- attr(pmf, "terms")
        }
        am <- TRUE
        if (!control$keepData) 
            rm(data)
        gsname <- if (is.list(formula)) 
            "gam.setup.list"
        else "gam.setup"
        G <- do.call(gsname, list(formula = gp, pterms = pterms, 
                                  data = mf, knots = knots, sp = sp, min.sp = NULL, H = H, 
                                  absorb.cons = TRUE, sparse.cons = 0, select = FALSE, 
                                  idLinksBases = control$idLinksBases, scale.penalty = control$scalePenalty, 
                                  paraPen = NULL, drop.intercept = TRUE), envir = environment(gam))
        G$family <- family
        if (ncol(G$X) > nrow(G$X)) 
            stop("Model has more coefficients than data")
        G$terms <- terms
        G$mf <- mf
        G$cl <- cl
        G$am <- am
        if (is.null(G$offset)) 
            G$offset <- rep(0, G$n)
        G$min.edf <- G$nsdf
        if (G$m) 
            for (i in 1:G$m) G$min.edf <- G$min.edf + G$smooth[[i]]$null.space.dim
        G$formula <- formula
        G$pred.formula <- gp$pred.formula
        environment(G$formula) <- environment(formula)
        G$conv.tol <- control$mgcv.tol
        G$max.half <- control$mgcv.half
        data <- pdata.frame(cbind(G$y, G$X, index_data), index = names(index_data))
        colnames(data)[1] <- "y"
        max <- length(colnames(data))
        colnames(data) <- c(colnames(data)[1], paste("V", colnames(data)[2:(max - 
                                                                                2)], sep = ""), colnames(data)[((max - 1):max)])
        names <- colnames(data)[-((max - 1):max)]
        max <- max - 2
        G$formula <- pFormula(as.formula(paste(c(paste(names[1], 
                                                       " ~ -1 +"), names[2:(max)]), collapse = "+")))
        G$X <- mod.matrix(pFormula(G$formula), data = data, effect = "individual", 
                          model = "fd")
        G$y <- do.call("pmodel.response", list(pFormula(G$formula), 
                                               data = data, effect = "individual", model = "fd"), envir = environment(plm))
        G$n <- length(G$y)
        ind_vector <- as.factor(paste(sapply(rownames(G$X), function(s) strsplit(s,                                                                        "-")[[1]][1])))
        if (gls) {
            smi <- list()
            indiv <- unique(ind_vector)
            for (j in 1:length(indiv)) {
                dim_matrix_j <- sum(ind_vector == indiv[j])
                W_inv <- diag(dim_matrix_j)
                W_inv[row(W_inv) - 1 == col(W_inv)] <- -0.5
                W_inv[row(W_inv) + 1 == col(W_inv)] <- -0.5
                smi[[j]] <- solve(t(chol(W_inv)))
            }
            smi <- (.bdiag(smi))
            smi <- as(smi, "CsparseMatrix")
            G$X <- as.matrix(smi %*% G$X)
            G$y <- as.numeric(smi %*% G$y)
        }
        
        if (length(corMatrix)>0) {
            corMatrix<-.bdiag(lapply(corMatrix,function (x) solve(t(chol(x)))))
            corMatrix <- as(corMatrix, "CsparseMatrix")
            G$X <- as.matrix(corMatrix %*% G$X)
            G$y <- as.numeric(corMatrix %*% G$y)
        }
        
        inddata <- index_data[, 1]
        ind <- which(inddata[-length(inddata)] != inddata[-1]) + 
            1
        ind <- c(1, ind)
        if (sum(G$w == rep(1, length(G$w))) != length(G$w)) 
            G$w <- G$w[-ind]
        else G$w <- rep(1, G$n)
        G$w <- G$w/mean(G$w)
        G$offset <- rep(0, G$n)
        object <- do.call("estimate.gam", list(G, method, optimizer, 
                                               control, in.out, scale, gamma, ...), envir = environment(gam))
        if (!is.null(G$L)) {
            object$full.sp <- as.numeric(exp(G$L %*% log(object$sp) + 
                                                 G$lsp0))
            names(object$full.sp) <- names(G$lsp0)
        }
        names(object$sp) <- names(G$sp)
        if (gls) {
            object$smi <- smi
        }
        object$paraPen <- G$pP
        object$formula <- formula
        object$model <- G$mf
        object$na.action <- attr(G$mf, "na.action")
        object$control <- control
        object$terms <- G$terms
        object$pterms <- G$pterms
        object$assign <- G$assign
        object$contrasts <- G$contrasts
        object$index_data <- index_data
        object$index_diffdata <- ind_vector
        object$weights <- G$w
        object$residuals[G$w == 0] <- 0
        if (!is.null(G$Xcentre)) 
            object$Xcentre <- G$Xcentre
        if (control$keepData) 
            object$data <- data
        object$n <- object$df.null <- nrow(G$X)
        object$df.residual <- nrow(G$X) - sum(object$edf)
        object$sig2 <- sum(object$residuals^2 * G$w)/object$df.residual
        object$min.edf <- G$min.edf
        object$optimizer <- optimizer
        object$formula <- formula_orig
        object$gls <- gls
        object$linear.predictors <- NULL
        object$deviance <- object$null.deviance <- NULL
        object$offset <- NULL
        object$xlevels <- NULL
        object$cmX <- NULL
        object$min.edf <- NULL
        object$pred.formula <- NULL
        object$working.weights <- NULL
        object$rV <- NULL
        object$dw.drho <- NULL
        class(object) <- c("pam", "gam", "glm", "lm")
        object
    }
