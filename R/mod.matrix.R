mod.matrix <-
function(object, data, model = c("within", "fd"), ...) {
        effect <- "individual"
        rhs <- 1
        theta <- NULL
        model <- match.arg(model)
        formula <- object
        has.intercept <- FALSE
        X <- model.matrix(do.call("as.Formula", list(formula), 
                                  envir = environment(plm)), rhs = rhs, data = data, ...)
        index <- attr(data, "index")
        id <- index[[1]]
        if (any(is.na(id))) {
            stop("NA in the individual index variable")
        }
        time <- index[[2]]
        pdim <- pdim(data)
        cond <- id
        result <- switch(model, within = Within(X, cond), fd = p.diff(X, cond))
        result
    }
