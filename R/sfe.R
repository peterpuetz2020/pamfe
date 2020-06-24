sfe <-
function(..., k = -1, m = NA, sp = NULL) {
        xt <- NULL
        by <- NA
        fx <- FALSE
        id <- NULL
        bs <- "ps"
        vars <- as.list(substitute(list(...)))[-1]
        d <- length(vars)
        by.var <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
        if (by.var == ".")
            stop("by=. not allowed")
        term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
        if (term[1] == ".")
            stop("s(.) not yet supported.")
        if (d > 1)
            for (i in 2:d) {
                term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
                if (term[i] == ".")
                    stop("s(.) not yet supported.")
            }
        for (i in 1:d) term[i] <- attr(terms(reformulate(term[i])), "term.labels")
        k.new <- round(k)
        if (all.equal(k.new, k) != TRUE) {
            warning("argument k of s() should be integer and has been rounded")
        }
        k <- k.new
        if (length(unique(term)) != d)
            stop("Repeated variables as arguments of a smooth are not permitted")
        full.call <- paste("s(", term[1], sep = "")
        if (d > 1)
            for (i in 2:d) full.call <- paste(full.call, ",", term[i], sep = "")
        label <- paste(full.call, ")", sep = "")
        if (!is.null(id)) {
            if (length(id) > 1) {
                id <- id[1]
                warning("only first element of `id' used")
            }
            id <- as.character(id)
        }
        ret <- list(term = term, bs.dim = k, dim = d,
                    label = label, sp = sp)
        class(ret) <- paste(bs, ".smooth.spec", sep = "")
        ret
    }
