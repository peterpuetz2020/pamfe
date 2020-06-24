p.diff <-
function(x, cond) {
        cond <- as.numeric(cond)
        n <- ifelse(is.matrix(x), nrow(x), length(x))
        cond <- c(NA, cond[2:n] - cond[1:(n - 1)])
        cond[cond != 0] <- NA
        result <- rbind(NA, x[2:n, , drop = FALSE] - x[1:(n - 1), , drop = FALSE])
        result[is.na(cond), ] <- NA
        result <- na.omit(result)
        attr(result, "na.action") <- NULL
        result
    }
