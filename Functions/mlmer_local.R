# Function from omics package
# Removed unnecessary calculus to obtain only residuals.

mlmer_local <- function(formula, data=NULL, vars, lrt=TRUE, save.residuals=FALSE,
                        save.ranks=TRUE, ncores = 2L) {
    if (missing(data)) {
        data <- NULL
    }
    
    Y <- get(omics:::response.name(formula, data=data), envir=environment(formula))
    
    tmp <- formula
    tmp[[2]] <- NULL
    lf <- lme4::lFormula(tmp, data, REML=FALSE, na.action=na.pass,
                         control=lme4::lmerControl(check.nobs.vs.nlev="ignore",
                                                   check.nobs.vs.nRE="ignore",
                                                   check.rankX="ignore",
                                                   check.scaleX="ignore",
                                                   check.formula.LHS="ignore"))
    labs <- labels(terms(lme4::nobars(formula), data=data))
    if (missing(vars)) {
        vars <- labs
    }
    
    mm <- lf$X
    idx <- which(attr(lf$X, "assign") %in% match(vars, labs))
    if (length(idx) == 0 & !save.residuals) {
        stop("No variables selected")
    }

    vars <- colnames(mm)[idx]
    colnames(mm) <- sprintf("V%d", 1:ncol(mm))

    formula <- as.formula(sprintf("y ~ %s - 1",
                                  paste0(c(
                                      sprintf("(%s)", lme4::findbars(formula)), colnames(mm)
                                  ), collapse=" + ")
    ))
    re.labs <- as.character(attr(terms(lf$fr), "predvars.random")[-1])
    model.data <- data.frame(mm, mget(re.labs, as.environment(lf$fr)))
    
    if (save.residuals) {
        residuals <- matrix(NA, nrow(Y), ncol(Y), dimnames=dimnames(Y))
    }

    res <- lapply( 1:ncol(Y), function(i) {
        model.data$y <- Y[,i]
        model <- try(
            lme4::lmer(formula, data=model.data, REML=FALSE, na.action=na.exclude),
            silent=TRUE )
        if (inherits(model, "try-error")) {
            return(NA)
        }
        
        return( list(residual = resid(model, type="response"),
                     name = colnames(Y)[i] ))
    })
    
    residuals = do.call(rbind,  lapply(res, "[[", "residual"))
    rownames(residuals) = do.call(rbind,  lapply(res, "[[", "name"))
    
    return( residuals)
}
