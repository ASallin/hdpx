# helpers SCRIPT ----------------------------------------------------------

## Helper functions

## Functions for improved error bounds
### Modified version of the code accompanying the paper:
###   Shah, R. D. and Samworth, R. J. (2013), Variable selection with error
###   control: Another look at Stability Selection, J. Roy. Statist. Soc., Ser.
###   B, 75, 55-80. DOI: 10.1111/j.1467-9868.2011.01034.x
###
### Original code available from
###   http://www.statslab.cam.ac.uk/~rds37/papers/r_concave_tail.R
### or
###   http://www.statslab.cam.ac.uk/~rjs57/r_concave_tail.R
D <- function(theta, which, B, r) {
  ## compute upper tail of r-concave distribution function
  ## If q = ceil{ B * 2 * theta} / B + 1/B,..., 1 return the tail probability.
  ## If q < ceil{ B * 2 * theta} / B return 1
  
  s <- 1/r
  thetaB <- theta * B
  k_start <- (ceiling(2 * thetaB) + 1)
  
  if (which < k_start)
    return(1)
  
  if(k_start > B)
    stop("theta to large")
  
  Find.a <- function(prev_a)
    uniroot(Calc.a, lower = 0.00001, upper = prev_a,
            tol = .Machine$double.eps^0.75)$root
  
  Calc.a <- function(a) {
    denom <- sum((a + 0:k)^s)
    num <- sum((0:k) * (a + 0:k)^s)
    num / denom - thetaB
  }
  
  OptimInt <- function(a, t, k, thetaB, s) {
    num <- (k + 1 - thetaB) * sum((a + 0:(t-1))^s)
    denom <- sum((k + 1 - (0:k)) * (a + 0:k)^s)
    1 - num / denom
  }
  
  ## initialize a
  a_vec <- rep(100000, B)
  
  ## compute a values
  for(k in k_start:B)
    a_vec[k] <- Find.a(a_vec[k-1])
  
  cur_optim <- rep(0, B)
  for (k in k_start:(B-1))
    cur_optim[k] <- optimize(f=OptimInt, lower = a_vec[k+1],
                             upper = a_vec[k],
                             t = which, k = k, thetaB = thetaB, s = s,
                             maximum  = TRUE)$objective
  return(max(cur_optim))
}

## minD function for error bound in case of r-concavity
minD <- function(q, p, pi, B, r = c(-1/2, -1/4)) {
  ## get the integer valued multiplier W of
  ##   pi = W * 1/(2 * B)
  which <- ceiling(signif(pi / (1/(2* B)), 10))
  maxQ <- maxQ(p, B)
  if (q > maxQ)
    stop(sQuote("q"), " must be <= ", maxQ,  call. = FALSE)
  min(c(1, D(q^2 / p^2, which - B, B, r[1]), D(q / p, which , 2*B, r[2])))
}


## Functions to compute the optimal cutoff and optimal q values

## function to find optimal cutoff in stabsel (when sampling.type = "SS")
optimal_cutoff <- function(p, q, PFER, B, assumption = "unimodal") {
  if (assumption == "unimodal") {
    ## cutoff values can only be multiples of 1/(2B)
    cutoffgrid <- 1/2 + (2:B)/(2*B)
    c_min <- min(0.5 + (q/p)^2, 0.5 + 1/(2*B) + 0.75 * (q/p)^2)
    cutoffgrid <- cutoffgrid[cutoffgrid > c_min]
    
    if (length(cutoffgrid) == 0) {
      maxQ <- max(floor(p * sqrt(0.5)),
                  floor(p * sqrt((0.5 - 1/(2*B)) * 4/3)))
      stop(sQuote("q"), " must be <= ", maxQ,  call. = FALSE)
    }
    
    upperbound <- rep(NA, length(cutoffgrid))
    for (i in 1:length(cutoffgrid))
      upperbound[i] <- q^2 / p / um_const(cutoffgrid[i], B, theta = q/p)
    cutoff <- cutoffgrid[upperbound < PFER][1]
    return(cutoff)
  } else {
    ## cutoff values can only be multiples of 1/(2B)
    cutoff <- (2*B):1/(2*B)
    cutoff <- cutoff[cutoff >= 0.5]
    for (i in 1:length(cutoff)) {
      if (minD(q, p, cutoff[i], B) * p > PFER) {
        if (i == 1)
          cutoff <- cutoff[i]
        else
          cutoff <- cutoff[i - 1]
        break
      }
    }
    return(tail(cutoff, 1))
  }
}

## function to find optimal q in stabsel (when sampling.type = "SS")
optimal_q <- function(p, cutoff, PFER, B, assumption = "unimodal") {
  if (assumption == "unimodal") {
    if (cutoff <= 0.75) {
      upper_q <- max(p * sqrt(cutoff - 0.5),
                     p * sqrt(4/3 * (cutoff - 0.5 - 1/(2*B))))
      ## q must be an integer < upper_q
      upper_q <- ceiling(upper_q - 1)
    } else {
      upper_q <- p
    }
    q <- uniroot(function(q)
      q^2 / p / um_const(cutoff, B, theta = q/p) - PFER,
      lower = 1, upper = upper_q)$root
    return(floor(q))
  } else {
    for (q in 1:maxQ(p, B)) {
      if (minD(q, p, cutoff, B) * p > PFER) {
        q <- q - 1
        break
      }
    }
    return(max(1, q))
  }
}

## obtain maximal value possible for q
maxQ <- function(p, B) {
  if(B <= 1)
    stop("B must be at least 2", call. = FALSE)
  
  fact_1 <- 4 * B / p
  tmpfct <- function(q)
    ceiling(q * fact_1) + 1 - 2 * B
  
  res <- tmpfct(1:p)
  length(res[res < 0])
}

## obtain constant for unimodal bound
um_const <- function(cutoff, B, theta) {
  if (cutoff <= 3/4) {
    if (cutoff < 1/2 + min(theta^2, 1 / (2*B) + 3/4 * theta^2))
      stop ("cutoff out of bounds", call. = FALSE)
    return( 2 * (2 * cutoff - 1 - 1/(2*B)) )
  } else {
    if (cutoff > 1)
      stop ("cutoff out of bounds", call. = FALSE)
    return( (1 + 1/B)/(4 * (1 - cutoff + 1 / (2*B))) )
  }
}


## Pre-processing functions for stabsel

## check if folds result from subsampling with p = 0.5.
check_folds <- function(folds, B, n, sampling.type) {
  if (!is.matrix(folds) || ncol(folds) != B || nrow(folds) != n ||
      !all(folds %in% c(0, 1)))
    stop(sQuote("folds"),
         " must be a binary or logical matrix with dimension nrow(x) times B",
         call. = FALSE)
  if (!all(colMeans(folds) %in% c(floor(n * 0.5) / n, ceiling(n * 0.5) / n)))
    warning("Subsamples are not of size n/2; results might be wrong",
            call. = FALSE)
  ## use complementary pairs?
  if (sampling.type == "SS") {
    folds <- cbind(folds, rep(1, n) - folds)
  }
  folds
}


# Subsamlpling method -----------------------------------------------------


## modified version of mboost's cv function
subsample <- function(weights, B = 100, strata = NULL) {
  n <- length(weights)
  
  if (is.null(strata)) strata <- gl(1, n)
  if (!is.factor(strata))
    stop(sQuote("strata"), " must be a factor")
  folds <- matrix(0, nrow = n, ncol = B)
  
  make_subsample <- function(n, B) {
    k <- floor(n * 0.5)
    indx <- rep(c(0, 1), c(n - k, k))
    replicate(B, sample(indx))[sample(1:n),, drop = FALSE]
  }
  
  ### <FIXME> handling of weights needs careful documentation </FIXME>
  for (s in levels(strata)) {
    indx <- which(strata == s)
    folds[indx,] <- make_subsample(length(indx), B = B) * weights[indx]
  }
  attr(folds, "type") <- paste(B, "-fold subsampling", sep = "")
  return(folds)
}



# fitfun SCRIPT -----------------------------------------------------------

glmnet.lasso <- function(x, y, q, type = c("conservative", "anticonservative"), ...) {
  if (!requireNamespace("glmnet", quietly=TRUE))
    stop("Package ", sQuote("glmnet"), " needed but not available")
  
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"),
            " is coerced to a model matrix without intercept")
    x <- model.matrix(~ . - 1, x)
  }
  
  if ("lambda" %in% names(list(...)))
    stop("It is not permitted to specify the penalty parameter ", sQuote("lambda"),
         " for lasso when used with stability selection.")
  
  ## fit model
  type <- match.arg(type)
  if (type == "conservative")
    fit <- suppressWarnings(glmnet::glmnet(x, y, pmax = q, ...))
  if (type == "anticonservative")
    fit <- glmnet::glmnet(x, y, dfmax = q - 1, ...)
  
  ## which coefficients are non-zero?
  selected <- predict(fit, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  ## compute selection paths
  cf <- fit$beta
  sequence <- as.matrix(cf != 0)
  ## return both
  return(list(selected = ret, path = sequence,
              coefs = cf)) # AS: return coefficient size
}

lars.lasso <- function(x, y, q, ...) {
  if (!requireNamespace("lars", quietly=TRUE))
    stop("Package ", sQuote("lars"), " needed but not available")
  
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"),
            " is coerced to a model matrix without intercept")
    x <- model.matrix(~ . - 1, x)
  }
  
  ## fit model
  fit <- lars::lars(x, y, max.steps = q, ...)
  
  ## which coefficients are non-zero?
  selected <- unlist(fit$actions)
  ## check if variables are removed again from the active set
  ## and remove these from selected
  if (any(selected < 0)) {
    idx <- which(selected < 0)
    idx <- c(idx, which(selected %in% abs(selected[idx])))
    selected <- selected[-idx]
  }
  
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  ## compute selection paths
  cf <- fit$beta
  sequence <- t(cf != 0)
  ## return both
  return(list(selected = ret, path = sequence,
              coefs = cf)) # AS: return coefficient size
}

lars.stepwise <- function(x, y, q, ...) {
  if (!requireNamespace("lars", quietly=TRUE))
    stop("Package ", sQuote("lars"), " needed but not available")
  
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"),
            " is coerced to a model matrix without intercept")
    x <- model.matrix(~ . - 1, x)
  }
  
  ## fit model
  fit <- lars::lars(x, y, max.steps = q, type = "stepwise", ...)
  
  ## which coefficients are non-zero?
  selected <- unlist(fit$actions)
  ## check if variables are removed again from the active set
  ## and remove these from selected
  if (any(selected < 0)) {
    idx <- which(selected < 0)
    idx <- c(idx, which(selected %in% abs(selected[idx])))
    selected <- selected[-idx]
  }
  
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  ## compute selection paths
  cf <- fit$beta
  sequence <- t(cf != 0)
  ## return both
  return(list(selected = ret, path = sequence,
              coefs = cf)) # AS: return coefficient size
}

glmnet.lasso_maxCoef <- function(x, y, q, ...) {
  if (!requireNamespace("glmnet", quietly=TRUE))
    stop("Package ", sQuote("glmnet"), " needed but not available")
  
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"),
            " is coerced to a model matrix without intercept")
    x <- model.matrix(~ . - 1, x)
  }
  
  args <- list(...)
  if (!("lambda" %in% names(args)) && length(args$lambda) != 1)
    stop("Please specify a fixed (!) value of ", sQuote("lambda"),
         ", which is small enough that at least ", sQuote("q"),
         " variables can be selected.")
  
  ## fit model
  fit <- glmnet::glmnet(x, y, ...)
  
  ## which coefficients are the q biggest
  selected <- order(coef(fit)[-1])[1:q]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  ## return selection
  return(list(selected = ret, path = NULL,
              coefs = cf)) # AS: return coefficient size
}






# Stabsel matrix ----------------------------------------------------------

## Generic implementation of stability selection
stabsel <- function(x, ...) {
  UseMethod("stabsel", x)
}


### TODO: parallelization ala cvrisk needed! Seems to work?
### TODO: Use same arguments for .mboost und .formula
### TODO: Should y be a matrix? Perhaps we need this for survival data which
### might be specified as a matrix?
stabsel.matrix <- function(x, y, fitfun = lars.lasso, args.fitfun = list(),
                           cutoff, q, PFER, 
                           folds = subsample(rep(1, nrow(x)), B = B),
                           B = ifelse(sampling.type == "MB", 100, 50),
                           assumption = c("unimodal", "r-concave", "none"),
                           sampling.type = c("SS", "MB"),
                           papply = mclapply, mc.preschedule = FALSE,
                           verbose = TRUE, FWER, eval = TRUE,
                           ...) {
  cll <- match.call()
  p <- ncol(x) ## TODO: what about intercept?
  
  ## use graphical model structure if fitfun is the right class
  graphical <- inherits(fitfun, "graphical_model")
  if (missing(y)) {
    if (!graphical) {
      ## probably meant to be a graphical model - here we issue warnings
      ## and set the function class if the function isn't a tagged as a 
      ## graphical model
      warning("No ", sQuote("y"), " supplied and ", sQuote("fitfun"), 
              " is not of class ", dQuote("graphical_model"),"\n",
              "however proceeding with graphical model analysis")
      graphical <- TRUE
      class(fitfun) <- c(class(fitfun), "graphical_model")
    }
    # set p and y for the graphical case
    y <- x
    p <- p * (p-1)/2
  } else {
    if (graphical) {
      stop("Both ", sQuote("y"), " and a graphical_model ", sQuote("fitfun"), " supplied")
    }
  }
  n <- nrow(x)
  
  if (is.null(colnames(x)))
    colnames(x) <- 1:ncol(x)
  
  ## needed here to make B and folds happy
  sampling.type <- match.arg(sampling.type)
  if (sampling.type == "MB")
    assumption <- "none"
  else
    assumption <- match.arg(assumption)
  
  ## make sure that y is a matrix
  if (!is.matrix(y))
    y <- matrix(y, ncol = 1)
  
  if (n != nrow(y))
    stop(sQuote("x"), " and ", sQuote("y"),
         " must have the same number of observations")
  
  ## define fitting function;
  ## the function implicitly knows x and y as it is defined in this environment
  fit_model <- function(i, folds, q, args.fitfun) {
    inbag <- as.logical(folds[, i])
    do.call(fitfun, c(list(x = x[inbag, ], y = y[inbag, ], q = q),
                      args.fitfun))
  }
  
  nms <- colnames(x)
  
  if (graphical) {
    ## graphical models need different names
    allnms <- outer(nms, nms, paste, sep=" : ")
    nms <- allnms[upper.tri(allnms)]
  }
  
  ret <- run_stabsel(fitter = fit_model, args.fitter = args.fitfun,
                     n = n, p = p, cutoff = cutoff, q = q, 
                     PFER = PFER, folds = folds, B = B, assumption = assumption,
                     sampling.type = sampling.type, papply = papply,
                     verbose = verbose, FWER = FWER, eval = eval, names = nms,
                     mc.preschedule = mc.preschedule, ...)
  ret$call <- cll
  ret$call[[1]] <- as.name("stabsel")
  return(ret)
}

stabsel.data.frame <- function(x, y, intercept = FALSE, ...) {
  if (intercept) {
    x <- model.matrix(~ ., x)
  } else {
    x <- model.matrix(~ . - 1, x)
  }
  stabsel(x, y, ...)
}



# stabsel_parameters ------------------------------------------------------

stabsel_parameters <- function(p, ...)
  UseMethod("stabsel_parameters")

stabsel_parameters.default <- function(p, cutoff, q, PFER, 
                                       B = ifelse(sampling.type == "MB", 100, 50),
                                       assumption = c("unimodal", "r-concave", "none"),
                                       sampling.type = c("SS", "MB"),
                                       verbose = FALSE, FWER, ...) {
  
  sampling.type <- match.arg(sampling.type)
  if (sampling.type == "MB")
    assumption <- "none"
  else
    assumption <- match.arg(assumption)
  
  
  ## only two of the four arguments can be specified
  if ((nmiss <- sum(missing(PFER), missing(cutoff),
                    missing(q), missing(FWER))) != 2) {
    if (nmiss > 2)
      stop("Two of the three argumnets ",
           sQuote("PFER"), ", ", sQuote("cutoff"), " and ", sQuote("q"),
           " must be specifed")
    if (nmiss < 2)
      stop("Only two of the three argumnets ",
           sQuote("PFER"), ", ", sQuote("cutoff"), " and ", sQuote("q"),
           " can be specifed at the same time")
  }
  
  if (!missing(FWER)) {
    if (!missing(PFER))
      stop(sQuote("FWER"), " and ", sQuote("PFER"),
           " cannot be spefified at the same time")
    PFER <- FWER
    warning(sQuote("FWER"), " is deprecated. Use ", sQuote("PFER"),
            " instead.")
  }
  
  if ((!missing(PFER) || !missing(FWER)) && PFER < 0)
    stop(sQuote("PFER"), " must be greater 0")
  
  if (!missing(cutoff) && (cutoff < 0.5 | cutoff > 1))
    stop(sQuote("cutoff"), " must be between 0.5 and 1")
  
  if (!missing(q)) {
    if (p < q)
      stop("Average number of selected effects ", sQuote("q"),
           " must be smaller \n  than the number of effects",
           " specified in the model ", sQuote("object"))
    if (q < 0)
      stop("Average number of selected effects ", sQuote("q"),
           " must be greater 0")
  }
  
  if (missing(cutoff)) {
    if (assumption == "none") {
      cutoff <- min(1, tmp <- (q^2 / (PFER * p) + 1) / 2)
      upperbound <- q^2 / p / (2 * cutoff - 1)
    } else {
      if (assumption == "unimodal") {
        cutoff <- tmp <- optimal_cutoff(p, q, PFER, B,
                                        assumption = assumption)
        upperbound <- q^2 / p / um_const(cutoff, B, theta = q/p)
      } else {
        cutoff <- tmp <- optimal_cutoff(p, q, PFER, B,
                                        assumption = assumption)
        upperbound <- minD(q, p, cutoff, B) * p
      }
    }
    upperbound <- signif(upperbound, 3)
    if (verbose && tmp > 0.9 && upperbound - PFER > PFER/2) {
      warning("Upper bound for PFER > ", PFER,
              " for the given value of ", sQuote("q"),
              " (true upper bound = ", round(upperbound, 2), ")")
    }
  }
  
  if (missing(q)) {
    if (assumption == "none") {
      q <- floor(sqrt(PFER * (2 * cutoff - 1) * p))
      upperbound <- q^2 / p / (2 * cutoff - 1)
    } else {
      if (assumption == "unimodal") {
        q <- optimal_q(p, cutoff, PFER, B, assumption = assumption)
        upperbound <- q^2 / p / um_const(cutoff, B, theta = q/p)
      } else {
        q <- optimal_q(p, cutoff, PFER, B, assumption = assumption)
        upperbound <- minD(q, p, cutoff, B) * p
      }
    }
    upperbound <- signif(upperbound, 3)
    if (verbose && upperbound - PFER > PFER/2)
      warning("Upper bound for PFER > ", PFER,
              " for the given value of ", sQuote("cutoff"),
              " (true upper bound = ", upperbound, ")")
  }
  
  if (missing(PFER)) {
    if (assumption == "none") {
      upperbound <- PFER <- q^2 / p / (2 * cutoff - 1)
    } else {
      if (assumption == "unimodal") {
        upperbound <- PFER <- q^2 / p / um_const(cutoff, B, theta = q/p)
      } else {
        upperbound <- PFER <- minD(q, p, cutoff, B) * p
      }
    }
    upperbound <- signif(upperbound, 3)
  }
  
  if (verbose && PFER >= p)
    warning("Upper bound for PFER larger than the number of effects.")
  
  res <- list(cutoff = cutoff, q = q, PFER = upperbound,
              specifiedPFER = PFER, p = p,
              B = B, sampling.type = sampling.type, assumption = assumption)
  class(res) <- "stabsel_parameters"
  res
}



# run_stabsel -------------------------------------------------------------

### the actual stability selection function (which is usually called by the
### generic stabsel function)
run_stabsel <- function(fitter, args.fitter, 
                        n, p, cutoff, q, PFER, folds, B, assumption,
                        sampling.type, papply, verbose, FWER, eval, names,
                        mc.preschedule = FALSE, ...) {
  
  folds <- check_folds(folds, B = B, n = n, sampling.type = sampling.type)
  # pars <- res
  pars <- stabsel_parameters(p = p, cutoff = cutoff, q = q,
                             PFER = PFER, B = B,
                             verbose = verbose, sampling.type = sampling.type,
                             assumption = assumption, FWER = FWER)
  cutoff <- pars$cutoff
  q <- pars$q
  PFER <- pars$PFER
  
  ## return parameter combination only if eval == FALSE
  if (!eval)
    return(pars)
  
  ## fit model on subsamples;
  ## Depending on papply, this is done sequentially or in parallel
  
  ## if mclappy is used consider mc.preschedule
  if (all.equal(papply, mclapply, check.environment = FALSE) == TRUE) {
    res <- suppressWarnings(
      papply(1:ncol(folds), function(...) try(fitter(...), silent = TRUE), folds = folds, q = q,
             args.fitfun = args.fitter, mc.preschedule = mc.preschedule, ...))
  } else {
    res <- papply(1:ncol(folds), function(...) try(fitter(...), silent = TRUE), folds = folds, q = q,
                  args.fitfun = args.fitter, ...)
  }
  
  ## if any errors occured remove results and issue a warning
  if (any(idx <- sapply(res, is.character))) {
    warning(sum(idx), " fold(s) encountered an error. ",
            "Results are based on ", ncol(folds) - sum(idx),
            " folds only.\n",
            "Original error message(s):\n",
            sapply(res[idx], function(x) x))
    res[idx] <- NULL
  }
  
  ## check results
  if (!is.list(res[[1]]) || !"selected" %in% names(res[[1]]))
    stop(sQuote("fitfun"), " must return a list with a named element called ",
         sQuote("selected"), ", and an optional second element called ", sQuote("path"))
  
  
  phat <- NULL
  if (!is.null(res[[1]]$path)) {
    ## extract selection paths
    paths <- lapply(res, function(x) x$path)
    # make path-matrices comparable
    steps <- sapply(paths, ncol)
    maxsteps <- max(steps)
    nms <- colnames(paths[[which.max(steps)]])
    paths <- lapply(paths, function(x) {
      if (ncol(x) < maxsteps) {
        x <- cbind(x, x[, rep(ncol(x), maxsteps - ncol(x))])
      }
      return(x)
    })
    phat <- paths[[1]]
    for (i in 2:length(paths))
      phat <- phat + paths[[i]]
    phat <- phat/length(paths)
    colnames(phat) <- nms
    rownames(phat) <- names
    
  }
  
  ## extract violations (only needed for boosting models)
  violations <- sapply(res, function(x)
    !is.null(attr(x, "violations")) && attr(x, "violations"))
  
  # AS: extract coefficients
  coefs <- lapply(res, function(x) x$coefs)
  
  ## extract selected variables
  res <- lapply(res, function(x) x$selected)
  res <- matrix(nrow = length(res), byrow = TRUE,
                unlist(res))
  colnames(res) <- names
  
  # AS: add coefficients
  ret <- list(phat = phat,
              selected = which(colMeans(res) >= cutoff),
              max = colMeans(res),
              coefs = coefs)
  ret <- c(ret, pars)
  
  ## return violations as attribute
  if (any(violations))
    attr(ret, "violations") <- violations
  
  class(ret) <- "stabsel"
  ret
}


## function to change PFER, cutoff or the assumption for a given stabsel object
stabsel.stabsel <- function(x, cutoff, PFER, assumption = x$assumption, ...) {
  
  assumption <- match.arg(assumption,
                          choices = c("unimodal", "r-concave", "none"))
  
  if (x$sampling.type == "MB" && assumption != "none")
    warning(sQuote('sampling.type == "MB"'), " but ",
            sQuote('assumption != "none"'))
  
  if (sum(missing(cutoff), missing(PFER)) == 0)
    stop("Only one of the two argumnets ",
         sQuote("PFER"), " and ", sQuote("cutoff"),
         " can be specifed")
  
  ## if nothing changed: nothing to do
  if (assumption == x$assumption) {
    if (sum(missing(cutoff), missing(PFER)) == 2)
      return(x)
    if (!missing(cutoff) && x$cutoff == cutoff)
      return(x)
    if (!missing(PFER) && x$PFER == PFER)
      return(x)
  }
  else {
    if (sum(missing(cutoff), missing(PFER)) == 2) {
      ## If both parameters PFER and cutoff were originally specified stop!
      if (!is.null(x$call[["PFER"]]) && !is.null(x$call[["cutoff"]]))
        stop("Specify one of ", sQuote("PFER"), " and ",
             sQuote("cutoff"))
      
      ## If originally only one of PFER and cutoff was specified use this
      ## parameter
      if (is.null(x$call[["PFER"]])) {
        cutoff <- x$cutoff
      }
      if (is.null(x$call[["cutoff"]])) {
        PFER <- x$specifiedPFER
      }
    }
  }
  if (!missing(cutoff)) {
    x$call[["cutoff"]] <- cutoff
    x$call[["q"]] <- x$q
    if (!is.null(x$call[["PFER"]]))
      x$call[["PFER"]] <- NULL
    x$specifiedPFER <- NULL
  }
  if (!missing(PFER)) {
    x$call[["PFER"]] <- PFER
    x$call[["q"]] <- x$q
    if (!is.null(x$call[["cutoff"]]))
      x$call[["cutoff"]] <- NULL
    x$specifiedPFER <- PFER
  }
  if (x$assumption != assumption)
    x$call[["assumption"]] <- assumption
  
  pars <- stabsel_parameters(p = x$p, cutoff = cutoff, q = x$q, PFER = PFER,
                             B = x$B, assumption = assumption,
                             sampling.type = x$sampling.type,
                             verbose = FALSE)
  
  cutoff <- pars$cutoff
  PFER <- pars$PFER
  
  ### now change results (by using new cutoff)
  x$selected <- which(x$max >= pars$cutoff)
  x$cutoff <- pars$cutoff
  x$PFER <- pars$PFER
  x$assumption <- assumption
  return(x)
}