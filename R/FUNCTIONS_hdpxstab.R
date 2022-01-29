# Method hdpxstab ---------------------------------------------------------

# Define method
hdpxstab <- function(x, ...) {
  UseMethod("hdpxstab", x)
}

# Function for prediction
run_hdpxstab <- function(fitfun = glinternet.lasso, 
                         args.fitfun, x, y,
                         p, cutoff, q, folds, B, 
                         parallel = ncol(folds)>10, 
                         coefficients.main = F,
                         ...){
  
  # Define variables
  n <- nrow(x)
  if(is.null(colnames(x))) {colnames(x) <- paste0("V", 1:ncol(x))}
  
  # Fit model function (from stabs package)
  fit_model <- function(i, folds, x, y, q, args.fitfun, coefficients.main, ...) {
    inbag <- as.logical(folds[, i])
    do.call(fitfun, c(list(x = x[inbag, ], y = y[inbag, ], q = q,
                           coefficients.main = coefficients.main),
                      list(...),
                      list(args.fitfun)))
  }
  
  # Fit function in parallel with future
  library("future.apply")
  
  if (isFALSE(parallel)){
    fits <- suppressWarnings(
      lapply(1:ncol(folds), try(fit_model, silent = T), 
             folds = folds, q = q, x = x, y = y, 
             coefficients.main = coefficients.main,
             args.fitfun = args.fitfun, ...)
    )
  } else {
    message("Preparing multisession for parallel processing.")
    plan(multisession)
    fits <- suppressWarnings(
      future_lapply(1:ncol(folds), try(fit_model, silent = T), 
                    folds = folds, q = q, x = x, y = y, 
                    args.fitfun = args.fitfun, 
                    coefficients.main = coefficients.main, ...)
    )
    plan(sequential) 
  }
  
  ## Issue a warning if errors
  if (any(idx <- sapply(fits, is.character))) {
    warning(sum(idx), " fold(s) encountered an error. ",
            "Results are based on ", ncol(folds) - sum(idx),
            " folds only.\n",
            "Original error message(s):\n",
            sapply(res[idx], function(x) x))
    fits[idx] <- NULL
  }
  
  ## Check results
  if (!is.list(fits[[1]]) || !"selected" %in% names(fits[[1]]))
    stop(sQuote("fitfun"), " must return a list with a named element called ",
         sQuote("selected"), ", and an optional second element called ", sQuote("path"))
  
  ## Extract paths
  phat <- NULL
  if (!is.null(fits[[1]]$path)) {
    paths <- lapply(fits, function(x) x$path)
    # Give all paths the same number of cols
    steps <- sapply(paths, ncol)
    maxsteps <- max(steps)
    nms <- colnames(paths[[which.max(steps)]])
    paths <- lapply(paths, function(x) {
      if (ncol(x) < maxsteps) {
        x <- cbind(x, x[, rep(ncol(x), maxsteps - ncol(x))])
      }
      return(x)
    })
    
    # My version: select the path with all possible selected rows (maxrows)
    
    # Get all rownames
    rownames.all <- sapply(paths, rownames)
    rownames.all <- Reduce(union,  rownames.all) 
    newphat <- matrix(NA, length(rownames.all), ncol(paths[[which.max(steps)]]),
                      dimnames=list(rownames.all))
    newphat <- ifelse(is.na(newphat), 0, newphat)
    phat <- newphat
    
    # Fill this matrix with the other matrices
    for (i in c(1:length(paths))){
      paths[[i]] <- ifelse(is.na(paths[[i]]), 0, paths[[i]])
      newphat[match(rownames(paths[[i]]), rownames(newphat)), ] <- paths[[i]]
      phat <- phat + newphat
    }
    
    phat <- phat/length(paths)
    phat <- phat[,2:ncol(phat)]
    # colnames(phat) <- nms
  }
  
  # Extract coefficients
  coefs <- lapply(fits, function(x) x$coefs)
  
  # Extract selected variables
  res <- lapply(fits, function(x) x$selected)
  res <- as.matrix(dplyr::bind_rows(res))
  res <- ifelse(is.na(res), FALSE, res)
  
  # AS: add coefficients
  ret <- list(phat = phat,
              selected = which(colMeans(res) >= cutoff),
              max = colMeans(res),
              coefs = coefs)
  pars <- list("cutoff" = cutoff, "q" = q, "n" = n)
  ret <- c(ret, pars)
  
  class(ret) <- "hdpxstab"
  ret
}



# Summary function for method
summary.hdpxstab <- function(x, decreasing = FALSE, print.all = TRUE, ...) {
  
  cat("Stability Selection")
  if (length(x$selected) > 0) {
    cat("\nSelected variables:\n")
    print(x$selected)
  } else {
    cat("\nNo variables selected\n")
  }
  cat("\nAverage selection probabilities:\n")
  if (print.all) {
    print(sort(x$max, decreasing = decreasing))
  } else {
    print(sort(x$max[x$max > 0], decreasing = decreasing))
  }
  
  cat("\n---\n")
  cat("\n")
  invisible(x)
}


# Plot function for method
plot.hdpxstab <- function(x, main = deparse(x$call), type = c("maxsel", "paths"),
                          maxprobsel = 0.5, ggplot = TRUE, xlab = NULL, ylab = NULL,
                          title = NULL, col = NULL, ymargin = 10, np = sum(x$max > 0),
                          labels = NULL, ...) {
  
  type <- match.arg(type)
  
  if(isFALSE(ggplot)){
    if (is.null(col))
      col <- hcl(h = 40, l = 50, c = x$max / max(x$max) * 490)
    if (type == "paths" && is.null(x$phat)) {
      warning("Stability paths ", sQuote("x$phat"), " are missing, ",
              "plot maximum selection frequency instead")
      type <- "maxsel"
    }
    if (type == "paths") {
      ## if par(mar) not set by user ahead of plotting
      if (all(par()[["mar"]] == c(5, 4, 4, 2) + 0.1))
        ..old.par <- par(mar = c(5, 4, 4, ymargin) + 0.1)
      h <- x$phat    
      h <- h[rowSums(h) > 0, , drop = FALSE]
      if (is.null(xlab)) {
        xlab <- "Step"
      }
      if (is.null(ylab)) {
        ylab <- "Selection probability"
      }
      matplot(t(h), type = "l", lty = 1,
              xlab = xlab, ylab = ylab,
              main = main, col = col[x$max > 0], ylim = c(0, 1), ...)
      abline(h = x$cutoff, lty = 1, col = "lightgray")
      if (is.null(labels))
        labels <- rownames(x$phat)
      axis(4, at = x$phat[rowSums(x$phat) > 0, ncol(x$phat)],
           labels = labels[rowSums(x$phat) > 0], las = 1, ...)
    } else {
      ## if par(mar) not set by user ahead of plotting
      if (all(par()[["mar"]] == c(5, 4, 4, 2) + 0.1))
        ..old.par <- par(mar = c(5, ymargin, 4, 2) + 0.1)
      if (np > length(x$max))
        stop(sQuote("np"), "is set too large")
      inc_freq <- x$max  ## inclusion frequency
      if (is.null(xlab)) {
        xlab <- expression(hat(pi))
      }
      if (is.null(ylab)) {
        ylab <- ""
      }
      plot(tail(sort(inc_freq), np), 1:np,
           type = "n", yaxt = "n", xlim = c(0, 1),
           ylab = ylab, xlab = xlab,
           main = main, ...)
      abline(h = 1:np, lty = "dotted", col = "grey")
      points(tail(sort(inc_freq), np), 1:np, pch = 19,
             col = col[tail(order(inc_freq), np)])
      if (is.null(labels))
        labels <- names(x$max)
      axis(2, at = 1:np, labels[tail(order(inc_freq), np)], las = 2, ...)
      ## add cutoff
      abline(v = x$cutoff, col = "grey")
    }
    if (exists("..old.par"))
      par(..old.par) # reset plotting settings
  }
  
  if (isTRUE(ggplot)){
    if (type == "paths") {
      ncol.phat <- ncol(x$phat)
      
      x$phat <- cbind(0, x$phat)
      colnames(x$phat) <- paste0("s", 0:ncol.phat)
      
      df <- as.data.frame(x$phat) %>% 
        rownames_to_column(var = "variable") %>% 
        pivot_longer(-variable, names_to = "path", values_to = "probability") %>% 
        mutate(path = factor(path, levels = colnames(x$phat)))
      
      graph <- ggplot(df, aes(x = path, y = probability, group = variable, color = variable)) 
      
      graph <- graph +
        geom_line(size = 1) + 
        geom_hline(yintercept = 0.75, size = 1, color = "gray") +
        scale_colour_hue(h = c(0, 180), l = 60, guide = 'none') +
        # scale_color_hue(l=40, c=35, guide = 'none') +
        # scale_color_brewer(palette="OrRd", guide = 'none') +
        scale_x_discrete(expand= expansion(mul = c(0, 0.25))) +
        coord_cartesian(ylim = c(maxprobsel, 1.05)) +
        geom_dl(aes(label = variable), method = list(dl.trans(x = x + 0.2), 
                                                     "last.bumpup", cex = 0.8)) +
        labs(x = ifelse(is.null(xlab), "Step", xlab),
             y = ifelse(is.null(ylab), "Selection probability", ylab),
             title = ifelse(is.null(title), "Stable selection path", title)) +
        theme_tufte() +
        theme(panel.grid.major.x = element_line(color = "gray90", size = 0.1),
              panel.grid.major.y = element_blank(),
              # plot.margin = ggplot2::margin(2, 1, 1, 1, "cm"),
              axis.line = element_line(colour = "black"))
    }
    
    if (type == "maxsel") {
      df <- x$max
      df <- sort(df, decreasing = F)
      df <- df[df > 0.45]
      
      df <- as.data.frame(df)
      colnames(df) <- "value"
      df$variable <- rownames(df)
      df$variable <- factor(df$variable, levels = df$variable)
      
      graph <- ggplot(df, aes(x = value, y = variable)) 
      graph <- graph +
        geom_point(color = "red") +
        scale_color_discrete(name="Coefficient: ") +
        labs(x = NULL, 
             y = NULL,
             title = "Selection probability") +
        theme_tufte() +
        theme(panel.grid.major.x = element_line(color = "gray90", size = 0.1),
              panel.grid.major.y = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.spacing = unit(2, "lines"),
              legend.position = "bottom")
    }
    
    graph
  }
}



# Internal functions ------------------------------------------------------

#' glinternet.lasso
#' Wrapper around glinternet from package `glinternet`
#' @keywords internal
glinternet.lasso <- function(x, y, q, args.fitfun, coefficients.main = F, ...) {
  
  if (!requireNamespace("glinternet", quietly=TRUE))
    stop("Package ", sQuote("glinternet"), " needed but not available")
  
  if (is.null(args.fitfun$numLevels))
    stop("numLevels must be specified in", sQuote("args.fitfun"),". ",
         "See the vignette of ", sQuote("glinternet"), "for more details.")
  
  numLevels <- args.fitfun$numLevels
  
  # if(isTRUE(coefficients.main)){
  #   message("Coefficients of main effects are also extracted for this model.")
  # } else {
  #   message("Only coefficients for selected predictors are extracted for this model.")
  # }
  
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"),
            " is coerced to a model matrix without intercept")
    x <- model.matrix(~ . - 1, x)
  }
  
  
  # ## fit model
  # fit <- glinternet::glinternet(x, y, numToFind = q, 
  #                               numLevels = numLevels, ...)
  # 
  # ## which coefficients are non-zero?
  # nfeat <- length(fit$activeSet)
  # 
  
  ### New version with q 
  ## fit model
  fit <- glinternet::glinternet(x, y, numLevels = numLevels, ...)
  
  ## which coefficients are non-zero?
  length.q <- lapply(fit$activeSet, function(x) sum(lengths(x)))
  nfeat <- length(which(length.q < q + 1))
  
  
  
  which.cat <- which(numLevels>1)
  which.cont <- which(numLevels==1)
  
  # Replace function
  replace.ind <- function(mat, which.numLevels) {
    if(is.matrix(mat)) {
      mat <- apply(mat,  c(1,2), function(x) which.numLevels[x])
    } else mat <- NULL
    mat
  }
  
  # Function to find the index for main variables and interactions
  index.fct <- function(fit.glinternet, which.numLevels.cat, 
                        which.numLevels.cont, x = x,
                        main = TRUE, lambda = NULL){
    
    if(is.null(lambda)) lambda <- length(fit.glinternet$activeSet)
    
    if(isTRUE(main)){
      cat <- which.numLevels.cat[fit.glinternet$activeSet[[lambda]]$cat]
      cont <- which.numLevels.cont[fit.glinternet$activeSet[[lambda]]$cont]
      
      ind <- c(cat, cont)
      
      cf <- coef(fit.glinternet, lambda)
      coefs <- unlist(c(lapply(cf[[1]]$mainEffectsCoef$cat[c(0:length(cat))], # the first n coefficients are the selected interactions, the rest are the main effects
                               function(l) l[[2]]), # take [[2]] for dummy == 1
                        cf[[1]]$mainEffectsCoef$cont[c(0:length(cont))]))
      
      coefs.names <- c()
      
      if(length(ind)>0){
        coefs.names <- c(colnames(x)[ind])
      } 
    }
    
    if(isFALSE(main)){
      ind <- 
        rbind(replace.ind(fit.glinternet$activeSet[[lambda]]$catcat, which.numLevels.cat),
              replace.ind(fit.glinternet$activeSet[[lambda]]$contcont, which.numLevels.cont),
              cbind(which.numLevels.cat[fit.glinternet$activeSet[[lambda]]$catcont[,1]],
                    which.numLevels.cont[fit.glinternet$activeSet[[lambda]]$catcont[,2]]))
      
      cf <- coef(fit.glinternet, lambda)
      coefs <- unlist(c(lapply(cf[[1]]$interactionsCoef$catcat, function(l) l[[1]]),
                        cf[[1]]$interactionsCoef$contcont,
                        lapply(cf[[1]]$interactionsCoef$catcont, function(l) l[[1]])))
      
      coefs.names <- c()
      
      if(nrow(ind)>0){
        for (kkk in 1:nrow(ind)){coefs.names <- c(coefs.names, 
                                                  paste0(colnames(x)[ind[kkk,1]], ":",
                                                         colnames(x)[ind[kkk,2]]))}  
      } 
    }
    
    list("ind" = ind,
         "coefs" = coefs,
         "coefs.names" = coefs.names)
    
  }
  
  main <- index.fct(fit, which.numLevels.cont = which.cont, x =x,
                    which.numLevels.cat = which.cat, main = TRUE, lambda = nfeat)
  int <- index.fct(fit, which.numLevels.cont = which.cont, x = x,
                   which.numLevels.cat = which.cat, main = FALSE, lambda = nfeat)
  
  
  # Create a dataset with the additional interactions
  if(nrow(int$ind)>0){
    names.var <- c(main$coefs.names, int$coefs.names)
  } else {
    names.var <- main$coefs.names
  }
  
  ret <- logical(length(names.var))
  ret[] <- TRUE
  names(ret) <- names.var
  
  
  ## Compute selection paths
  # For coefficient: coef is different than activeSet because coefficients for 
  # main effects are included.
  
  ## ACTIVESET
  # Coefs for the whole lambda path
  if(isFALSE(coefficients.main)){
    coefs <- unlist(c(main$coefs, int$coefs))
    coefs <- cbind(coefs)
    rownames(coefs) <- c(main$coefs.names, int$coefs.names)
    
    coefs.full <- coefs
    for (iii in 2:nfeat){
      cfi.main <- index.fct(fit, which.numLevels.cont = which.cont, 
                            which.numLevels.cat = which.cat, main = TRUE, 
                            x=x, lambda = iii)
      cfi.int <- index.fct(fit, which.numLevels.cont = which.cont, 
                           which.numLevels.cat = which.cat, main = F, 
                           x=x, lambda = iii)
      cfi <- unlist(c(cfi.main$coefs, cfi.int$coefs))
      cfi <- cbind(cfi)
      rownames(cfi) <- c(cfi.main$coefs.names, cfi.int$coefs.names)
      
      coefs.full <- cbind(coefs.full, cfi[match(rownames(coefs.full), rownames(cfi))])
    }
    
    sequence <- as.matrix(coefs.full != 0)
  }
  
  
  ## COEFFICIENT
  ### Start with the last element to have the full matrix
  if(isTRUE(coefficients.main)){
    
    # Function to extract coefficients and names from coef.glinternet
    index.fct.coef <- function(fit.glinternet, which.numLevels.cat,
                               which.numLevels.cont, x = x,
                               main = TRUE, lambda = NULL){
      
      if(is.null(lambda)) lambda <- length(fit.glinternet$activeSet)
      
      cf <- coef(fit.glinternet, lambdaIndex = lambda)
      
      # Coefs for nfeat
      if(isTRUE(main)){
        res <- list("var.coefs" = c(unlist(lapply(cf[[1]]$mainEffectsCoef$cat, function(l) l[[1]])),
                                    unlist(cf[[1]]$mainEffectsCoef$cont)),
                    "var.names" = c(colnames(x[, which.numLevels.cat[unlist(cf[[1]]$mainEffects$cat)], drop = F]),
                                    colnames(x[, which.numLevels.cont[unlist(cf[[1]]$mainEffects$cont)], drop = F])))
      }
      
      if(isFALSE(main)){
        cf.int.names <- rbind(replace.ind(unlist(cf[[1]]$interactions$catcat), which.numLevels.cat),
                              replace.ind(unlist(cf[[1]]$interactions$contcont), which.numLevels.cont),
                              cbind(replace.ind(unlist(cf[[1]]$interactions$catcont[,1, drop = F]), which.numLevels.cat),
                                    replace.ind(unlist(cf[[1]]$interactions$catcont[,2, drop = F]), which.numLevels.cont)))
        int.names <- c()
        
        if(!is.null(cf.int.names)){
          for (kkk in 1:nrow(cf.int.names)){int.names <- c(int.names, paste0(colnames(x)[cf.int.names[kkk,1]], ":",
                                                                             colnames(x)[cf.int.names[kkk,2]]))}
        } else {int.names <- NULL}      
        
        res <- list("var.coefs" = c(unlist(lapply(cf[[1]]$interactionsCoef$catcat, function(l) l[[1]])),
                                    unlist(lapply(cf[[1]]$interactionsCoef$contcont, function(l) l[[1]])),
                                    unlist(lapply(cf[[1]]$interactionsCoef$catcont, function(l) l[[1]]))),
                    "var.names" = int.names)
        rm(int.names, cf.int.names)
      }
      
      res
    }
    
    cf.main <- index.fct.coef(fit, which.numLevels.cat = which.cat,
                              which.numLevels.cont = which.cont, 
                              x = x,
                              main = TRUE, lambda = nfeat)
    cf.int <- index.fct.coef(fit, which.numLevels.cat = which.cat,
                             which.numLevels.cont = which.cont, 
                             x = x,
                             main = F, lambda = nfeat)
    
    # Coefs for the whole lambda path
    coefs <- unlist(c(cf.main$var.coefs, cf.int$var.coefs))
    coefs <- cbind(coefs)
    rownames(coefs) <- c(cf.main$var.names, cf.int$var.names)
    
    coefs.full <- coefs
    for (iii in 2:nfeat){
      cfi.main <- index.fct.coef(fit, which.numLevels.cat = which.cat,
                                 which.numLevels.cont = which.cont, 
                                 x = x,
                                 main = TRUE, lambda = iii)
      cfi.int <- index.fct.coef(fit, which.numLevels.cat = which.cat,
                                which.numLevels.cont = which.cont, 
                                x = x,
                                main = F, lambda = iii)
      cfi <- unlist(c(cfi.main$var.coefs, cfi.int$var.coefs))
      cfi <- cbind(cfi)
      rownames(cfi) <- c(cfi.main$var.names, cfi.int$var.names)
      
      coefs.full <- cbind(coefs.full, cfi[match(rownames(coefs.full), rownames(cfi))])
    }
    coefs.full <- coefs.full[,-1] 
    colnames(coefs.full) <- paste0("s", 1:(nfeat-1))
    rownames(coefs.full) <- c(cf.main$var.names, cf.int$var.names)
    
    sequence <- as.matrix(coefs.full != 0)
  }
  
  
  ## return both
  return(list(selected = ret, 
              path = sequence,
              coefs = coefs.full)) 
}


#' subsample
#' Prepares subsampling for the stable selection. Function taken from package 
#' "stabs".
# Subsample function from stabsel
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
  
  # Here, subsampling of 50% of observations within strata
  for (s in levels(strata)) {
    indx <- which(strata == s)
    folds[indx,] <- make_subsample(length(indx), B = B) * weights[indx]
  }
  
  folds
}