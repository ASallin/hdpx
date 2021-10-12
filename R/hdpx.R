
# Method hdpx ---------------------------------------------------------

# Define method
hdpx <- function(x, ...) {
  UseMethod("hdpx", x)
}


# Main function: hdpx
fit.hdpx <- function(object.hdpxtrans = NULL, folds, list.methods, args.glinternet,
                     parallel = TRUE, B = NULL, K = NULL, strata = NULL, 
                     x = NULL, y = NULL){
  
  # Preparation
  if(!inherits(object.hdpxtrans, "hdpxtrans")){
    if(is.null(x)) stop("Give a x matrix.")
    if(is.null(y)) stop("Give a y vector.")
  }
  
  x <- object.hdpxtrans$x
  y <- object.hdpxtrans$y
  
  if(is.null(folds)) {
    if (is.null(B) | is.null(K)) stop("folds is NULL. Provide K and B.")
    message(paste0("computing folds (B=", B, ", K = ", K, ")"))
    folds <- subsample.cv(weights = rep(1, nrow(x)), B = B,
                          strata = strata,
                          K = K)
  }
  
  K <- folds$K
  B <- folds$B
  folds.cv <- folds$folds
  
  # here, use do.call(cbind, folds.cv) to cbind all matrices from folds.cv
  if(length(which(apply(folds.cv,1,function(x) length(which(x==2)))==0)) > 0){
    stop("The folds.cv matrix has rows that are not in the test sample.")
  }
  
  # Option if folds is NULL
  fit_model <- function(i, folds, x, y, list.methods, 
                        args.glinternet = NULL, args.ranger = NULL){
    x.train.i <- x[which(folds[,i]==1), ]
    y.train.i <- y[which(folds[,i]==1)]
    
    x.test.i <- x[which(folds[,i]==2), ]
    y.test.i <- y[which(folds[,i]==2)]
    
    predictions <- list()
    
    for (mm in list.methods){
      predictions[[mm]] <- 
        do.call(wrapper_hdpx,
                c(list(x.train = x.train.i, x.test = x.test.i,
                       y.train = y.train.i, y.test = y.test.i,
                       method = mm, args.glinternet, args.ranger)))
    }
    
    predictions
  }
  
  library("future.apply")
  
  if (isFALSE(parallel)){
    fits <- suppressWarnings(
      # lapply(1:ncol(folds.cv), try(fit_model, silent = T),
      lapply(1:length(folds.cv), try(fit_model, silent = T),
             folds = folds.cv, x = x, y = y, 
             list.methods = list.methods, 
             args.glinternet = args.glinternet)
    )
  } else {
    message("Preparing multisession for parallel processing.")
    plan(multisession)
    fits <- suppressWarnings(
      future_lapply(1:ncol(folds.cv), try(fit_model, silent = T), 
                    folds = folds.cv, x = x, y = y,
                    list.methods = list.methods, 
                    args.glinternet = args.glinternet)
    )
    plan(sequential) 
  }
  
  # Create a matrix of y for each folds.cv==2
  y.folds <- apply(folds.cv, 2, function(x) y[x==2])
  
  # bind cols within each list
  fits <- lapply(fits, function(x) do.call(cbind, x)) 
  
  # Return a list of prediction (one list per row of "folds")
  names(fits) <- c(1:ncol(folds.cv))
  
  ret <- list(yhat = fits,
              folds = folds,
              y = y, 
              x = x,
              data = object.hdpxtrans$data,
              vars = object.hdpxtrans$vars)
  
  call <- list("methods" = list.methods)
  ret <- c(ret, call)
  
  class(ret) <- "hdpx"
  ret
}



# Predict fitted values from fitted object
predict.hdpx <- function(object, residuals = TRUE,
                         method = c("mean", "lowest.rmse", "superlearner"),
                         prediction.matrix = NULL, y = NULL, B = NULL, K = NULL,
                         methods = NULL, folds.m = NULL){
  
  if (inherits(object, "hdpx")){
    B <- object$folds$B
    K <- object$folds$K
    prediction.matrix <- object$yhat
    n.k <- ncol(prediction.matrix[[1]]) #here, I changed from length(object$methods)
    y <- object$y
    folds.m <- object$folds$folds
    y.folds <- apply(folds.m, 2, function(x) y[x==2]) # Create a list of y per col of folds
  } else {
    warning("Object is not of class", sQuote("hdpx"), ". Provide prediction matrix.")
    if(is.null(prediction.matrix)) stop("Provide prediction matrix.")
    if(is.null(y)) stop("Provide y matrix.")
    if(is.null(folds)) stop("Provide folds matrix.")
  }
  
  if(is.null(method)) method <- "mean"
  method <- match.arg(method)

  if(is.null(object$vars$level.demeaning) & isTRUE(residuals)){
    warning(sQuote("residuals"), "is TRUE but no demeaning. Setting ", sQuote("residuals"), " to FALSE.")
    residuals <- FALSE
  }
    
  # Compute RMSE for each folds
  rmse <- function(yhat, y){sqrt(mean((yhat-y)^2, na.rm = T)) }
  rmse.list <- mapply(prediction.matrix, y.folds, FUN = function(x,y) apply(x, 2, rmse, y = y),
                      SIMPLIFY = F) # CHECK
  
  # Variant 1: take the row mean
  if(method == "mean"){
    yhat <- lapply(prediction.matrix, function(x) apply(x, 1, mean, na.rm = T)) 
    yhat <- lapply(yhat, as.matrix)
  }
  
  # Variant 2: take the observation from the fold and method with the lowest RMSE
  if(method == "lowest.rmse"){
    
    lowest.rmse.fct <- function(rmse.list, prediction.matrix){
      rmse.vector <- unlist(rmse.list)
      
      index.min.rmse <- prediction.matrix != 0
      index.min.rmse <- t(t(index.min.rmse) * rmse.vector) # elementwise multiplication
      
      ind <- apply(index.min.rmse,1,function(x) which.min(x))
      ind[lengths(ind)==0] <- NA
      ind <- unlist(ind)
      
      yhat <- c()
      for(i in 1:nrow(prediction.matrix)){
        yhat[i] <- prediction.matrix[i, ind[i]]  ## CHECK THIS
      }
      yhat <- as.matrix(yhat)
      rm(index.min.rmse, ind)
      yhat
    }
    
    yhat <- pmap(list(rmse.list, prediction.matrix), lowest.rmse.fct)
    
  }
  
  # Variant 3: use a superlearner
  if(method == "superlearner"){
    
    ensembles <- mapply(function(x,y) ensemble(predictors = x, y = y, k = ncol(x)),
                        prediction.matrix, y.folds, SIMPLIFY = F)
    # str(ensembles)
    yhat <- lapply(ensembles, function(x) x[["weighted.yhat"]])
  }
  
  
  # Place all yhat together:
  # 1) repeat an empty matrix in a list with B elements
  yhat.final <- lapply(seq_len(B), function(x) matrix(NA, nrow = nrow(folds.m), ncol = 1))
  # 2) Sequence: list with B elements of K cols
  list.folds <- split(1:length(yhat), rep(1:B, each = K))
  
  for (jj in seq_len(B)){
    for (kk in list.folds[[jj]]){
      yhat.final[[jj]][which(folds.m[,kk]==2),] <- yhat[[kk]]
    }
  }
  
  if(length(which(lapply(yhat.final, nrow)==nrow(folds.m))) != B) stop("Problem with yhat.")
  if (length(which(lapply(yhat.final, function(x) which(is.na(x))) > 0)) != 0) stop("Wrong final matrix.")
  
  
  # Remove the demeaning:
  if(isFALSE(residuals)){
    mean.y <- cbind("outcome" = object$data[[paste0(object$vars$outcome)]],
                    "cell.demeaning" = object$vars$vector.demeaning)
    mean.y <- data.table(mean.y)
    mean.y <- mean.y[, mean.y := mean(outcome), by = cell.demeaning]$mean.y

    yhat.final <- lapply(yhat.final, function(x) x + mean.y)
  }
  
  
  ret <- list(yhat = yhat.final,
              y = y,
              folds = folds.m,
              method = method)
  ret
}






# Internal functions ------------------------------------------------------

# combine.hdpx
# Combine together hdpx objects from fit.hdpx in order to predict them all together.
combine.hdpx <- function(hdpx.1, hdpx.2){
  
  if (!inherits(hdpx.1, "hdpx") | !inherits(hdpx.2, "hdpx")) stop("Arguments must be hdpx objects.")
  
  if(!setequal(hdpx.1$y, hdpx.2$y)) stop("Arguments must share same y.")
  if(!setequal(hdpx.1$folds$folds, hdpx.2$folds$folds)) stop("Arguments must share same folds.")
  
  ret <- list(yhat = mapply(cbind, hdpx.1$yhat, hdpx.2$yhat, SIMPLIFY=FALSE),
              y = hdpx.1$y,
              folds = hdpx.1$folds,
              data = hdpx.1$data,
              vars = hdpx.1$vars)
  
  call <- list("methods" = c(hdpx.1$methods, hdpx.2$methods))
  ret <- c(ret, 
           call)
  class(ret) <- "hdpx"
  ret
}


# subsample.cv
# Sampling prediction function
subsample.cv <- function(weights, B = 1, strata = NULL, 
                         K = 10, sample.train.ratio = 0.8) {
  
  if (is.null(strata)) strata <- gl(1, n)
  if (!is.factor(strata))
    stop(sQuote("strata"), " must be a factor")
  
  
  # Subfunctions
  ## IMPROVE FOR SPEED: get rid of merge
  cv.sampling <- function(n, K, strata, sample.train.ratio) {
    
    
    if(length(levels(strata))>1){
      matrix.folds <- matrix(as.numeric(levels(strata)))
      matrix.folds <- cbind(matrix.folds, runif(length(levels(strata))))
      matrix.folds <- cbind(matrix.folds, as.numeric(cut(matrix.folds[,2], 
                                                         quantile(matrix.folds[,2], probs = seq(0, 1, 1/K)),
                                                         include.lowest = TRUE)))
      matrix.folds <- matrix.folds[, c(1,3)]
      colnames(matrix.folds) <- c("strata", "fold")
      
      matrix.strata <- matrix(as.numeric(strata))
      colnames(matrix.strata) <- "strata"
      
      matrix.cv <- merge(matrix.folds, matrix.strata, by = "strata")
      matrix.cv[, "fold"] <- as.factor(matrix.cv[, "fold"])
      matrix.cv <- model.matrix(~-1 + strata + fold, matrix.cv)  
      
      # 2 means test, 1 means train
      matrix.cv[, c(1:K+1)] <- matrix.cv[, c(1:K+1)] + 1
      dim(matrix.cv)
      
      # Randomly draw X observations in the train sample
      matrix.cv <- 
        apply(matrix.cv[, c(1:K+1)], 2, 
              function(x) {x[sample(which(x == 1), floor((1-sample.train.ratio)*sum(x == 1)))] <- 0
              return(x)})
      
    }
    
    if(length(levels(strata))==1){
      split   <- runif(length(weights))
      folds   <- as.numeric(cut(split, quantile(split, probs = seq(0, 1, 1/K)),
                                include.lowest = TRUE)) 
      rm(split)
      
      matrix.cv <- as.matrix(as.factor(folds))
      colnames(matrix.cv) <- "fold"
      matrix.cv <- model.matrix(~-1 + fold, as.data.frame(matrix.cv) )
      matrix.cv <- apply(matrix.cv, 2, as.numeric) + 1
      
    }
    
    colnames(matrix.cv) <- NULL
    matrix.cv
  }
  
  # Diagnostics (to cancel)
  # test <- cv.sampling(n, K = 10, strata = strata, sample.train.ratio=0.01)
  # count(as_tibble(test), fold1);  count(as_tibble(test), fold2)
  # summary(apply(test, 1, function(x) sum(x == 2)))
  
  # Settings
  n <- length(weights)
  folds <- replicate(n = B, cv.sampling(n, K = K, strata = strata, sample.train.ratio),
                     simplify = F)
  
  folds <- do.call(cbind, folds)
  folds <- apply(folds, 2, function(x) x * weights)
  
  return(list(folds = folds,
              K = K, 
              B=B))
}




#'wrapper_hdpx
#' This function is a wrapper around ML predicting functions.
#' 
#' @param x.train Matrix of covariates for the training dataset. The matrix must 
#' contain the treatment variable if the goal is to estimate the outcome nuisance
#' parameter. Name of treatment is given in "treatment".
#' @param x.test Matrix of covariates for the test dataset. The matrix must 
#' contain the treatment variable if the goal is to estimate the outcome nuisance
#' parameter. Name of treatment is given in "treatment".
#' @param y.train Vector of outcomes
#' @param y.test Vector of outcomes for the test dataset
#' 
#' @return A vector of predicted outcomes for the test dataset
#' 
#' @keywords internal

wrapper_hdpx <- function(x.train, x.test, y.train, y.test, 
                         method = c("ranger.rf", "glinternet.lasso"), 
                         args.glinternet = NULL, args.ranger = NULL,
                         ...){
  
  # Prepare
  method <- match.arg(method)
  
  x.train <- data.matrix(x.train)
  x.test <- data.matrix(x.test)
  y.train <- data.matrix(y.train)
  y.test <- data.matrix(y.test)
  
  # Arguments
  args.glinternet.lasso <- 
    c(list(nFolds = 10),
      args.glinternet)
  
  args.ranger.rf <- 
    c(list(num.trees = 1000,
           classification = FALSE,
           importance = "impurity", 
           verbose = T),
      args.ranger)
  
  
  # Train on training sample
  if(method == "glinternet.lasso") {fit <- glinternet_fit(x.train, y.train, args = c(args.glinternet.lasso))}
  if(method == "ranger.rf")   {fit <- ranger_fit(x.train, y.train, args = args.ranger.rf)}
  
  # Predict on test sample
  if(method == "glinternet.lasso") {y.hat.test <- predict.glinternet_fit(fit, x.train, xnew = x.test)}
  if(method == "ranger.rf")    {y.hat.test <- predict.ranger_fit(fit, x.train, xnew = x.test)}
  
  return(y.hat.test)
}



#' RF wrapper around ranger
#' This function estimates random forest lasso regression based on the \code{\link{ranger}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes (no df format)
#' @param args List of arguments passed to  \code{\link{ranger}}
#' @import ranger
#'
#' @return An object with S3 class "ranger"
#'
#' @keywords internal

ranger_fit = function(x, y, args=list()) {
  
  data.ranger <- as.data.frame(cbind(y, x)) 
  names(data.ranger)[1] <- "Y"
  
  message("Training random forest...")
  rf <- do.call(ranger,
                c(list(paste("Y ~ ."),
                       data = data.ranger), 
                  args))
  
  return(rf)
}



#' Prediction based on ranger
#' @param lasso_fit Output of \code{\link{ranger}} or \code{\link{ranger_fit}}
#' @param x Covariate matrix of training sample
#' @param xnew Covariate matrix of test sample
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#'
#' @keywords internal

predict.ranger_fit = function(ranger_fit, x, xnew = NULL) {
  
  if (is.null(xnew)) xnew <- x
  
  data.ranger <- as.data.frame(xnew)
  
  fit_ranger <- predict(ranger_fit,
                        data = data.ranger,
                        type = "response")$predictions
  
  return(as.matrix(fit_ranger))
  
} 


#' Wrapper around glinternet
#' This function estimates cross-validated group lasso regression based on the \code{\link{glinternet}} package
#'
#' @param x Matrix of covariates (number of observations times number of covariates matrix)
#' @param y vector of outcomes (no df format)
#' @param args List of arguments passed to  \code{\link{glinternet}}
#' @import biglasso
#'
#' @return An object with S3 class "glinternet.cv "
#'
#' @keywords internal

glinternet_fit = function(x, y, args=list()) {
  
  if(is.null(args$numLevels)){stop("numLevels must be provided.")}
  
  message(paste0("Training glinternet..."))
  lasso <- do.call(glinternet.cv,
                   c(list(X = x,
                          Y = y, 
                          verbose = F),
                     args))
  
  return(lasso)
}



#' Prediction based on glinternet
#' @param glinternet_fit Output of \code{\link{glinternet}}
#' @param x Covariate matrix of training sample
#' @param xnew Covariate matrix of test sample
#'
#' @return Returns list containing:
#' \item{prediction}{vector of predictions for xnew}
#'
#' @keywords internal

predict.glinternet_fit = function(glinternet_fit, x, xnew = NULL) {
  
  if (is.null(xnew)) xnew <- x
  
  message(paste0("Predicting glinternet on test set..."))
  
  fit <- predict(glinternet_fit, 
                 X = xnew, 
                 lambdaType = "lambdaHat")
  
  return(as.matrix(fit))
}


#' superlearner ensemble learner function
#' computes a grid search of all possible weights that minimizes the RMSE within
#' each K folds. Programming inspired by Chernozhukov et al. for implementation. 
#' @param k is the number of predictors that should be included in the learning 
#' computation. 

ensemble <- function(predictors, k, y){
  
  error <- function(yhat, y){
    err         <- sqrt(mean((yhat-y)^2, na.rm = T))
    mis         <- sum(abs(as.numeric(yhat > .5)-(as.numeric(y))))/length(y)   
    
    return(list(err = err, mis = mis));
  }
  
  # Lots of NaN with Ridge. Change them to 0 in the matrix.
  predictors <- round(predictors, 10) 
  # predictors[is.nan(predictors)] = Inf
  
  # If no colnames, give colnames
  colnames(predictors) <- paste0("method.", 1:ncol(predictors))
  
  # Compute MSE for each column of methods (for information)
  MSE.method <- apply(predictors, 2, error, y = y) %>% 
    map(., pluck("err")) %>% 
    do.call(rbind, .)
  MSE.method <- cbind(rownames(MSE.method), MSE.method) 
  rownames(MSE.method) <- NULL
  
  # Grid search
  # If there is NaN values, take k-number of NaN
  # if (which(is.nan(as.numeric(MSE.method[, 2]))) > 0){
  # k <- k - length(which(is.nan(as.numeric(MSE.method[, 2]))))
  # }
  
  # Take only the k methods with lowest RMSE
  predictors <- predictors[, order(MSE.method[, 2], decreasing = F)[1:k]]
  
  # Prepare grid search
  if(k < 4)  {lst     <- lapply(numeric(k), function(x) as.vector(seq(0, 1, 0.01))) }
  if(k == 4) {lst     <- lapply(numeric(k), function(x) as.vector(seq(0, 1, 0.04))) }
  if(k == 5) {lst     <- lapply(numeric(k), function(x) as.vector(seq(0, 1, 0.05))) }
  if(k > 5)  {lst     <- lapply(numeric(k), function(x) as.vector(seq(0, 1, 0.10))) }
  
  gr                  <- as.matrix(expand.grid(lst))
  weight              <- gr[rowSums(gr)==1,]
  rm(gr)
  
  # Weight predictors
  weighted.predictors <- predictors %*% t(weight)
  
  # Compute MSE for each column of the weighted predictions
  MSE <- apply(weighted.predictors, MARGIN = 2, error, y = y) %>% 
    map(., pluck("err")) 
  
  # Check which weight combination is the winning one
  RMSE.weights        <- weight[which.min(MSE), 1:k]
  names(RMSE.weights) <- colnames(predictors)
  
  # Select column with smallest MSE
  Yhat  <- weighted.predictors[, which.min(MSE)] %>% 
    as.matrix()
  
  # Return
  return(list("weighted.yhat" = Yhat,
              "RMSE.methods" = MSE.method,
              "RMSE.weights" = RMSE.weights))
}


