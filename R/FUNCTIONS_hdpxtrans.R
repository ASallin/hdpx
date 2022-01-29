# Define method
hdpxtrans <- function(x, ...) {
  UseMethod("hdpxtrans", x)
}

# Print
print.hdpxtrans <- function(x){
  cat("Data prepared for high dimensional peer-effects analysis")
  if (x$state == "types"){
    cat(" with all types.\n")
    cat("\nDimension of new dataset:")
    print(dim(x$pe))
  }
  if (x$state == "loo"){
    cat(" with all types and leave-one-out variables.\n")
    cat("\nDimension of new dataset:")
    print(dim(x$loo))
  }
  cat("\n---\n")
  cat("\n")
  invisible(x)
}


#' types.hdpxtrans
#' Takes a dataframe that contains peer effects and computes types and interactions
#' between types.
#' @param df a dataframe
#' @param types a character vector containing the variables that are used as types
#' @param level.interactions an integer giving the "depth" of interactions between 
#' types. If level.interactions = 1, no interactions are computed.
#' @param full.types is F by default. If TRUE, the most aggregated "depth" of 
#' interactions is computed.
#' @param remove.collinear removes highly collinear variables. 
types.hdpxtrans <- function(df, types, level.interactions = 1, 
                            full.types = FALSE, remove.collinear = FALSE){
  
  # Necessary packages loading
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(data.table, fastDummies)
  
  # Warnings
  if (is.null(types)) stop("Types must be given.")
  
  # Select which peer effects we want
  pe <- df[, names(df) %in% types] 
  
  # Create dummies (full rank!)
  if (level.interactions == 1){
    pe <- dummy_cols(pe, select_columns = c(names(pe)), 
                     remove_first_dummy = T,
                     remove_selected_columns = T)
  } else {
    pe <- dummy_cols(pe, select_columns = c(names(pe)), remove_first_dummy = F,
                     remove_selected_columns = F)
    pe <- pe[, !(colnames(pe) %in% c(grep("_1", colnames(pe), value=T)))]
  }
  
  # Create a matrix of interactions: either TYPES or interactions with a given 
  # number of levels
  if (isTRUE(full.types)) {level.interactions <- seq_len(length(types))}
  if (isFALSE(full.types)) {level.interactions <- seq_len(level.interactions)}
  
  if (length(level.interactions) > 1){
    names.interactions <- c()
    for(iii in level.interactions){
      names.interactions <- c(names.interactions, combn(colnames(pe), iii, FUN=paste, collapse=':'))
    }
  } else {
    names.interactions <- combn(colnames(pe), level.interactions, FUN=paste, collapse=':')
  }
  
  # Create columns
  pe <- transformations(pe, poly = 1, full = T)
  pe <- pe[, names.interactions]
  dim(pe)
  
  # Remove the dummy if it is always zero
  pe <- pe[, colSums(pe) > 0]
  
  # Remove the dummy if it is always one
  pe <- pe[, colSums(pe) < nrow(pe)]
  
  # Remove the main types that are collinear (female-male)
  library(corrplot); library(caret)
  # corrplot(cor(pe), tl.pos='n')
  
  if (isTRUE(remove.collinear)) {
    drops <- findLinearCombos(pe)
    if(!is.null(drops$remove)) pe <- pe[, -c(drops$remove)]
  }
  
  # Remove the type if it is smaller than X
  colnames(pe[, colSums(pe) < 0.001*nrow(pe)])
  pe <- pe[, colSums(pe) > 0.001*nrow(pe)]
  dim(pe)
  
  # Names.pe
  names.pe <- data.frame(names = colnames(pe),
                         pe = paste0("type_", 1:ncol(pe)))
  colnames(pe) <- paste0("type_", 1:ncol(pe))
  
  names.pe$Type <- str_remove(names.pe$pe, "type_")
  names.pe$names <- str_replace_all(names.pe$names, "SN.before", "SNbefore")
  names.pe$names <- str_replace_all(names.pe$names, ":", ".")
  
  ret <- list("pe" = pe,
              "names.pe" = names.pe,
              "data" = data.table(df),
              "state" = "types")
  class(ret) <- "hdpxtrans"
  
  ret
}


# Create LoO
loo.hdpxtrans <- function(types, interaction.types.loo = TRUE, FE_size = NULL,
                          interactions.FE_size = TRUE, peer.level){
  
  if (!inherits(types, "hdpxtrans")) stop("This function applies to object of class ",
                                       sQuote("hdpxtrans"), ".")
  
  pe <- data.table(cbind(types$data, types$pe))
  cols.type <- colnames(types$pe)[startsWith(colnames(types$pe), "type")]
  
  # Generate LoO
  LoO <- pe[, (paste0("LoO_", cols.type)):= lapply(.SD, function(x){(sum(x)-x)/(length(x)-1)}), 
            .SDcols = cols.type, 
            by = list(pe[[peer.level]])]
  
  LoO <- LoO[, grep("type|LoO", names(LoO), value=T), with=FALSE]
  # LoO[, lapply(.SD, mean), .SDcols =  grep("LoO", names(LoO), value=T)]
  
  
  # Interact the LoO with the types
  if (isTRUE(interaction.types.loo)){
    LoO <- as.data.frame(LoO)
    
    cols.LoO <- combn(colnames(LoO), 2, FUN=paste, collapse=':')
    cols.LoO <- cols.LoO[!(cols.LoO %in% 
                             c(grep("^type_\\d+:type_\\d+$", cols.LoO, value=T),
                               grep("^LoO_type_\\d+:LoO_type_\\d+$", cols.LoO, value=T)))]
    
    LoO.inter <- hdinteractions(names.cols = cols.LoO, df = LoO, K = 7)
    LoO <- data.table(cbind(LoO, LoO.inter))
  }
  
  
  # Interact all LoO and interactions with variables of cell size
  ### ! To extend to other potential variables
  if (!is.null(FE_size) & isTRUE(interactions.FE_size)){
    LoO_other <- LoO[, grep("^LoO", colnames(LoO), value=T), with = FALSE] 
    LoO_other <- bind_cols("FE_size" = pe[[FE_size]], LoO_other)
    
    LoO_other <- data.table(model.matrix(~ -1 + FE_size *., data = LoO_other) )
    LoO_other <- LoO_other[, c(grep("^LoO", colnames(LoO), value=T)) := NULL] 
    LoO <- cbind(LoO, LoO_other)
  }
  
  if (!is.null(FE_size) & isFALSE(interactions.FE_size)){
    LoO <- bind_cols(LoO, "FE_size" = pe[[FE_size]])
  }
  
  ret <- list("loo" = LoO,
              "names.pe" = types$names.pe,
              "data" = types$data,
              "state" = "loo")
  class(ret) <- "hdpxtrans"
  
  ret 
}


# Finalize LoO and types into final dataset and demean
final.hdpxtrans <- function(loo, outcome, demean = c("all", "cont", "none"),
                            cell.demeaning = NULL, trend.factor = NULL, 
                            list.vars = NULL, list.vars.notdemeaned = NULL){
  
  demean <- match.arg(demean)
  
  if (!inherits(loo, "hdpxtrans")) stop("This function applies to object of class ",
                                          sQuote("hdpxtrans"), ".")
  
  if (loo$state == "types") message("Dataset does not contain leave-one-out variables: 
                                    prepare types only.")
  
  # Necessary packages loading
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(data.table, lfe)
  
  if (demean != "none" & is.null(cell.demeaning)) stop("Provide the demeaning level.")
  if (!is.character(outcome)) stop("The name of the outcome variable must be provided.")
  if (demean != "none" & !is.character(cell.demeaning)) stop("The name of the demeaning strata must be provided.")
  if (!is.null(list.vars) & !is.character(list.vars)) stop(sQuote("list.vars"), 
                                                           "must be a character vector.")
  if (!is.null(list.vars.notdemeaned) & !is.character(list.vars.notdemeaned)) 
    stop(sQuote("list.vars"), "must be a character vector.")
  if(!is.null(trend.factor) & !is.character(trend.factor)) stop(sQuote("trend.factor"),
                                                                "must be a character vector.")
    
  
  # Check variable names
  if(!is.null(list.vars)) list.vars <- make.names(list.vars)
  if(!is.null(list.vars.notdemeaned)) list.vars.notdemeaned < make.names(list.vars.notdemeaned)
  colnames(loo$data) <- make.names(colnames(loo$data))
  colnames(loo$loo) <- make.names(colnames(loo$loo))
  
  # Work with df
  df <- cbind(loo$data, loo$loo)
  df <- data.table(df)
  
  if(demean == "none") cols.select <- c(outcome, colnames(loo$loo))
  if(demean != "none") cols.select <- c(outcome, cell.demeaning, colnames(loo$loo))
  if(!is.null(list.vars)) cols.select <- c(cols.select, list.vars)
  if(!is.null(trend.factor)) cols.select <- c(cols.select, trend.factor)
  
  df <- df[, ..cols.select]
  setnames(df, old = outcome, new = "outcome")
  if(demean != "none") setnames(df, old = cell.demeaning, new = "cell.demeaning")
  if(!is.null(trend.factor)) setnames(df, old = trend.factor, new = "trend.factor")
  
  # Prepare variables to be demeaned
  if(demean == "cont") names_vars <- colnames(df)[!grepl("^type_[0-9]+$", colnames(df))]
  if(demean == "all") names_vars <- colnames(df)
  if(demean == "none")  names_vars <- NULL
  
  names_vars <- names_vars[!names_vars %in% "cell.demeaning"]
  names_vars <- names_vars[!names_vars %in% "trend.factor"]
  
  if (!is.null(list.vars.notdemeaned)) {
    names_vars <- names_vars[!names_vars %in% list.vars.notdemeaned]
  }
  
  # Demean
  if(demean != "none"){
    # Use lfe model for demeaning
    if(!is.null(trend.factor)) {
      list_factor <- list(as.factor(df[["cell.demeaning"]]),
                          as.factor(df[["trend.factor"]]))
    } else {
      list_factor <- list(as.factor(df[["cell.demeaning"]]))  
      }
    
    df.demeaned <- lfe::demeanlist(df[, ..names_vars],
                                   fl = list_factor)
      
  }
  
  if(demean == "cont") df <- cbind(df[, !..names_vars], 
                                   df.demeaned)
  if(demean == "all") df <- df.demeaned
  
  if (!is.null(list.vars.notdemeaned)) {
    df <- cbind(df, cbind(loo$loo, loo$data)[, ..list.vars.notdemeaned])
  }
  
  df <- df[, .SD, .SDcols = unique(names(df))]
  
  # Prepare matrices for estimation
  if(demean == "none") x <- as.matrix(df[, !c("outcome")]) 
  if(demean != "none") x <- as.matrix(df[, !c("outcome", "cell.demeaning")]) 
  if(demean != "none") {vector.demeaning <- df[, cell.demeaning] } else {vector.demeaning <- NULL}
  y <- df[, outcome] 
  
  ret <- list("finaldf" = df,
              "names.pe" = loo$names.pe,
              "data" = loo$data,
              "state" = demean,
              "x" = x, 
              "y" = y,
              vars = list("vector.demeaning" = vector.demeaning,
                          "level.demeaning" = cell.demeaning,
                          "outcome" = outcome))
  
  class(ret) <- "hdpxtrans"
  
  ret 
}
  
