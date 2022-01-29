# Function for fractionalization
fractionalize <- function(group.vars, data, leave_one_out = F){
  
  # Creation of the types variables as a factor variables
  if (leave_one_out == T){
    
    number.types <- length(group.vars)
    
    for (i in 1:number.types){
      data[, paste0("LoO_", i)] <- data[, group.vars[i]] 
    }
    }
  
  data <- data %>% 
    group_by_at(vars(all_of(group.vars))) %>% 
    mutate(types = cur_group_id()) %>% 
    ungroup()
    
  # Create dummies from the type variables
  data <- dummy_cols(data, select_columns = "types")
    
  number.types <- max(data$types)
  
  cat(paste0("The number of types is ", number.types, "."))
  
  
  # Check for missing types
  which(is.na(data$types))
  
  # Compute the total of students with type X per cohort 
  data <- data %>% 
    group_by(schooltrackyear) %>% 
    mutate(across(.cols = starts_with("types_"),
                  .fns = list(total = sum))) %>% 
    ungroup()
  
  
  # Create leave-one-out mean (with base R for simplicity)
  for (i in seq_len(number.types)) {
    data[, paste0("leaveoneoutTYPE_", i)] <- 
      (data[, paste0("types_", i, "_total")] - data[, paste0("types_", i)]) /
      (data[, "schooltrackyear_size"] -1)
  }
  
  
  # Create fractionalization score
  data <- data %>% 
    mutate(across(.cols = starts_with("leaveoneoutTYPE_"),
                  .fns = list("pw2" = ~.x^2),
                  .names = "{.col}_{.fn}")) %>% 
    mutate(total_s = rowSums(select(., ends_with("_pw2"))),
           FRAC = 1-total_s)
  
  data %>% 
    select(schooltrackyear, starts_with("types_"), 
           schooltrackyear_size, starts_with("leaveoneoutTYPE_"), FRAC) %>% 
    arrange(schooltrackyear)  
  
  # Standardize the fractionalization measure
  data <- data %>% 
    mutate(FRAC_std = FRAC/max(FRAC))
  # mutate(FRAC_std = (FRAC-mean(FRAC))/sd(FRAC))
  

  # Check Fractionalization measure
  summary(data$FRAC);  summary(data$FRAC_std)
  
  # Return
  return(list("Frac" = as.matrix(data$FRAC),
              "Frac_std" = as.matrix(data$FRAC_std),
              "Types" = as.matrix(data$types),
              "Data" = data))
}
# Function to create factoralization index --------------------------------

fractionalization <- function(fe_types, fe_level, data){
  
  # 1) Create group indices
  data <- data %>% 
    group_by_at(vars(all_of(fe_types))) %>% 
    mutate(types = cur_group_id()) %>% 
    ungroup()
  
  # 2) Create types indicators
  data <- dummy_cols(data, select_columns = "types")
  
  number.types <- max(data$types)
  
  cat(paste0("The number of types is ", number.types, "."))
  
  
  # 3) Check for missing types
  which(is.na(data$types))
  
  
  # 4) Compute the total of students with type X per cohort 
  data <- data %>% 
    group_by_at(vars(all_of(fe_level))) %>% 
    mutate(across(.cols = starts_with("types_"),
                  .fns = list(total = sum))) %>% 
    ungroup()
  
  
  # 5) Create leave-one-out mean (with base R for simplicity)
  for (i in seq_len(number.types)) {
    data[, paste0("leaveoneoutTYPE_", i)] <- 
      (data[, paste0("types_", i, "_total")] - data[, paste0("types_", i)]) /
      (data[, "schooltrackyear_size"] -1)
  }
  
  
  # 6) Create fractionalization score
  data <- data %>% 
    mutate(across(.cols = starts_with("leaveoneoutTYPE_"),
                  .fns = list("pw2" = ~.x^2),
                  .names = "{.col}_{.fn}")) %>% 
    mutate(total_s = rowSums(select(., ends_with("_pw2"))),
           FRAC = 1-total_s)
  
  # data %>% 
    # select(schooltrackyear, starts_with("types_"), 
           # schooltrackyear_size, starts_with("leaveoneoutTYPE_"), FRAC) %>% 
    # arrange(schooltrackyear)  
  
  
  # 7) Standardize the fractionalization measure
  data <- data %>% 
    mutate(FRAC_std = FRAC/max(FRAC))
  # mutate(FRAC_std = (FRAC-mean(FRAC))/sd(FRAC))
  
  
  # Check Fractionalization measure
  summary(data$FRAC);  summary(data$FRAC_std)
  
  # Return
  return(list("Frac" = as.matrix(data$FRAC),
              "Frac_std" = as.matrix(data$FRAC_std)))
}




# Function for transforming datasets  -------------------------------------

transformations <- function(dataset, poly, full = FALSE) {
  
  # Create vector of binaries
  vector_binary <- names(
    dataset[
      which(apply(dataset, 2, max, na.rm=T)==1)]
  )
  
  
  # Create polynomials
  if (poly == 1) {
    x <- dataset
  }
  
  else if (poly == 2) {
    x <- dataset %>% 
      mutate_at(vars(-vector_binary), funs(two = .^2))
  }
  
  else if (poly == 3) {
    x <- dataset %>% 
      mutate_at(vars(-vector_binary), funs(two = .^2)) %>% 
      mutate_at(vars(-vector_binary, -matches("two")), funs(three = .^3))
  }
  
  # Interact all variables
  if (isTRUE(full)){
    formula_interaction <- as.formula(paste0("~ .^", ncol(x)))
  } else {
    formula_interaction <- as.formula(~ .*.)
  }
  # previous_na_action <- options('na.action')
  # options(na.action = "na.pass")
  x <- model.matrix(formula_interaction, x)[, c(-1)]
  # options(na.action=previous_na_action$na.action)
  
  return(x)
}
# Extract variable importance for lasso -----------------------------------

biglasso.importance <- function(biglasso.object, X.train){
  
  coefficients <- coef(biglasso.object, s = "lambda.min") 
  
  coefficients <- as.matrix(coefficients, rownames(coefficients))
  
  importance <- as.matrix(coefficients[which(coefficients != 0),])
  
  
  # For ranking of variable importance, I follow 
  # https://stats.stackexchange.com/questions/14853/variable-importance-from-glmnet/211396#211396
  # The non-standardized coefficient is the amount of change in the result for 
  # a unit change in the predictor. Normalized coeffcifient is the change in 
  # the result for a 1-standard-deviation change in the predictor; 
  # to get this, you must multiply by the SD
  sds <- apply(X.train, 2, sd)
  std_coefs <- coefficients[-1, 1] * sds
  
  rank.importance <- as.matrix(std_coefs[which(std_coefs != 0)])
  
  rank.importance <- rank.importance[order(abs(rank.importance), decreasing = T), , drop = F]
  
  list.return <- list("importance" = importance,
                      "rank.importance" = rank.importance)
  return(list.return)
  
}




# Cosmetics function for names of variables -------------------------------

# Cosmetics: change the name of variables into the types for better readability
# Add some regex syntax to the vector of names
cosmetics_name <- function(names.pe, character.vector.to.change){
  
  if(isFALSE(exists("names.pe"))){stop("Need the object 'names.pe' for names.")}
  
  # Define pattern and replacement
  pattern <- as.character(names.pe[, 2])
  replacement <- as.character(names.pe[, 1])
  
  # Replace values in character vector of interest
  names.new <- 
    stri_replace_all_regex(character.vector.to.change, 
                           paste0(pattern, "$"), 
                           replacement,
                           vectorize_all = FALSE) %>% 
    stri_replace_all_regex(., 
                           paste0(pattern, ".LoO"), 
                           paste0(replacement, ".LoO"),
                           vectorize_all = FALSE) %>% 
    stri_replace_all_regex(., 
                           paste0(pattern, ":"), 
                           paste0(replacement, ":"),
                           vectorize_all = FALSE) 
  
  # Change factor levels
  levels <- c(as.character(names.pe[, 1]),
              as.character(names.new[which(!names.new %in% as.character(names.pe[, 1]))])) 
  levels <-  unique(levels)
  names.new <- factor(names.new, levels = levels)
  
  return(names.new)
}



# Simulation functions ----------------------------------------------------


  # Simulation Carrell et al ------------------------------------------------
simulation.Carrell <- function(DT, n.sims, level.randomization, level.treatment,
                               names.pe){
  
  
  # Rename the variables in data.table for speed
  setnames(DT, old = level.randomization, new = "level.rand")
  setnames(DT, old = level.treatment, new = "level.treat")
  
  draws <- list()
  
  # Simulation
  for (jj in seq_len(n.sims)){
    
    # Generate random year variable
    DT <- DT[, rand.int := level.treat[sample(.N, replace = TRUE)],
             by = list(level.rand)]
    
    # Compute mean of covariates per random year
    cols <- c(names.pe$pe)
    
    draws[[jj]] <- 
      DT[ , lapply(.SD, mean), by = list(level.rand, rand.int), .SDcols = cols][, simcount := jj]
    
    cat(paste0("Draw", jj, "... "))
  }
  
  cat("Random draws generated.")
  
  # Bind all elements of the list
  random.carrell <- do.call(rbind, draws) %>% 
    arrange(level.rand, rand.int)
  
  # Merge with actual distribution
  cols <- c(names.pe$pe)
  original.sd <- DT[ , lapply(.SD, mean), 
                               by = list(level.rand, 
                                         level.treat), 
                               .SDcols = cols] 
  original.sd <- original.sd[order(level.rand, 
                                   level.treat)]
  
  random.carrell.final <- left_join(random.carrell, original.sd,
                                    by = c("level.rand", "rand.int" = "level.treat"))
  
  
  # Conduct tests per variable
  cat("\n Conducting tests: ")
  test.results <- list()
  
  for (kk in seq_len(length(names.pe$pe))) {
    
    tmpx  <- paste0("type_", kk, ".x")
    tmpy  <- paste0("type_", kk, ".y")
    
    # Calculate share of coefficients that are larger than the observed ones
    TEST.dt <- random.carrell.final[, type := get(tmpx)-get(tmpy)]
    TEST.dt <- TEST.dt[, .(level.rand, rand.int, type, simcount)]
    TEST.dt <- TEST.dt[, above := ifelse(type > 0, 1, 0)]
    TEST.dt <- TEST.dt[, prop := sum(above), by = list(level.rand, rand.int)]
    TEST.dt <- TEST.dt[, sum := .N, by = list(level.rand, rand.int)]
    TEST.dt <- TEST.dt[, share_above := prop/sum][, sum := NULL]
    
    # Take only one observation per group
    TEST.dt <- TEST.dt[ ,`:=`(simcount = NULL, above = NULL,
                              type = NULL)]
    # Take first observation per group
    TEST.dt <- TEST.dt[, .SD[1], list(level.rand, rand.int)]
    
    # Split the data in multiple groups into a list to conduct test
    df.list <- split(TEST.dt, list(TEST.dt$level.rand), drop=TRUE) 
    
    # Conduct a KS test against uniform distribution for each group
    test.list <- 
      map(df.list, ~pluck(., "share_above") %>% 
            as.matrix() %>% 
            ks.test(., "punif"))
    
    # Extract the p-values for each test and count the number of times it is not rejected
    matrix.test <- 
      map(test.list, pluck, "p.value") %>% 
      as.matrix()  
    
    n.nonrejections <- length(matrix.test[which(matrix.test[, 1] > 0.05), ])
    
    cat(paste0(kk,"..."))
    test.results[[kk]] <- list("type" = paste0("type_", kk), 
                               "nonrejections" = n.nonrejections)
  }
  
  return(test.results)
}




# Simulation Bifulco et al ------------------------------------------------
# DT <- data_full.dt
# n.sims = 50
# demeaning <- NULL
# trend.factor <- "swjahr"
# rm(list_factor)


simulation.Bifulco <- function(DT, n.sims, level.randomization, level.treatment,
                               names.pe, demeaning = NULL, trend.factor = NULL){
  
  
  # Rename the variables in data.table for speed
  setnames(DT, old = level.randomization, new = "level.rand", skip_absent=TRUE)
  setnames(DT, old = level.treatment, new = "level.treat", skip_absent=TRUE)
  
  cols <- c(names.pe$pe)
  
  if(!is.null(trend.factor) & !is.null(demeaning)) {
    setnames(DT, old = trend.factor, new = "trend.factor", skip_absent=TRUE)
    list_factor <- list(as.factor(DT[["level.rand"]]),
                        as.factor(DT[["trend.factor"]]))} 
  if (is.null(trend.factor) & !is.null(demeaning)) {
    list_factor <- list(as.factor(DT[["level.rand"]])) 
  }
  
  draws <- list(); draws.sd <- list()
  
  for (jj in seq_len(n.sims)){
    
    # Simulation step
    DT.loop <- DT[, rand.int := level.treat[sample(.N, replace = TRUE)], 
                       by = list(level.rand)]
    
   # If demeaning
    if(is.null(demeaning)){
      # Per draw, take share of type
      draws[[jj]] <- 
        DT.loop[ , lapply(.SD, mean), by = list(rand.int, level.rand), .SDcols = cols] [order(-rand.int, -level.rand)]
      
      draws.sd[[jj]] <- draws[[jj]][ , lapply(.SD, sd), .SDcols = cols]
    } else {
      
      # For demeaning
      demeaned.mean.matrix <- 
        lfe::demeanlist(DT.loop[, ..cols],
                        fl = list_factor)
      
      demeaned.mean.matrix <- as.data.table(cbind(demeaned.mean.matrix, 
                                                  "level.rand" = DT.loop[["level.rand"]],
                                                  "rand.int" = DT.loop[["rand.int"]]))
      draws[[jj]] <- demeaned.mean.matrix[ , lapply(.SD, mean), 
                                           by = list(level.rand, rand.int), .SDcols = cols] [order(-level.rand, rand.int)]
      
      draws.sd[[jj]] <- draws[[jj]][ , lapply(.SD, sd), .SDcols = cols]
    }    
    
    # Message
    cat(paste0("."))
  }
  
  # Bind all elements
  random <- bind_rows(draws) %>% 
    arrange(level.rand, rand.int)
  # random <- draws[[1]]
  
  random.sd <- bind_rows(draws.sd) 
  # random.sd <- draws.sd[[1]]
  
  # Compute the mean and mean sd per type per draw
  mean.matrix <- data.table(random) %>% 
    .[ , lapply(.SD, mean), 
       by = list(level.rand, rand.int),
       .SDcols = cols] 
  sd.matrix <- data.table(random.sd) %>% 
    .[ , lapply(.SD, mean), .SDcols = names(random.sd)] 
  
  
  # Compute share in original data
  if(is.null(demeaning)){
    # Per draw, take share of type
    original.mean <- DT[ , lapply(.SD, mean), by = list(level.rand, level.treat), .SDcols = cols]
    # original.mean <- DT[ , lapply(.SD, mean), .SDcols = cols] 
  } else {
    demeaned.data <- 
      lfe::demeanlist(DT[, ..cols],
                      fl = list_factor)
    demeaned.data <- as.data.table(cbind(demeaned.data, 
                                         "level.rand" = DT[["level.rand"]],
                                         "level.treat" = DT[["level.treat"]]))
    original.mean <- demeaned.data[ , lapply(.SD, mean), 
                                    by = list(level.rand, level.treat),
                                    .SDcols = cols] 
  }
  
  
  return(list("original.mean" = original.mean,
              "mean.matrix" = mean.matrix,
              "sd.matrix"  = sd.matrix))
         
}       



# High-dimensional interactions -------------------------------------------

hd_internal <- function(names.cols, df){
  require(data.table)
  
  # Split the colnames before and after the ":"
  cols <- c(stringr::str_extract(names.cols, ".+(?=\\:)"), 
            stringr::str_extract(names.cols, "(?<=\\:).+"))
  
  df <- df[, ..cols]
  newdf <- df[,1]*df[,2]
  setnames(newdf, names.cols)
  
  return(newdf)
}

hdinteractions <- function(names.cols, df, parallel = FALSE, 
                           sequential = FALSE, K = 100){
  
  require(parallel)
  
  df <- data.table(df)
  
  # if(isTRUE(parallel)){
  #   no_cores <- detectCores()
  #   clust <- makeCluster(no_cores) 
  #   
  #   list.cols <- parLapply(clust, names.cols, hd_internal, df = df)
  #   stopCluster(clust)
  # } else {
  #   
  #   list.cols <- lapply(names.cols, hd_internal, df = df)
  # }
  # 
  # # Select variables with colSums > 0
  # list.cols <- list.cols[which(lapply(list.cols, colSums) > 0)]
  # 
  # # Bind into df
  # x <- bind_cols(list.cols)
  
  
  # Sequential
  split   <- runif(length(names.cols))
  folds   <- as.numeric(cut(split, quantile(split, probs = seq(0, 1, 1/K)),
                            include.lowest = TRUE)) 
  
  list.dfs <- list()
  for (jj in 1:K){
    list.cols <- lapply(names.cols[folds == jj], hd_internal, df = df)
    list.cols <- list.cols[which(lapply(list.cols, colSums) > 0)]
    
    list.dfs[[jj]] <- bind_cols(list.cols)
    rm(list.cols)
    cat(paste0(jj, "..."))
  }
  x <- bind_cols(list.dfs)
  
  
  return(x)
  
}




# Helper function to create felm coeffs from glinternet output. -----------
# path.excel <- "Script/csv/vars.selected.postfelm.xlsx"
# sheet.excel <- paste0("interaction.", var.clustering.level)
# stab.object <- stab.glinternet.interaction
# trans.felm <- final.felm

# Prepare graphs and create variables from excel file
post.felm <- function(sheet.excel, path.excel, stab.object, var.clustering.level, 
                      trans.felm){
  
  # EXCEL FILE:
  #' - selected.var is the name of the selected var as obtained by glinternet$selected
  #' - label is created label. Must be filled for variables that require split sample. NA otherwise
  #' - regression.var is the var that we want to capture from the regression output. Same as selected
  #'   var, but differs only if split sample
  #' - formula: formula passed for felm estimation
  #' - split.sample: gives the type the sample is split on for interactions types x LoO
  #' - split.value: value for the split sample. Either 1 or 0
  #' - create.var: indicator if a variable must be created. Creates interaction + LoO
  #' - create.var1: first term of interaction
  #' - create.var2: second term of interaction
  #' - create.var.int: name for the new LoO of the created interaction.
  
  # Manually:
  xls_tibble <- read_excel(path.excel, sheet = sheet.excel)%>% 
    as_tibble()
  
  if(isFALSE(setequal(colnames(xls_tibble), c("selected.var", "label", "regression.var", "formula",
                                   "split.sample", "split.value", "create.var",
                                   "create.var1", "create.var2", "create.var.int")))){
    stop("Colnames of .xls must be specified differently.")
  }
  
  dat <- trans.felm$finaldf
  
  # Create variables
  create.tab <- xls_tibble[xls_tibble$create.var == 1, ]
  if(nrow(create.tab)>0){
    
    cols.new.data <- list()
    
    for (jj in 1:nrow(create.tab)){
      name_1 <- as.character(create.tab[jj, "create.var1"] )
      name_2 <- as.character(create.tab[jj, "create.var2"] )
      name_new <- as.character(create.tab[jj, "create.var.int"] )
      
      cols <- c(name_1, name_2, var.clustering.level)
      dfloop <- dat[, ..cols]
      colnames(dfloop) <- c("int1", "int2", "cluster")
      dfloop <- dfloop[, interaction := int1 * int2]
      dfloop <- dfloop[, LoO := lapply(.SD, function(x){(sum(x)-x)/(length(x)-1)}), 
                       .SDcols = "interaction", 
                       by = list(cluster)]
      dfloop <- setnames(dfloop, "LoO", name_new)
      
      cols.new.data[[jj]] <- dfloop[, ..name_new]
    }
    
    new.cols <- do.call(cbind, cols.new.data)
    dat <- cbind(dat, new.cols)
  }
  
  # Create tibble for estimation with split sample
  felm.df <- mutate(xls_tibble, data = list(dat))
  felm.df <- mutate(felm.df, data = pmap(list(split.sample, data, split.value), 
                                         ~if(is.na(..1)){..2} else{..2[which(..2[[paste0(..1)]] == ..3), ]}))
  
  # Compute FELM
  if(var.clustering.level == "schooltrackyear"){
    felm.coeff <- pmap(list(felm.df$formula, felm.df$data), 
                       ~felm(as.formula(paste0("outcome ~ ", ..1, "|schooltrackyear| 0 |classfix")),
                             data = ..2))
  } else {
    felm.coeff <- pmap(list(felm.df$formula, felm.df$data), 
                       ~felm(as.formula(paste0("outcome ~ ", ..1, "| swjahr + schooltrack | 0 | schooltrackyear")),
                             data = ..2))
  }
  
  felm.coeff <- lapply(felm.coeff, function(x) tidy(x,  conf.int = T, conf.level = 0.95))
  
  felm.coeff.list <- list()
  for(ii in 1:length(felm.coeff)){
    felm.coeff.list[[ii]] <- felm.coeff[[ii]][
      which(as.matrix(felm.coeff[[ii]][,1]) %in% as.matrix(felm.df[ii, 3])), ]
  }
  
  felm.coeff.list <- do.call(rbind, felm.coeff.list)
  felm.coeff <- cbind(felm.df[, c(1,2)], felm.coeff.list)
  
  return(felm.coeff)
}
