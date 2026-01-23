### --- Bin (numeric columns) and compute relative frequencies tables ----------
bin_and_table <- function(syn_col, real_train_col, real_test_col = NULL, n_bins = NULL, table = TRUE) {
  
  # Ensure x and y are numeric vectors
  if (is.numeric(real_train_col) & (uniqueN(na.omit(real_train_col)) > 2)) {
    
    range_real_train_col <- diff(range(real_train_col, na.rm = T))
    
    # Determine bin_width using real_train (no pretty, only Freedman-Diaconis or fixed n_bins)
    if (is.null(n_bins)) {
      FD_bin_width <- 2 * IQR(real_train_col, na.rm = T) / (length(na.omit(real_train_col))^(1/3))
      n_bins <- round(range_real_train_col / FD_bin_width)
    }
    
    bin_width <- range_real_train_col / n_bins
    
    # extend to cover full range of combined data
    combined <- c(syn_col, real_train_col, real_test_col)
    n_extra_bins_lo <-ceiling((min(real_train_col, na.rm = T) - min(combined, na.rm = T)) / bin_width)
    n_extra_bins_hi <- ceiling((max(combined, na.rm = T) - max(real_train_col, na.rm = T)) / bin_width)
    
    breaks <- seq(min(real_train_col, na.rm = TRUE) - n_extra_bins_lo * bin_width,
                  max(real_train_col, na.rm = TRUE) + n_extra_bins_hi * bin_width,
                  length.out = n_bins + n_extra_bins_lo + n_extra_bins_hi + 1)
    
    # Bin train, test and syn using same breaks
    real_train_col_binned <- cut(real_train_col, breaks = breaks, include.lowest = TRUE)
    syn_col_binned <- cut(syn_col, breaks = breaks, include.lowest = TRUE)
    
    if (!is.null(real_test_col)) {
      real_test_col_binned <- cut(real_test_col, breaks = breaks, include.lowest = TRUE)
    } else {
      real_test_col_binned <- factor(NULL)
    }
    
  } else {
    # For categorical data, convert to factors
    real_train_col_binned <- factor(real_train_col)
    syn_col_binned <- factor(syn_col)
    if (!is.null(real_test_col)) {
      real_test_col_binned <- factor(real_test_col)
    } else {
      real_test_col_binned <- NULL
    }
  }
  
  # Align factor levels and compute relative frequencies
  all_levels <- union(union(levels(syn_col_binned), levels(real_train_col_binned)), levels(real_test_col_binned))
  real_train_col_binned <- factor(real_train_col_binned, levels = all_levels)
  syn_col_binned <- factor(syn_col_binned, levels = all_levels)
  if (!is.null(real_test_col)) {
    real_test_col_binned <- factor(real_test_col_binned, levels = all_levels)
  } else {
    real_test_col_binned <- NULL
  }
  
  if (!table) {
    list(real_train_col_binned = real_train_col_binned, real_test_col_binned = real_test_col_binned, syn_col_binned = syn_col_binned)
  } else {
    p_real_train_col <- table(real_train_col_binned) / length(na.omit(real_train_col))
    p_syn_col <- table(syn_col_binned) / length(na.omit(syn_col))
    
    if (!is.null(real_test_col)) {
      p_real_test_col <- table(real_test_col_binned) / length(na.omit(real_test_col))
    } else {
      p_real_test_col <- NULL
    }
    
    list(p_real_train_col = p_real_train_col, p_real_test_col = p_real_test_col, p_syn_col = p_syn_col)
  }
}


### --- Mixed-type encoding: min-max scaling (numeric) + one-hot (categorical) -
encode_mixed_dt <- function(syn, real_train, real_test = NULL, base_normalization = real_train, base_OH = rbind(real_train, real_test), scale_OH = 1, output = "data.table") {
  
  # Identify numeric and categorical columns
  num_cols <- names(real_train)[sapply(real_train, is.numeric)]
  cat_cols <- setdiff(names(real_train), num_cols)
  
  # One-hot encode categorical columns using caret dummyVars on real_train data
  if (length(cat_cols) != 0) {
    
    if (!is.null(base_OH)) {
      
      dummy <- dummyVars(~ ., data = base_OH[, ..cat_cols], fullRank = FALSE)
      real_train_cat_oh <- predict(dummy, newdata = real_train[, ..cat_cols]) * scale_OH
      syn_cat_oh <- predict(dummy, newdata = syn[, ..cat_cols]) * scale_OH
      
      if (!is.null(real_test)) {
        real_test_cat_oh <- predict(dummy, newdata = real_test[, ..cat_cols]) * scale_OH
      } else {
        real_test_cat_oh <- NULL
      }
      
    } else {
      
      real_train_cat_oh <- real_train[, ..cat_cols]
      syn_cat_oh <- syn[, ..cat_cols]
      
      if (!is.null(real_test)) {
        real_test_cat_oh <- real_test[, ..cat_cols]
      } else {
        real_test_cat_oh <- NULL
      }
      
    }
  } else {
    real_train_cat_oh <- NULL
    syn_cat_oh <- NULL
    real_test_cat_oh <- NULL
  }
  
  if (length(num_cols) != 0) {
    
    # min-max scaling on numeric features
    if (!is.null(base_normalization)) {
      mins <- apply(base_normalization[, ..num_cols], 2, min, na.rm = T)
      maxs <- apply(base_normalization[, ..num_cols], 2, max, na.rm = T)
    } else {
      mins <- 0
      maxs <- 1
    }
    scale_minmax <- function(x) sweep(sweep(x, 2, mins, "-"), 2, maxs - mins, "/")
    real_train_num_scaled <- scale_minmax(real_train[, ..num_cols])
    syn_num_scaled <- scale_minmax(syn[, ..num_cols])
    
    if (!is.null(real_test)) {
      real_test_num_scaled <- scale_minmax(real_test[, ..num_cols])
    } else {
      real_test_num_scaled <- NULL
    }
    
  } else {
    real_train_num_scaled <- NULL
    syn_num_scaled <- NULL
    real_test_num_scaled <- NULL
  }
  
  # Combine numeric + categorical
  real_train_final <- data.table(real_train_num_scaled, real_train_cat_oh)
  syn_final <- data.table(syn_num_scaled, syn_cat_oh)
  real_test_final <- if(!(is.null(real_test_num_scaled) & is.null(real_test_cat_oh))) data.table(real_test_num_scaled, real_test_cat_oh) else NULL
  
  if (output == "data.table") {
    list(real_train = real_train_final, real_test = real_test_final, syn = syn_final)
  } else if (output == "matrix") {
    list(real_train = as.matrix(real_train_final), real_test = if (!is.null(real_test_final)) as.matrix(real_test_final) else NULL, syn = as.matrix(syn_final))
  } else {
    stop("Invalid output type. Choose either 'matrix', 'data.table', or 'data.dfa'.")
  }
}


### --- Mixed-type Gower-style distance with flexible normalization ranges -----
mixed_distance <- function(x, y, base_normalization = NULL, top_n = NULL, pair_x = NULL, pair_y = NULL, eps = 1e-08, weights = NULL, ignore_case = FALSE,
                           n_chunks_x = 1, parallel_chunks_x = FALSE) {
  
  on.exit(gc(), add = TRUE)
  
  gower_flexible_ranges <-function (x, y, pair_x, pair_y, n = NULL, eps = 1e-08, ranges, weights = NULL, ignore_case= FALSE) 
  {
    stopifnot(is.numeric(eps), eps > 0, nrow(x) > 0, nrow(y) > 
                0, ncol(x) > 0, ncol(y) > 0)
    nthread <- 1L # parallelize in outer loop only
    if (is.null(pair_x) & is.null(pair_y)) {
      xnames <- if (ignore_case) 
        toupper(names(x))
      else names(x)
      ynames <- if (ignore_case) 
        toupper(names(y))
      else names(y)
      pair <- match(xnames, ynames, nomatch = 0L)
    } else if (is.null(pair_x)) {
      pair <- pair_y
    } else {
      if (is.character(pair_x) & is.character(pair_y)) {
        m <- match(names(x), pair_x, nomatch = 0)
        pair_x <- pair_x[m]
        pair_y <- pair_y[m]
      }
      pair <- numeric(ncol(x))
      pair[pair_x] <- pair_y
    }
    if (!any(pair > 0)) {
      message("Nothing to compare")
      return(if (is.null(n)) {
        invisible(numeric(0))
      } else {
        invisible(list(distance = matrix(0)[0, 0], index = matrix(0L)[0, 
                                                                      0]))
      })
    }
    if (!is.null(weights) && (any(weights < 0) || !all(is.finite(weights)))) {
      stop("At least one element of 'weights' is not a finite nonnegative number", 
           call. = FALSE)
    }
    if (!is.null(weights) && length(weights) < length(pair)) {
      msg <- sprintf("%d weights specified, expected %d", length(weights), 
                     length(pair))
      stop(msg, call. = FALSE)
    }
    if (is.null(weights)) 
      weights <- rep(1, ncol(x))
    
    factor_x <- sapply(x, is.factor)
    factor_y <- sapply(y, is.factor)
    for (i in seq_along(pair)) {
      if (pair[i] == 0) 
        next
      iy <- pair[i]
      if (!factor_x[i] & !factor_y[iy]) 
        next
      if (factor_x[i] && !factor_y[iy]) {
        stop("Column ", i, " of x is of class factor while matching column ", 
             pair[i], " of y is of class ", class(y[[iy]]))
      }
      if (!factor_x[i] && factor_y[iy]) {
        stop("Column ", i, " of x is of class ", class(x[[i]]), 
             " while matching column ", pair[i], " of y is of class factor")
      }
      if (factor_x[i] && factor_y[iy]) {
        if (!isTRUE(all.equal(levels(x[[i]]), levels(y[[iy]])))) {
          stop("Levels in column ", i, " of x do  not match those of column ", 
               pair[i], " in y.")
        }
      }
    }
    factor_pair <- as.integer(factor_x)
    eps <- as.double(eps)
    pair <- as.integer(pair - 1L)
    if (is.null(n)) {
      .Call("R_gower", x, y, ranges, pair, factor_pair, eps, 
            weights, nthread, PACKAGE = "gower")
    } else {
      L <- .Call("R_gower_topn", x, y, ranges, pair, factor_pair, 
                 as.integer(n), eps, weights, nthread, PACKAGE = "gower")
      names(L) <- c("index", "distance")
      dim(L$index) <- c(n, nrow(x))
      dim(L$distance) <- dim(L$index)
      dimnames(L$index) <- list(topn = NULL, row = NULL)
      dimnames(L$distance) <- dimnames(L$index)
      L
    }
  }
  
  if (is.null(base_normalization)) {
    base_normalization <- rbind(x, y)
  }
  if (is.data.table(base_normalization)) {
    ranges <- base_normalization[, sapply(.SD, \(col) if (is.numeric(col)) diff(range(col, na.rm = T)) else 1)]
  } else if (is.data.frame(base_normalization)) {
    ranges <- sapply(base_normalization, \(col) if(is.numeric(col)) diff(range(col, na.rm = T))[1] else 1)
  } else if (is.numeric(base_normalization)) {
    ranges <- base_normalization
  } else {
    stop("base_normalization must be NULL (defaults to gower ranges), a data.frame, a data.table, or a numeric vector")
  }
  
  chunk_indices_x <- split(x[,.I], ceiling(x[,.I] / (x[, .N] / n_chunks_x)))
  
  if (is.null(top_n)) {
    
    idy <- y[, .I]
    
    if (parallel_chunks_x && n_chunks_x > 1) {
      result <- foreach (i = seq_len(n_chunks_x), .combine = "c") %dopar% {
        idx_i <- chunk_indices_x[[i]]
        idx <- CJ(x_row = idx_i, y_row = idy)
        gower_flexible_ranges(x[idx$x_row], y[idx$y_row], pair_x = pair_x, pair_y = pair_y, eps = eps, n = top_n, ranges = ranges, weights = weights, ignore_case = ignore_case)
      }
    } else {
      result <- foreach (i = seq_len(n_chunks_x), .combine = "c") %do% {
        idx_i <- chunk_indices_x[[i]]
        idx <- CJ(x_row = idx_i, y_row = idy)
        gower_flexible_ranges(x[idx$x_row], y[idx$y_row], pair_x = pair_x, pair_y = pair_y, eps = eps, n = top_n, ranges = ranges, weights = weights, ignore_case = ignore_case)
      }
    }
  } else {
    if (parallel_chunks_x & n_chunks_x > 1) {
      result <- foreach (i = seq_len(n_chunks_x), .combine = \(a,b) Map(\(x,y) {`dimnames<-`(cbind(x,y), dimnames(x))}, a, b)) %dopar% {
        gower_flexible_ranges(x[chunk_indices_x[[i]]], y, pair_x = pair_x, pair_y = pair_y, eps = eps, n = top_n, ranges = ranges, weights = weights, ignore_case = ignore_case)
      }
    } else {
      result <- foreach (i = seq_len(n_chunks_x), .combine = \(a,b) Map(\(x,y) {`dimnames<-`(cbind(x,y), dimnames(x))}, a, b)) %do% {
        gower_flexible_ranges(x[chunk_indices_x[[i]]], y, pair_x = pair_x, pair_y = pair_y, eps = eps, n = top_n, ranges = ranges, weights = weights, ignore_case = ignore_case)
      }
    }
  }
  
  result
}


### --- Mixed-type correlation: Pearson, Cramer's V eta ------------------------

# Calculate correlation for pairs given variable types
cor_mixed_pair <- function(x, y, type_x, type_y) {
  if (all(c(type_x, type_y) %in% c("numeric", "binary"))) {
    if (type_x == "binary") {
      x <- as.numeric(as.factor(x)) - 1
    }
    if (type_y == "binary") {
      y <- as.numeric(as.factor(y)) - 1
    }
    cor(x, y, use = "pairwise.complete.obs") # Pearson
  } else if (all(c(type_x, type_y) %in% c("binary", "categorical"))) {
    # drop unused factor levels after omitting NAs
    pair <- na.omit(data.table(x = x, y = y))
    x <- factor(pair$x)
    y <- factor(pair$y)
    vcd::assocstats(table(x,y))$cramer # Cramer's V
  } else if (all(c("numeric", "categorical") %in% c(type_x, type_y))) {
    anova <- ifelse(type_x == "numeric", summary(aov(x ~ y)), summary(aov(y ~ x)))
    SS_between <- anova[[1]]$`Sum Sq`[1]
    SS_total <- sum(anova[[1]]$`Sum Sq`)
    sqrt(SS_between / SS_total) # eta
  } else {
    NA  # Fallback for unexpected types
  }
}

# Convert table output to correlation matrix
cor_mixed_2matrix <- function(cor_mixed_table) {
  # output as correlation matrix with variables as row and column names
  cor_full <- rbind(cor_mixed_table,
                    cor_mixed_table[, .(variable1 = variable2,
                                        variable2 = variable1,
                                        type.x = type.y,
                                        type.y = type.x, cor_type, cor)])
  cor_matrix <- dcast(cor_full, variable1 ~ variable2, value.var = "cor")
  cor_matrix <- as.matrix(cor_matrix[,-1], rownames = cor_matrix[[1]])
  diag(cor_matrix) <- 1
  cor_matrix
}

# Mixed-type correlation function for datasets
cor_mixed <- function(data_, matrix = F) {
  data <- as.data.table(copy(data_))
  types <- data[, .(variable = factor(names(data), levels = names(data)), type = sapply(.SD, \(x) {
    if (uniqueN(x[!is.na(x)]) == 1) "degenerate"
    else if (uniqueN(x[!is.na(x)]) == 2) "binary"
    else if (is.logical(x)) "binary"
    else if (is.numeric(x)) "numeric"
    else "categorical"
  }))]
  # get all combinations of variables in two columns
  type_combinations <- CJ(variable1 = types$variable, variable2 = types$variable, sorted = F)[as.numeric(variable1)  < as.numeric(variable2)]
  type_combinations <- merge(type_combinations, types, by.x = "variable1", by.y = "variable", sort = F)
  type_combinations <- merge(type_combinations, types, by.x = "variable2", by.y = "variable", sort = F)
  setcolorder(type_combinations, "variable1", "variable2")
  
  type_combinations[, c("cor_type", "cor") := .(NA_character_, NA_real_)]
  
  type_combinations[type.x == "degenerate" | type.y == "degenerate",
                    cor_type := "degenerate"]
  
  type_combinations[is.na(cor_type), cor_type := apply(.SD, 1, \(x) {
    pair <- na.omit(data[, .SD, .SDcols = c(x[[1]], x[[2]])])
    type1 <- x[[3]]
    type2 <- x[[4]]
    if (any(sapply(pair, uniqueN) == 1)) "degenerate"
    else if (all(c(type1, type2) %in% c("numeric", "binary"))) "Pearson"
    else if (all(c(type1, type2) %in% c("binary", "categorical"))) "Cramer's V"
    else if (all(c("numeric", "categorical") %in% c(type1, type2))) "eta"
    else NA_character_
  })]
  
  # use vectorised Pearson function cor
  if (type_combinations[, any(cor_type == "Pearson")]) {
    pearson_vars <- as.character(types[type %in% c("numeric", "binary"), variable])
    suppressWarnings({
      pearson <- cor(data[, lapply(.SD, \(col) {
        if (is.numeric(col)) col
        else as.numeric(as.factor(col)) - 1
      }), .SDcols = pearson_vars], use = "pairwise.complete.obs")
    })
    
    type_combinations[all(c(type.x, type.y) %in% c("numeric", "binary")), cor := pearson[lower.tri(pearson)]]
  }
  
  type_combinations[cor_type == "degenerate", cor := 0]
  
  type_combinations[is.na(cor), cor := apply(.SD, 1,\(x) cor_mixed_pair(
    data[, .SD, .SDcols = x[[1]]][[1]], 
    data[, .SD, .SDcols = x[[2]]][[1]],
    x[[3]],
    x[[4]]
  ))]
  
  if (!matrix) {
    return(type_combinations[])
  } else {
    return(cor_mixed_2matrix(type_combinations[]))
  }
}
