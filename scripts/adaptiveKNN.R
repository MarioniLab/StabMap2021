getArgMax = function(M, return_colnames = FALSE) {
  # For each row in a matrix calculate the first index
  # which gives the maximum value
  # keep the rownames
  # if return_colnames then extract the column name,
  # otherwise just the index
  m = max.col(M, ties.method = "first")
  
  if (return_colnames) {
    m <- colnames(M)[m]
  }
  
  names(m) <- rownames(M)
  return(m)
}

getAdaptiveK = function(A,
                        labels = NULL,
                        local = NULL,
                        outputPerCell = TRUE,
                        ...) {
  
  # adaptive k selection for KNN classification
  # Given an accuracy matrix A, with rows corresponding to cells
  # and columns corresponding to candidate k values, with values
  # themselves corresponding to accuracy values (either binary 
  # for single classification, or continuous after multiple 
  # classification)
  # and given an optional factor labelling/grouping of cells
  # identify the k that maximises the accuracy for cells belonging
  # to that label/group
  # if no labelling given, expect a cell-cell similarity network
  # to identify the k that maximises the accuracy for cells within
  # that neighbourhood
  # if neither are given, simply treat all cells as if they have
  # the same labelling/grouping.
  
  # ... includes return_colnames, whether to give the
  # colnames of the best selected, or just the index, 
  # which is default FALSE
  
  # if outputPerCell then return a vector of adaptive k
  # values for each cell, not just for each label type
  # (used for when labels is given)
  
  # if both labels and local given, labels will be 
  # prioritised
  
  # local is a neighbourhood index representation
  # as typically output using BiocNeighbors::findKNN()

  # example data generation
  # data = matrix(rpois(10*20, 10), 10, 20) # 10 genes, 20 cells
  # local = BiocNeighbors::findKNN(t(data), k = 5, get.distance = FALSE)$index
  # A = matrix(runif(100),20,5)
  # colnames(A) <- paste0("K_", 1:5)
  # labels = factor(rep(letters[1:2], each = 10))
  
  require(Matrix)
  require(SpatialUtils)
  
  if (is.null(labels) & is.null(local)) {
    labels = factor(rep("All", nrow(A)))
  }
  
  if (!is.null(labels)) {
    if (class(labels) != "factor") {
      labels <- factor(labels)
    }
    L = fac2sparse(labels)
    
    LA = L %*% A
    
    k_best = getArgMax(LA, ...)
    
    if (outputPerCell) {
      k_best <- k_best[labels]
      names(k_best) <- rownames(A)
    }
    
    return(k_best) 
  }

  # if function still running, then use the neighbours in local
  # ensure that self is also included
  local_self = cbind(seq_len(nrow(A)), local)
  
  LA = apply(A, 2, function(a) rowSums(vectorSubset(a, local_self)))
  
  k_best = getArgMax(LA, ...)
  names(k_best) <- rownames(A)
  
  return(k_best)
}

getModeFirst = function(x, first) {
  # identify the mode of x among the first values
  # x is a character or a factor
  # first is an integer
  # x = knn_class[1,]
  # first = query_best_k[i]
  names(which.max(table(x[1:first])[unique(x[1:first])]))
}

adaptiveKNN = function(knn,
                       class,
                       k_local) {
  
  require(SpatialUtils)
  
  # knn is a k-nearest neighbour matrix, giving the 
  # indices of the training set that the query is 
  # closest to. Rows are the query cells, columns
  # are the NNs, should be a large value. Typically
  # output using BiocNeighbors::queryKNN(,,k = max(k_local))
  
  # class is the labels associated with the training
  # set
  
  # k_local is an integer vector length of the training
  # set, giving the local k to use
  # if k_local is given as a single integer, then
  # that value is used as k for all observations
  
  # example data
  # data = matrix(rpois(10*20, 10), 10, 20) # 10 genes, 20 cells
  # local = BiocNeighbors::findKNN(t(data), k = 5, get.distance = FALSE)$index
  # A = matrix(runif(100),20,5)
  # colnames(A) <- paste0("K_", 1:5)
  # labels = factor(rep(letters[1:2], each = 10))
  # k_local = getAdaptiveK(A, labels = labels)
  # data_2 = matrix(rpois(10*30, 10), 10, 30) # 10 genes, 30 cells
  # knn = BiocNeighbors::queryKNN(t(data), t(data_2), k = 5, get.distance = FALSE)$index
  # class = labels
  
  if (length(k_local) == 1) {
    k_local <- rep(k_local, max(knn))
  }
  
  # first use 1NN to identify the local best k value
  query_best_k = k_local[knn[,1]]
  
  if (any(query_best_k > ncol(knn))) {
    warning("k is larger than nearest neighbours provided, taking all neighbours given")
  }
  
  # convert the KNN names to the class labels
  knn_class = vectorSubset(class, knn)
  
  # extract the most frequent among the nearest local best k neighbours
  # new_class = sapply(1:nrow(knn), function(i) getModeFirst(knn_class[i,], query_best_k[i]))
  # same as:
  new_class = mapply(getModeFirst, split(knn_class, seq_len(nrow(knn_class))), query_best_k,
                     USE.NAMES = FALSE)
  
  return(new_class)
}

isEqual = function(x,y) {
  # returns an integer vector
  1*(as.character(x) == as.character(y))
}

getBinaryAccuracy = function(knn,
                             k_values,
                             class_train,
                             class_true) {
  
  # output is a binary accuracy matrix A
  
  # knn is a k-nearest neighbour matrix, giving the 
  # indices of the training set that the query is 
  # closest to. Rows are the query cells, columns
  # are the NNs, should be a large value. Typically
  # output using BiocNeighbors::queryKNN(,,k = max(k_values))
  
  # k_values is an integer vector of the values of k
  # to consider for extracting accuracy
  # if k_values has names then pass these to 
  # colnames of A
  
  # class_train is a factor or character vector of
  # classes that corresponds to the indices given
  # within knn
  
  # class_true is a factor or character vector that
  # corresponds to the rows of knn
  # if class_true has names then pass these to rownames
  # of A
  
  # example data
  # data = matrix(rpois(10*20, 10), 10, 20) # 10 genes, 20 cells
  # local = BiocNeighbors::findKNN(t(data), k = 5, get.distance = FALSE)$index
  # A = matrix(runif(100),20,5)
  # colnames(A) <- paste0("K_", 1:5)
  # labels = factor(rep(letters[1:2], each = 10))
  # k_local = getAdaptiveK(A, labels = labels)
  # data_2 = matrix(rpois(10*30, 10), 10, 30) # 10 genes, 30 cells
  # knn = BiocNeighbors::queryKNN(t(data), t(data_2), k = 5, get.distance = FALSE)$index
  # class = labels
  # class_train = labels
  # class_true = rep(letters[1], 30)
  # k_values = c(1,3,5)
  
  if (max(k_values) > ncol(knn)) {
    stop("largest k value is larger than neighbours provided in knn,
         select a smaller k value or provide more neighbours")
  }
  
  # class_pred = adaptiveKNN(knn, class = class_train, k_local = k_values[1])
  # isEqual(class_pred, class_true)
  
  class_pred = mapply(adaptiveKNN, k_values, MoreArgs = list(class = class_train, knn = knn))
  A = apply(class_pred, 2, isEqual, y = class_true)
  
  rownames(A) <- names(class_true)
  colnames(A) <- names(k_values)
  
  return(A)
}

allEqual = function(x) {
  all(x == x[1])
}

combineBinaryAccuracies = function(A_list) {
  # A_list is a list containing matrices
  # each matrix must have the same number of columns
  # and contain rownames
  
  # example data
  # A_1 = matrix(rbinom(50, 1, 0.5), 10, 5)
  # rownames(A_1) <- letters[1:10]
  # A_2 = matrix(rbinom(60, 1, 0.5), 12, 5)
  # rownames(A_2) <- letters[5:16]
  # A_list = list(A_1, A_2)

  if (!allEqual(unlist(lapply(A_list, ncol)))) {
    stop("each matrix in A_list must have the same number of columns")
  }
  
  if (any(unlist(lapply(A_list, function(x) is.null(rownames(x)))))) {
    stop("each matrix in A_list must have associated rownames")
  }
  
  all_rows = unlist(lapply(A_list, rownames))
  
  A_exp = do.call(rbind, A_list)
  A_split = split.data.frame(A_exp, all_rows)
  A_means = lapply(A_split, colMeans, na.rm = TRUE)
  
  A = do.call(rbind, A_means)
  
  return(A)  
}