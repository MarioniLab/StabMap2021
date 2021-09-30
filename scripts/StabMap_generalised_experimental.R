############################
# simulate some data
if (FALSE) {
  full_expr = matrix(rnorm(300*150), nrow = 300, ncol = 150,
                     dimnames = list(paste0("gene_", 1:300),
                                     paste0("cell_", 1:150)))
  
  assay_list = list(
    D_R = full_expr[1:150, 1:50],
    D_j = full_expr[76:215, 51:100],
    D_i = full_expr[151:300, 101:150]
  )
  
  # labels list
  labels_list = list(
    D_R = rep(letters[1:5], each = 10)
  )
  
  # whether the data should be treated as reference
  reference_list = list(
    D_R = TRUE,
    D_j = FALSE,
    D_i = TRUE
  )
}
############################

reWeightEmbedding = function(embedding, weights = NULL, factor = 1e6) {
  
  # embedding is a cells x dimensions matrix
  # weights is an optional named list (names correspond to all before the underscore)
  # for weighting each embedding component
  # factor is a multiplicative value to avoid tiny numbers
  cols = as.character(interaction(gsub("_.+", "", colnames(embedding)),
                                  ifelse(grepl("_PC", colnames(embedding)), "PC", "LD"), sep = "_"))
  
  cols_split = split(colnames(embedding), cols)
  
  if (is.null(weights)) {
    weights = lapply(cols_split, function(x) 1)
  }
  
  norms = lapply(cols_split, function(cols) {
    sum(abs(embedding[,cols]))
  })
  
  norms_long = unsplit(norms, cols)
  weights_long = factor*unlist(weights)[cols]
  
  embedding_norm = t( (t(embedding) / norms_long) * weights_long)
  if (FALSE) {
    barplot(colSums(embedding_norm^2))
    boxplot(embedding_norm)
  }
  return(embedding_norm)
}


# feature overlaps between datasets
plotFeatureOverlaps = function(assay_list) {
  require(UpSetR)
  g = upset(as.data.frame(1*do.call(cbind, lapply(assay_list, function(x) Reduce(union, lapply(assay_list, rownames)) %in% rownames(x)))))
  print(g)
  return(g)
}

# generate a network of the datasets
featureNetwork = function(assay_list) {
  require(igraph)
  # given a list of assays, generate a 
  # network relating the datasets to each other
  # in terms of number of shared features
  # (rownames)
  
  datasets = names(assay_list)
  
  pairs = t(combn(datasets, 2))
  
  edge_weights = apply(pairs, 1, function(x) {
    length(Reduce(intersect, lapply(assay_list[x], rownames)))
  })
  
  pairs_overlapping = pairs[edge_weights != 0,, drop = FALSE]
  edge_weights_overlapping = edge_weights[edge_weights != 0]
  
  g = graph.edgelist(pairs_overlapping, directed = FALSE)
  E(g)$weight <- edge_weights_overlapping
  
  return(g)
}

# the start of the generalised StabMap function

# input: assay_list
# optional input: labels_list
# optional input: reference_list

# labels_list is a list of sample labels,
# assumed to be the same order as the colnames of the
# assay_list entry of the same name

# reference_list is a list of logical values
# as to whether this embedding should be included

# treating one of the datasets as a reference, identify the 
# paths for each of the other datasets

"%**%" <- function(X,Y) {
  # multiply two matrices but first reorder rows of the
  # second matrix
  
  # if Y is a list, then include a scaling subtraction
  # assumed given as second item
  
  if (is.list(Y)) {
    rmeans = Y[[2]]
    Y <- Y[[1]]
  } else {
    rmeans = rep(0, nrow(Y))
    names(rmeans) <- rownames(Y)
  }
  
  features = intersect(colnames(X), rownames(Y))
  XY = (X[,features] - rmeans[features]) %*% Y[features,]
  return(XY)
}

# "%pred%" <- function(a, b) {
#   # a is a matrix
#   # b is a model to pass through predict
#   
#   if (class(b) == "lda") {
#     am = MASS:::predict.lda(b, newdata = a)$x
#   }
#   if (class(b) == "svm") {
#     am = attr(predict(b, newdata = a, decision.values = TRUE), "decision.values")
#   }
#   
#   return(cbind(am))
#   
# }

"%projpred%" <- function(a, b) {
  # a is a matrix
  # b is a list of a 1. projection matrix
  # and 2. a model to pass through predict
  # if b is not a list, then just do normal
  # matrix multiplication
  # alternatively, if b is already an lda 
  # object then just perform the prediction
  
  if (class(b)[1] == "lda") {
    features = rownames(b$scaling)
    am = MASS:::predict.lda(b, newdata = a[,features])$x
    return(am)
  }
  
  if (!is.list(b)) {
    return(a %*% b)
  }
  
  ab = a %*% b[[1]]
  if (class(b[[2]]) == "lda") {
    am = MASS:::predict.lda(b[[2]], newdata = a)$x
  }
  if (class(b[[2]]) == "svm") {
    am = attr(predict(b[[2]], newdata = a, decision.values = TRUE), "decision.values")
  }
  
  return(cbind(ab,am))
  
}

"%*1%" <- function(a, b) {
  if (is.list(a)) {
    a <- a[[1]]
  }
  # if (is.list(b)) {
  #   b <- b[[1]]
  # }
  cbind(intercept = 1, a) %**% b
}


calculatePCA_irlba = function(data, nPC = 50, ntop = 500) {
  
  require(irlba)
  
  # get the relevant genes via calculatPCA
  features = rownames(attr(calculatePCA(data, ncomponents = nPC, ntop = ntop, scale = FALSE),
                           "rotation"))
  
  out_pca = prcomp_irlba(t(data[features,]), n = nPC,
                         center = FALSE, scale. = FALSE)
  out_loadings = out_pca$rotation
  rownames(out_loadings) <- features
  out_scores = out_pca$x
  rownames(out_scores) <- colnames(data)
  attr(out_scores, "rotation") <- out_loadings
  return(out_scores)
}

runOps = function(obj, ops, leftToRight = TRUE) {
  if (leftToRight) {
    out = obj[[1]]
    for (i in 1:length(ops)) {
      out <- get(ops[[i]])(out, obj[[i+1]])
    }
  } else {
    out = obj[[length(obj)]]
    for (i in length(ops):1) {
      out <- get(ops[[i]])(obj[[i]], out)
    }
  }
  return(out)
}


stabMapGeneralised = function(assay_list,
                              labels_list = NULL,
                              reference_list = sapply(names(assay_list), function(x) TRUE, simplify = FALSE),
                              reference_features_list = lapply(assay_list, rownames),
                              ncomponentsReference = 50,
                              ncomponentsSubset = 50,
                              suppressMessages = TRUE,
                              projectAll = FALSE,
                              maxFeatures = 1000,
                              plot = TRUE,
                              scale.center = FALSE,
                              scale.scale = FALSE) {
  
  # require packages
  require(igraph)
  require(scater)
  
  # defensive programming, check various things
  # the columns of each assay_list (cells) should all have different names
  # each of the assays should have rownames and colnames
  
  # assay_list should have names
  
  # if labels_list given the entries should match the ncol(assay_list)
  
  # remove messages from calculatePCA
  if (suppressMessages) {
    sm = function(expr) {
      suppressWarnings(suppressMessages(eval(expr)))
    }
  } else {
    sm = function(expr) {
      eval(expr)
    }
  }
  
  # reference_list should have names
  # if it's not a list, then should be a character vector of the names
  # of assay_list
  if (is.character(reference_list)) {
    reference_list <- sapply(names(assay_list), function(x) x %in% reference_list, simplify = FALSE)
  }
  
  if (plot) plotFeatureOverlaps(assay_list)
  
  assay_network = featureNetwork(assay_list)
  if (plot) plot(assay_network)
  
  # check whether the network is a connected component
  # the number of components should be 1:
  if (components(assay_network)$no != 1) {
    stop("feature network is not connected, features must overlap in some way")
  }
  
  # if needed, scale the data
  if (any(c(scale.center, scale.scale))) {
    assay_list <- lapply(
      assay_list,
      function(x) t(scale(t(x), center = scale.center, scale = scale.scale))
    )
  }
  
  
  all_embeddings_list = list()
  
  for (reference_dataset in names(assay_list)) {
    # reference_dataset = names(assay_list)[1]
    
    # if not a reference just go to the next one
    if (!reference_list[[reference_dataset]]) next
    
    message(paste0("treating \"", reference_dataset, "\" as reference"))
    
    # when the graph has a weight, then by default it will
    # use them
    to_nodes = names(sort(distances(assay_network, to = reference_dataset)[,reference_dataset]))
    all_paths = lapply(all_shortest_paths(assay_network, 
                                          from = reference_dataset)$res, names)
    names(all_paths) <- unlist(lapply(all_paths, function(x) rev(x)[1]))
    all_paths <- all_paths[to_nodes]
    
    # the PC space of the reference dataset
    nPC = min(ncomponentsReference, nrow(assay_list[[reference_dataset]]))
    
    for (projectionType in c("PC", "LD")) {
      # projectionType = "LD"
      
      if (projectionType %in% "PC") {
        
        
        
        if (TRUE) {
          reference_scores = sm(calculatePCA(assay_list[[reference_dataset]][reference_features_list[[reference_dataset]],], ncomponents = nPC,
                                             scale = FALSE))
          attr(reference_scores, "loadings") <- list(attr(reference_scores, "loadings"),
                                                     rowMeans(assay_list[[reference_dataset]]))
        } 
        
        # if labels are given, identify the rotation of PCs that represent LDs
        d_nPC = diag(nPC)
        colnames(d_nPC) <- paste0(reference_dataset, "_", colnames(reference_scores))
        
        P_0 = d_nPC
        
      }
      
      
      if (projectionType %in% "LD") {
        
        if (is.null(labels_list[[reference_dataset]])) next
        
        require(MASS)
        
        message(paste0("labels provided for \"", reference_dataset, "\", adding LD components"))
        
        # features = Reduce(intersect, lapply(assay_list[path_current[1:2]], rownames))
        features = Reduce(intersect, lapply(assay_list, rownames))
        # features = rownames(assay_list[[reference_dataset]])
        if (length(features) == 0) {
          features = reference_features_list[[reference_dataset]]
        }
        
        if (length(features) > maxFeatures) {
          message("more input features than maxFeatures, subsetting features using variance ranking")
          require(scran)
          genevars = modelGeneVar(assay_list[[reference_dataset]][features,])
          genevars_sorted = genevars[order(genevars$bio, decreasing = TRUE),]
          features <- rownames(genevars_sorted)[seq_len(maxFeatures)]
        }
        
        labels_train = labels_list[[reference_dataset]]
        
        require(Matrix)
        # experimental: remove features with zero variance for LDA
        vars = rowMaxs(apply(fac2sparse(labels_train), 1, function(x)
          rowWeightedVars(assay_list[[reference_dataset]][features,],x)), na.rm = TRUE)
        if (any(vars == 0)) message("removing genes with zero intra-class variance")
        features <- features[vars > 0]
        
        
        data_train = t(assay_list[[reference_dataset]][features,])
        
        lda.fit = sm(lda(data_train[!is.na(labels_train),],
                         grouping = labels_train[!is.na(labels_train)]))
        colnames(lda.fit$scaling) <- paste0(reference_dataset, "_", colnames(lda.fit$scaling))
        # lda.scaling = lda.fit$scaling
        # colnames(lda.scaling) <- paste0(reference_dataset, "_", colnames(lda.scaling))
        
        if (FALSE) {
          P_0 = lda.fit
        }
        # reference_scores = t(assay_list[[reference_dataset]][features,])
        reference_scores = t(assay_list[[reference_dataset]][features,]) %projpred% lda.fit
        
        d_nLD = diag(ncol(reference_scores))
        colnames(d_nLD) <- colnames(reference_scores)
        
        P_0 = d_nLD
        
        
        
        # if (!is.null(labels_list[[reference_dataset]])) {
        if (FALSE) {
          
          require(MASS)
          
          message(paste0("labels provided for \"", reference_dataset, "\", adding LD components"))
          
          labels_train = labels_list[[reference_dataset]]
          data_train = reference_scores
          
          
          if (TRUE) {
            lda.fit = lda(data_train[!is.na(labels_train),],
                          grouping = labels_train[!is.na(labels_train)])
            # lda.scaling = lda.fit$scaling
            # colnames(lda.scaling) <- paste0(reference_dataset, "_", colnames(lda.scaling))
            
            # experimental:
            # colnames(lda.fit$scaling) <- paste0(reference_dataset, "_", colnames(lda.fit$scaling))
          } else {
            # experimental:
            require(e1071)
            lda.fit = svm(x = data_train[!is.na(labels_train),],
                          y = factor(labels_train[!is.na(labels_train)]),
                          kernel = "linear")
            message("fitted labels model")
          }
          
          if (FALSE) {
            # P_0 = cbind(d_nPC, lda.scaling)
          } # experimental:
          P_0 = list(d_nPC, lda.fit)
        } else {
          # because the reference scores has nPC components
          # P_0 = d_nPC
        }
        
        
      }
      
      
      embedding_list = list()
      
      for (path in all_paths) {
        # path = all_paths[[2]]
        
        message(paste0("generating embedding for path with reference \"",
                       reference_dataset,
                       "\": ",
                       paste0(rev(paste0("\"", path, "\"")), collapse = " -> ")))
        
        if (identical(path, reference_dataset)) {
          if (FALSE) {
            embedding_list[[reference_dataset]] <- reference_scores %*% P_0
          }
          # experimental:
          embedding_list[[reference_dataset]] <- reference_scores %projpred% P_0
          next
        }
        
        path_current = path
        
        P = P_0
        
        # experimental:
        obj = list(P)
        ops = list("%projpred%")
        
        while(length(path_current) > 1) {
          
          features_current = Reduce(intersect, lapply(assay_list[path_current[1:2]], rownames))
          if (path_current[1] == reference_dataset) {
            current_scores = as.matrix(reference_scores)
          } else {
            current_obj = obj[[length(obj)]]
            if (is.list(current_obj)) current_obj <- current_obj[[1]]
            scores_features = setdiff(rownames(current_obj), "intercept")
            # if (FALSE) {
            #   current_scores = as.matrix(t(assay_list[[path_current[1]]][rownames(P),]))
            # }
            current_scores = as.matrix(t(assay_list[[path_current[1]]][scores_features,]))
          }
          
          
          if (length(path_current) > 2) {
            nPC_sub = min(ncomponentsSubset, length(features_current))
            
            # if (TRUE) {
            dimred_current = sm(calculatePCA(assay_list[[path_current[1]]][features_current,],
                                             ncomponents = nPC_sub, scale = FALSE))
            attr(dimred_current, "rotation") <- list(attr(dimred_current, "rotation"),
                                                     rowMeans(assay_list[[path_current[1]]][features_current,]))
            # } else {
            #   # experimental: replace calculatePCA with irlba
            #   dimred_current = sm(calculatePCA_irlba(assay_list[[path_current[1]]][features_current,],
            #                                          nPC = nPC_sub))
            # }
            loadings_current = attr(dimred_current, "rotation")
            
            # if (length(path_current) > 2) {
            # if (FALSE) {
            coef = lm.fit(cbind(intercept = 1, dimred_current), current_scores)$coefficients
            coef <- na.omit(coef)
            
            # even more experimental:
            # obj <- c(coef, obj)
            obj[[length(obj) + 1]] <- coef
            
            # if (TRUE) {
            # if (length(path_current) == 2) {
            # obj <- c(loadings_current, obj)
            # if (TRUE) {
            obj[[length(obj) + 1]] <- loadings_current
            # }
            # ops <- c("%projpred%",ops)
            ops <- c("%*1%", ops)
            ops <- c("%**%", ops)
            # }
            
            
          } else {
            # if (FALSE) {
            # coef = lm.fit(cbind(intercept = 1, as.matrix(t(assay_list[[path_current[1]]][features_current,]))), current_scores)$coefficients
            # }
            
            
            # experimental: 
            # if there are more than 5000 in features_current,
            # then restrict to HVGs
            # (set as a parameter maxFeatures ?)
            
            # maxFeatures = 5000
            if (length(features_current) > maxFeatures) {
              message("more input features than maxFeatures, subsetting features using variance ranking")
              require(scran)
              genevars = modelGeneVar(assay_list[[path_current[1]]][features_current,])
              genevars_sorted = genevars[order(genevars$bio, decreasing = TRUE),]
              features_current <- rownames(genevars_sorted)[seq_len(maxFeatures)]
            }
            
            # print("fitting linear models")
            
            coef = lm.fit(
              cbind(
                intercept = 1, 
                as.matrix(t(assay_list[[path_current[1]]][features_current,]))
              ),
              current_scores
            )$coefficients
            coef <- na.omit(coef)
            
            # print("fitted linear models")
            
            # even more experimental:
            # obj <- c(coef, obj)
            obj[[length(obj) + 1]] <- coef
            
            # ops <- c("%*1%",ops)
            ops <- c("%*1%", ops)
            
            
          }
          # if length(path_current) == 2 then this is the last step,
          # i.e. can use all genes in the remaining assay for coefficient estimation
          
          
          # experimental: checking for numerical overflow?
          # matList[[length(matList) + 1]] <- coef
          # matList[[length(matList) + 1]] <- cbind(1, loadings_current)
          
          # # projection matrix
          # if (FALSE) {
          #   P <- cbind(interecept = 1, loadings_current) %*% coef %*% P
          # }
          # # experimental:
          # if (FALSE) {
          #   P <- cbind(intercept = 1, loadings_current) %*% coef %projpred% P
          # }
          
          ## experimental:
          path_previous <- path_current[1]
          
          path_current <- path_current[-1]
        }
        
        # now that the path is just length 1, perform the projection
        if (FALSE) {
          embedding_list[[path_current]] <- t(assay_list[[path_current]]) %**% P
        } else {
          # experimental:
          # obj <- c(t(assay_list[[path_current]]), obj)
          obj[[length(obj) + 1]] <- t(assay_list[[path_current]])
          if (length(obj) - length(ops) != 1) {
            ops <- c("%**%",ops)
          }
          embedding_list[[path_current]] <- runOps(rev(obj), ops, leftToRight = FALSE)
        }
        
        ## experimental: also project the prior data too
        # embedding_list[[path_previous]] <- t(assay_list[[path_previous]]) %**% P
        if (projectAll) {
          obj[[length(obj)]] <- t(assay_list[[path_previous]])
          # ops <- c("%**%",ops)
          embedding_list[[path_previous]] <- runOps(rev(obj), ops, leftToRight = FALSE)
          
          for (path_previous_previous in path) {
            
            if (path_previous_previous == path_previous) break
            
            # path_previous_previous = path[which(path == path_previous) - 1]
            features_previous_previous = rownames(assay_list[[path_previous_previous]])
            
            previous_previous_ind = max(which(unlist(lapply(obj, function(x){
              if (is.list(x)) {
                x <- x[[1]]
              }
              any(features_previous_previous %in% rownames(x))
            }))))
            
            obj_previous_previous = obj[1:previous_previous_ind]
            ops_previous_previous = rev(rev(ops)[1:(length(obj_previous_previous) - 1)])
            
            obj_previous_previous[[length(obj_previous_previous) + 1]] <- t(assay_list[[path_previous_previous]])
            ops_previous_previous <- c("%**%",ops_previous_previous)
            
            embedding_list[[path_previous_previous]] <- runOps(rev(obj_previous_previous), ops_previous_previous, leftToRight = FALSE)
          }
          
        }
      }
      
      embedding = as.matrix(do.call(rbind, embedding_list))
      
      # recenter the embedding:
      embedding <- t(t(embedding) - colMeans(embedding))
      
      all_embeddings_list[[paste0(reference_dataset, "_", projectionType)]] <- embedding
    }
    
  }
  
  # experimental: make sure the rownames are matching when cbinding
  all_cells = rownames(all_embeddings_list[[1]])
  
  # all_embeddings = do.call(cbind, all_embeddings_list)
  all_embeddings = do.call(cbind, lapply(all_embeddings_list, "[", all_cells, ))
  
  return(all_embeddings)
}

getInteractionFactor = function(sce, levels) {
  # sce is a SingleCellExperiment object
  # levels is a character vector of the colData name to extract a factor for
  
  do.call(interaction, 
          sapply(levels, function(x) 
            colData(sce)[,x], simplify = FALSE))
}



###### more simulation

if (FALSE) {
  # first examine the feature relationships:
  plotFeatureOverlaps(assay_list)
  plot(featureNetwork(assay_list))
  out = stabMapGeneralised(assay_list, 
                           labels_list = labels_list,
                           reference_list = c("D_R"),
                           ncomponentsReference = 20, ncomponentsSubset = 20)
  
  # testing imputeEmbedding
  # imp = imputeEmbedding(assay_list, embedding = out, reference = colnames(assay_list[[1]]), query = colnames(assay_list[[2]]))
  # 
}

