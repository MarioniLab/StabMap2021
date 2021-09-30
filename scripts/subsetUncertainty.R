# functions for calculating uncdertainty under feature subsetting
# using distributional changes of jaccard distances among
# neighbours

embeddingUncertainty = function(embedding_reference,
                                embedding_joint,
                                k = 50) {
  
  require(igraph)
  require(bluster)
  
  # k = 50
  graph = makeKNNGraph(embedding_reference, k = k)
  V(graph)$name <- rownames(embedding_reference)
  message("made KNN graph using reference embedding")
  
  graph_sim = igraph::similarity(graph, method = "jaccard")
  rownames(graph_sim) <- V(graph)$name
  colnames(graph_sim) <- V(graph)$name
  
  full_sim = graph_sim
  message("extracted similarities using reference embedding")
  
  # identify nearest neighbours
  tmp_r = embedding_joint[rownames(embedding_reference),]
  
  
  graph = makeKNNGraph(tmp_r, k = k)
  V(graph)$name <- rownames(tmp_r)
  message("made KNN graph using joint embedding")
  
  graph_sim = igraph::similarity(graph, method = "jaccard")
  rownames(graph_sim) <- V(graph)$name
  colnames(graph_sim) <- V(graph)$name
  
  subset_sim = graph_sim
  message("extracted similarities using joint embedding")
  
  
  ref_knn_name = queryNamedKNN(tmp_r,
                               embedding_joint,
                               k = k)
  # ref_knn_name = apply(ref_knn, 2, function(x) rownames(tmp_r)[x])
  # rownames(ref_knn_name) <- rownames(jointPCs)
  message("identified nearest neighbours between reference and joint embeddings")
  
  # extract cell-specific uncertainty score
  uncertainty_scores = sapply(rownames(ref_knn_name), function(i) {
    # if (verbose) print(i)
    ref_sim_nn = full_sim[ref_knn_name[i,], ref_knn_name[i,]]
    ref_sim_sub_nn = subset_sim[ref_knn_name[i,], ref_knn_name[i,]]
    
    stat = suppressWarnings({ks.test(c(ref_sim_nn[lower.tri(ref_sim_nn)]),
                                     c(ref_sim_sub_nn[lower.tri(ref_sim_sub_nn)]))$stat})
    names(stat) <- NULL
    return(stat)
  })
  message("extracted uncertainty scores for all cells in joint embedding")
  
  return(uncertainty_scores)
}

generateSimilarity = function(SCE, k = 50, batchFactor = NULL, HVGs = NULL) {
  # SCE is a single cell experiment object containing the gene expression
  # in "logcounts" slot, otherwise a genes x cells matrix of logcounts
  # k is the number of nearest neighbours in the estimated KNN network
  # batchFactor is a factor matching columns of SCE specifying batches for MNN correction
  # after PCA
  # HVGs is optional set of genes to calculate similarity
  
  require(scran)
  require(SingleCellExperiment)
  require(bluster)
  require(igraph)
  require(scater)
  require(batchelor)
  
  if (!"logcounts" %in% names(assays(SCE))) {
    require(scuttle)
    SCE <- logNormCounts(SCE)
  }
  
  if (is.null(HVGs)) {
    fit = modelGeneVar(logcounts(SCE))
    HVGs = getTopHVGs(fit)
  }
  
  SCE <- runPCA(SCE, subset_row = HVGs)
  
  if (!is.null(batchFactor)) {
    SCE_corrected <- fastMNN(SCE, batch = batchFactor)
    PCs = reducedDim(SCE_corrected, "corrected")
  } else {
    PCs = reducedDim(SCE, "PCA")
  }
  
  graph = makeKNNGraph(PCs, k = k)
  V(graph)$name <- colnames(SCE)
  
  graph_sim = igraph::similarity(graph, method = "jaccard")
  rownames(graph_sim) <- V(graph)$name
  colnames(graph_sim) <- V(graph)$name
  
  return(graph_sim)
}

getSubsetUncertainty = function(SCE,
                                querySCE = NULL,
                                subsetGenes = NULL,
                                k = 50, 
                                full_sim = NULL,
                                plot = FALSE,
                                plotAdditional = NULL,
                                verbose = FALSE,
                                jointBatchFactor = NULL,
                                returnAdditional = NULL,
                                correctSubset = FALSE,
                                ...) {
  
  # output is a named numeric vector of uncertainty values for each cell
  # if querySCE is provided, uncertainty values will include these
  # cells too
  
  # SCE is a SingleCellExperiment object of the reference dataset
  # querySCE is a SingleCellExperiment object of the query dataset,
  # if subsetGenes is NULL then the rownames of these are given as the 
  # subset
  # subsetGenes is a character vector of genes to subset with, this can
  # be NULL if querySCE is provided
  # k integer is the number of nearest neighbours
  # full_sim is a square matrix assumed to be the similarity of the 
  # reference data given as SCE, which can be generated a priori 
  # using generateSimilarity(SCE)
  # jointBatchFactor is a named factor that should have values for SCE and
  # querySCE if provided, which will be included as a batch via interaction
  # with the Reference and Query batch
  # returnAdditional is a character vector of any additional objects to be 
  # returned along with the uncertainty score, e.g. to also extract the 
  # joint PCs set returnAdditional = "jointPCs", or "g" for the plot

  # correctSubset is a logical whether the similarity of the reference
  # data using the subset of genes should be calculated with batch correction
  # applied or not. default FALSE
  
  require(igraph)
  require(BiocNeighbors)
  
  if (is.null(full_sim)) {
    full_sim = generateSimilarity(SCE, ...)
  }
  
  # combining SCE objects is nontrivial in general
  if (is.null(querySCE)) {
    if (is.null(subsetGenes)) stop("Either querySCE or subsetGenes needs to be provided")
    jointSCE = SCE[subsetGenes,]
    batchFactor = rep(c("Reference"), times = c(ncol(SCE)))
  } else {
    if (is.null(subsetGenes)) {
      jointSCE = cbind(SCE, querySCE)[rownames(querySCE),]
      batchFactor = rep(c("Reference", "Query"), times = c(ncol(SCE), ncol(querySCE)))
    } else {
      jointSCE = cbind(SCE, querySCE)[subsetGenes,]
      batchFactor = rep(c("Reference"), times = c(ncol(SCE)))
    }
  }
  
  # add the additional batch factor if given
  if (!is.null(jointBatchFactor)) {
    batchFactor <- interaction(batchFactor, jointBatchFactor[colnames(jointSCE)])
  }
  
  # extract similarity of the subsetted genes, ... will pass the batchFactor defined in the
  # start of the function, despite being redefined above
  dots <- list(...)
  dots[["HVGs"]] <- rownames(jointSCE)
  dots[["SCE"]] <- SCE
  if (!correctSubset) {
  dots[["batchFactor"]] <- NULL
  }
  subset_sim = do.call(generateSimilarity, dots)
  
  # concatenate and batch correct the reference and query datasets (if applicable)
  jointSCE <- logNormCounts(jointSCE)
  jointSCE <- runPCA(jointSCE)
  if (length(unique(batchFactor)) != 1) {
    jointSCE_corrected <- fastMNN(jointSCE, batch = batchFactor)
    jointPCs = reducedDim(jointSCE_corrected, "corrected")
  } else {
    jointPCs = reducedDim(jointSCE, "PCA")
  }
  
  # identify nearest neighbours
  tmp_r = jointPCs[colnames(SCE),]
  
  ref_knn = queryKNN(tmp_r,
                     query = jointPCs,
                     k = k)$index
  ref_knn_name = apply(ref_knn, 2, function(x) rownames(tmp_r)[x])
  rownames(ref_knn_name) <- rownames(jointPCs)
  
  # extract cell-specific uncertainty score
  uncertainty_scores = sapply(rownames(ref_knn_name), function(i) {
    if (verbose) print(i)
    ref_sim_nn = full_sim[ref_knn_name[i,], ref_knn_name[i,]]
    ref_sim_sub_nn = subset_sim[ref_knn_name[i,], ref_knn_name[i,]]
    
    stat = suppressWarnings({ks.test(c(ref_sim_nn[lower.tri(ref_sim_nn)]),
                                     c(ref_sim_sub_nn[lower.tri(ref_sim_sub_nn)]))$stat})
    names(stat) <- NULL
    return(stat)
  })
  
  # generate UMAP for plotting
  if (plot) {
    jointSCE$uncertainty = uncertainty_scores
    reducedDim(jointSCE, "UMAP") <- calculateUMAP(t(jointPCs))
    g = plotUMAP(jointSCE, colour_by = "uncertainty")
    if (!is.null(plotAdditional)) {
      # e.g. plotAdditional = list("celltype", scale_colour_manual(values = celltype_colours))
      require(patchwork)
      gAdditional = plotUMAP(jointSCE, colour_by = plotAdditional[[1]])
      gAll = gAdditional + plotAdditional[[2]] + labs(colour = plotAdditional[[1]]) + g
      print(gAll)
    } else {
      print(g)
    }
  }
  
  if (!is.null(returnAdditional)) {
    out = mget(c("uncertainty_scores", intersect(ls(), returnAdditional)))
    return(out)
  } else {
    return(uncertainty_scores)
  }
}