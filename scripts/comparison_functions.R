mapPCA = function(
  # SCE_list,
  # genes = Reduce(intersect, lapply(SCE_list, rownames)),
  # assayNames = rep("logcounts", length(SCE_list)),
  assay_list,
  features = Reduce(intersect, lapply(assay_list, rownames)),
  nPCs = 50,
  multiBatchPCA = TRUE) {
  
  # Wrapper function to perform PCA on overlapping features
  
  # assay_list is a list of assay matrices (rows are features columns are features)
  # features is the features to include for PCA
  # nPCs is the number of PCs to use
  # multiBatchPCA is a logical as to whether multiBatchPCA should be run,
  # otherwise irlba PCA
  
  # previous input was a list of SCE objects
  # all_assay = do.call(cbind, lapply(mapply(assay, SCE_list, assayNames), "[", genes, ))
  
  all_assay = do.call(cbind, lapply(assay_list, "[", features, ))
  
  batchFac = rep(seq_len(length(assay_list)), times = unlist(lapply(assay_list, ncol)))
  
  if (multiBatchPCA) {
    require(batchelor)
    pca_all <- batchelor::multiBatchPCA(all_assay,
                                        d = nPCs,
                                        batch = batchFac,
                                        preserve.single = TRUE)[[1]]
  } else {
    require(irlba)
    pca_all <- irlba::prcomp_irlba(t(all_assay), n = nPCs, scale. = TRUE)$x
    rownames(pca_all) <- colnames(all_assay)
  }
  
  return(pca_all)
}

UINMF_wrap = function(
  SCE_list = NULL,
  counts_list = NULL,
  ncomponentsSubset = 50) {
  
  # wrapper function to perform the UINMF method
  
  # SCE_list is a list of SingleCellExperiment objects
  # counts_list is a list of read counts (not logcounts!)
  # if SCE_list is given, then counts_list is ignored
  
  # ensure package version is the one with UINMF method included
  # following vignette downloaded 29 april 2021
  # http://htmlpreview.github.io/?https://github.com/welch-lab/liger/blob/master/vignettes/UINMF_vignette.html
  # installed using 
  # library(devtools); install_github("welch-lab/liger", ref = "U_algorithm", force = TRUE)
  require(liger)
  # require(SingleCellExperiment)
  
  # SCE_list is a list of single cell experiment objects
  # that must have a counts assay slot
  # SCE_list = list(X = SCE_X, Y = SCE_Y)
  # counts_list = list(X = readRDS("/Users/ghazan01/Downloads/AADik2b2-Qo3os2QSWXdIAbna/OSMFISH.vin.RDS"),
  #                 Y = readRDS("/Users/ghazan01/Downloads/AADik2b2-Qo3os2QSWXdIAbna/DROPVIZ.vin.RDS"))
  
  # the method takes counts then performs the normalisation
  
  if (!is.null(SCE_list)) {
    counts_list = lapply(SCE_list, assay, "counts")
  }
  
  # saveRDS(counts_list, file = "/Users/ghazan01/Dropbox/Backup/StabMap/StabMap2021/data/counts_list.Rds")
  # saveRDS(counts_list, file = "/Users/ghazan01/Dropbox/Backup/StabMap/StabMap2021/data/counts_list_2.Rds")
  
  # reorder counts by number of features
  counts_list <- counts_list[order(unlist(lapply(counts_list,nrow)))]
  
  # if there are no names, give the list some names
  if (is.null(names(counts_list)[1])) {
    names(counts_list) <- paste0("data_", seq_len(length(counts_list)))
  }
  
  # create liger object
  data_liger = createLiger(counts_list)
  
  # normalise the data
  data_liger <- liger::normalize(data_liger)
  
  # data_liger <- selectGenes(data_liger, unshared = TRUE,
  #                           unshared.datasets = as.list(seq_len(length(counts_list))),
  #                           unshared.thresh = 0.4)
  
  # select features from the unshared set, from all data modalities
  # include flag for which modalities actually have unshared features to begin with
  shared_features = Reduce(intersect, lapply(counts_list, rownames))
  has_unshared = unlist(lapply(counts_list, function(x) any(!rownames(x) %in% shared_features)))
  
  # only use unshared mode if there are unshared features
  unshared = any(has_unshared)
  if (!unshared) {
    message("no unshared features, running liger without UINMF")
  }
  
  # using parameter from vignette
  data_liger <- selectGenes(data_liger, unshared = unshared,
                            unshared.datasets = as.list(which(has_unshared)),
                            unshared.thresh = 0.4)
  
  # rescale the data
  data_liger <- scaleNotCenter(data_liger)
  
  # set the number of dimensions to be bounded by the number of scaled data
  ncomponentsSubset <- min(ncomponentsSubset, ncol(data_liger@scale.data[[1]]))
  
  # perform the UINMF, using parameters from vignette
  data_liger <- optimizeALS(data_liger, k=ncomponentsSubset, use.unshared = unshared)
  
  # quantile normalise
  # use the default of the larger dataset
  data_liger <- quantile_norm(data_liger)
  
  # extract out the H matrix (cells by k (k = number of reduced dimensions))
  H = data_liger@H.norm
  
  return(H)
}

make_an = function(assay) {
  
  # helper function to use in MultiMAP_wrap
  # assay is a matrix-like object with rows being features
  # and columns being cells
  
  require(anndata)
  require(Matrix) # if assay is sparse this is needed
  
  ad <- AnnData(
    X = t(assay),
    var = data.frame(group = rep("", nrow(assay)), row.names = rownames(assay)),
    obs = data.frame(type = rep("", ncol(assay)), row.names = colnames(assay))
  )
  
  return(ad)
}

MultiMAP_wrap = function(assay_list, verbose = FALSE) {
  
  # wrapper function for MultiMAP
  # assay_list is a list of assays with rows being features and columns cells
  # currently restricted to first two entries of assay_list
  # todo - use paste or some other way to allow all list objects to pass through
  # to MultiMAP
  
  # requires Python 3.8 (for a compatible version of numba)
  require(reticulate)
  
  # to view version:
  # py_run_string("from importlib.metadata import version")
  # py_run_string("print(version('MultiMAP'))")

  py_run_string("import scanpy as sc")
  py_run_string("import anndata")
  py_run_string("import MultiMAP")
  
  if (verbose) message("imported python modules")
  
  # generate two random names for global assignment
  if ("ad_1" %in% ls(envir = .GlobalEnv)) {
    ad_1_name = paste0(c("ad_1_", sample(letters, 20, replace = TRUE)), collapse = "")
  } else {
    ad_1_name = "ad_1"
  }
  
  # generate two random names for global assignment
  if ("ad_2" %in% ls(envir = .GlobalEnv)) {
    ad_2_name = paste0(c("ad_2_", sample(letters, 20, replace = TRUE)), collapse = "")
  } else {
    ad_2_name = "ad_2"
  }
  
  # assign these globally (YUCK)
  assign(ad_1_name, make_an(assay_list[[1]]), envir = .GlobalEnv)
  assign(ad_2_name, make_an(assay_list[[2]]), envir = .GlobalEnv)
  
  # ad_1 = make_an(assay_list[[1]])
  # ad_2 = make_an(assay_list[[2]])
  
  if (verbose) message("created anndata objects, assigned to R global environment")
  
  # add the anndata objects into python
  py_run_string(paste0("ad_1 = r.", ad_1_name, ".copy()"))
  py_run_string(paste0("ad_2 = r.", ad_2_name, ".copy()"))
  
  # remove those from the global environment
  rm(list = c(ad_1_name, ad_2_name), envir = .GlobalEnv)
  if (verbose) message("removed anndata objects from R global environment")
  
  # scale and PCA for each object
  py_run_string("sc.pp.scale(ad_1)")
  py_run_string("sc.pp.pca(ad_1)")
  
  # dim(py_eval("ad_1.obsm['X_pca']"))
  
  py_run_string("sc.pp.scale(ad_2)")
  py_run_string("sc.pp.pca(ad_2)")
  
  # dim(py_eval("ad_2.obsm['X_pca']"))
  
  if (verbose) message("performed PCA on anndata objects")
  
  # Perform the MultiMAP:
  py_run_string("adata = MultiMAP.MultiMAP_Integration([ad_1, ad_2], ['X_pca', 'X_pca'])")
  if (verbose) message("successfully performed MultiMAP")
  
  out = py_eval("adata.obsm['X_multimap']")
  obs = py_to_r(py_eval("adata.obs"))
  
  rownames(out) <- rownames(obs)
  
  return(out)
}


# jaccard related functions

embeddingJaccard = function(embedding1, embedding2, k1 = 50, k2 = 50) {
  # given two embeddings, with the same cells
  # calculate the jaccard of the neighbourhoods
  # ... passes to queryNamedKNN()
  
  NN1 = queryNamedKNN(embedding1, embedding1, k = k1)
  NN2 = queryNamedKNN(embedding2, embedding2, k = k2)
  
  jac = neighbourhoodJaccard(NN1, NN2)
  
  return(jac)
}


jaccard = function(x,y) {
  union = length(union(x,y))
  int = length(intersect(x,y))
  return(int/union)
}

neighbourhoodJaccard = function(NN1, NN2) {
  
  # NN1 and NN2 are cells x n1 or n2 matrices, where n is the number of neighbours
  # to consider (they need not be the same)
  # entries are character of the cell within that neighbourhood
  # output is jaccard similarity of these sets
  # typically NN1 and NN2 are output from queryNamedKNN()
  # rownames of NN1 are prioritised
  
  if (nrow(NN1) != nrow(NN2)) stop("NN should have same number of rows")
  
  NN1_list = split.data.frame(NN1, seq_len(nrow(NN1)))
  NN2_list = split.data.frame(NN2, seq_len(nrow(NN2)))
  
  jac = mapply(jaccard, NN1_list, NN2_list)
  names(jac) <- rownames(NN1)
  
  return(jac)
}

