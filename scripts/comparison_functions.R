Seurat_wrap = function(
  SCE_list = NULL,
  counts_list = NULL) {
  
  require(Seurat)
  
  # wrapper function to perform Seurat integration
  
  if (!is.null(SCE_list)) {
    counts_list = lapply(SCE_list, assay, "counts")
  }
  
  seurat.list = sapply(names(counts_list), function(x){
    CreateSeuratObject(counts_list[[x]], project = x)
  }, simplify = FALSE)
  
  # normalize and identify variable features for each dataset independently
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = Inf)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = Inf)
  
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, reference = c(1,2), reduction = "rpca")
  seurat.combined <- IntegrateData(anchorset = seurat.anchors)
  
  seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
  seurat.combined <- RunPCA(seurat.combined, verbose = FALSE)
  
  seurat.PC = seurat.combined@reductions$pca@cell.embeddings
  
  return(seurat.PC)
}


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
  # UPDATE: since UINMF moved to main branch, install using just
  # library(devtools); install_github('welch-lab/liger')
  # for which package will be installed under name "rliger" to not be confused
  # with liger: Lightweight Iterative Geneset Enrichment
  require(rliger)
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
  data_liger <- rliger::normalize(data_liger)
  
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

make_an_MultiVI = function(assay, modalityValue) {
  
  # helper function to use in MultiMAP_wrap
  # assay is a matrix-like object with rows being features
  # and columns being cells
  
  require(anndata)
  require(Matrix) # if assay is sparse this is needed
  
  ad <- AnnData(
    X = t(assay),
    var = data.frame(group = rep("", nrow(assay)), row.names = rownames(assay),
                     modality = modalityValue),
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

cobolt_wrap = function(assay_list, verbose = FALSE) {
  
  # wrapper function for Cobolt
  # python module code in ../../cobolt
  # downloaded from github after trying to install from pip
  # assay_list is a list of assays with rows being features and columns cells
  # currently restricted to first two entries of assay_list
  # assay_list should contain counts, not logcounts
  
  # cobolt wants data matrices to be split by cells and by omic layer
  # so it is easier to split this according to the assay mosaic data 
  # topology than in cobolt itself
  
  # the main comparison is for the current figure 3, where we perform
  # disjoint mosaic data integration with no shared features among
  # the RNA and ATAC modalities. therefore write this wrapper function 
  # to deal with this situation
  
  # first identify which is the multiomic dataset
  g = mosaicDataTopology(assay_list)
  multiomic_data = names(which(degree(g) == 2))
  if (length(multiomic_data) != 1) stop("wrapper function can only handle a single bridge data")
  singleomic_data = names(which(degree(g) == 1))
  
  # then split the multiomic dataset into two
  multiomic_split = sapply(singleomic_data, function(nm){
    assay_list[[multiomic_data]][intersect(rownames(assay_list[[nm]]), rownames(assay_list[[multiomic_data]])),]
  }, simplify = FALSE)
  
  # requires Python 3.8
  require(reticulate)
  
  # to view version:
  # py_run_string("from importlib.metadata import version")
  # py_run_string("print(version('Cobolt'))")
  # for some reason I have to add the local path
  py_run_string("import sys")
  py_run_string("sys.path.append(\"/Users/ghazan01/Dropbox/Backup/StabMap/cobolt/cobolt\")")
  
  py_run_string("from cobolt.utils import SingleData, MultiomicDataset")
  py_run_string("from cobolt.model import Cobolt")
  py_run_string("import os")
  py_run_string("from scipy import io")
  py_run_string("from scipy import sparse")
  py_run_string("import pandas as pd")
  py_run_string("import numpy")
  
  if (verbose) message("imported python modules")
  
  # create the single datas in Cobolt input format
  sapply(singleomic_data, function(nm) {
    rnames_cob <<- rownames(assay_list[[nm]])
    rn_cob <<- "rnames_cob"
    py_run_string(paste0("rnames = r.", rn_cob, ".copy()"))
    py_run_string(paste0("rnames = numpy.asarray(rnames)"))
    
    cnames_cob <<- colnames(assay_list[[nm]])
    cn_cob <<- "cnames_cob"
    py_run_string(paste0("cnames = r.", cn_cob, ".copy()"))
    py_run_string(paste0("cnames = numpy.asarray(cnames)"))
    
    counts_cob <<- t(assay_list[[nm]])
    Cn_cob <<- "counts_cob"
    py_run_string(paste0("counts = r.", Cn_cob, ".copy()"))
    py_run_string(paste0("Scounts = sparse.csr_matrix(counts)"))
    
    py_run_string(paste0(nm, " = SingleData(\"",nm,"\", \"",nm,"\", rnames, Scounts, cnames)"))
  }, simplify = FALSE)
  
  # create the multiomic data in Cobolt input format
  # i.e. two separate omics
  
  sapply(names(multiomic_split), function(nm) {
    rnames_cob <<- rownames(multiomic_split[[nm]])
    rn_cob <<- "rnames_cob"
    py_run_string(paste0("rnames = r.", rn_cob, ".copy()"))
    py_run_string(paste0("rnames = numpy.asarray(rnames)"))
    
    cnames_cob <<- colnames(multiomic_split[[nm]])
    cn_cob <<- "cnames_cob"
    py_run_string(paste0("cnames = r.", cn_cob, ".copy()"))
    py_run_string(paste0("cnames = numpy.asarray(cnames)"))
    
    counts_cob <<- t(multiomic_split[[nm]])
    Cn_cob <<- "counts_cob"
    py_run_string(paste0("counts = r.", Cn_cob, ".copy()"))
    py_run_string(paste0("Scounts = sparse.csr_matrix(counts)"))
    
    py_run_string(paste0(multiomic_data, "_", nm, " = SingleData(\"",nm,"\", \"",multiomic_data,"\", rnames, Scounts, cnames)"))
  }, simplify = FALSE)
  
  # now put together multiomic data in cobolt format
  
  if (verbose) message("Single omic data objects created")

  py_run_string(paste0("multi_dt = MultiomicDataset.from_singledata(",
                       paste0(singleomic_data, collapse = ", "),
                       ", ",
                       paste0(multiomic_data, "_", names(multiomic_split), collapse = ", "),
                       ")"))
  py_run_string("print(multi_dt)")
  
  if (verbose) message("Multiomic data object created")
  
  # now run Cobolt
  
  # following the vignette instructions
  
  py_run_string("model = Cobolt(dataset=multi_dt, lr=0.005, n_latent=10)")
  py_run_string("model.train(num_epochs=100)")
  
  # model finished!
  if (verbose) message("Cobolt model fit completed")
  
  # calculate all corrected variables
  py_run_string("model.calc_all_latent()")
  
  # access the corrected latent variables
  py_run_string("latent = model.get_all_latent()")
  
  # no need to get all the latent variables
  # py_run_string("latent_raw = model.get_all_latent(correction=False)")
  
  if (verbose) message("successfully performed Cobolt")
  
  out = py_eval("latent")
  
  cobolt = out[[1]]
  
  cellnames = gsub(".*~", "", out[[2]])
  
  rownames(cobolt) <- cellnames
  
  return(cobolt)
}


cobolt2_wrap = function(assay_list, verbose = FALSE) {
  
  # wrapper function for Cobolt with two input datasets
  # python module code in ../../cobolt
  # downloaded from github after trying to install from pip
  # assay_list is a list of assays with rows being features and columns cells
  # currently restricted to first two entries of assay_list
  # assay_list should contain counts, not logcounts
  
  # multiomic_data is a character name of the assay that is multiomic
  
  # cobolt wants data matrices to be split by cells and by omic layer
  # so it is easier to split this according to the assay mosaic data 
  # topology than in cobolt itself
  
  # the main comparison is for the current figure 3, where we perform
  # disjoint mosaic data integration with no shared features among
  # the RNA and ATAC modalities. therefore write this wrapper function 
  # to deal with this situation
  
  # if (length(assay_list) != 2) stop("wrapper function can only handle two input assays")
  
  # first identify which is the multiomic dataset
  g = mosaicDataTopology(assay_list)
  # multiomic_data = names(which(degree(g) == 2))
  # multiomic_data = "Multiome"
  
  # take every feature and assign it to the combination of data that it is in
  allfeatures = Reduce(union, lapply(assay_list, rownames))
  mem = do.call(cbind,lapply(lapply(assay_list, rownames), function(x) allfeatures %in% x))
  memp = apply(mem, 1, function(x) paste0(colnames(mem)[x], collapse = "_"))
  names(memp) <- allfeatures
  
  # split all the assays based on memp
  assay_list_split = lapply(assay_list, function(assay){
    split.data.frame(assay, memp[rownames(assay)])
  })
  
  # requires Python 3.8
  require(reticulate)
  
  # to view version:
  # py_run_string("from importlib.metadata import version")
  # py_run_string("print(version('Cobolt'))")
  # for some reason I have to add the local path
  py_run_string("import sys")
  py_run_string("sys.path.append(\"/Users/ghazan01/Dropbox/Backup/StabMap/cobolt/cobolt\")")
  
  py_run_string("from cobolt.utils import SingleData, MultiomicDataset")
  py_run_string("from cobolt.model import Cobolt")
  py_run_string("import os")
  py_run_string("from scipy import io")
  py_run_string("from scipy import sparse")
  py_run_string("import pandas as pd")
  py_run_string("import numpy")
  
  if (verbose) message("imported python modules")
  
  # create the single datas in Cobolt input format
  
  singleomicnames_raw = sapply(names(assay_list_split), function(dataset) {
    
    sapply(names(assay_list_split[[dataset]]), function(nm) {
      
      rnames_cob <<- rownames(assay_list_split[[dataset]][[nm]])
      rn_cob <<- "rnames_cob"
      py_run_string(paste0("rnames = r.", rn_cob, ".copy()"))
      py_run_string(paste0("rnames = numpy.asarray(rnames)"))
      
      cnames_cob <<- colnames(assay_list_split[[dataset]][[nm]])
      cn_cob <<- "cnames_cob"
      py_run_string(paste0("cnames = r.", cn_cob, ".copy()"))
      py_run_string(paste0("cnames = numpy.asarray(cnames)"))
      
      counts_cob <<- t(assay_list_split[[dataset]][[nm]])
      Cn_cob <<- "counts_cob"
      py_run_string(paste0("counts = r.", Cn_cob, ".copy()"))
      py_run_string(paste0("Scounts = sparse.csr_matrix(counts)"))
      
      print(paste0(dataset, "_", nm))
      
      py_run_string(paste0(dataset, "_", nm, " = SingleData(\"",nm,"\", \"",dataset,"\", rnames, Scounts, cnames)"))
      
      return(paste0(dataset, "_", nm))
    }, simplify = FALSE)
    
  }, simplify = FALSE)
  singleomicnames <- unlist(singleomicnames_raw, recursive = TRUE)
  
  # now put together multiomic data in cobolt format
  
  if (verbose) message("Single omic data objects created")
  
  py_run_string(paste0("multi_dt = MultiomicDataset.from_singledata(",
                       paste0(singleomicnames, collapse = ", "),
                       ")"))
  py_run_string("print(multi_dt)")
  
  if (verbose) message("Multiomic data object created")
  
  # now run Cobolt
  
  # following the vignette instructions
  py_run_string("model = Cobolt(dataset=multi_dt, lr=0.005, n_latent=10)")
  py_run_string("model.train(num_epochs=100)")
  
  # model finished!
  if (verbose) message("Cobolt model fit completed")
  
  # calculate all corrected variables
  py_run_string("model.calc_all_latent()")
  
  # access the corrected latent variables
  py_run_string("latent = model.get_all_latent()")
  
  # no need to get all the latent variables
  # py_run_string("latent_raw = model.get_all_latent(correction=False)")
  
  if (verbose) message("successfully performed Cobolt")
  
  out = py_eval("latent")
  
  cobolt = out[[1]]
  
  cellnames = gsub(".*~", "", out[[2]])
  
  rownames(cobolt) <- cellnames
  
  return(cobolt)
}


multiVI_wrap = function(assay_list, verbose = FALSE) {
  
  # MultiVI wants count data
  # wrapper function for MultiVI
  
  # MultiVI seems to be extremely specific to the task of integrating Multiome
  # with ATAC and/or RNA-seq
  # unsure what you'd need to do if not in this analysis situation
  # Very important!!!!!!
  # especially needs to be in the order Multiome, RNA and chromatin
  # else arguments are to be specified by name
  
  # assay_list is a list of assays with rows being features and columns cells
  # currently restricted to first two entries of assay_list
  
  # requires Python 3.8 (for a compatible version of numba)
  require(reticulate)
  
  message("MultiVI requires data input as Multiome, RNA and Chromatin, in that order!")
  
  # to view version:
  py_run_string("import scvi")
  py_run_string("import anndata")
  py_run_string("import numpy as np")
  py_run_string("import scanpy as sc")
  
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
  
  # generate two random names for global assignment
  if ("ad_3" %in% ls(envir = .GlobalEnv)) {
    ad_3_name = paste0(c("ad_3_", sample(letters, 20, replace = TRUE)), collapse = "")
  } else {
    ad_3_name = "ad_3"
  }
  
  # assign these globally (YUCK)
  assign(ad_1_name, make_an_MultiVI(assay_list[[1]], ifelse(rownames(assay_list[[1]]) %in% rownames(assay_list[[2]]),
                                                            "Gene Expression", "Peaks")), envir = .GlobalEnv)
  assign(ad_2_name, make_an_MultiVI(assay_list[[2]], "Gene Expression"), envir = .GlobalEnv)
  assign(ad_3_name, make_an_MultiVI(assay_list[[3]], "Peaks"), envir = .GlobalEnv)
  
  # ad_1 = make_an(assay_list[[1]])
  # ad_2 = make_an(assay_list[[2]])
  
  if (verbose) message("created anndata objects, assigned to R global environment")
  
  # add the anndata objects into python
  py_run_string(paste0("ad_1 = r.", ad_1_name, ".copy()"))
  py_run_string(paste0("ad_2 = r.", ad_2_name, ".copy()"))
  py_run_string(paste0("ad_3 = r.", ad_3_name, ".copy()"))
  
  # remove those from the global environment
  rm(list = c(ad_1_name, ad_2_name, ad_3_name), envir = .GlobalEnv)
  if (verbose) message("removed anndata objects from R global environment")
  
  # concatenate annadata objects
  # We can now use the organizing method from scvi to concatenate these anndata
  py_run_string("adata_mvi = scvi.data.organize_multiome_anndatas(ad_1, ad_2, ad_3)")
  
  # apparently MultiVI requires features to be ordered such that genes appear
  # before genomic regions
  py_run_string("adata_mvi = adata_mvi[:, adata_mvi.var[\"modality\"].argsort()].copy()")
  # py_run_string("print(adata_mvi.var)")
  
  # filter fewer than 1%
  py_run_string("sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))")
  # py_run_string("print(adata_mvi.shape)")
  
  # make modality the batch annotation
  py_run_string("scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key=\"modality\")")

  # determine architecture for each modality 
  # sensitive to indentation!
  py_run_string("mvi = scvi.model.MULTIVI(
    adata_mvi,
    n_genes=(adata_mvi.var[\"modality\"]==\"Gene Expression\").sum(),
    n_regions=(adata_mvi.var[\"modality\"]==\"Peaks\").sum(),
)
mvi.view_anndata_setup()")
  
  # train the MultiVI model
  py_run_string("mvi.train()")
  
  # extract the latent space
  py_run_string("adata_mvi.obsm[\"MultiVI_latent\"] = mvi.get_latent_representation()")
  # py_run_string("sc.pp.neighbors(adata_mvi, use_rep=\"MultiVI_latent\")")
  # py_run_string("sc.tl.umap(adata_mvi, min_dist=0.2)")
  # py_run_string("sc.pl.umap(adata_mvi, color='modality')")
  
  out = py_eval("adata_mvi.obsm[\"MultiVI_latent\"]")
  obs = py_to_r(py_eval("adata_mvi.obs"))
  
  # MultiVI appends a data type to the rownames, so use regex to remove
  rownames(out) <- gsub("_.*", "", rownames(obs))
  
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

