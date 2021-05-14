stabMapNonOverlapping = function(M_R,
                                 M_A,
                                 R,
                                 A,
                                 l_GP = 50,
                                 l_G = 50,
                                 l_P = 50) {
    
    # to-do: incorporate labels via LD components
    # to-do: loosen restriction on dataset sizes
    # to-do: think about labels on the R or A modalities
    
    # perform stabmap when there is some data that is completely non-overlapping (or disjoint)
    # current set-up, there are three datasets
    # with the prime example in mind:
    # M 10X multiome (RNA + ATAC) data with two submatrices M_R and M_A, of size GxN_m and PxN_m respectively
    # R scRNA-seq data of size GxN_r
    # A ATAC-seq data of size PxN_a
    # for now, assume that the dataset sizes are matched between the modalities
    
    # generate some example data
    # G = 500 # number of genes
    # P = 1000 # number of atac peak
    # N_m = 500 # number of multimodal cells
    # N_r = 100 # number of RNA cells
    # N_a = 200 # number of ATAC cells
    # 
    # M_R = matrix(rnorm(N_m*G), G, N_m)
    # M_A = matrix(rnorm(N_m*P), P, N_m)
    # colnames(M_R) <- colnames(M_A) <- paste0("M_cell_", 1:ncol(M_R))
    # rownames(M_R) <- paste0("RNA_", 1:nrow(M_R))
    # rownames(M_A) <- paste0("ATAC_", 1:nrow(M_A))
    # 
    # R = matrix(rnorm(N_r*G), G, N_r)
    # rownames(R) <- rownames(M_R)
    # colnames(R) <- paste0("R_cell_", 1:ncol(R))
    # 
    # A = matrix(rnorm(N_a*P), P, N_a)
    # rownames(A) <- rownames(M_A)
    # colnames(A) <- paste0("A_cell_", 1:ncol(A))
    
    M = rbind(M_R, M_A)
    
    # generate three PC loadings from the data
    # PC of the full data
    PC_GP = calculatePCA(M, ncomponents = l_GP)
    L_GP = attr(PC_GP, "rotation")
    
    # PC of the RNA
    PC_G = calculatePCA(M_R, ncomponents = l_G)
    L_G = attr(PC_G, "rotation")
    
    # PC of the ATAC
    PC_P = calculatePCA(M_A, ncomponents = l_P)
    L_P = attr(PC_P, "rotation")
    
    # generate two sets of regression coefficients
    # RNA PCs onto the joint PCs
    B_R = lm.fit(PC_G, PC_GP)[["coefficients"]]
    
    # ATAC PCs onto the joint PCs
    B_A = lm.fit(PC_P, PC_GP)[["coefficients"]]
    
    # project the RNA data into the joint space
    W_r = t(R[rownames(L_G),]) %*% L_G %*% B_R
    
    # project the ATAC data into the joint space
    W_a = t(A[rownames(L_P),]) %*% L_P %*% B_A
    
    # project the multimodal data into the joint space
    # this is just the same as extracting the scores
    # W_m = t(M[rownames(L_GP),]) %*% L_GP
    W_m = PC_GP
    
    # bind these cells together
    W = rbind(W_r, W_m, W_a)
    
    return(W)
}


stabMapSubsetResidual = function(
    SCE,
    genes,
    assayName = "logcounts",
    ncomponentsFull = 50,
    ncomponentsSubset = 50,
    sparse = FALSE) {
    
    # given a dataset and a subset of genes,
    # calculate the error matrix of PC values between 
    # the two sets of features
    
    # TODO: get the residual from LDA grouping (does this make sense?)
    
    PCs = list(calculatePCA(assay(SCE, assayName), ncomponents = ncomponentsFull))
    PCs_genes = list(calculatePCA(assay(SCE, assayName), ncomponents = ncomponentsSubset, subset_row = genes))
    
    PCs_genes_loadings = attr(PCs_genes[[1]], "rotation")
    
    if (sparse) {
        # Fit lasso models with lambda selected with CV to estimate the best 
        # linear combination of PCs among the subset features
        sparseCoef = function(PCs_X, PCs_X_genes) {
            coefList_X = list()
            for (i in seq_len(ncol(PCs_X))) {
                # print(i)
                fit.cv = cv.glmnet(x = PCs_X_genes, y = PCs_X[,i], intercept = FALSE)
                fit = glmnet(x = PCs_X_genes, y = PCs_X[,i],
                             lambda = fit.cv["lambda.min"][[1]], intercept = FALSE)
                coefList_X[[i]] <- coef(fit)[-1]
            }
            coef_X = do.call(cbind,coefList_X)
            return(coef_X)
        }
        
        coefs = mapply(sparseCoef, PCs, PCs_genes, SIMPLIFY = FALSE)
    } else {
        coefs = lapply(mapply(lm.fit, PCs_genes, PCs, SIMPLIFY = FALSE), "[[", "coefficients")
    }
    
    all_assay = assay(SCE, assayName)[genes,]
    
    matMult = function(m1_rnames, m2, m3, m1) {
        t(m1[m1_rnames,]) %*% m2 %*% m3 
    }
    
    emb = matMult(m1_rnames = rownames(PCs_genes_loadings),
                  m2 = PCs_genes_loadings,
                  m3 = coefs[[1]],
                  m1 = all_assay)
    
    # fitted minus true values
    res = emb - PCs[[1]]
    
    return(res)
}


stabMapComparative = function(
    # SCE_X,
    # SCE_Y,
    # SCE_list,
    assay_list,
    # genes = intersect(rownames(SCE_X), rownames(SCE_Y)),
    features = Reduce(intersect, lapply(assay_list, rownames)),
    # assayNames = rep("logcounts", length(assay_list)),
    stabilise = rep(TRUE, length(assay_list)),
    # assayNameX = "logcounts",
    # assayNameY = "logcounts",
    ncomponentsFull = 50,
    ncomponentsSubset = 50,
    sparse = FALSE) {
    
    # assay_list is a list of assays (rows features columns cells)
    # features are by default the intersection across all of the assay_list objects
    # stabilise is a logical vector indicating if stabilised components should
    # be output for each of the assay_list objects, e.g. if assay_list corresponds
    # to a reference and a query i.e. list(reference_assay, query_assay), then 
    # stabilise should be c(TRUE, FALSE)
    
    # TODO: lookup table for features that should be considered joint
    # between the various data (cross-species context)
    # TODO: set maximum PCs based on features numbers
    
    # Given two (or more) datasets without a set reference vs query structure
    # (i.e. comparative scenario), perform stabMap
    
    require(scater)
    if (sparse) require(glmnet)
    
    # if some object with a dimension (usually a assay) is given, 
    # place it into a list with a single object
    if (!is.null(dim(assay_list))) {
        assay_list <- list(assay_list)
    }
    
    # PCs = lapply(mapply(assay, SCE_list, assayNames), calculatePCA, ncomponents = ncomponentsFull)
    # PCs_genes = lapply(mapply(assay, SCE_list, assayNames), calculatePCA, ncomponents = ncomponentsSubset, subset_row = genes)
    
    # PCs from full feature space per dataset
    PCs = lapply(assay_list, calculatePCA, ncomponents = ncomponentsFull)
    PCs_genes = lapply(assay_list, calculatePCA, ncomponents = ncomponentsSubset, subset_row = features)
    
    PCs_genes_loadings = lapply(PCs_genes, attr, "rotation")
    
    if (sparse) {
        # Fit lasso models with lambda selected with CV to estimate the best 
        # linear combination of PCs among the subset features
        
        coefs = mapply(sparseCoef, PCs, PCs_genes, SIMPLIFY = FALSE)
        
        # coefList_Y = list()
        # for (i in seq_len(ncol(PCs_Y))) {
        #     print(i)
        #     fit.cv = cv.glmnet(x = PCs_Y_genes, y = PCs_Y[,i], intercept = FALSE)
        #     fit = glmnet(x = PCs_Y_genes, y = PCs_Y[,i],
        #                  lambda = fit.cv["lambda.min"][[1]], intercept = FALSE)
        #     coefList_Y[[i]] <- coef(fit)[-1]
        # }
        # coef_Y = do.call(cbind,coefList_Y)
    } else {
        
        coefs = lapply(mapply(lm.fit, PCs_genes, PCs, SIMPLIFY = FALSE), "[[", "coefficients")
        
        # fit = lm.fit(x = PCs_X_genes, y = PCs_X)
        # coef_X = fit$coefficients
        # 
        # fit = lm.fit(x = PCs_Y_genes, y = PCs_Y)
        # coef_Y = fit$coefficients
    }
    
    all_assay = do.call(cbind, lapply(assay_list, "[", features, ))
    # all_assay = do.call(cbind, lapply(mapply(assay, SCE_list, assayNames), "[", features, ))
    
    # all_assay = cbind(assay(SCE_X, assayNameX)[genes,],
    #                   assay(SCE_Y, assayNameY)[genes,])
    
    # all_assay_sub = lapply(lapply(PCs_genes_loadings, rownames))
    # all_assay_sub = mapply("[", all_assay, lapply(PCs_genes_loadings, rownames), MoreArgs = list())
    
    # all_assay_sub_X = all_assay[rownames(PCs_X_genes_loadings), ]
    # all_assay_sub_Y = all_assay[rownames(PCs_Y_genes_loadings), ]
    
    # d_sub = d[rownames(referencePCs_genes_loadings),]
    
    matMult = function(m1_rnames, m2, m3, m1) {
        t(m1[m1_rnames,]) %*% m2 %*% m3 
    }
    
    emb_list = mapply(matMult,
                      m1_rnames = lapply(PCs_genes_loadings, rownames),
                      m2 = PCs_genes_loadings,
                      m3 = coefs,
                      MoreArgs = list(m1 = all_assay))
    
    # emb_X = t(all_assay_sub_X) %*% PCs_X_genes_loadings %*% coef_X
    # emb_Y = t(all_assay_sub_Y) %*% PCs_Y_genes_loadings %*% coef_Y
    
    # emb = cbind(emb_X, emb_Y)
    
    emb = do.call(cbind, emb_list[stabilise])
    
    return(emb)
    
}

sparseCoef = function(PCs_X, PCs_X_genes) {
    require(glmnet)
    coefList_X = list()
    for (i in seq_len(ncol(PCs_X))) {
        # print(i)
        fit.cv = cv.glmnet(x = PCs_X_genes, y = PCs_X[,i], intercept = FALSE)
        fit = glmnet(x = PCs_X_genes, y = PCs_X[,i],
                     lambda = fit.cv["lambda.min"][[1]], intercept = FALSE)
        coefList_X[[i]] <- coef(fit)[-1]
    }
    coef_X = do.call(cbind,coefList_X)
    return(coef_X)
}

stabMapContinuous = function(referenceSCE,
                             querySCE,
                             genes = intersect(rownames(querySCE), rownames(referenceSCE)),
                             assayNameReference = "logcounts",
                             assayNameQuery = "logcounts",
                             ncomponentsFull = 50,
                             ncomponentsSubset = 50,
                             sparse = TRUE) {
    
    # stabMap using continuous information between reference and query data
    
    if (sparse) require(glmnet)
    
    referencePCs = calculatePCA(assay(referenceSCE, assayNameReference), ncomponents = ncomponentsFull)
    referencePCs_genes = calculatePCA(assay(referenceSCE, assayNameReference)[genes,],
                                      ncomponents = ncomponentsSubset)
    
    # loadings are saved in the attributes
    referencePCs_genes_loadings = attr(referencePCs_genes, "rotation")
    
    if (sparse) {
        # Fit lasso models with lambda selected with CV to estimate the best 
        # linear combination of PCs among the subset features
        coefList = list()
        # i = 1
        for (i in seq_len(ncol(referencePCs))) {
            print(i)
            fit.cv = cv.glmnet(x = referencePCs_genes, y = referencePCs[,i], intercept = FALSE)
            fit = glmnet(x = referencePCs_genes, y = referencePCs[,i],
                         lambda = fit.cv["lambda.min"][[1]], intercept = FALSE)
            coefList[[i]] <- coef(fit)[-1]
        }
        coef = do.call(cbind,coefList)
    } else {
        # fit = lm.fit(x = referencePCs_genes, y = referencePCs[,i])
        # coefList[[i]] <- fit$coefficients
        fit = lm.fit(x = referencePCs_genes, y = referencePCs)
        coef = fit$coefficients
    }
    
    if (is.null(querySCE)) {
        all_assay = cbind(assay(referenceSCE, assayNameReference)[genes,])
    } else {
        all_assay = cbind(assay(referenceSCE, assayNameReference)[genes,],
                          assay(querySCE, assayNameQuery)[genes,])
    }
    
    all_assay_sub = all_assay[rownames(referencePCs_genes_loadings), ]
    
    # d_sub = d[rownames(referencePCs_genes_loadings),]
    
    emb = t(all_assay_sub) %*% referencePCs_genes_loadings %*% coef
    
    return(emb)
    
}


stabMapLabelled = function(
    referenceSCE,
    # querySCE,
    # reference_assay,
    query_assay,
    features = intersect(rownames(referenceSCE), rownames(query_assay)),
    labels = "celltype",
    assayNameReference = "logcounts",
    # assayNameQuery = "logcounts",
    prop_explained = 1,
    selectionLFC = ifelse(length(features) >= 1000, 0.1, 0),
    selectBestLDA = TRUE,
    ...) {
    
    # TODO: allow a list of querySCE objects?
    # TODO: allow sparse matrix in selectionLFC
    
    # referenceSCE is a SingleCellExperiment object containing 
    # an assay of assayNameReference, and a colData factor of grouping
    # query_assay is a features x cells query matrix
    # features are those which should be considered
    # labels is the labels to identify LDs
    # assayNameReference is the assay name to extract from referenceSCE
    # prop_explained is the proportion of explained variation from the LDs
    # selectionLFC is the minimum absolute LFC to include genes for inclusion
    # in LDs
    # selectBestLDA is a logical whether a cross validation step should
    # be run to avoid overfitting
    # ... passed to lda_best
    
    require(Matrix)
    require(MASS)
    
    # selectionLFC is a log-fold change to use for selecting features among genes first
    
    if (selectionLFC != 0) {
        
        print("performing fold-change feature selection")
        
        sce = SingleCellExperiment(assays = list("expr" = assay(referenceSCE, assayNameReference)[features,]))
        colData(sce) <- colData(referenceSCE)
        
        # make assay dense (to fix)
        assay(sce, "expr") <- as.matrix(assay(sce, "expr"))
        
        # only compare against groupings that have at least 10 cells
        sub = names(which(table(colData(sce)[,labels]) >= 10))
        
        sel = getMarkerArray(sce,
                             assayName = "expr",
                             group = labels,
                             subset = colData(sce)[,labels] %in% sub,
                             verbose = FALSE)
        
        # only keep genes with at least selectionLFC log fold change
        topLFC = rowMax(abs(sel[,,"LFC"]))
        
        features_sel = features[topLFC >= selectionLFC]
        
        features <- features_sel
        
        message(length(features), " features with large enough LFC")
        
    }
    
    # prefilter cells and genes if needed
    # v = fac2sparse(colData(referenceSCE)[,grouping]) %*% t(assay(referenceSCE, assayNameReference)[genes,])
    # vars = apply(assay(referenceSCE, assayNameReference)[genes,], 1, function(x)
    # max(tapply(x, colData(referenceSCE)[,grouping], var), na.rm = TRUE))
    # genes <- genes[vars > 0]
    
    vars = rowMaxs(apply(fac2sparse(colData(referenceSCE)[,labels]), 1, function(x)
        rowWeightedVars(assay(referenceSCE, assayNameReference)[features,],
                        x)), na.rm = TRUE)
    if (any(vars == 0)) warning ("removing genes with zero intra-class variance")
    features <- features[vars > 0]
    
    # build LDA model for genes
    data_train = t(assay(referenceSCE, assayNameReference)[features,])
    labels_train = colData(referenceSCE)[,labels]
    
    # remove cells with missing labels
    if (selectBestLDA) {
        fit = lda_best(data = data_train[!is.na(labels_train),],
                       labels = labels_train[!is.na(labels_train)], ...)
    } else {
        fit = lda(data_train[!is.na(labels_train),],
                  grouping = labels_train[!is.na(labels_train)])
    }
    
    # proportion of between group variation explained by the linear
    # discriminants
    prop_rel = fit$svd^2/sum(fit$svd^2)
    n_prop_rel = which(cumsum(prop_rel) >= prop_explained)[1]
    
    if (is.null(query_assay)) {
        all_assay = assay(referenceSCE, assayNameReference)[features,]
    } else {
        all_assay = cbind(assay(referenceSCE, assayNameReference)[features,],
                          query_assay[features,])
    }
    
    resub = predict(object = fit, newdata = t(all_assay[rownames(fit$scaling),]))
    
    if (prop_explained == 1) {
        resub_out = resub$x
    } else {
        resub_out = resub$x[,seq_len(n_prop_rel)]  
    }
    return(resub_out)
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
    
    # TODO: add stop when features not appearing in data etc
    
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

reducedMNN_batchFactor = function(embedding,
                                  batchFactor,
                                  ...) {
    # batch correct within this embedding, wrapper around reducedMNN
    # batchFactor is a named vector that is matched
    # ... passed to reducedMNN
    
    require(batchelor)
    
    # for some reason reducedMNN only takes dense matrices
    embedding <- as.matrix(embedding)
    
    batchFactor_used = batchFactor[rownames(embedding)]
    
    out = reducedMNN(embedding, batch = batchFactor_used, ...)
    resub_corrected = out$corrected
    
    return(resub_corrected)
}

ComBat_batchFactor = function(embedding,
                              batchFactor,
                              ...) {
    # batch correct within this embedding, wrapper around sva's ComBat
    # batchFactor is a named vector that is matched
    # ... passed to reducedMNN
    
    require(sva)
    
    batchFactor_used = batchFactor[rownames(embedding)]
    
    out = ComBat(t(embedding), batch = batchFactor_used, ...)
    resub_corrected = t(out)
    
    return(resub_corrected)
}

Harmony_batchFactor = function(embedding,
                              batchFactor, ...) {
    # batch correct within this embedding, wrapper around sva's ComBat
    # batchFactor is a named vector that is matched
    # ... passed to reducedMNN
    
    require(harmony)
    
    batchFactor_used = batchFactor[rownames(embedding)]
    
    out = HarmonyMatrix(t(embedding), batchFactor_used, do_pca = FALSE, ...)
    resub_corrected = t(out)
    
    return(resub_corrected)
}

batchCorrect_batchFactor = function(embedding,
                                    batchFactor,
                                    PARAM) {
    # batch correct within this embedding, wrapper around sva's ComBat
    # batchFactor is a named vector that is matched
    # PARAM passed to batchCorrect
    
    require(batchelor)
    
    batchFactor_used = batchFactor[rownames(embedding)]
    
    out = batchCorrect(t(embedding), batch = batchFactor_used, PARAM = PARAM)
    resub_corrected = t(assays(out)[[1]])
    
    return(resub_corrected)
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

embeddingJaccard = function(embedding1, embedding2, k1 = 50, k2 = 50) {
    # given two embeddings, with the same cells
    # calculate the jaccard of the neighbourhoods
    # ... passes to queryNamedKNN()
    
    NN1 = queryNamedKNN(embedding1, embedding1, k = k1)
    NN2 = queryNamedKNN(embedding2, embedding2, k = k2)
    
    jac = neighbourhoodJaccard(NN1, NN2)
    
    return(jac)
}








vectorSubset = function(vec, mat){
    # vec is a named vector
    # mat is a matrix containing the names or indices for which you want
    # to get the entries of vec
    
    vmat = c(mat)
    vvec = vec[vmat]
    
    vecmat = matrix(vvec, nrow = nrow(mat), ncol = ncol(mat))
    colnames(vecmat) <- colnames(mat)
    rownames(vecmat) <- rownames(mat)
    
    return(vecmat)
}

vectorMatch = function(vec, mat, vecnames){
    # vec is an unnamed vector
    # vecnames is the names of vec
    # mat is a matrix containing the names or indices for which you want
    # to get the entries of vec, matching vecnames
    
    vmat = c(mat)
    
    vecind = match(vmat,vecnames)
    
    vvec = vec[vecind]
    
    vecmat = matrix(vvec, nrow = nrow(mat), ncol = ncol(mat))
    colnames(vecmat) <- colnames(mat)
    rownames(vecmat) <- rownames(mat)
    
    return(vecmat)
}


calculateUMAP_rnames = function(embedding, ...) {
    # pass to calculateUMAP from scater
    # but keep rownames
    require(scater)
    message("calculating UMAP...")
    embedding_UMAP = calculateUMAP(t(embedding), ...)
    rownames(embedding_UMAP) <- rownames(embedding)
    return(embedding_UMAP)
}

# not used currently
get_umap_graph <- function(umap_output) {
    require(Matrix)
    require(igraph)
    # extracts the graph information and outputs a weighted igraph
    dists.knn <- umap_output$knn[["distances"]]
    indx.knn <- umap_output$knn[["indexes"]]
    m.adj <- Matrix(0, nrow=nrow(indx.knn), ncol=nrow(indx.knn), sparse=TRUE) 
    rownames(m.adj) <- colnames(m.adj) <- rownames(indx.knn)
    for (i in seq_len(nrow(m.adj))) {
        m.adj[i,rownames(indx.knn)[indx.knn[i,]]] <- dists.knn[i,] 
    }
    igr <- graph_from_adjacency_matrix(m.adj, weighted=TRUE)
    return(igr)
}


getMarkerArray = function(sce, 
                          group, 
                          assayName = "logcounts", 
                          block = NULL, 
                          subset = NULL,
                          verbose = TRUE, 
                          pseudobulk = FALSE) {
    require(scran)
    
    # sce is a SingleCellExperiment object
    # group is a character for the group to split by
    # block is a character vector of the block to use
    # subset is a logical of the same length as ncol(sce) for which cells to 
    # use for testing
    # verbose if TRUE 
    # pseudobulk whether you perform pseudobulk by averaging values first
    
    if (!is.null(subset)) {
        sce_sub = sce[,subset]
    } else {
        sce_sub <- sce
    }
    
    genes = rownames(sce_sub)
    groups = sort(unique(colData(sce_sub)[, group]))
    
    if (!is.null(block)) {
        blockFactor = getInteractionFactor(sce_sub, levels = block)
    } else {
        blockFactor = NULL
    }
    
    marker_array = array(NA, dim = c(nrow(sce_sub), length(groups), 5),
                         dimnames = list(
                             rownames(sce_sub),
                             groups,
                             NULL
                         ))
    
    # per group level, output a matrix of genes and test outputs
    
    exprs = assay(sce_sub, assayName)
    
    for (grouplevel in groups) {
        
        if (verbose) print(grouplevel)
        
        grps = colData(sce_sub)[, group] == grouplevel
        
        if (pseudobulk) {
            
            colData(sce_sub)[, "grouplevel"] <- grps
            # then replace SCE with the pseudobulk one
            sce_pseudo <- sumCountsAcrossCells(assay(sce_sub, assayName),
                                               ids = colData(sce_sub)[,c("grouplevel", "grouplevel", block)],
                                               average = TRUE)
            # assayName <- "sum"
            blockFactor <- NULL
            
            exprs <- assay(sce_pseudo, "sum")
            grps = colData(sce_pseudo)[,"grouplevel"]
        }
        
        out = findMarkers(exprs, 
                          groups = grps,
                          block = blockFactor,
                          test.type = "t",
                          direction = "any",
                          log.p = FALSE,
                          sorted = FALSE)
        out_TRUE <- out[["TRUE"]]
        
        meanIn = rowMeans(exprs[, grps, drop = FALSE])
        meanOut = rowMeans(exprs[, !grps, drop = FALSE])
        
        out_mat = cbind("p.value" = out_TRUE[, "p.value"],
                        "FDR" = out_TRUE[, "FDR"],
                        "LFC" = out_TRUE[, "logFC.FALSE"],
                        "meanIn" = meanIn,
                        "meanOut" = meanOut)
        
        # markerList[[grouplevel]] <- 
        marker_array[,grouplevel,] <- out_mat
        
    }
    
    dimnames(marker_array)[[3]] <- colnames(out_mat)
    
    return(marker_array)
    
    # e.g. marker_array[,"Gut",] to see all genes for this grouplevel
    # e.g. marker_array["Ttn",,] to see all p-values for this gene
}

markerArrayPlot = function(marker_array,
                           grouplevel,
                           otherGroupLabel = NULL,
                           LFC = 0.5,
                           p.value = 0.05,
                           FDR = NULL,
                           onlyLabelTopRanked = NULL) {
    # marker_array is the output from getMarkerArray
    # grouplevel is the 2nd dimension name to graph
    # onlyLabelTopRanked can be NULL, FALSE, or a number
    
    require(ggplot2)
    require(ggrepel)
    
    mk = as.data.frame(marker_array[,grouplevel,])
    if (is.null(FDR)) {
        mk$sig <- mk$p.value < p.value & abs(mk$LFC) > LFC
    } else {
        mk$sig <- mk$FDR < FDR & abs(mk$LFC) > LFC
    }
    mk$p.valueLFC <- mk$p.value
    mk$p.valueLFC[abs(mk$LFC) < LFC] <- 1
    mk$rank <- rank(mk$p.valueLFC)
    mk$lab <- rownames(mk)
    mk$gene <- rownames(mk)
    if (is.numeric(onlyLabelTopRanked)) {
        mk$lab[mk$rank > onlyLabelTopRanked] <- NA
    } else {
        mk$lab[!mk$sig] <- NA
    }
    
    if (!is.null(otherGroupLabel)) {
        gtitle = paste0("High LFC indicates higher expression in ", grouplevel,
                        ", compared to ", otherGroupLabel)
    } else {
        gtitle = paste0("High LFC indicates higher expression in ", grouplevel)
    }
    
    if (is.null(FDR)) {
        ltitle = paste0("Significantly DE (P-value < ", p.value,")")
    } else {
        ltitle = paste0("Significantly DE (FDR-adjusted P-value < ", FDR,")")    
    }
    
    g = ggplot(mk, aes(x = LFC, y = -log10(p.value), label = gene)) + 
        geom_point(aes(colour = sig)) + 
        theme_classic() +
        theme(legend.position = "bottom") +
        scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey")) +
        geom_vline(xintercept = 0, linetype = "dotted") +
        ylab("-log10(P-value)") +
        xlab("Log Fold Change") +
        ggtitle(gtitle) +
        guides(colour = guide_legend(title = ltitle,
                                     override.aes = list(size = 5))) +
        NULL
    
    if (!isFALSE(onlyLabelTopRanked)) {
        g = g + geom_text_repel(aes(label = lab))
    }
    
    if (is.null(FDR)) {
        g = g + geom_hline(yintercept = -log10(p.value), linetype = "dotted")
    }
    
    return(g)
}

lda_best = function(data, labels, nFold = 10, testPct = 0.2) {
    
    require(MASS)
    
    # with many features LDA can overfit to the training data
    # perform LDA repeatedly, and calculate the cross
    # validation error rate
    # return the lda object with the best accuracy
    
    fitList = replicate(nFold, {
        train = sample(c(TRUE,FALSE), length(labels), prob = c(1-testPct,testPct), replace = TRUE)
        
        data_train = data[train, ]
        labels_train = labels[train]
        
        vars = rowMaxs(apply(fac2sparse(labels_train), 1, function(x)
            rowWeightedVars(t(data_train), x)), na.rm = TRUE)
        # if (any(vars == 0)) warning ("removing genes with zero intra-class variance")
        features <- colnames(data_train)[vars > 0]
        
        fit = lda(data_train[,features], grouping = labels_train)
        pred = predict(fit, data[!train, rownames(fit$scaling)])$class
        true = labels[!train]
        acc = mean(isEqual(pred, true))
        return(list(acc = acc, fit = fit))
    }, simplify = FALSE)
    
    best = which.max(unlist(lapply(fitList, "[[", "acc")))
    
    return(fitList[[best]][["fit"]])
    
}
