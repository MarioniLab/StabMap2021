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

getSVMModels = function(X, Ys) {
    # X is a cells x features matrix
    # of explanatory variables
    # like a reducedDim of a SCE
    # Ys is another different cells x 
    # features matrix
    # the cells in X and Y should be matched
    # output is a list of model objects,
    # of the same length as nrow(Y)
    fitList = apply(Ys, 2, function(y) {
        message("fitting...")
        fit = svm(x = X, y = y, kernel = "polynomial", scale = FALSE)
        return(fit)
    })
    return(fitList)
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
        
        if (FALSE) {
        coefs = mapply(sparseCoef, PCs, PCs_genes, SIMPLIFY = FALSE)
        }
        
        # test with SVM
        require(e1071)
        coefs = mapply(getSVMModels, PCs_genes, PCs, SIMPLIFY = FALSE)
        
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
    
    # test with SVM:
    if (FALSE) {
    matMult = function(m1_rnames, m2, m3, m1) {
        t(m1[m1_rnames,]) %*% m2 %pred% m3 
    }
    
    emb_list = mapply(matMult,
                      m1_rnames = lapply(PCs_genes_loadings, rownames),
                      m2 = PCs_genes_loadings,
                      m3 = coefs,
                      MoreArgs = list(m1 = all_assay))
    }
    
    if (TRUE) {
    matMult = function(m1_rnames, m2, m3, m1) {
        t(m1[m1_rnames,]) %*% m2 %*% m3 
    }
    
    emb_list = mapply(matMult,
                      m1_rnames = lapply(PCs_genes_loadings, rownames),
                      m2 = PCs_genes_loadings,
                      m3 = coefs,
                      MoreArgs = list(m1 = all_assay))
    }
    
    # emb_X = t(all_assay_sub_X) %*% PCs_X_genes_loadings %*% coef_X
    # emb_Y = t(all_assay_sub_Y) %*% PCs_Y_genes_loadings %*% coef_Y
    
    # emb = cbind(emb_X, emb_Y)
    
    emb = do.call(cbind, emb_list[stabilise])
    
    return(emb)
    
}

stabMapContinuous = function(referenceSCE,
                             querySCE,
                             genes = intersect(rownames(querySCE), rownames(referenceSCE)),
                             assayNameReference = "logcounts",
                             assayNameQuery = "logcounts",
                             ncomponentsFull = 50,
                             ncomponentsSubset = 50,
                             sparse = TRUE) {
    
    # superceded by stabMapComparative!
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
    # require(matrixStats)
    
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
    
    assay(referenceSCE, assayNameReference) <- as.matrix(assay(referenceSCE, assayNameReference))
    
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

