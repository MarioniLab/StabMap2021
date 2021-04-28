stabMapComparative = function(SCE_X,
                              SCE_Y,
                              genes = intersect(rownames(SCE_X), rownames(SCE_Y)),
                              assayNameX = "logcounts",
                              assayNameY = "logcounts",
                              ncomponentsFull = 50,
                              ncomponentsSubset = 50,
                              sparse = FALSE) {
    
    # TODO: use list or ... to expand to more than just two datasets
    # TODO: lookup table for features that should be considered joint
    # between the various data
    
    # Given two (or more) datasets without a set reference vs query structure
    # (i.e. comparative scenario), perform stabMap
    
    if (sparse) require(glmnet)
    
    # calculate PCs for each dataset
    PCs_X = calculatePCA(assay(SCE_X, assayNameX), ncomponents = ncomponentsFull)
    PCs_X_genes = calculatePCA(assay(SCE_X, assayNameX)[genes,],
                               ncomponents = ncomponentsSubset)
    
    PCs_Y = calculatePCA(assay(SCE_Y, assayNameY), ncomponents = ncomponentsFull)
    PCs_Y_genes = calculatePCA(assay(SCE_Y, assayNameY)[genes,],
                               ncomponents = ncomponentsSubset)
    
    # loadings are saved in the attributes
    PCs_X_genes_loadings = attr(PCs_X_genes, "rotation")
    PCs_Y_genes_loadings = attr(PCs_Y_genes, "rotation")
    
    if (sparse) {
        # Fit lasso models with lambda selected with CV to estimate the best 
        # linear combination of PCs among the subset features
        coefList_X = list()
        coefList_Y = list()
        # i = 1
        for (i in seq_len(ncol(PCs_X))) {
            print(i)
            fit.cv = cv.glmnet(x = PCs_X_genes, y = PCs_X[,i], intercept = FALSE)
            fit = glmnet(x = PCs_X_genes, y = PCs_X[,i],
                         lambda = fit.cv["lambda.min"][[1]], intercept = FALSE)
            coefList_X[[i]] <- coef(fit)[-1]
        }
        coef_X = do.call(cbind,coefList_X)
        for (i in seq_len(ncol(PCs_Y))) {
            print(i)
            fit.cv = cv.glmnet(x = PCs_Y_genes, y = PCs_Y[,i], intercept = FALSE)
            fit = glmnet(x = PCs_Y_genes, y = PCs_Y[,i],
                         lambda = fit.cv["lambda.min"][[1]], intercept = FALSE)
            coefList_Y[[i]] <- coef(fit)[-1]
        }
        coef_Y = do.call(cbind,coefList_Y)
    } else {
        # fit = lm.fit(x = referencePCs_genes, y = referencePCs[,i])
        # coefList[[i]] <- fit$coefficients
        fit = lm.fit(x = PCs_X_genes, y = PCs_X)
        coef_X = fit$coefficients
        
        fit = lm.fit(x = PCs_Y_genes, y = PCs_Y)
        coef_Y = fit$coefficients
    }
    
    
    all_assay = cbind(assay(SCE_X, assayNameX)[genes,],
                      assay(SCE_Y, assayNameY)[genes,])
    
    all_assay_sub_X = all_assay[rownames(PCs_X_genes_loadings), ]
    all_assay_sub_Y = all_assay[rownames(PCs_Y_genes_loadings), ]
    
    # d_sub = d[rownames(referencePCs_genes_loadings),]
    
    emb_X = t(all_assay_sub_X) %*% PCs_X_genes_loadings %*% coef_X
    emb_Y = t(all_assay_sub_Y) %*% PCs_Y_genes_loadings %*% coef_Y
    
    emb = cbind(emb_X, emb_Y)
    
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


stabMapLabelled = function(referenceSCE,
                           querySCE,
                           genes = intersect(rownames(querySCE), rownames(referenceSCE)),
                           grouping = "celltype_parsed_sub",
                           assayNameReference = "logcounts",
                           assayNameQuery = "logcounts",
                           prop_explained = 1,
                           selectionLFC = ifelse(length(genes) >= 1000, 0.1, 0)) {
    
    # stabmap using LDA embedding
    
    # selectionLFC is a log-fold change to use for selecting features among genes first
    
    if (selectionLFC != 0) {
        
        print("performing fold-change feature selection")
        
        # make assay dense (not ideal)
        assay(referenceSCE, assayNameReference) <- as.matrix(assay(referenceSCE, assayNameReference))
        
        # only compare against groupings that have at least 10 cells
        sub = names(which(table(colData(referenceSCE)[,grouping]) >= 10))
        
        sel = SpatialUtils::getMarkerArray(referenceSCE,
                                           assayName = assayNameReference,
                                           group = grouping,
                                           subset = colData(referenceSCE)[,grouping] %in% sub,
                                           verbose = FALSE)
        
        # only keep genes with at least selectionLFC log fold change
        topLFC = rowMax(abs(sel[genes,,"LFC"]))
        
        genes_sel = genes[topLFC >= selectionLFC]
        
        genes <- genes_sel
        
        print(length(genes))
        
    }
    
    require(MASS)
    
    # prefilter cells and genes if needed
    # v = fac2sparse(colData(referenceSCE)[,grouping]) %*% t(assay(referenceSCE, assayNameReference)[genes,])
    # vars = apply(assay(referenceSCE, assayNameReference)[genes,], 1, function(x)
    # max(tapply(x, colData(referenceSCE)[,grouping], var), na.rm = TRUE))
    # genes <- genes[vars > 0]
    
    vars = rowMaxs(apply(fac2sparse(colData(referenceSCE)[,grouping]), 1, function(x)
        rowWeightedVars(assay(referenceSCE, assayNameReference)[genes,],
                        x)), na.rm = TRUE)
    if (any(vars == 0)) warning ("removing genes with zero intra-class variance")
    genes <- genes[vars > 0]
    
    # build LDA model for genes
    fit = lda(t(assay(referenceSCE, assayNameReference)[genes,]),
              grouping = colData(referenceSCE)[,grouping])
    
    # proportion of between group variation explained by the linear
    # discriminants
    prop_rel = fit$svd^2/sum(fit$svd^2)
    n_prop_rel = which(cumsum(prop_rel) >= prop_explained)[1]
    
    if (is.null(querySCE)) {
        all_assay = cbind(assay(referenceSCE, assayNameReference)[genes,])
    } else {
        all_assay = cbind(assay(referenceSCE, assayNameReference)[genes,],
                          assay(querySCE, assayNameQuery)[genes,])
    }
    
    resub = predict(object = fit, newdata = t(all_assay))
    
    if (prop_explained == 1) {
        resub_out = resub$x
    } else {
        resub_out = resub$x[,seq_len(n_prop_rel)]  
    }
    return(resub_out)
}



mapPCA = function(referenceSCE,
                  querySCE,
                  genes = intersect(rownames(querySCE), rownames(referenceSCE)),
                  assayNameReference = "logcounts",
                  assayNameQuery = "logcounts",
                  nPCs = 50,
                  multiBatchPCA = TRUE) {
    
    # get PCA embedding - there is no stability aspect here
    # and is primarily used for comparison
    
    if (is.null(querySCE)) {
        all_assay = cbind(assay(referenceSCE, assayNameReference)[genes,])
    } else {
        all_assay = cbind(assay(referenceSCE, assayNameReference)[genes,],
                          assay(querySCE, assayNameQuery)[genes,])
    }
    
    if (multiBatchPCA) {
        require(batchelor)
        pca_all <- batchelor::multiBatchPCA(all_assay,
                                            d = nPCs,
                                            batch = rep(c("reference","query"),
                                                        times = c(ncol(referenceSCE), ncol(querySCE))),
                                            preserve.single = TRUE)[[1]]
    } else {
        require(irlba)
        pca_all <- irlba::prcomp_irlba(t(all_assay), n = nPCs, scale. = TRUE)$x
        rownames(pca_all) <- colnames(all_assay)
    }
    
    return(pca_all)
}

reducedMNN_batchFactor = function(LD_embedding,
                                  batchFactor) {
    # batch correct within this embedding, wrapper around reducedMNN
    # batchFactor is a named vector that is matched
    
    
    require(batchelor)
    
    batchFactor_used = batchFactor[rownames(LD_embedding)]
    
    out = reducedMNN(LD_embedding, batch = batchFactor_used)
    resub_corrected = out$corrected
    
    return(resub_corrected)
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
    
    # extract similarity of the subsetted genes
    subset_sim = generateSimilarity(SCE, HVGs = rownames(jointSCE))
    
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
