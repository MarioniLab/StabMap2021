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
  
  out = HarmonyMatrix(as.matrix(t(embedding)), batchFactor_used, do_pca = FALSE, ...)
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
