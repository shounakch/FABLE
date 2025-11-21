cleanData <- function(rawData, needs.processing = TRUE, var.prop = 0.10) {
  
  ##Return centered data matrix after cleaning (if required).
  
  if(needs.processing) {
    
    ##M3Drop processing
        
    # rawData = as.data.frame(rawData)
    # rownames(rawData) = rawData[,1]
    # rawData = rawData[,-1]
    # 
    # ## Use M3Drop to drop irrelevant variables
    # 
    # YNorm = M3Drop::M3DropConvertData(rawData, is.counts = FALSE, is.log=TRUE) #requires p*n matrix, not n*p
    # M3Drop_genes = M3Drop::M3DropFeatureSelection(YNorm, mt_method = "bon",
    #                                               mt_threshold = 0.05)
    # 
    # YNorm = as.data.frame(YNorm)
    # YNormImp = YNorm %>% filter(row.names(YNorm) %in% M3Drop_genes$Gene)
    # 
    # YLog = log2(YNormImp + 1)
    # Y_uncen = t(YLog) #uncentered, center for LFM use
    
    ## Variance processing of data (choose top 5%)
    
    rawData = as.data.frame(rawData)
    rownames(rawData) = rawData[,1]
    rawData = rawData[,-1]
    
    rawData = log2(rawData + 1)
    rawData = t(rawData)
    meanRawData = colMeans(rawData)
    centData = sweep(rawData, 2, meanRawData)
    colVars = apply(centData, 2, var)
    highVarIndices = which(colVars >= quantile(colVars, 1 - var.prop))
    
    centData = centData[,highVarIndices]
    
  }else {
    
    rawData = t(rawData)
    
    meanRawData = colMeans(rawData)
    centData = sweep(rawData, 2, meanRawData)
    
    colVars = apply(centData, 2, var)
    highVarIndices = which(colVars >= quantile(colVars, 1 - var.prop))
    
    centData = centData[,highVarIndices]
    
  }
  
  return(centData)
  
}