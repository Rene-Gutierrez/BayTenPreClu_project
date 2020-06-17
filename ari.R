ari <- function(clu1, clu2){
  #  Build the Contingency Table
  ## Number of Partitions for the Predicted Clustering
  preCluLen <- length(clu1)
  ## Number of Partitions for the Real Clustering
  reaCluLen <- length(clu2)
  ## Initializes the COntingency Table
  conTab <- matrix(NA,
                   nrow = preCluLen,
                   ncol = reaCluLen)
  ## Populates the Matrix
  for(i in 1:preCluLen){
    for(j in 1:reaCluLen){
      conTab[i,j] <- sum(clu1[[i]] %in% clu2[[j]])
    }
  }
  
  ## Obtains the Contingency Table Margins
  colMar <- rowSums(conTab)
  rowMar <- colSums(conTab)
  
  # Obtians the Diferent Combinations
  ## Single Entries
  comn <- conTab * (conTab - 1) / 2
  ## Column Margins Combinations 
  comC <- colMar * (colMar - 1) / 2
  ## Row Margins Combinations
  comR <- rowMar * (rowMar - 1) / 2
  ## Total COmbinations
  comT <- n * (n - 1) / 2
  
  # Computes the ARI
  res <- (sum(comn) - sum(comR) * sum(comC) / comT) /(0.5*(sum(comC) + sum(comR)) - sum(comR) * sum(comC) / comT)
  
  # Returns the Value
  return(res)
}