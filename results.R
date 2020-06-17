# Result Plots

## General Parameters
the       <- theta * 10
reaCluNum <- length(l)
tenSiz    <- "big"

### Adjusted Rand Index

samClu <- out[101:1100]
reaClu <- list()
for(i in 1:3){
    reaClu[[i]] <- seq(1:n)[zTC == i]
}

## proMat

# Open a pdf file
filNam <- paste0("./out/proMat_", reaCluNum ,"_", the, "_", tenSiz, ".pdf")
pdf(file   = filNam,
    width  = 5,
    height = 5) 
# 2. Create a plot
pheatmap(mat            = pM / nmcmc,
         border_color   = NA,
         treeheight_row = 0,
         treeheight_col = 0)
# Close the pdf file
dev.off()

# Open a pdf file
filNam <- paste0("./out/cluNum_", reaCluNum ,"_", the, "_", tenSiz, ".pdf")
pdf(file   = filNam,
    width  = 5,
    height = 5) 
# 2. Create a plot
hist(cluNum,
     breaks = 10,
     freq = FALSE,
     xlab = "Number Of Clusters",
     main = "")
abline(v   = median(cluNum),
       col = "red",
       lwd = 3)
# Close the pdf file
dev.off() 

# Cluster Distribution
## Cluster Number Initialization
cluNum <- numeric(length = nmcmc)
## Obtains the Number of Clusters
for(i in (burnin + 1):(burnin + nmcmc)){
    cluNum[i] <- length(out[[i]])
}

# Open a pdf file
filNam <- paste0("./out/ari_", reaCluNum ,"_", the, "_", tenSiz, ".pdf")
pdf(file   = filNam,
    width  = 5,
    height = 5) 
# 2. Create a plot
hist(ariSam,
     breaks = 10,
     freq = FALSE,
     xlab = "ARI",
     main = "")
abline(v   = median(ariSam),
       col = "red",
       lwd = 3)
# Close the pdf file
dev.off() 

# Writes the ARI for the Estimated number of Clusters
# Partition Point Estimate
## Gets the Partition Estimate
estPar <- salso::salso(psm = pM / nmcmc)$estimate
## Builds the Clustering List
estClu <- list()
for(i in 1:max(estPar)){
    estClu[[i]] <- seq(1:n)[estPar == i]
}
## ARI of the Best Estimate
ariEst <- ari(clu1 = reaClu, clu2 = estClu)
filNam <- paste0("./out/estAri_", reaCluNum ,"_", the, "_", tenSiz, ".txt")
write.table(ariEst, file = filNam)