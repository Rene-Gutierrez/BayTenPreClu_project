###############################################################################
###
### Simulated Data Generation for Mixture Tensor Graphical Models and
### Clustering
###
###############################################################################

### Set-up
set.seed(240221)
library(BayTenPreClu)
library(BayTenGraMod)
library(Matrix)
library(Tlasso)
library(pheatmap)

### Parameter Settings
ii <- 1                     # Replication Number
p <- c(40, 50, 60)          # Matrix size for each dimension
d <- length(p)              # Number of tensor dimensions
r <- c(0.15, 0.10, 0.05)    # Sparity level for Random Matrices Only
l <- c(1/4, 1/4, 1/4, 1/4)       # Mixture Weights
b <- 103                    # Degrees of Freedom for Precision Matrices
n <- 100                    # Number of Observations
nmcmc  <-1000               # Number of MCMC samples
burnin <- 100               # Burn-In Period
theta  <- 1                 # Dirichlet Process Prior Parameter
### General Suffix
suffix <- paste0(p[1], "-", p[2], "-", p[3], "-",
                 r[1], "-", r[2], "-", r[3], "-",
                 l[1], "-", l[2], "-", l[3], "-",
                 n, "-", nmcmc, "-", burnin, "-", ii)


### Generates the Precision, Covariance and Adjacency matrix for each
### component
## Number of Components
H   <- length(l)
K   <- length(p)
PSE <- list()
for(h in 1:H){
  PSE[[h]] <- spaPreMatGen(p          = p,
                           type       = "R",
                           sparsity   = r,
                           covariance = FALSE)
}

### Identity list
I <- list()
for(i in 1:K){
  I[[i]] <- diag(p[i])
}

### Asigns each Matrix List
O <- list()
S <- list()
E <- list()
for(h in 1: H){
  O[[h]] <- PSE[[h]]$P
  S[[h]] <- PSE[[h]]$S
  E[[h]] <- PSE[[h]]$E
}

# Generates the Samples
print("Creating Samples")
## Initialization of Cample Tensor List
lT <- list()
## Component mebership
zT  <- matrix(data = NA, nrow = n, ncol = H)
zF  <- matrix(data = NA, nrow = n, ncol = H)
zF2  <- matrix(data = NA, nrow = n, ncol = 2)
zF5  <- matrix(data = NA, nrow = n, ncol = 5)
## Component membership by Category
zTC <- numeric(n)
zFC <- numeric(n)
zFC2 <- numeric(n)
zFC5 <- numeric(n)
for(i in 1:n){
  ## Component Membership
  zT[i,] <- rmultinom(n    = 1,
                      size = 1,
                      prob = l)
  ## Component Membership (Fake)
  zF[i,] <- rmultinom(n    = 1,
                      size = 1,
                      prob = l)
  zF2[i,] <- rmultinom(n    = 1,
                       size = 1,
                       prob = rep(1/2, 2))
  zF5[i,] <- rmultinom(n    = 1,
                       size = 1,
                       prob = rep(1/5, 5))
  ## By Category
  zTC[i] <- which(zT[i,] == 1)
  zFC[i] <- which(zF[i,] == 1)
  zFC2[i] <- which(zF2[i,] == 1)
  zFC5[i] <- which(zF5[i,] == 1)
  ## Samples According to the Category
  lT[[i]] <- TNormalSampler(n         = 1,
                            SigmaList = S[[zTC[i]]])[[1]]
}

### Computes the K-mode Matricization of each tensor and takes its covariance matrix for each mode and tensor
kCov <- list()
### For every mode
for(k in 1:K){
  kM <- matrix(data = NA, nrow = p[k] * p[k], ncol = n)
  ### For Every tensor
  for(i in 1:n){
    kM[,i] <- as.vector(cov(kModMat(tensor = lT[[i]], mode = k)))
  }
  kCov[[k]] <- kM
}


### Runs the Stochastic Search Clustering
iP <- list()
iP[[1]] <- seq(1,n)
out <- PMBTC(iP = iP,
             theta    = theta,
             kPseCov  = kCov,
             burnin   = burnin,
             nmcmc    = nmcmc,
             progress = TRUE)

### Computes the Probability Matrix
pM <- proMat(sP = out[(burnin + 1):(nmcmc + burnin)],
             n  = n)

## Cluster Distribution
### Cluster Number Initialization
cluNum <- numeric(length = nmcmc)
### Obtains the Number of Clusters
for(i in (burnin + 1):(burnin + nmcmc)){
  cluNum[i] <- length(out[[i]])
}

# ARI Distribution
## Clusters
samClu <- out[101:1100]
reaClu <- list()
for(i in 1:3){
  reaClu[[i]] <- seq(1:n)[zTC == i]
}

## ARI Samples
ariSam <- numeric(length = nmcmc)
## Computes the ARI for Each Sample
for(i in 1:nmcmc){
  ariSam[i] <- ari(clu1 = samClu[[i]],
                   clu2 = reaClu)
}

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