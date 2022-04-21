# Helpers
buildERMsequenceFast <- function(eps = 1e-16){
  ERMsequence <- matrix(nrow = d, ncol = 3, data = 0)
  for(ii in c(1,2,3))   ## change here for the fast version
  {
    if(ii==1) {setPoints <- 1:(ntrain/2)}        ## ERM on first half
    if(ii==2) {setPoints <- (ntrain/2+1):ntrain} ## ERM on second half
    if(ii==3) {setPoints <- 1:ntrain}            ## ERM on full sample
    
    ## Definining the objective and the corresponding gradient
    fn <- function(theta){
      X <- Xtrain[setPoints,]
      Y <- Ytrain[setPoints]
      Z <-  as.numeric(X%*%theta)
      return(-mean(Y*log(sigmoid(Z))+(1-Y)*log(1-sigmoid(Z))) +lambda*dot(theta,theta)/2)
    }
    gr <- function(theta){
      X <- Xtrain[setPoints,]
      Y <- Ytrain[setPoints]
      Z <- as.numeric(X%*%theta)
      return(-as.numeric((Y-sigmoid(Z))%*%X/length(setPoints)) + lambda*theta)
    }
    par <- numeric(d)
    ERMsequence[,ii] <- optim(par, fn, gr = gr, method = "BFGS", hessian = TRUE)$par
  }
  return(ERMsequence)
}

computeExpectation <- function(ERMfull,NMC, sigma2){
  #theta_samples <- get_sample(type = distribution, mean=ERMfull, variance2=sigma2, n_samples=NMC)
  result <- loss(Ytrain,predictor(Xtrain,theta_samplesTS))
  return(mean(result))
}

OptimEtaEmpV <- function(etaGrid, Vn, KL, n, Ndelta){
  etaGridSize <- length(etaGrid)
  result <- numeric(etaGridSize)
  for(etaInd in 1:etaGridSize){
    eta <- etaGrid[etaInd]
    result[etaInd] <- Vn*(-log(1-eta*b)/(eta*b*b) - 1/b) +
      (KL  + log(etaGridSize/Ndelta))/(eta*n)
  }
  argmin <- which.min(result)
  val <- result[argmin]
  return(list(val=val,etaOpt=argmin))
}

# PAC-Bayes un-expected Bernstein bound

## vanilla
MGG <- function(NMC, sigma2){
  Ndelta <- delta
  # compute empirical loss & the variance
  Ln <- mean(loss(Ytrain,predictor(Xtrain,theta_samplesTS)))
  Vn <- mean(loss(Ytrain,predictor(Xtrain,theta_samplesTS))^2)
  
  # Grid of eta  
  etaGridSize <- ceil(log(0.5*sqrt(ntrain/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))  
  
  # Compute the KL 
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  
  # compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = Vn,
                      KL = KL,
                      n = ntrain,
                      Ndelta = Ndelta)
  etaOpt <- tmp$etaOpt
  val1 <- tmp$val
  
  # compute the bound
  val <- Ln + val1
  
  return(list(val=val,KL=KL,etaOpt=c(etaOpt,Inf)))
}

## Forward
MGG_FW <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta
  
  # compute empirical loss & the variance
  Ln <- mean(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],theta_samplesTS)))
  Vn <- mean(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],theta_samplesTS))^2)
  
  # Grid of eta  
  etaGridSize <- ceil(log(0.5*sqrt(nhalf/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  # Computing the KL (KL(rho, pi_S1))
  KL <-  KLGauss(ERMs[,3],ERMs[,1], sigma2)
  
  # compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = Vn,
                      KL = KL + log(sigma2GridSize),
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOpt <- tmp$etaOpt
  val1 <- tmp$val
  
  # compute the bound
  val <- Ln + val1
  
  return(list(val=val,KL=KL,etaOpt=c(etaOpt,Inf)))
}

## Backward
MGG_BW <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta
  
  # compute empirical loss & the variance
  Ln <- mean(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samplesTS)))
  Vn <- mean(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samplesTS))^2)
  
  # Grid of eta  
  etaGridSize <- ceil(log(0.5*sqrt(nhalf/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  # Computing the KL (KL(rho, pi_S2))
  KL <-  KLGauss(ERMs[,3],ERMs[,2], sigma2)
  
  # compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = Vn,
                      KL = KL + log(sigma2GridSize),
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOpt <- tmp$etaOpt
  val1 <- tmp$val
  
  # compute the bound
  val <- Ln + val1
  
  return(list(val=val,KL=KL,etaOpt=c(etaOpt,Inf)))
}

## Forward + Excess
MGG_FWEL <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/2 # for Term1 and RefTerm1

  # Grid of eta
  etaGridSize <- ceil(log(0.5*sqrt(ntrain/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  # RefTerm1
  ## compute loss & variance
  losses <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERMs[,1]))
  Ln2 <- mean(losses)
  Vn2 <- mean(losses^2)
  ## compute KL
  KL2 <- 0
  ## compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = Vn2,
                      KL = KL2 + log(sigma2GridSize),
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOpt2 <- tmp$etaOpt
  val2 <- tmp$val
  
  RefTerm1 <- Ln2 + val2
  
  # Term1
  ## compute reference loss
  losses <- t(matrix(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERMs[,1])), nrow=NMC, ncol=nhalf, byrow=TRUE))
  ## compute excess loss
  Diffloss <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samplesTS))-losses
  ## compute loss & variance
  Ln1 <- mean(Diffloss)
  Vn1 <- mean(Diffloss^2)
  ## compute KL
  KL1 <- KLGauss(ERMs[,3],ERMs[,1], sigma2)

  # compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = Vn1,
                      KL = KL1 + log(sigma2GridSize),
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOpt1 <- tmp$etaOpt
  val1 <- tmp$val
  
  Term1 <- Ln1 + val1
  
  # compute the bound
  val <- Term1 + RefTerm1
  return(list(val=val,KL=KL1,etaOpt=c(etaOpt1,etaOpt2)))
}

## Backward + Excess
MGG_BWEL <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/2 # for Term2 and RefTerm2
  
  # Grid of eta
  etaGridSize <- ceil(log(0.5*sqrt(ntrain/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  # RefTerm2
  ## compute loss & variance
  losses <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERMs[,2]))
  Ln2 <- mean(losses)
  Vn2 <- mean(losses^2)
  ## compute KL
  KL2 <- 0
  ## compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                       Vn = Vn2,
                       KL = KL2 + log(sigma2GridSize),
                       n = nhalf,
                       Ndelta = Ndelta)
  etaOpt2 <- tmp$etaOpt
  val2 <- tmp$val
  
  RefTerm2 <- Ln2 + val2
  
  # Term2
  ## compute reference loss
  losses <- t(matrix(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERMs[,2])), nrow=NMC, ncol=nhalf, byrow=TRUE))
  ## compute excess loss
  Diffloss <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),], theta_samplesTS))-losses
  ## compute loss & variance
  Ln1 <- mean(Diffloss)
  Vn1 <- mean(Diffloss^2)
  ## compute KL
  KL1 <- KLGauss(ERMs[,3],ERMs[,2], sigma2)
  
  # compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = Vn1,
                      KL = KL1 + log(sigma2GridSize),
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOpt1 <- tmp$etaOpt
  val1 <- tmp$val
  
  Term2 <- Ln1 + val1
  
  # compute the bound
  val <- Term2 + RefTerm2
  
  return(list(val=val,KL=KL1,etaOpt=c(etaOpt1,etaOpt2)))
}

## For Average + Excess
COMP <- function(ERMfull,ERM1,ERM2,sigma2){
  #val1 <- sum(square_diff(ERMfull,ERM1))
  #val2 <- sum(square_diff(ERMfull,ERM2))
  #result <-   (val1+val2)/(2*sigma2)
  val1 <- KLGauss(ERMfull, ERM1, sigma2)
  val2 <- KLGauss(ERMfull, ERM2, sigma2)
  result <- val1+val2
  return(result)
}

VnTerm <- function(ERMfull,ERM1,ERM2,NMC,sigma2){
  #theta_samples <- get_sample(type = distribution, mean=ERMfull, variance2=sigma2, n_samples=NMC)
  result  <- matrix(nrow = ntrain, ncol = NMC, data = NA)
  loss1 <- t(matrix(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERM1)), nrow=NMC, ncol=ntrain/2, byrow=TRUE))
  loss2 <- t(matrix(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERM2)), nrow=NMC, ncol=ntrain/2, byrow=TRUE))
  result[1:(ntrain/2),] <- square_diff(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samplesTS)),loss1)
  result[(ntrain/2+1):ntrain,]<-square_diff(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samplesTS)),loss2)
  return(mean(result))
}

VnPrimeTerm <- function(ERM1,ERM2,NMC){
  res <- numeric(ntrain)
  res[1:(ntrain/2)] <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERM1))^2
  res[(ntrain/2+1):ntrain] <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERM2))^2
  return(mean(res))
}

## Average + Excess
MGG_Avg <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/2 # for (Term1,2) and (RefTerm1,2)
  Ln <- mean(loss(Ytrain,predictor(Xtrain,theta_samplesTS)))
  
  # Grid of eta
  etaGridSize <- ceil(log(0.5*sqrt(ntrain/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  # Term1,2
  Vn1 <- VnTerm(ERMfull = ERMs[,3],
                   ERM1 = ERMs[,2],
                   ERM2 = ERMs[,1],
                   NMC = NMC, sigma2=sigma2)
  KL1 <- 0.5*COMP(ERMfull = ERMs[,3],
             ERM1 = ERMs[,2],
             ERM2 = ERMs[,1],
             sigma2 = sigma2)
  tmp1 <- OptimEtaEmpV(etaGrid = etaGrid,
                       Vn = Vn1,
                       KL = KL1 + log(sigma2GridSize),
                       n = nhalf,
                       Ndelta = Ndelta)
  etaOpt1 <- tmp1$etaOpt
  val1 <- tmp1$val

  # RefTerm1,2
  Vn2 <- VnPrimeTerm(ERM1 = ERMs[,2],
                            ERM2 = ERMs[,1],
                            NMC = NMC)
  KL2 <- 0
  tmp2 <- OptimEtaEmpV(etaGrid = etaGrid,
                       Vn = Vn2,
                       KL = KL2 + log(sigma2GridSize),
                       n = nhalf,
                       Ndelta = Ndelta/2) # union bound of \hat L(hS1,S2) & \hat L(hS2,S1)
  etaOpt2 <- tmp2$etaOpt
  val2 <- tmp2$val
  
  val <- Ln + val1 + val2
  return(list(val=val,KL=KL1,etaOpt=c(etaOpt1,etaOpt2)))
}

## NOT USING!!
MGG_IF <- function(NMC,sigma2){
  Ln <- mean(loss(Ytrain,predictor(Xtrain,theta_samplesTS)))
  
  # Grid of eta
  etaGridSize <- ceil(log(0.5*sqrt(ntrain/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  vnTerm <- VnTerm(ERMfull = ERMs[,3],
                          ERM1 = ERMs[,2],
                          ERM2 = ERMs[,1],
                          NMC = NMC, sigma2=sigma2)
  KL <- COMP(ERMfull = ERMs[,3],
                   ERM1 = ERMs[,2],
                   ERM2 = ERMs[,1],
                   sigma2 = sigma2)
  tmp1 <- OptimEtaVn(etaGrid,
                     vnTerm = vnTerm,
                     compTerm = KL)
  etaOpt1 <- tmp1$etaOpt
  val1 <- tmp1$val
  
  vnTermPrim <- VnPrimeTerm(ERM1 = ERMs[,2],
                                ERM2 = ERMs[,1],
                                NMC = NMC)
  tmp2 <- OptimEtaVnPrime(etaGrid,
                          vnTermPrim = vnTermPrim)
  etaOpt2 <- tmp2$etaOpt
  val2 <- tmp2$val
  
  val <- Ln + val1 + val2
  return(list(val=val,val1=val1,val2=val2,vnTerm=vnTerm,vnTermPrim=vnTermPrim,KL=KL,etaOpt=c(etaOpt1,etaOpt2)))
}


