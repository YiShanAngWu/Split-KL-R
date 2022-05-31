# PAC-Bayes un-expected Bernstein bound

# Helpers
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
  # for stats
  VnTerm <- Vn*(-log(1-etaGrid[argmin]*b)/(etaGrid[argmin]*b*b) - 1/b)
  LinTerm <- (KL  + log(etaGridSize/Ndelta))/(etaGrid[argmin]*n)
  return(list(val=val,VnTerm=VnTerm,LinTerm=LinTerm,etaOpt=argmin))
}

## vanilla
MGG <- function(NMC, sigma2){
  b <- 1
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
  val <- min(c(1, val))
  return(list(val=val,KL=KL, Term1=val, Term2=0, ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=Ln, L2=0, ExL1=0, ExL2=0, RefL1=0, RefL2=0, ExL1_0rate=0, ExL2_0rate=0,etaOpt1=etaOpt,etaOpt2=Inf))
}

## Forward Informed Prior
MGG_FW <- function(NMC, sigma2){
  b <- 1
  nhalf <- ntrain/2
  Ndelta <- delta
  
  # compute empirical loss & the variance
  L1 <- mean(Loss1)
  V1 <- mean(Loss1^2)
  
  # Grid of eta  
  etaGridSize <- ceil(log(0.5*sqrt(nhalf/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  # Computing the KL (KL(rho, pi_S1))
  KL <-  KLGauss(ERMs[,3],ERMs[,1], sigma2)
  
  # compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = V1,
                      KL = KL + log(sigma2GridSize),
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOpt <- tmp$etaOpt
  val1 <- tmp$val
  
  # compute the bound
  val <- L1 + val1
  val <- min(c(1, val))
  return(list(val=val,KL=KL, Term1=val, Term2=0, ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=L1, L2=0, ExL1=0, ExL2=0, RefL1=0, RefL2=0, ExL1_0rate=0, ExL2_0rate=0,etaOpt1=etaOpt,etaOpt2=Inf))
}

## Backward Informed Prior
MGG_BW <- function(NMC, sigma2){
  b <- 1
  nhalf <- ntrain/2
  Ndelta <- delta
  
  # compute empirical loss & the variance
  L2 <- mean(Loss2)
  V2 <- mean(Loss2^2)
  
  # Grid of eta  
  etaGridSize <- ceil(log(0.5*sqrt(nhalf/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  # Computing the KL (KL(rho, pi_S2))
  KL <-  KLGauss(ERMs[,3],ERMs[,2], sigma2)
  
  # compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = V2,
                      KL = KL + log(sigma2GridSize),
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOpt <- tmp$etaOpt
  val1 <- tmp$val
  
  # compute the bound
  val <- L2 + val1
  val <- min(c(1, val))
  return(list(val=val,KL=KL, Term1=0, Term2=val,  ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=0, L2=L2, ExL1=0, ExL2=0, RefL1=0, RefL2=0, ExL1_0rate=0, ExL2_0rate=0,etaOpt1=etaOpt,etaOpt2=Inf))
}

## Forward Informed Prior + Excess
MGG_FWEL <- function(NMC, sigma2){
  b <- 1
  nhalf <- ntrain/2
  Ndelta <- delta/2 # for Term1 and RefTerm1

  # Grid of eta
  etaGridSize <- ceil(log(0.5*sqrt(ntrain/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  # ExTerm1
  ## compute reference loss
  ExL1 <- mean(Excessloss1)
  ExV1 <- mean(Excessloss1^2)
  # compute the rate of Excess Loss = 0
  ExL1_0rate <- as.vector(table(Excessloss1)['0'])/nhalf/NMC
  
  # compute KL
  ## 1. Informed Prior
  KL1 <- KLGauss(ERMs[,3],ERMs[,1], sigma2)
  ## 2. NO Informed Prior
  #ratio <- initsigma2/(sigma2)
  #KL1 <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])

  # compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = ExV1,
                      KL = KL1 + log(sigma2GridSize), # 1. informed prior
                      #KL = KL1, # 2. NO informed prior
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOptExL1 <- tmp$etaOpt
  valExL1 <- tmp$val
  
  ExTerm1 <- ExL1 + valExL1
  ExTerm1 <- min(c(b, ExTerm1))
  
  # RefTerm1
  RefL1 <- mean(Refloss1)
  RefTerm1 <- bin_inv_sup(nhalf, nhalf*RefL1, Ndelta)
  
  # compute the bound
  val <- ExTerm1 + RefTerm1
  return(list(val=val,KL=KL1, Term1=0, Term2=0, ExTerm1=ExTerm1, ExTerm2=0, RefTerm1=RefTerm1, RefTerm2=0,
              L1=0, L2=0, ExL1=ExL1, ExL2=0, RefL1=RefL1, RefL2=0, ExL1_0rate=ExL1_0rate, ExL2_0rate=0,
              etaOpt1=etaOptExL1,etaOpt2=0))
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
  
  # ExTerm2
  ## compute loss & variance
  ExL2 <- mean(Excessloss2)
  ExV2 <- mean(Excessloss2^2)
  # compute the rate of Excess Loss = 0
  ExL2_0rate <- as.vector(table(Excessloss2)['0'])/nhalf/NMC
  ## compute KL
  KL2 <- KLGauss(ERMs[,3],ERMs[,2], sigma2)
  
  # compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = ExV2,
                      KL = KL2 + log(sigma2GridSize),
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOptExL2 <- tmp$etaOpt
  valExL2 <- tmp$val
  
  ExTerm2 <- ExL2 + valExL2
  ExTerm2 <- min(c(b, ExTerm2))
  
  # RefTerm1
  RefL2 <- mean(Refloss2)
  RefTerm2 <- bin_inv_sup(nhalf, nhalf*RefL2, Ndelta)
  
  # compute the bound
  val <- ExTerm2 + RefTerm2
  
  return(list(val=val,KL=KL2, Term1=0, Term2=0, ExTerm1=0, ExTerm2=ExTerm2, RefTerm1=0, RefTerm2=RefTerm2,
              L1=0, L2=0, ExL1=0, ExL2=ExL2, RefL1=0, RefL2=RefL2, ExL1_0rate=0, ExL2_0rate=ExL2_0rate,
              etaOpt1=etaOptExL1,etaOpt2=0))
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
MGG_AvgEx <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/2 # for (ExTerm1,2) and (RefTerm1,2)

  # Grid of eta
  etaGridSize <- ceil(log(0.5*sqrt(ntrain/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  # compute excess loss
  ExL1 <- mean(Excessloss1)
  ExL2 <- mean(Excessloss2)
  ExL1_0rate <- as.vector(table(Excessloss1)['0'])/nhalf/NMC
  ExL2_0rate <- as.vector(table(Excessloss2)['0'])/nhalf/NMC
  ExV1 <- mean(Excessloss1^2)
  ExV2 <- mean(Excessloss2^2)
  ExVn <- 0.5*ExV1 + 0.5*ExV2
  
  # ExTerm1,2
  # compute KL
  ## 1. Informed Prior
  #KL <- 0.5*(KLGauss(ERMs[,3],ERMs[,1], sigma2) + KLGauss(ERMs[,3],ERMs[,2], sigma2))
  ## 2. NO Informed Prior
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  
  tmp1 <- OptimEtaEmpV(etaGrid = etaGrid,
                       Vn = ExVn,
                       KL = KL + log(sigma2GridSize), # informed prior
                       #KL = KL, # NO informed prior
                       n = nhalf,
                       Ndelta = Ndelta)
  etaOptExL <- tmp1$etaOpt
  valExL <- tmp1$val
  VnTermExL <- tmp1$VnTerm
  LinTermExL <- tmp1$LinTerm
  
  ExTerm <- 0.5*ExL1 + 0.5*ExL2 + valExL
  ExTerm <- min(c(b, ExTerm))

  # compute reference terms
  RefL1 <- mean(Refloss1)
  RefTerm1 <- bin_inv_sup(nhalf, nhalf*RefL1, Ndelta)
  RefL2 <- mean(Refloss2)
  RefTerm2 <- bin_inv_sup(nhalf, nhalf*RefL2, Ndelta)
  
  val <- ExTerm + 0.5*RefTerm1 + 0.5*RefTerm2
  return(list(val=val,KL=KL, Term1=0, Term2=0, ExTerm1=ExTerm, ExTerm2=0, RefTerm1=RefTerm1, RefTerm2=RefTerm2,
              L1=0, L2=0, ExL1=ExL1, ExL2=ExL2, RefL1=RefL1, RefL2=RefL2, ExL1_0rate=ExL1_0rate, ExL2_0rate=ExL2_0rate,
              VnTermExL=VnTermExL, LinTermExL=LinTermExL,
              etaOpt1=etaOptExL,etaOpt2=0))
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


