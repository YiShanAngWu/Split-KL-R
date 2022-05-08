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
  
  return(list(val=val,KL=KL, Term1=val, Term2=0, ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=L1, L2=0, ExL1=0, ExL2=0, RefL1=0, RefL2=0,etaOpt=c(etaOpt,Inf)))
}

## Backward
MGG_BW <- function(NMC, sigma2){
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
  
  return(list(val=val,KL=KL, Term1=0, Term2=val,  ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=0, L2=L2, ExL1=0, ExL2=0, RefL1=0, RefL2=0,etaOpt=c(etaOpt,Inf)))
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
  RefL1 <- mean(Refloss1)
  RefV1 <- mean(Refloss1^2)
  ## compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = RefV1,
                      KL = 0 + log(sigma2GridSize),
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOptRefL1 <- tmp$etaOpt
  valRefL1 <- tmp$val
  
  RefTerm1 <- RefL1 + valRefL1
  
  # ExTerm1
  ## compute reference loss
  ExL1 <- mean(Excessloss1)
  ExV1 <- mean(Excessloss1^2)
  # compute the rate of Excess Loss = 0
  ExL1_0rate <- as.vector(table(Excessloss1)['0'])/nhalf/NMC
  ## compute KL
  KL1 <- KLGauss(ERMs[,3],ERMs[,1], sigma2)

  # compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                      Vn = ExV1,
                      KL = KL1 + log(sigma2GridSize),
                      n = nhalf,
                      Ndelta = Ndelta)
  etaOptExL1 <- tmp$etaOpt
  valExL1 <- tmp$val
  
  ExTerm1 <- ExL1 + valExL1
  
  # compute the bound
  val <- ExTerm1 + RefTerm1
  return(list(val=val,KL=KL1, Term1=0, Term2=0, ExTerm1=ExTerm1, ExTerm2=0, RefTerm1=RefTerm1, RefTerm2=0,
              L1=0, L2=0, ExL1=ExL1, ExL2=0, RefL1=RefL1, RefL2=0, ExL1_0rate=ExL1_0rate, ExL2_0rate=0,
              etaOpt=c(etaOptExL1,etaOptRefL1)))
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
  RefL2 <- mean(Refloss2)
  RefV2 <- mean(Refloss2^2)
  ## compute eta-related terms
  tmp <- OptimEtaEmpV(etaGrid = etaGrid,
                       Vn = RefV2,
                       KL = 0 + log(sigma2GridSize),
                       n = nhalf,
                       Ndelta = Ndelta)
  etaOptRefL2 <- tmp$etaOpt
  valRefL2 <- tmp$val
  
  RefTerm2 <- RefL2 + valRefL2
  
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
  
  # compute the bound
  val <- ExTerm2 + RefTerm2
  
  return(list(val=val,KL=KL2, Term1=0, Term2=0, ExTerm1=0, ExTerm2=ExTerm2, RefTerm1=0, RefTerm2=RefTerm2,
              L1=0, L2=0, ExL1=0, ExL2=ExL2, RefL1=0, RefL2=RefL2, ExL1_0rate=0, ExL2_0rate=ExL2_0rate,
              etaOpt=c(etaOptExL2,etaOptRefL2)))
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
  
  # compute reference loss
  loss2 <- t(matrix(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERMs[,2])), nrow=NMC, ncol=nhalf, byrow=TRUE))
  loss1 <- t(matrix(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERMs[,1])), nrow=NMC, ncol=nhalf, byrow=TRUE))
  Vn2 <- 0.5*(mean(loss1^2)+mean(loss2^2))
  
  # compute excess loss
  Diffloss2 <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samplesTS))-loss2
  Diffloss1 <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samplesTS))-loss1
  Vn1 <- 0.5*(mean(Diffloss1^2)+mean(Diffloss2^2))
  
  # Grid of eta
  etaGridSize <- ceil(log(0.5*sqrt(ntrain/log(1/delta)))/log(rho))
  etaGrid <- numeric(etaGridSize)
  for(jj in 1:(etaGridSize))
    etaGrid[jj] <- 1/(b*rho^(jj))
  
  # Term1,2
  #Vn1 <- VnTerm(ERMfull = ERMs[,3],
  #                 ERM1 = ERMs[,2],
  #                 ERM2 = ERMs[,1],
  #                 NMC = NMC, sigma2=sigma2)
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
  # For statistics
  Term1 <- 0.5*(mean(Diffloss1)+mean(Diffloss2)) + val1

  # RefTerm1,2
  #Vn2 <- VnPrimeTerm(ERM1 = ERMs[,2],
  #                          ERM2 = ERMs[,1],
  #                          NMC = NMC)
  KL2 <- 0
  tmp2 <- OptimEtaEmpV(etaGrid = etaGrid,
                       Vn = Vn2,
                       KL = KL2 + log(sigma2GridSize),
                       n = nhalf,
                       Ndelta = Ndelta/2) # union bound of \hat L(hS1,S2) & \hat L(hS2,S1)
  etaOpt2 <- tmp2$etaOpt
  val2 <- tmp2$val
  # For statistics
  RefTerm1 <- 0.5*(mean(loss1^2)+mean(loss2^2)) + val2
  
  val <- Ln + val1 + val2
  return(list(val=val,KL=KL1, Term1=Term1, Term2=Term1, RefTerm1=RefTerm1, RefTerm2=RefTerm1,
              ExL1=mean(Diffloss1), ExL2=mean(Diffloss2), RefL1=mean(loss1), RefL2=mean(loss2),
              etaOpt=c(etaOpt1,etaOpt2)))
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


