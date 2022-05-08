# Maurer's bound (PBkl)

## Vanilla
PBkl <- function(NMC, sigma2){
  Ndelta <- delta
  Ln <- mean(loss(Ytrain,predictor(Xtrain,theta_samplesTS)))
  # Computing the KL 
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  # Computing the complexity
  RHS <- (KL + log(2*sqrt(ntrain)/Ndelta))/ntrain
  # Computing the bound
  val <- kl_inv_sup(Ln, RHS)
  
  return(list(val=val,KL=KL))
}

## Forward
PBkl_FW <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta
  L1 <- mean(Loss1)
  # Computing the KL (KL(rho, pi_S1))
  KL <-  KLGauss(ERMs[,3],ERMs[,1], sigma2)
  
  RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  val <- kl_inv_sup(L1, RHS)
  
  return(list(val=val,KL=KL, Term1=val, Term2=0, ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=L1, L2=0, ExL1=0, ExL2=0, RefL1=0, RefL2=0))
}

## Backward
PBkl_BW <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta

  L2 <- mean(Loss2)
  # Computing the KL (KL(rho, pi_S2))
  KL <-  KLGauss(ERMs[,3],ERMs[,2], sigma2)
  
  RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  val <- kl_inv_sup(L2, RHS)
  
  return(list(val=val,KL=KL, Term1=0, Term2=val,  ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=0, L2=L2, ExL1=0, ExL2=0, RefL1=0, RefL2=0))
}

## Forward + Excess
PBkl_FWEL <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/2
  b <- 1
  a <- -1

  # compute ExTerm1
  KL1 <- KLGauss(ERMs[,3], ERMs[,1], sigma2)
  RHS1 <- (KL1 + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  ExL1 <- mean(Excessloss1)
  LHS1 <- (ExL1-a)/(b-a)
  ExTerm1 <- (b-a)*kl_inv_sup(LHS1, RHS1)+a
  # compute the rate of Excess Loss = 0
  ExL1_0rate <- as.vector(table(Excessloss1)['0'])/nhalf/NMC
  
  # compute RefTerm1
  RefL1 <- mean(Refloss1)
  RefTerm1 <- bin_inv_sup(nhalf, nhalf*RefL1, Ndelta)
  
  # compute the bound
  val <- ExTerm1 + RefTerm1
  return(list(val=val,KL=KL1, Term1=0, Term2=0, ExTerm1=ExTerm1, ExTerm2=0, RefTerm1=RefTerm1, RefTerm2=0,
              L1=0, L2=0, ExL1=ExL1, ExL2=0, RefL1=RefL1, RefL2=0, ExL1_0rate=ExL1_0rate, ExL2_0rate=0))
}

## Backward + Excess
PBkl_BWEL <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/2
  b <- 1
  a <- -1
  
  # compute ExTerm2
  KL2 <- KLGauss(ERMs[,3], ERMs[,2], sigma2)
  RHS2 <- (KL2 + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  ExL2 <- mean(Excessloss2)
  LHS2 <- (ExL2-a)/(b-a)
  ExTerm2 <- (b-a)*kl_inv_sup(LHS2, RHS2)+a
  # compute the rate of Excess Loss = 0
  ExL2_0rate <- as.vector(table(Excessloss2)['0'])/nhalf/NMC
  
  # compute RefTerm2
  RefL2 <- mean(Refloss2)
  RefTerm2 <- bin_inv_sup(nhalf, nhalf*RefL2, Ndelta)
  
  # compute the bound
  val <- ExTerm2 + RefTerm2
  return(list(val=val,KL=KL2, Term1=0, Term2=0, ExTerm1=0, ExTerm2=ExTerm2, RefTerm1=0, RefTerm2=RefTerm2,
              L1=0, L2=0, ExL1=0, ExL2=ExL2, RefL1=0, RefL2=RefL2, ExL1_0rate=0, ExL2_0rate=ExL2_0rate))
}

## Average + Excess
PBkl_Avg <- function(NMC,sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/4
  b <- 1
  a <- -1
  
  # compute reference loss
  loss2 <- t(matrix(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERMs[,2])), nrow=NMC, ncol=nhalf, byrow=TRUE))
  loss1 <- t(matrix(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERMs[,1])), nrow=NMC, ncol=nhalf, byrow=TRUE))
  
  # compute excess loss
  Diffloss2 <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samplesTS))-loss2
  Diffloss1 <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samplesTS))-loss1
  
  # compute Delta(h,h_S2) term
  KL2 <- KLGauss(ERMs[,3], ERMs[,2], sigma2)
  RHS2 <- (KL2 + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  LHS2 <- (mean(Diffloss2)-a)/(b-a)
  Term2 <- (b-a)*kl_inv_sup(LHS2, RHS2)+a
  
  # compute Delta(h,h_S1) term
  KL1 <- KLGauss(ERMs[,3], ERMs[,1], sigma2)
  RHS1 <- (KL1 + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  LHS1 <- (mean(Diffloss1)-a)/(b-a)
  Term1 <- (b-a)*kl_inv_sup(LHS1, RHS1)+a
  
  # compute reference term
  RefTerm1 <- bin_inv_sup(nhalf, nhalf*mean(loss1), Ndelta)
  RefTerm2 <- bin_inv_sup(nhalf, nhalf*mean(loss2), Ndelta)
  
  val <- 0.5*(Term1 + Term2 + RefTerm1 + RefTerm2)
  return(list(val=val,KL=0.5*(KL1+KL2), 
              Term1=Term1, Term2=Term2, RefTerm1=RefTerm1, RefTerm2=RefTerm2,
              ExL1=mean(Diffloss1), ExL2=mean(Diffloss2), RefL1=mean(loss1), RefL2=mean(loss2)))
}