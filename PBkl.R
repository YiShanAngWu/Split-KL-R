# Maurer's bound (PBkl)

## Vanilla
PBkl <- function(NMC, sigma2){
  Ndelta <- delta

  # Computing the KL 
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  # Computing the complexity
  RHS <- (KL + log(2*sqrt(ntrain)/Ndelta))/ntrain
  # Computing the bound
  val <- kl_inv_sup(Ln, RHS)
  return(list(val=val,KL=KL, Term1=val, Term2=0, ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=Ln, L2=0, ExL1=0, ExL2=0, RefL1=0, RefL2=0, ExL1_0rate=0, ExL2_0rate=0))
}

## Forward
PBkl_FW <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta
  L1 <- mean(Loss1)
  
  # Computing the KL
  ## 1. Informed Prior with 2sqrt(n)
  #KL <-  KLGauss(ERMs[,3],ERMs[,1], sigma2)
  #RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  ## 2. Informed Prior without 2sqrt(n)
  KL <-  KLGauss(ERMs[,3],ERMs[,1], sigma2)
  RHS <- (KL + log(sigma2GridSize/Ndelta))/nhalf
  ## 3. NO Informed Prior with 2sqrt(n)
  #ratio <- initsigma2/(sigma2)
  #KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  #RHS <- (KL + log(2*sqrt(nhalf)/Ndelta))/nhalf
  ## 4. NO Informed Prior without 2sqrt(n)
  #ratio <- initsigma2/(sigma2)
  #KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  #RHS <- (KL + log(1/Ndelta))/nhalf
  
  val <- kl_inv_sup(L1, RHS)
  return(list(val=val,KL=KL, Term1=val, Term2=0, ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=L1, L2=0, ExL1=0, ExL2=0, RefL1=0, RefL2=0, ExL1_0rate=0, ExL2_0rate=0))
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
              L1=0, L2=L2, ExL1=0, ExL2=0, RefL1=0, RefL2=0, ExL1_0rate=0, ExL2_0rate=0))
}

## Average
PBkl_Avg <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta
  
  L1 <- mean(Loss1)
  L2 <- mean(Loss2)

  # Computing the KL
  ## 1. Informed Prior with 2sqrt(n)
  #KL <-  0.5*(KLGauss(ERMs[,3],ERMs[,1], sigma2) + KLGauss(ERMs[,3],ERMs[,2], sigma2))
  #RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  ## 2. Informed Prior without 2sqrt(n)
  #KL <-  0.5*(KLGauss(ERMs[,3],ERMs[,1], sigma2) + KLGauss(ERMs[,3],ERMs[,2], sigma2))
  #RHS <- (KL + log(sigma2GridSize/Ndelta))/nhalf
  ## 3. NO Informed Prior with 2sqrt(n)
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  RHS <- (KL + log(2*sqrt(nhalf)/Ndelta))/nhalf
  ## 4. NO Informed Prior without 2sqrt(n)
  #ratio <- initsigma2/(sigma2)
  #KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  #RHS <- (KL + log(1/Ndelta))/nhalf

  val <- kl_inv_sup(0.5*(L1+L2), RHS)
  
  return(list(val=val,KL=KL, Term1=val, Term2=0,  ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=L1, L2=L2, ExL1=0, ExL2=0, RefL1=0, RefL2=0, ExL1_0rate=0, ExL2_0rate=0))
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
PBkl_AvgEx <- function(NMC,sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/4
  b <- 1
  a <- -1

  # compute ExTerm2
  KL2 <- KLGauss(ERMs[,3], ERMs[,2], sigma2)
  RHS2 <- (KL2 + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  ExL2 <- mean(Excessloss2)
  LHS2 <- (ExL2-a)/(b-a)
  ExTerm2 <- (b-a)*kl_inv_sup(LHS2, RHS2)+a
  ExL2_0rate <- as.vector(table(Excessloss2)['0'])/nhalf/NMC
  
  # compute ExTerm1
  KL1 <- KLGauss(ERMs[,3], ERMs[,1], sigma2)
  RHS1 <- (KL1 + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  ExL1 <- mean(Excessloss1)
  LHS1 <- (ExL1-a)/(b-a)
  ExTerm1 <- (b-a)*kl_inv_sup(LHS1, RHS1)+a
  ExL1_0rate <- as.vector(table(Excessloss1)['0'])/nhalf/NMC
  
  # compute reference term
  RefL1 <- mean(Refloss1)
  RefTerm1 <- bin_inv_sup(nhalf, nhalf*RefL1, Ndelta)
  RefL2 <- mean(Refloss2)
  RefTerm2 <- bin_inv_sup(nhalf, nhalf*RefL2, Ndelta)
  
  val <- 0.5*(ExTerm1 + ExTerm2 + RefTerm1 + RefTerm2)
  
  return(list(val=val,KL=(KL1+KL2)/2, Term1=0, Term2=0, ExTerm1=ExTerm1, ExTerm2=ExTerm2, RefTerm1=RefTerm1, RefTerm2=RefTerm2,
              L1=0, L2=0, ExL1=ExL1, ExL2=ExL2, RefL1=RefL1, RefL2=RefL2, ExL1_0rate=ExL1_0rate, ExL2_0rate=ExL2_0rate))
}