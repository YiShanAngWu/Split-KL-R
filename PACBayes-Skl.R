# PAC-Bayes Split-kl bound

## Vanilla
### r.v.\in[0,1]. Doesn't make sense to set \mu=1/2, but if \mu=0, the result is the same as boundPBKL
PBSkl <- function(NMC, sigma2){
  # compute loss
  Lnloss <- loss(Ytrain,predictor(Xtrain,theta_samplesTS))
  
  # split loss
  mu <- 0
  LnlossP <- apply(Lnloss, 1:2, function(x) max(c(0, x-mu)))
  LnlossM <- apply(Lnloss, 1:2, function(x) max(c(0, mu-x)))
  
  # compute complexity term
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  
  if(mu==0){
    Ndelta <- delta
    RHS <- (KL + log(2*sqrt(ntrain)/Ndelta))/ntrain
    PlusLHS <- mean(LnlossP)
    MinusLHS <- 0
  }else{
    Ndelta <- delta/2
    RHS <- (KL + log(4*sqrt(ntrain)/Ndelta))/ntrain
    PlusLHS <- mean(LnlossP)/(1-mu)
    MinusLHS <- mean(LnlossM)/(mu-0)
  }
  # compute Plus Term
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  
  # compute Minus Term
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  val <- mu + (1-mu)*PlusTerm - (mu-0)*MinusTerm
  return(list(val=val, KL=KL, Term1=val, Term2=0,  ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=mean(Lnloss), L2=0, ExL1=0, ExL2=0, RefL1=0, RefL2=0, ExL1_0rate=0, ExL2_0rate=0,
              L1P=mean(LnlossP), L1M=mean(LnlossM), L2P=0, L2M=0, ExL1P=0, ExL1M=0, ExL2P=0, ExL2M=0, muOpt=mu))
}

## Forward
### r.v.\in[0,1]. Doesn't make sense to set \mu=1/2, but if \mu=0, the result is the same as boundPBKL_half
PBSkl_FW <- function(NMC, sigma2){
  nhalf <- ntrain/2
  
  # split loss
  mu <- 0
  Loss1P <- apply(Loss1, 1:2, function(x) max(c(0, x-mu)))
  Loss1M <- apply(Loss1, 1:2, function(x) max(c(0, mu-x)))
  
  # compute complexity term
  KL <-  KLGauss(ERMs[,3],ERMs[,1], sigma2)
  
  if(mu==0){
    Ndelta <- delta
    RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
    PlusLHS <- mean(Loss1P)
    MinusLHS <- 0
  }else{
    Ndelta <- delta/2
    RHS <- (KL + log(4*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
    PlusLHS <- mean(Loss1P)/(1-mu)
    MinusLHS <- mean(Loss1M)/(mu-0)
  }
  
  # compute Plus Term
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  
  # compute Minus Term
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  val <- mu + (1-mu)*PlusTerm - (mu-0)*MinusTerm
  return(list(val=val, KL=KL, Term1=val, Term2=0,  ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=mean(Loss1), L2=0, ExL1=0, ExL2=0, RefL1=0, RefL2=0, ExL1_0rate=0, ExL2_0rate=0,
              L1P=mean(Loss1P), L1M=mean(Loss1M), L2P=0, L2M=0, ExL1P=0, ExL1M=0, ExL2P=0, ExL2M=0, muOpt=mu))
}

## Backward
### r.v.\in[0,1]. Doesn't make sense to set \mu=1/2, but if \mu=0, the result is the same as boundPBKL_half
PBSkl_BW <- function(NMC, sigma2){
  nhalf <- ntrain/2
  
  # split loss
  mu <- 0
  Loss2P <- apply(Loss2, 1:2, function(x) max(c(0, x-mu)))
  Loss2M <- apply(Loss2, 1:2, function(x) max(c(0, mu-x)))
  
  # compute complexity term
  KL <-  KLGauss(ERMs[,3],ERMs[,2], sigma2)
  
  if(mu==0){
    Ndelta <- delta
    RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
    PlusLHS <- mean(Loss2P)
    MinusLHS <- 0
  }else{
    Ndelta <- delta/2
    RHS <- (KL + log(4*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
    PlusLHS <- mean(Loss2P)/(1-mu)
    MinusLHS <- mean(Loss2M)/(mu-0)
  }
  
  # compute Plus Term
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  
  # compute Minus Term
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  val <- mu + (1-mu)*PlusTerm - (mu-0)*MinusTerm
  return(list(val=val, KL=KL, Term1=0, Term2=val,  ExTerm1=0, ExTerm2=0, RefTerm1=0, RefTerm2=0,
              L1=0, L2=mean(Loss2), ExL1=0, ExL2=0, RefL1=0, RefL2=0, ExL1_0rate=0, ExL2_0rate=0,
              L1P=0, L1M=0, L2P=mean(Loss2P), L2M=mean(Loss2M), ExL1P=0, ExL1M=0, ExL2P=0, ExL2M=0, muOpt=mu))
}

## Forward + Excess
PBSkl_FWEL <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/3
  b <- 1
  a <- -1
  
  # split excess loss
  mu <- 0
  Excessloss1P <- apply(Excessloss1, 1:2, function(x) max(c(0, x-mu)))
  Excessloss1M <- apply(Excessloss1, 1:2, function(x) max(c(0, mu-x)))
  # only for stats
  ExL1 <- mean(Excessloss1)
  ExL1_0rate <- as.vector(table(Excessloss1)['0'])/nhalf/NMC
  
  # compute ExTerm1
  ## compute complexity
  KL <- KLGauss(ERMs[,3], ERMs[,1], sigma2)
  # Computing the KL 
  #ratio <- initsigma2/(sigma2)
  #KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  
  #RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  RHS <- (KL + log(sigma2GridSize/Ndelta))/nhalf # conjecture
  #RHS <- (KL + log(2*sqrt(nhalf)/Ndelta))/nhalf # no informed prior
  #RHS <- (KL + log(1/Ndelta))/nhalf # no informed prior, conjecture
  
  ## compute Plus Term
  PlusLHS <- mean(Excessloss1P)/(b-mu)
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  ## compute Minus Term
  MinusLHS <- mean(Excessloss1M)/(mu-a)
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  ExTerm1 <- mu + (b-mu)*PlusTerm - (mu-a)*MinusTerm
  ExTerm1 <- min(c(b, ExTerm1))

  # RefTerm1
  RefL1 <- mean(Refloss1)
  RefTerm1 <- bin_inv_sup(nhalf, nhalf*RefL1, Ndelta)

  val <- ExTerm1 + RefTerm1
  
  return(list(val=val,KL=KL, Term1=0, Term2=0, ExTerm1=ExTerm1, ExTerm2=0, RefTerm1=RefTerm1, RefTerm2=0,
              L1=0, L2=0, ExL1=ExL1, ExL2=0, RefL1=RefL1, RefL2=0, ExL1_0rate=ExL1_0rate, ExL2_0rate=0,
              L1P=0, L1M=0, L2P=0, L2M=0, ExL1P=mean(Excessloss1P), ExL1M=mean(Excessloss1M), ExL2P=0, ExL2M=0, muOpt=mu))
}

## Backward + Excess
PBSkl_BWEL <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/3
  b <- 1
  a <- -1
  
  # split excess loss
  mu <- 0
  Excessloss2P <- apply(Excessloss2, 1:2, function(x) max(c(0, x-mu)))
  Excessloss2M <- apply(Excessloss2, 1:2, function(x) max(c(0, mu-x)))
  # only for stats
  ExL2 <- mean(Excessloss2)
  ExL2_0rate <- as.vector(table(Excessloss2)['0'])/nhalf/NMC
  
  # compute ExTerm2
  ## compute complexity
  KL <- KLGauss(ERMs[,3], ERMs[,2], sigma2)
  RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  ## compute Plus Term
  PlusLHS <- mean(Excessloss2P)/(b-mu)
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  ## compute Minus Term
  MinusLHS <- mean(Excessloss2M)/(mu-a)
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  ExTerm2 <- mu + (b-mu)*PlusTerm - (mu-a)*MinusTerm
  ExTerm2 <- min(c(b, ExTerm2))
  
  # RefTerm2
  RefL2 <- mean(Refloss2)
  RefTerm2 <- bin_inv_sup(nhalf, nhalf*RefL2, Ndelta)

  val <- ExTerm2 + RefTerm2
  return(list(val=val,KL=KL, Term1=0, Term2=0, ExTerm1=0, ExTerm2=ExTerm2, RefTerm1=0, RefTerm2=RefTerm2,
              L1=0, L2=0, ExL1=0, ExL2=ExL2, RefL1=0, RefL2=RefL2, ExL1_0rate=0, ExL2_0rate=ExL2_0rate,
              L1P=0, L1M=0, L2P=0, L2M=0, ExL1P=0, ExL1M=0, ExL2P=mean(Excessloss2P), ExL2M=mean(Excessloss2M), muOpt=mu))
}

## Average + Excess
PBSkl_AvgEx <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/4
  b <- 1
  a <- -1

  # split excess loss
  mu <- 0
  Excessloss1P <- apply(Excessloss1, 1:2, function(x) max(c(0, x-mu)))
  Excessloss1M <- apply(Excessloss1, 1:2, function(x) max(c(0, mu-x)))
  Excessloss2P <- apply(Excessloss2, 1:2, function(x) max(c(0, x-mu)))
  Excessloss2M <- apply(Excessloss2, 1:2, function(x) max(c(0, mu-x)))
  # only for stats
  ExL1 <- mean(Excessloss1)
  ExL2 <- mean(Excessloss2)
  ExL1_0rate <- as.vector(table(Excessloss1)['0'])/nhalf/NMC
  ExL2_0rate <- as.vector(table(Excessloss2)['0'])/nhalf/NMC
  
  # compute complexity term
  ## 1. Informed Prior with 2sqrt(n)
  #KL <- 0.5*(KLGauss(ERMs[,3],ERMs[,1], sigma2) + KLGauss(ERMs[,3],ERMs[,2], sigma2))
  #RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  ## 2. Informed Prior without 2sqrt(n)
  #KL <- 0.5*(KLGauss(ERMs[,3],ERMs[,1], sigma2) + KLGauss(ERMs[,3],ERMs[,2], sigma2))
  #RHS <- (KL + log(sigma2GridSize/Ndelta))/nhalf
  ## 3. NO Informed Prior with 2sqrt(n)
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  RHS <- (KL + log(2*sqrt(nhalf)/Ndelta))/nhalf
  ## 4. NO Informed Prior without 2sqrt(n)
  #ratio <- initsigma2/(sigma2)
  #KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  #RHS <- (KL + log(1/Ndelta))/nhalf
  
  # compute Plus Term
  ExL1P <- mean(Excessloss1P)
  ExL2P <- mean(Excessloss2P)
  PlusLHS <- (0.5*ExL1P+0.5*ExL2P)/(b-mu)
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  
  # compute Minus Term
  ExL1M <- mean(Excessloss1M)
  ExL2M <- mean(Excessloss2M)
  MinusLHS <- (0.5*ExL1M+0.5*ExL2M)/(mu-a)
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  # compute reference term
  RefL1 <- mean(Refloss1)
  RefL2 <- mean(Refloss2)
  RefTerm1 <- bin_inv_sup(nhalf, nhalf*RefL1, Ndelta)
  RefTerm2 <- bin_inv_sup(nhalf, nhalf*RefL2, Ndelta)
  
  ExTerm <- mu + (b-mu)*PlusTerm - (mu-a)*MinusTerm 
  ExTerm <- min(c(b, ExTerm))
  val <- ExTerm + 0.5*(RefTerm1+RefTerm2)
  
  return(list(val=val,KL=KL, Term1=0, Term2=0, ExTerm1=PlusTerm, ExTerm2=MinusTerm, RefTerm1=RefTerm1, RefTerm2=RefTerm2,
              L1=0, L2=0, ExL1=ExL1, ExL2=ExL2, RefL1=RefL1, RefL2=RefL2, ExL1_0rate=ExL1_0rate, ExL2_0rate=ExL2_0rate,
              L1P=0, L1M=0, L2P=0, L2M=0, ExL1P=ExL1P, ExL1M=ExL1M, ExL2P=ExL2P, ExL2M=ExL2M, muOpt=mu))
}