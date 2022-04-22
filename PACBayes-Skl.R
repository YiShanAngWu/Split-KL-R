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
  return(list(val=val, KL=KL, PlusTerm=PlusTerm, MinusTerm=MinusTerm, RefTerm1=0, RefTerm2=0))
}

## Forward
### r.v.\in[0,1]. Doesn't make sense to set \mu=1/2, but if \mu=0, the result is the same as boundPBKL_half
PBSkl_FW <- function(NMC, sigma2){
  nhalf <- ntrain/2
  
  # compute loss
  Lnloss <- loss(Ytrain[(ntrain/2+1):ntrain], predictor(Xtrain[(ntrain/2+1):ntrain,],theta_samplesTS))
  
  # split loss
  mu <- 0
  LnlossP <- apply(Lnloss, 1:2, function(x) max(c(0, x-mu)))
  LnlossM <- apply(Lnloss, 1:2, function(x) max(c(0, mu-x)))
  
  # compute complexity term
  KL <-  KLGauss(ERMs[,3],ERMs[,1], sigma2)
  
  if(mu==0){
    Ndelta <- delta
    RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
    PlusLHS <- mean(LnlossP)
    MinusLHS <- 0
  }else{
    Ndelta <- delta/2
    RHS <- (KL + log(4*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
    PlusLHS <- mean(LnlossP)/(1-mu)
    MinusLHS <- mean(LnlossM)/(mu-0)
  }
  
  # compute Plus Term
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  
  # compute Minus Term
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  val <- mu + (1-mu)*PlusTerm - (mu-0)*MinusTerm
  return(list(val=val, KL=KL, PlusTerm=PlusTerm, MinusTerm=MinusTerm, RefTerm1=0, RefTerm2=0))
}

## Backward
### r.v.\in[0,1]. Doesn't make sense to set \mu=1/2, but if \mu=0, the result is the same as boundPBKL_half
PBSkl_BW <- function(NMC, sigma2){
  nhalf <- ntrain/2
  
  # compute loss
  Lnloss <- loss(Ytrain[1:(ntrain/2)], predictor(Xtrain[1:(ntrain/2),],theta_samplesTS))
  
  # split loss
  mu <- 0
  LnlossP <- apply(Lnloss, 1:2, function(x) max(c(0, x-mu)))
  LnlossM <- apply(Lnloss, 1:2, function(x) max(c(0, mu-x)))
  
  # compute complexity term
  KL <-  KLGauss(ERMs[,3],ERMs[,2], sigma2)
  
  if(mu==0){
    Ndelta <- delta
    RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
    PlusLHS <- mean(LnlossP)
    MinusLHS <- 0
  }else{
    Ndelta <- delta/2
    RHS <- (KL + log(4*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
    PlusLHS <- mean(LnlossP)/(1-mu)
    MinusLHS <- mean(LnlossM)/(mu-0)
  }
  
  # compute Plus Term
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  
  # compute Minus Term
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  val <- mu + (1-mu)*PlusTerm - (mu-0)*MinusTerm
  return(list(val=val, KL=KL, PlusTerm=PlusTerm, MinusTerm=MinusTerm, RefTerm1=0, RefTerm2=0))
}

## Forward + Excess
PBSkl_FWEL <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/3
  b <- 1
  a <- -1
  
  # compute reference loss
  losses <- t(matrix(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERMs[,1])), nrow=NMC, ncol=nhalf, byrow=TRUE))
  # compute excess loss
  Diffloss <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samplesTS))-losses
  
  # split excess loss
  mu <- 0
  DifflossP <- apply(Diffloss, 1:2, function(x) max(c(0, x-mu)))
  DifflossM <- apply(Diffloss, 1:2, function(x) max(c(0, mu-x)))
  
  # compute Delta(h,h_S1) term
  ## compute complexity
  KL <- KLGauss(ERMs[,3], ERMs[,1], sigma2)
  RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  ## compute Plus Term
  PlusLHS <- mean(DifflossP)/(b-mu)
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  ## compute Minus Term
  MinusLHS <- mean(DifflossM)/(mu-a)
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  Term1 <- mu + (b-mu)*PlusTerm + (mu-a)*MinusTerm
  
  # RefTerm h_S1
  RefTerm1 <- bin_inv_sup(nhalf, nhalf*mean(losses), Ndelta)
  
  val <- Term1 + RefTerm1
  
  return(list(val=val, KL=KL, Term1=Term1, Term2=0, RefTerm1=RefTerm1, RefTerm2=0,
              ExL1=mean(Diffloss), ExL2=0, RefL1=mean(losses), RefL2=0,
              ExL1P=mean(DifflossP), ExL1M=mean(DifflossM), ExL2P=0, ExL2M=0))
}

## Backward + Excess
PBSkl_BWEL <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/3
  b <- 1
  a <- -1
  
  # compute reference loss
  losses <- t(matrix(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERMs[,2])), nrow=NMC, ncol=nhalf, byrow=TRUE))
  # compute excess loss
  Diffloss <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),], theta_samplesTS))-losses

  # split excess loss
  mu <- 0
  DifflossP <- apply(Diffloss, 1:2, function(x) max(c(0, x-mu)))
  DifflossM <- apply(Diffloss, 1:2, function(x) max(c(0, mu-x)))
  
  # compute Delta(h,h_S2) term
  ## compute complexity
  KL <- KLGauss(ERMs[,3], ERMs[,2], sigma2)
  RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  ## compute Plus Term
  PlusLHS <- mean(DifflossP)/(b-mu)
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  ## compute Minus Term
  MinusLHS <- mean(DifflossM)/(mu-a)
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  Term2 <- mu + (b-mu)*PlusTerm + (mu-a)*MinusTerm
  
  # RefTerm h_S2
  RefTerm2 <- bin_inv_sup(nhalf, nhalf*mean(losses), Ndelta)
  
  val <- Term2 + RefTerm2
  
  return(list(val=val, KL=KL, Term1=0, Term2=Term2, RefTerm1=0, RefTerm2=RefTerm2,
              ExL1=0, ExL2=mean(Diffloss), RefL1=0, RefL2=mean(losses),
              ExL1P=0, ExL1M=0, ExL2P=mean(DifflossP), ExL2M=mean(DifflossM)))
}

## Average + Excess
PBSkl_Avg <- function(NMC, sigma2){
  nhalf <- ntrain/2
  Ndelta <- delta/4
  b <- 1
  a <- -1

  # compute reference loss
  loss2 <- t(matrix(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERMs[,2])), nrow=NMC, ncol=ntrain/2, byrow=TRUE))
  loss1 <- t(matrix(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERMs[,1])), nrow=NMC, ncol=ntrain/2, byrow=TRUE))
  
  # compute excess loss
  Diffloss <- matrix(nrow = ntrain, ncol = NMC, data = NA)
  Diffloss[1:(ntrain/2),] <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samplesTS))-loss2
  Diffloss[(ntrain/2+1):ntrain,] <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samplesTS))-loss1
  
  # split excess loss
  mu <- 0
  DifflossP <- apply(Diffloss, 1:2, function(x) max(c(0, x-mu)))
  DifflossM <- apply(Diffloss, 1:2, function(x) max(c(0, mu-x)))
  
  # compute complexity term
  KL <- 0.5*COMP(ERMfull = ERMs[,3],
                   ERM1 = ERMs[,2],
                   ERM2 = ERMs[,1],
                   sigma2 = sigma2)
  RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/Ndelta))/nhalf
  
  # compute Plus Term
  PlusLHS <- mean(DifflossP)/(b-mu)
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  
  # compute Minus Term
  MinusLHS <- mean(DifflossM)/(mu-a)
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  # compute reference term
  RefTerm1 <- bin_inv_sup(nhalf, nhalf*mean(loss1), Ndelta)
  RefTerm2 <- bin_inv_sup(nhalf, nhalf*mean(loss2), Ndelta)
  
  val <- mu + (1-mu)*PlusTerm + (mu+1)*MinusTerm + 0.5*(RefTerm1+RefTerm2)
  return(list(val=val, KL=KL, Term1=PlusTerm, Term2=MinusTerm, RefTerm1=RefTerm1, RefTerm2=RefTerm2))
}