# PAC-Bayes Split-kl bound
## PAC-Bayes split-kl bound without informed prior
### r.v.\in[0,1]. Doesn't make sense to set \mu=1/2, but if \mu=0, the result is the same as boundPBKL
boundSkl <- function(NMC, sigma2){
  ERMfull <- ERMs[,3]
  
  # compute loss
  Lnloss <- loss(Ytrain,predictor(Xtrain,theta_samplesTS))
  
  # split loss
  mu <- 0
  LnlossP <- apply(Lnloss, 1:2, function(x) max(c(0, x-mu)))
  LnlossM <- apply(Lnloss, 1:2, function(x) max(c(0, mu-x)))
  
  # compute complexity term
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMfull,ERMfull)
  
  if(mu==0){
    RHS <- (KL + log(2*sqrt(ntrain)/delta))/ntrain
    PlusLHS <- mean(LnlossP)
    MinusLHS <- 0
  }else{
    RHS <- (KL + log(4*sqrt(ntrain)/delta))/ntrain
    PlusLHS <- mean(LnlossP)/(1-mu)
    MinusLHS <- mean(LnlossM)/(mu-0)
  }
  # compute Plus Term
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  
  # compute Minus Term
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  val <- mu + (1-mu)*PlusTerm - (mu-0)*MinusTerm
  return(list(val=val, KL=KL, PlusTerm=PlusTerm, MinusTerm=MinusTerm))
}

## PAC-Bayes split-kl bound on half the data  (trained on S1, bound on S2)
### r.v.\in[0,1]. Doesn't make sense to set \mu=1/2, but if \mu=0, the result is the same as boundPBKL_half
boundSkl_half <- function(NMC, sigma2){

  # compute loss
  nhalf <- ntrain/2
  Lnloss <- loss(Ytrain[(ntrain/2+1):ntrain], predictor(Xtrain[(ntrain/2+1):ntrain,],theta_samplesTS))
  
  # split loss
  mu <- 0
  LnlossP <- apply(Lnloss, 1:2, function(x) max(c(0, x-mu)))
  LnlossM <- apply(Lnloss, 1:2, function(x) max(c(0, mu-x)))
  
  # compute complexity term
  KL <-  KLGauss(ERMs[,3],ERMs[,1], sigma2)
  
  if(mu==0){
    RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/delta))/nhalf
    PlusLHS <- mean(LnlossP)
    MinusLHS <- 0
  }else{
    RHS <- (KL + log(4*sigma2GridSize*sqrt(nhalf)/delta))/nhalf
    PlusLHS <- mean(LnlossP)/(1-mu)
    MinusLHS <- mean(LnlossM)/(mu-0)
  }
  
  # compute Plus Term
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  
  # compute Minus Term
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  val <- mu + (1-mu)*PlusTerm - (mu-0)*MinusTerm
  return(list(val=val, KL=KL, PlusTerm=PlusTerm, MinusTerm=MinusTerm))
}


## PAC-Bayes split-kl bound with informed prior
boundSkl_IF <- function(NMC, sigma2){
  ERMfull <- ERMs[,3]
  ERM1 <- ERMs[,2]
  ERM2 <- ERMs[,1]
  #theta_samplesTS <- get_sample(type = distribution, mean=ERMfull, variance2=sigma2, n_samples=NMC)
  # compute reference loss
  loss1 <- t(matrix(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERM1)), nrow=NMC, ncol=ntrain/2, byrow=TRUE))
  loss2 <- t(matrix(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERM2)), nrow=NMC, ncol=ntrain/2, byrow=TRUE))
  
  # compute excess loss
  Diffloss <- matrix(nrow = ntrain, ncol = NMC, data = NA)
  Diffloss[1:(ntrain/2),] <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samplesTS))-loss1
  Diffloss[(ntrain/2+1):ntrain,] <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samplesTS))-loss2
  
  # split excess loss
  mu <- 0
  DifflossP <- apply(Diffloss, 1:2, function(x) max(c(0, x-mu)))
  DifflossM <- apply(Diffloss, 1:2, function(x) max(c(0, mu-x)))
  
  # compute complexity term
  compTerm <- COMP(ERMfull = ERMs[,3],
                   ERM1 = ERMs[,2],
                   ERM2 = ERMs[,1],
                   sigma2 = sigma2)
  RHS <- (compTerm + 2*log(8*sigma2GridSize*sqrt(ntrain/2)/delta))/ntrain
  
  # compute Plus Term
  PlusLHS <- mean(DifflossP)/(1-mu)
  PlusTerm <- kl_inv_sup(PlusLHS, RHS)
  
  # compute Minus Term
  MinusLHS <- mean(DifflossM)/(mu+1)
  MinusTerm <- kl_inv_inf(MinusLHS, RHS)
  
  # compute reference term
  RefTerm <- bin_inv_sup(ntrain/2, ntrain/2*mean(loss1), delta/4) +  bin_inv_sup(ntrain/2, ntrain/2*mean(loss2), delta/4)
  
  val <- mu + (1-mu)*PlusTerm + (mu+1)*MinusTerm + 0.5*RefTerm
  return(list(val=val, compTerm=compTerm, PlusLHS=PlusLHS, MinusLHS=MinusLHS, RefTerm=RefTerm))
}