### Tolstikhin-Seldin bound
boundPBEB <- function(NMC, sigma2){
  c1 <- c2 <- 1.15
  # Computing the KL
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3]) + log(sigma2GridSize)
    
  v1 <- ceiling((1/log(c1))*log(sqrt((exp(1)-2)*ntrain/(4*log(2/delta))))+1)
  v2 <- ceiling((1/log(c2))*log(.5*sqrt((ntrain-1)/(log(2/delta))+1)+.5))
  
  theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,3], variance2=sigma2, NMC)
  Loss <- loss(Ytrain,predictor(Xtrain,theta_samplesTS))
  Ln <-  mean(Loss)
  VarTS <- mean(apply(X = Loss, MARGIN = 2, FUN = var))
  
  Vn <- VarTS + (1+c2)*sqrt(VarTS*(KL+log(2*v2/delta))/(2*(ntrain-1)))+2*c2*(KL+log(2*v2/delta))/(ntrain-1)
  Vnbar <- min(Vn,1/4)
  
  if(sqrt((KL+log(2*v1/delta))/((exp(1)-2)*Vnbar)) <= sqrt(ntrain)){
    val <- (1+c1)*sqrt((exp(1)-2)*Vnbar*(KL+log(2*v1/delta))/ntrain)
  }else{
    val <- 2*(KL+log(2*v1/delta))/ntrain
  }
  return(list(val=val,VarTS=VarTS,KL=KL))
}

### Maurer's bound (PBkl)
boundPBKL <- function(Ln,sigma2){
  # Computing the KL 
  ratio <- initsigma2/(sigma2)
  #KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3]) + log(sigma2GridSize)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3])
  RHS <- (KL + log(2*sqrt(ntrain)/delta))/ntrain
  val <- kl_inv_sup(Ln, RHS)
  
  return(list(val=val,KL=KL))
}

boundCatoni <- function(Ln,sigma2){
  alpha <- rho
  # Computing the KL 
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3]) + log(sigma2GridSize)
  
  D <- KL - log(delta)
  sqrtTerm <- sqrt(2*alpha * D/(ntrain*Ln*(1-Ln)))
  Num <-  1 - exp(-Ln*sqrtTerm - 
                    alpha/ntrain * (D + 2 * log(log(alpha^2*ntrain*sqrtTerm/log(alpha)))))
  Den <- 1 - exp(-sqrtTerm)
  val <- Num/Den
  return(list(val=val,KL=KL))
}


## Bounds on half the data (trained on S1, bound on S2)
### Tolstikhin-Seldin bound
boundPBEB_half <- function(NMC,sigma2){
  c1 <- c2 <- 1.15
  
  nhalf <- ntrain/2
  # Computing the KL 
  KL <-  (1/(2*sigma2))*sum(square_diff(ERMs[,1],ERMs[,2]))+ log(sigma2GridSize)

  v1 <- ceiling((1/log(c1))*log(sqrt((exp(1)-2)*nhalf/(4*log(2/delta))))+1)
  v2 <- ceiling((1/log(c2))*log(.5*sqrt((nhalf-1)/(log(2/delta))+1)+.5))
  
  theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,2],variance2=sigma2, NMC)
  Loss <- loss(Ytrain[1:nhalf],predictor(Xtrain[1:nhalf,],theta_samplesTS))
  Ln <-  mean(Loss)
  VarTS <- mean(apply(X = Loss, MARGIN = 2, FUN = var))
  
  Vn <- VarTS + (1+c2)*sqrt(VarTS*(KL+log(2*v2/delta))/(2*(nhalf-1)))+2*c2*(KL+log(2*v2/delta))/(ntrain-1)
  Vnbar <- min(Vn,1/4)
  
  if(sqrt((KL+log(2*v1/delta))/((exp(1)-2)*Vnbar)) <= sqrt(nhalf)){
    val <- (1+c1)*sqrt((exp(1)-2)*Vnbar*(KL+log(2*v1/delta))/nhalf)
  }else{
    val <- 2*(KL+log(2*v1/delta))/nhalf
  }
  return(list(val=val,VarTS=VarTS,KL=KL))
}

### Maurer's bound
boundPBKL_half <- function(NMC,sigma2){
  theta_samplesTS <- get_sample(type = distribution, variance2=sigma2, mean=ERMs[,3], NMC)
  nhalf <- ntrain/2
  Ln <- mean(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],theta_samplesTS)))
  # Computing the KL (KL(rho, pi_S1))
  KL <-  KLGauss(ERMs[,3],ERMs[,1], sigma2)
  
  RHS <- (KL + log(2*sigma2GridSize*sqrt(nhalf)/delta))/nhalf
  val <- kl_inv_sup(Ln, RHS)
  
  return(list(val=val,KL=KL))
}

boundCatoni_half <- function(NMC,sigma2){
  alpha <- 2
  theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,2],variance2=sigma2, NMC)
  nhalf <- ntrain/2
  Ln <- mean(loss(Ytrain[1:nhalf],predictor(Xtrain[1:nhalf,],theta_samplesTS)))
  
  # Computing the KL 
  KL <-(1/(2*sigma2))*sum(square_diff(ERMs[,1],ERMs[,2]))+ log(sigma2GridSize)
  D <- KL - log(delta)
  
  sqrtTerm <- sqrt(2*alpha * D/(nhalf*Ln*(1-Ln)))
  Num <-  1 - exp(-Ln*sqrtTerm - alpha/nhalf * (D + 2 * log(log(alpha^2*nhalf*sqrtTerm)/log(alpha))))
  Den <- 1 - exp(-sqrtTerm)
  val <- Num/Den
  return(list(val=val,KL=KL))
}

## Bounds with average informed prior

### Maurer's bound
boundPBKL_IF <- function(NMC,sigma2){
  nhalf <- ntrain/2
  ERMfull <- ERMs[,3]
  ERM1 <- ERMs[,2]
  ERM2 <- ERMs[,1]
  theta_samples <- get_sample(type = distribution, mean=ERMfull, variance2=sigma2, n_samples=NMC)
  
  # compute reference loss
  loss1 <- t(matrix(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERM1)), nrow=NMC, ncol=nhalf, byrow=TRUE))
  loss2 <- t(matrix(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERM2)), nrow=NMC, ncol=nhalf, byrow=TRUE))
  
  # compute excess loss
  Diffloss1 <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samples))-loss1
  Diffloss2 <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samples))-loss2

  # compute Delta(h,h_S2) term
  KL1 <- KLGauss(ERMfull, ERM1, sigma2)
  RHS1 <- (KL1 + log(8*sigma2GridSize*sqrt(nhalf)/delta))/nhalf
  LHS1 <- (mean(Diffloss1)+1)/2
  Term1 <- 2*kl_inv_sup(LHS1, RHS1)-1

  # compute Delta(h,h_S1) term
  KL2 <- KLGauss(ERMfull, ERM2, sigma2)
  RHS2 <- (KL2 + log(8*sigma2GridSize*sqrt(nhalf)/delta))/nhalf
  LHS2 <- (mean(Diffloss2)+1)/2
  Term2 <- 2*kl_inv_sup(LHS2, RHS2)-1  

  # compute reference term
  RefTerm <- bin_inv_sup(nhalf, nhalf*mean(loss1), delta/4) +  bin_inv_sup(nhalf, nhalf*mean(loss2), delta/4)
  
  val <- 0.5*(Term1 + Term2 + RefTerm)
  return(list(val=val,KL=(KL1+KL2), Term1=Term1, Term2=Term2, RefTerm=RefTerm))
}