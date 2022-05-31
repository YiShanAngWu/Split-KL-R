######################################################################
## R script for computing bounds                                    ##
## Paper: Split-kl and PAC-Bayes-split-kl Inequality                ##
######################################################################

rm(list=ls())
library(MASS)
library(tictoc)
library(lamW)
library(Matrix)
library(pracma)
library(dplyr)
library(tidyr)
library(stringr)
library(caret)
library(R.matlab)
par(mfrow=c(1,1))
path <- "C:\\Users\\lrc379\\OneDrive - University of Copenhagen\\Desktop\\Projects\\PAC-Bayes\\PAC-Bayes-ShiftedKL\\Split-KL-R"

## Experimental setup
set.seed(123)
data_option = "haberman"          # Options are: "sigmoid-synthetic", "haberman", "breast-cancer", 
                                  # "tictactoe", "bank-notes", "kr-vs-kp", "spam", "svmguide1", "mushroom", "adults"
problem_type = "classification"   # No other problem type is supported currently  
distribution <- "gaussian"        # No other distribution is supported currently 
lambda <- 0.01                    # Regularization for the logistic regression
initsigma2 <- 0.5                 # Prior variance for the Gaussian-ERM distribution 

## Global constants
rho <- 2            # This is for the grid-eta discretization
delta <- 0.05       # The probability Threshold for the bound
NMC <- 100          # how many Monte Carlo iterations for the Gaussian expectation?

## End of Experimental Setup
source(paste(path, "gen-data.R", sep="/"))
source(paste(path, "utils.R", sep="/"))
source(paste(path, "PBkl.R", sep="/"))
source(paste(path, "PBUB.R", sep="/"))
source(paste(path, "PBSkl.R", sep="/"))
str <- paste(c(data_option),collapse='-')

## Initializing 
nbRepet <- 2
bound <- array(dim = c(nbRepet,3), data = Inf)
Lntrain <- Lntest <- bestSigma2 <- KL <- array(dim = c(nbRepet,3), data = NA) #posterior
LnERMtrain <- LnERMtest <- array(dim = c(nbRepet), data = NA) # center
Term1 <- Term2 <- ExTerm1 <- ExTerm2 <- RefTerm1 <- RefTerm2 <- array(dim = c(nbRepet,3), data = NA)
L1 <- L2 <- ExL1 <- ExL2 <- RefL1 <- RefL2 <- ExL1_0rate <- ExL2_0rate <- array(dim = c(nbRepet,3), data = NA)
VnTermExL <- LinTermExL <- array(dim = c(nbRepet), data = NA) # for PBUB
L1P <- L1M <- L2P <- L2M <- ExL1P <- ExL1M <- ExL2P <- ExL2M <- muOpt <- etaOpt1 <- etaOpt2 <- array(dim = c(nbRepet), data = NA) # for Skl

## Experiment
for(irepet in 1:nbRepet){
  ## Generate data
  if(irepet%%5 ==1){
    tmp <- gendata(option = data_option)
    Xfull <- tmp$X
    Yfull <- tmp$Y
    d <- tmp$d
    nfull <- length(Yfull)
    if(irepet==1)
      print(c("Size of the full data=",nfull))
  }
  # 5-fold train-test split 
  ifelse(irepet %% 5==0,ind <- 4,ind <- irepet %% 5 - 1)
  ifelse((nfull-floor(nfull/5))%% 2==0, ntest<-floor(nfull/5), ntest <-ceil(nfull/5))
  testIND <- c((ind*ntest+1):min((ind+1)*ntest,nfull))
  Xtrain <- Xfull[-testIND,]
  Ytrain <- Yfull[-testIND]
  Xtest <- Xfull[testIND,]
  Ytest <- Yfull[testIND]

  ntrain <- length(Ytrain)
  if(irepet==1)
    print(c("Size of the training data=",ntrain))
  
  # Build grid of posterior variance for the Gaussian-ERM distribution 
  sigma2GridSize <- ceil(log(ntrain)/log(2))
  sigma2Grid <- numeric(sigma2GridSize)
  for(jj in 1:sigma2GridSize)
    sigma2Grid[jj] <- 1/(2^(jj))
  
  ## Compute ERMs
  ERMs <- buildERMsequenceFast()
  
  ## Compute the train/test error of ERMs[,3]
  LnERMtrain[irepet] <- mean(loss(Ytrain,predictor(Xtrain,ERMs[,3])))
  LnERMtest[irepet] <- mean(loss(Ytest,predictor(Xtest,ERMs[,3])))

  ## Loop over the grid of sigma2 and pick the best one for each method
  for(sigma2 in sigma2Grid){
    ## Compute the shared empirical error
    theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,3],variance2=sigma2, NMC)
    Ln <- mean(loss(Ytrain,predictor(Xtrain,theta_samplesTS)))
    
    ## Compute Excess Losses & Reference Losses
    ### TESTING: For FW/BW without Excess Loss
    Loss1 <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samplesTS))
    Loss2 <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samplesTS))
    ### For Excess Loss
    Refloss1 <- t(matrix(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERMs[,1])), nrow=NMC, ncol=ntrain/2, byrow=TRUE))
    Refloss2 <- t(matrix(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERMs[,2])), nrow=NMC, ncol=ntrain/2, byrow=TRUE))
    Excessloss1 <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samplesTS))-Refloss1
    Excessloss2 <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samplesTS))-Refloss2

    ## PBkl bound
    tmpBPBkl <- PBkl_Avg(NMC, sigma2)
    if(tmpBPBkl$val < bound[irepet,1]){
      bound[irepet,1] <-  tmpBPBkl$val
      Term1[irepet,1] <-  tmpBPBkl$Term1
      Term2[irepet,1] <-  tmpBPBkl$Term2
      ExTerm1[irepet,1] <-  tmpBPBkl$ExTerm1
      ExTerm2[irepet,1] <-  tmpBPBkl$ExTerm2
      RefTerm1[irepet,1] <-  tmpBPBkl$RefTerm1
      RefTerm2[irepet,1] <-  tmpBPBkl$RefTerm2
      L1[irepet,1] <-  tmpBPBkl$L1
      L2[irepet,1] <-  tmpBPBkl$L2
      ExL1[irepet,1] <-  tmpBPBkl$ExL1
      ExL2[irepet,1] <-  tmpBPBkl$ExL2
      ExL1_0rate[irepet,1] <- tmpBPBkl$ExL1_0rate
      ExL2_0rate[irepet,1] <- tmpBPBkl$ExL2_0rate
      RefL1[irepet,1] <-  tmpBPBkl$RefL1
      RefL2[irepet,1] <-  tmpBPBkl$RefL2
      KL[irepet,1] <- tmpBPBkl$KL
      bestSigma2[irepet,1] <- sigma2
    }
    
    ## PBUB bound 
    tmpBPBUB <- PBUB_AvgEx(NMC,sigma2)
    if(tmpBPBUB$val < bound[irepet,2]){
      bound[irepet,2] <- tmpBPBUB$val
      Term1[irepet,2] <-  tmpBPBUB$Term1
      Term2[irepet,2] <-  tmpBPBUB$Term2
      ExTerm1[irepet,2] <-  tmpBPBUB$ExTerm1
      ExTerm2[irepet,2] <-  tmpBPBUB$ExTerm2
      RefTerm1[irepet,2] <-  tmpBPBUB$RefTerm1
      RefTerm2[irepet,2] <-  tmpBPBUB$RefTerm2
      L1[irepet,2] <-  tmpBPBUB$L1
      L2[irepet,2] <-  tmpBPBUB$L2
      ExL1[irepet,2] <-  tmpBPBUB$ExL1
      ExL2[irepet,2] <-  tmpBPBUB$ExL2
      ExL1_0rate[irepet,2] <- tmpBPBUB$ExL1_0rate
      ExL2_0rate[irepet,2] <- tmpBPBUB$ExL2_0rate
      RefL1[irepet,2] <-  tmpBPBUB$RefL1
      RefL2[irepet,2] <-  tmpBPBUB$RefL2
      KL[irepet,2] <- tmpBPBUB$KL
      VnTermExL[irepet] <- tmpBPBUB$VnTermExL
      LinTermExL[irepet] <- tmpBPBUB$LinTermExL
      bestSigma2[irepet,2] <- sigma2
      etaOpt1[irepet] <- tmpBPBUB$etaOpt1
      etaOpt2[irepet] <- tmpBPBUB$etaOpt2
    }
    
    ## PBSkl bound
    tmpBPBSkl <- PBSkl_AvgEx(NMC, sigma2)
    if(tmpBPBSkl$val < bound[irepet,3]){
      bound[irepet,3] <-  tmpBPBSkl$val
      Term1[irepet,3] <-  tmpBPBSkl$Term1
      Term2[irepet,3] <-  tmpBPBSkl$Term2
      ExTerm1[irepet,3] <-  tmpBPBSkl$ExTerm1
      ExTerm2[irepet,3] <-  tmpBPBSkl$ExTerm2
      RefTerm1[irepet,3] <-  tmpBPBSkl$RefTerm1
      RefTerm2[irepet,3] <-  tmpBPBSkl$RefTerm2
      L1[irepet,3] <-  tmpBPBSkl$L1
      L2[irepet,3] <-  tmpBPBSkl$L2
      ExL1[irepet,3] <-  tmpBPBSkl$ExL1
      ExL2[irepet,3] <-  tmpBPBSkl$ExL2
      ExL1_0rate[irepet,3] <- tmpBPBSkl$ExL1_0rate
      ExL2_0rate[irepet,3] <- tmpBPBSkl$ExL2_0rate
      RefL1[irepet,3] <-  tmpBPBSkl$RefL1
      RefL2[irepet,3] <-  tmpBPBSkl$RefL2
      L1P[irepet] <- tmpBPBSkl$L1P
      L1M[irepet] <- tmpBPBSkl$L1M
      L2P[irepet] <- tmpBPBSkl$L2P
      L2M[irepet] <- tmpBPBSkl$L2M
      ExL1P[irepet] <- tmpBPBSkl$ExL1P
      ExL1M[irepet] <- tmpBPBSkl$ExL1M
      ExL2P[irepet] <- tmpBPBSkl$ExL2P
      ExL2M[irepet] <- tmpBPBSkl$ExL2M
      KL[irepet,3] <- tmpBPBSkl$KL
      bestSigma2[irepet,3] <- sigma2
      muOpt[irepet] <- tmpBPBSkl$muOpt
    }
  }
  for(mtd in 1:3){
    theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,3],variance2=bestSigma2[irepet,mtd], NMC)
    Lntrain[irepet,mtd] <- mean(loss(Ytrain,predictor(Xtrain,theta_samplesTS)))
    Lntest[irepet,mtd] <- mean(loss(Ytest,predictor(Xtest,theta_samplesTS)))
  }
}

## Prepare the table to write
Table <- data.frame(LnERMtrain,LnERMtest,Lntrain,Lntest,bestSigma2,bound,KL,
                    L1, Term1, L2, Term2,
                    ExL1,ExTerm1,ExL2,ExTerm2,ExL1_0rate,ExL2_0rate,
                    RefL1,RefTerm1,RefL2,RefTerm2,
                    L1P,L1M,L2P,L2M,
                    ExL1P,ExL1M,ExL2P,ExL2M,VnTermExL,LinTermExL,
                    muOpt,etaOpt1,etaOpt2)
colnames(Table) <- c("LnERMtrain","LnERMtest",
                     "Lntrain_PBkl","Lntrain_PBUB","Lntrain_PBSkl",
                     "Lntest_PBkl","Lntest_PBUB","Lntest_PBSkl",
                     "bestSigma2_PBkl","bestSigma2_PBUB","bestSigma2_PBSkl",
                     "bound_PBkl","bound_PBUB","bound_PBSkl",
                     "KL_PBkl","KL_PBUB","KL_PBSkl",
                     "L1_PBkl","L1_PBUB","L1_PBSkl",
                     "Term1_PBkl","Term1_PBUB","Term1_PBSkl",
                     "L2_PBkl","L2_PBUB","L2_PBSkl",
                     "Term2_PBkl","Term2_PBUB","Term2_SPBkl",
                     "ExL1_PBkl","ExL1_PBUB","ExL1_PBSkl",
                     "ExTerm1_PBkl","ExTerm1_PBUB","ExTerm1_PBSkl",
                     "ExL2_PBkl","ExL2_PBUB","ExL2_PBSkl",
                     "ExTerm2_PBkl","ExTerm2_PBUB","ExTerm2_PBSkl",
                     "ExL1_0rate_PBkl","ExL1_0rate_PBUB","ExL1_0rate_PBSkl",
                     "ExL2_0rate_PBkl","ExL2_0rate_PBUB","ExL2_0rate_PBSkl",
                     "RefL1_PBkl","RefL1_PBUB","RefL1_PBSkl",
                     "RefTerm1_PBkl","RefTerm1_PBUB","RefTerm1_PBSkl",
                     "RefL2_PBkl","RefL2_PBUB","RefL2_PBSkl",
                     "RefTerm2_PBkl","RefTerm2_PBUB","RefTerm2_SPBkl",
                     "L1P","L1M","L2P","L2M",
                     "ExL1P","ExL1M","ExL2P","ExL2M",
                     "VnTermExL","LinTermExL",
                     "muOpt","etaOpt1","etaOpt2")

## write
if(!dir.exists("out")){dir.create(file.path(path,"out"), showWarnings=F)}
outpath <- paste(c(path,"\\out\\",data_option,"-Avg.csv"),collapse="")
write.csv(Table,outpath, row.names=FALSE)


## Print
if(TRUE){
  meansbound <- apply(bound, 2, mean)
  varsbound <- apply(bound, 2, var)
  meanstest <- apply(Lntest, 2, mean)
  varstest <- apply(Lntest, 2, var)
  meanssigma <- apply(bestSigma2, 2, mean)
  print(paste(c(data_option, ". ERM test error=", round(mean(LnERMtest),3), " (", round(var(LnERMtest),3), " )"),collapse=""))
  print(paste(c("PBkl bound=", round(meansbound[1],3), " (", round(varsbound[1],3), ") ",
                ", PBUB Bound=",  round(meansbound[2],3), " (", round(varsbound[2],3), ") ",
             ", PBSkl bound=", round(meansbound[3],3), " (", round(varsbound[3],3), ") "
             ),collapse = ""))
  print(paste(c("PBkl test=", round(meanstest[1],3), " (", round(varstest[1],3), ") ",
                ", PBUB test=",  round(meanstest[2],3), " (", round(varstest[2],3), ") ",
                ", PBSkl test=", round(meanstest[3],3), " (", round(varstest[3],3), ") "
                ),collapse = ""))
  print(paste(c("PBkl sigma=", round(meanssigma[1],3),
                ", PBUB sigma=",  round(meanssigma[2],3),
                ", PBSkl sigma=", round(meanssigma[3],3)
                ),collapse = ""))
}




  

