######################################################################
## R script for computing bounds                                    ##
## Preprint: PAC-Bayes Un-Expected Bernstein Inequality             ##
## URL: https://arxiv.org/abs/1905.13367		            ##
## Authors (alphabetical order): Benjamin Guedj & Zakaria Mhammedi  ## 
######################################################################

rm(list=ls())
library(MASS)
library(tictoc)
library(lamW)
library(Matrix)
library(pracma)
library(stringr)
library(caret)
library(R.matlab)

par(mfrow=c(1,1))
path <- "C:\\Users\\lrc379\\OneDrive - University of Copenhagen\\Desktop\\Projects\\PAC-Bayes\\PAC-Bayes-ShiftedKL\\Split-KL-R"

## Experimental setup
set.seed(123)
data_option = "adults"          # Options are: "sigmoid-synthetic", "haberman", "breast-cancer", 
                                  # "tictactoe", "bank-notes", "kr-vs-kp", "spam", "mushroom", "adults"
problem_type = "classification"   # No other problem type is supported currently  
distribution <- "gaussian"        # No other distribution is supported currently 
lambda <- 0.01                    # Regularization for the logistic regression
                      
initsigma2 <- 0.5                 # Prior variance for the Gaussian-ERM distribution 
#half <-  F                        # Set to true so that half the data is used to build a prior for PACBayes-Others (i.e. not our method)
#IF <-  F                          # Set to true so that Informed Prior is used for all methods

## following block is only relevant for synthetic data
nb.seq <- 10
#nmin <- 800
#nmax <- 8000
#nb.grid <- 2*round(seq(nmin,nmax,length.out = nb.seq)/2)
##theta_star <-  c(3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6,2,6,4,3,3,8,3,2,7,9,5,0,2,8,8,4,1,9,7,1,6,9,3,9,9,3,7,5,1)
#theta_star <- c(3,1,4,1,5,9,2,6,5,3)
#d <- length(theta_star)
#noise <- T
#proportion <- 0.1

## Global constants
rho <- 2            # This is for the grid-eta discretization
delta <- 0.05       # The probability Threshold for the bound
NMC <- 100          # how many Monte Carlo iterations for the Gaussian expectation?

## Definining the loss and the predictor type
loss <- function(a,b) abs(a-b)
b <- 1    # upper bound on the loss
predictor <- function(x,theta) round(1/(1+exp(-x%*%theta)))

## End of Experimental Setup
source(paste(path, "gen-data.R", sep="/"))
source(paste(path, "utils.R", sep="/"))
source(paste(path, "PACBayes-kl.R", sep="/"))
source(paste(path, "PACBayes-MGG.R", sep="/"))
source(paste(path, "PACBayes-Skl.R", sep="/"))
#str <- paste(c(data_option,"imformedPrior",IF),collapse='-')
str <- paste(c(data_option),collapse='-')

# Initializing 
nbRepet <- 5
bound <- array(dim = c(nbRepet,3,nb.seq), data = Inf)
Lntrain <- Lntest <- bestSigma2 <- array(dim = c(nbRepet,3,nb.seq), data = NA) #posterior
LnERMtrain <- LnERMtest <- array(dim = c(nbRepet, nb.seq), data = NA) # center
Term1 <- Term2 <- RefTerm1 <- RefTerm2 <- array(dim = c(nbRepet,3,nb.seq), data = NA)
ExL1 <- ExL2 <- RefL1 <- RefL2 <- array(dim = c(nbRepet,3,nb.seq), data = NA)
ExL1P <- ExL1M <- ExL2P <- ExL2M <- array(dim = c(nbRepet, nb.seq), data = NA) # for Skl

pb <- txtProgressBar(min = 0, max = nb.seq, style = 3)
for(inb in 1:nb.seq){
  for(irepet in 1:nbRepet){
    ## Generate data
    if(irepet%%5 ==1 || grepl("synthetic", data_option, fixed=TRUE)){
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
    LnERMtrain[irepet,inb] <- mean(loss(Ytrain,predictor(Xtrain,ERMs[,3])))
    LnERMtest[irepet,inb] <- mean(loss(Ytest,predictor(Xtest,ERMs[,3])))

    ## Loop over the grid of sigma2 and pick the best one for each method
    for(sigma2 in sigma2Grid){
      ## Compute the shared empirical error
      theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,3],variance2=sigma2, NMC)
      Ln <- mean(loss(Ytrain,predictor(Xtrain,theta_samplesTS)))

      ## Maurer bound
      tmpBkl <- PBkl_FWEL(NMC, sigma2)
      #ifelse(IF,tmpBKL<- boundPBKL_half(NMC, sigma2),tmpBKL<-boundPBKL(NMC, sigma2))
      if(tmpBkl$val < bound[irepet,1,inb]){
        bound[irepet,1,inb] <-  tmpBkl$val
        Term1[irepet,1,inb] <-  tmpBkl$Term1
        Term2[irepet,1,inb] <-  tmpBkl$Term2
        RefTerm1[irepet,1,inb] <-  tmpBkl$RefTerm1
        RefTerm2[irepet,1,inb] <-  tmpBkl$RefTerm2
        ExL1[irepet,1,inb] <-  tmpBkl$ExL1
        ExL2[irepet,1,inb] <-  tmpBkl$ExL2
        RefL1[irepet,1,inb] <-  tmpBkl$RefL1
        RefL2[irepet,1,inb] <-  tmpBkl$RefL2
        bestSigma2[irepet,1,inb] <- sigma2
      }
      
      ## MGG bound 
      #ifelse(IF,tmpBMGG<- boundMGG_half(NMC,sigma2),tmpBMGG<-boundMGG(NMC,sigma2))
      tmpBMGG <- MGG_FWEL(NMC,sigma2)
      if(tmpBMGG$val < bound[irepet,2,inb]){
        bound[irepet,2,inb] <- tmpBMGG$val
        Term1[irepet,2,inb] <-  tmpBMGG$Term1
        Term2[irepet,2,inb] <-  tmpBMGG$Term2
        RefTerm1[irepet,2,inb] <-  tmpBMGG$RefTerm1
        RefTerm2[irepet,2,inb] <-  tmpBMGG$RefTerm2
        ExL1[irepet,2,inb] <-  tmpBMGG$ExL1
        ExL2[irepet,2,inb] <-  tmpBMGG$ExL2
        RefL1[irepet,2,inb] <-  tmpBMGG$RefL1
        RefL2[irepet,2,inb] <-  tmpBMGG$RefL2
        bestSigma2[irepet,2,inb] <- sigma2
      }
      
      ## Split-kl bound
      #ifelse(IF, tmpBSkl <-boundSkl_IF(NMC, sigma2), tmpBSkl <-boundSkl(NMC, sigma2))
      tmpBSkl <- PBSkl_FWEL(NMC, sigma2)
      if(tmpBSkl$val < bound[irepet,3,inb]){
        bound[irepet,3,inb] <-  tmpBSkl$val
        Term1[irepet,3,inb] <-  tmpBSkl$Term1
        Term2[irepet,3,inb] <-  tmpBSkl$Term2
        RefTerm1[irepet,3,inb] <-  tmpBSkl$RefTerm1
        RefTerm2[irepet,3,inb] <-  tmpBSkl$RefTerm2
        ExL1[irepet,3,inb] <-  tmpBSkl$ExL1
        ExL2[irepet,3,inb] <-  tmpBSkl$ExL2
        RefL1[irepet,3,inb] <-  tmpBSkl$RefL1
        RefL2[irepet,3,inb] <-  tmpBSkl$RefL2
        ExL1P[irepet,inb] <- tmpBSkl$ExL1P
        ExL1M[irepet,inb] <- tmpBSkl$ExL1M
        ExL2P[irepet,inb] <- tmpBSkl$ExL2P
        ExL2M[irepet,inb] <- tmpBSkl$ExL2M
        bestSigma2[irepet,3,inb] <- sigma2
      }
    }
    for(mtd in 1:3){
      theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,3],variance2=bestSigma2[irepet,mtd,inb], NMC)
      Lntrain[irepet,mtd,inb] <- mean(loss(Ytrain,predictor(Xtrain,theta_samplesTS)))
      Lntest[irepet,mtd,inb] <- mean(loss(Ytest,predictor(Xtest,theta_samplesTS)))
    }
  }
  #setTxtProgressBar(pb, inb)
  if(!grepl("synthetic",data_option, fixed=TRUE)){
    bound <- bound[,,1]
    Term1 <- Term1[,,1]
    Term2 <- Term2[,,1]
    RefTerm1 <- RefTerm1[,,1]
    RefTerm2 <- RefTerm2[,,1]
    ExL1 <- ExL1[,,1]
    ExL2 <- ExL2[,,1]
    RefL1 <- RefL1[,,1]
    RefL2 <- RefL2[,,1]
    ExL1P <- ExL1P[,1]
    ExL1M <- ExL1M[,1]
    ExL2P <- ExL2P[,1]
    ExL2M <- ExL2M[,1]
    bestSigma2 <- bestSigma2[,,1]
    LnERMtrain <- LnERMtrain[,1]
    LnERMtest <- LnERMtest[,1]
    Lntrain <- Lntrain[,,1]
    Lntest <- Lntest[,,1]
    break
  }
}
str <- paste(c(str,d),collapse="-")

if(!grepl("synthetic",data_option, fixed=TRUE)){
  Table <- data.frame(LnERMtrain,LnERMtest,Lntrain,Lntest,bestSigma2,bound,
                      ExL1,Term1,ExL2,Term2,RefL1,RefTerm1,RefL2,RefTerm2,
                      ExL1P,ExL1M,ExL2P,ExL2M)
  colnames(Table) <- c("LnERMtrain","LnERMtest",
                       "Lntrain_PBkl","Lntrain_MGG","Lntrain_Skl",
                       "Lntest_PBkl","Lntest_MGG","Lntest_Skl",
                       "bestSigma2_PBkl","bestSigma2_MGG","bestSigma2_Skl",
                       "bound_PBkl","bound_MGG","bound_Skl",
                       "ExL1_PBkl","ExL1_MGG","ExL1_Skl",
                       "Term1_PBkl","Term1_MGG","Term1_Skl",
                       "ExL2_PBkl","ExL2_MGG","ExL2_Skl",
                       "Term2_PBkl","Term2_MGG","Term2_Skl",
                       "RefL1_PBkl","RefL1_MGG","RefL1_Skl",
                       "RefTerm1_PBkl","RefTerm2_MGG","RefTerm2_Skl",
                       "RefL2_PBkl","RefL2_MGG","RefL2_Skl",
                       "RefTerm2_PBkl","RefTerm2_MGG","RefTerm2_Skl",
                       "ExL1P","ExL1M","ExL2P","ExL2M")
  print(Table)
  
  # write
  if(!dir.exists("out")){dir.create(file.path(path,"out"), showWarnings=F)}
  outpath <- paste(c(path,"\\out\\",data_option,".csv"),collapse="")
  write.csv(Table,outpath, row.names=FALSE)
}

if(!grepl("synthetic",data_option, fixed=TRUE)){
  meansbound <- apply(bound, 2, mean)
  varsbound <- apply(bound, 2, var)
  meanstest <- apply(Lntest, 2, mean)
  varstest <- apply(Lntest, 2, var)
  meanssigma <- apply(bestSigma2, 2, mean)
  print(paste(c(data_option, ". ERM test error=", round(mean(LnERMtest),3), " (", round(var(LnERMtest),3), " )"),collapse=""))
  print(paste(c("Maurer bound=", round(meansbound[1],3), " (", round(varsbound[1],3), ") ",
                ", MGG Bound=",  round(meansbound[2],3), " (", round(varsbound[2],3), ") ",
             ", SplitKL bound=", round(meansbound[3],3), " (", round(varsbound[3],3), ") "
             ),collapse = ""))
  print(paste(c("Maurer test=", round(meanstest[1],3), " (", round(varstest[1],3), ") ",
                ", MGG test=",  round(meanstest[2],3), " (", round(varstest[2],3), ") ",
                ", SplitKL test=", round(meanstest[3],3), " (", round(varstest[3],3), ") "
                ),collapse = ""))
  print(paste(c("Maurer sigma=", round(meanssigma[1],3),
                ", MGG sigma=",  round(meanssigma[2],3),
                ", SplitKL sigma=", round(meanssigma[3],3)
                ),collapse = ""))
}else{
  MeanBKL <- apply(X = bound[,1,], MARGIN = 2, FUN = mean)
  MeanBProb <- apply(X = bound[,2,], MARGIN = 2, FUN = mean)
  MeanBSkl <- apply(X = bound[,3,], MARGIN = 2, FUN = mean)
  
  ## Saving bound 
  writeMat(paste(c(path, "/save/results_",str,".mat"), collapse = ""), labpcexport = bound)
  writeMat(paste(c(path, "/save/LnERMtrain_",str,".mat"), collapse = ""), labpcexport = LnERMtrain)
}




  
