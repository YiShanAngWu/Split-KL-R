gendata <- function(option){
  # Classification datastes
  if(strcmpi(option,"sigmoid-synthetic")){
    # Train data
    ndata <- round(5/4 *nb.grid[inb]) # 1/4 of the data will be used as the test set 
    X <- rnorm(n = d * ndata , mean = rep(0,d), sd = rep(1,d))
    dim(X) <- c(ndata,d)
    Y <- round(as.numeric(1/(1+exp(-X%*%theta_star))))
    if(noise){
      set <- sample(x = 1:ndata,size = round(proportion*ndata),replace = F)
      Y[set] <- (1-Y[set])
    }
    return(list("X"=X,"Y"=Y,"d"=d))
  }
  else if(strcmpi(data_option,"breast-cancer")){
    d <- 9 
    df <- read.table(paste(path,"Datasets/breast-cancer-wisconsin.data",sep="/"), header = FALSE, sep=",")
    XY <- data.matrix(df)
    # Dicard the first ID column
    XY <- XY[,-1]
    print(c("Effective dimension=",d))
  }
  else if(strcmpi(data_option,"bank-notes")){
    d <- 4
    df <- read.table(paste(path,"Datasets/data_banknote_authentication.txt",sep="/"), header = FALSE, sep=",")
    XY <- data.matrix(df)
    print(c("Effective dimension=",d))
  }
  else if(strcmpi(data_option,"haberman")){
    d <- 3
    df <- read.table(paste(path,"Datasets/haberman.data",sep="/"), header = FALSE, sep=",")
    XY <- data.matrix(df)
    print(c("Effective dimension=",d))
  }
  else if(strcmpi(data_option,"kr-vs-kp")){
    d <- 36  
    df <- read.table(paste(path,"Datasets/kr-vs-kp.data",sep="/"), header = FALSE, sep=",")
    dmy <-  dummyVars(" ~ .", data = df[,-(d+1)])
    encoded <- data.frame(predict(dmy, newdata = df[,-(d+1)]))
    XY <- cbind(encoded,df[,d+1])
    XY <- as.matrix(XY)
    d <- dim(XY)[2] - 1 # Changing the effective dimension
    print(c("Effective dimension=",d))
    XY[,d+1] <- str_replace_all(XY[,d+1],c("won" = "1", "nowin" = "0"))
    XY<-apply(XY, 2, as.numeric)
  }
  else if(strcmpi(data_option,"mushroom")){
    d <- 22  
    df <- read.table(paste(path,"Datasets/agaricus-lepiota.data",sep="/"), header = FALSE, sep=",")
    dmy <-  dummyVars(" ~ .", data = df[,-c(1,17)]) 
    encoded <- data.frame(predict(dmy, newdata = df[,-c(1,17)])) # attribute 17 is the "veil type". We remove it because it is the same for all 
    XY <- cbind(encoded,df[,1])
    XY <- as.matrix(XY)
    d <- dim(XY)[2] - 1 # Changing the effective dimension
    print(c("Effective dimension=",d))
    XY[,d+1] <- str_replace_all(XY[,d+1],c("e" = "1", "p" = "0"))
    XY <- apply(XY, 2, as.numeric)
  }
  else if(strcmpi(data_option,"adults")){
    d <- 14
    cont <- c(1,3,5,11,12,13,d+1) # Columns with continuous data plus the label column
    df <- read.table(paste(path,"Datasets/adult-nospaces.data",sep="/"), header = FALSE, sep=",")
    ind <- which(df == "?", arr.ind=TRUE)
    df <- df[-ind[,1],]
    ind <- which(df == "?", arr.ind=TRUE)
    dmy <-  dummyVars(" ~ .", data = df[,-cont]) 
    encoded <- data.frame(predict(dmy, newdata = df[,-cont])) 
    XY <- cbind(encoded,df[,cont])
    XY <- as.matrix(XY)
    d <- dim(XY)[2] - 1 # Changing the effective dimension
    print(c("Effective dimension=",d))
    XY[,d+1] <- str_replace_all(XY[,d+1],c("<=50K" = "0", ">50K" = "1"))
    XY <- apply(XY, 2, as.numeric)
  }
  else if(strcmpi(data_option, "spam")){
    d <- 57
    df <- read.table(paste(path,"Datasets/spambase.data",sep="/"), header = FALSE, sep=",")
    XY <- as.matrix(df)
    print(c("Effective dimension=",d))
  }
  else if(strcmpi(data_option,"tictactoe")){
    d <- 9 
    df <- read.table(paste(path,"Datasets/tic-tac-toe.data",sep="/"), header = FALSE, sep=",")
    dmy <-  dummyVars(" ~ .", data = df[,-(d+1)]) 
    encoded <- data.frame(predict(dmy, newdata = df[,-(d+1)]))
    XY <- cbind(encoded,df[,d+1])
    XY <- as.matrix(XY)
    d <- dim(XY)[2] - 1 # Changing the effective dimension
    print(c("Effective dimension=",d))
    XY[,d+1] <- str_replace_all(XY[,d+1],c("positive" = "1", "negative" = "0"))
    XY <- apply(XY, 2, as.numeric)
  }
  else if(strcmpi(data_option,"svmguide1")){
    d <- 4
    df <- read.table(paste(path,"Datasets/svmguide1.data",sep="/"), header = FALSE, sep="")
    df <- df %>% relocate((1), .after = (d+1))
    for(i in 1:d){
      df <- separate(df, (i), c(NA, toString(i)), sep=":")
    }
    names(df) <- NULL
    XY <- as.matrix(df)
    XY <- apply(XY, 2, as.numeric)
    print(c("Effective dimension=",d))
  }
  else if(strcmpi(data_option,"splice")){
    d <- 60
    df <- read.table(paste(path,"Datasets/splice.data",sep="/"), header = FALSE, sep="")
    df <- df %>% relocate((1), .after = (d+1))
    for(i in 1:d){
      df <- separate(df, (i), c(NA, toString(i)), sep=":")
    }
    names(df) <- NULL
    XY <- as.matrix(df)
    XY[,d+1] <- str_replace_all(XY[,d+1],c("+1" = "1", "-1" = "0"))
    XY <- apply(XY, 2, as.numeric)
    print(c("Effective dimension=",d))
  }
  else{
    stop("Error: data option unsupported")
  }
  
  XY <- XY[sample(nrow(XY)),]  # Shuffle data
  X <- XY[,1:d]
  Y <- XY[,d+1]
  ndata <- length(Y)
  
  ## Common pre-processing depending on the problem type
  # Make the problem 0-1
  Ymin <- min(Y)
  Ymax <- max(Y)
  Y <- (Y - Ymin)/(Ymax-Ymin)
  Xmin <- numeric(d)
  Xmax <- numeric(d)
  
  # Feature scaling
  for(j in 1:d){
    Xmin[j] <- min(X[,j])
    Xmax[j] <- max(X[,j])
    diff <- (Xmax[j]-Xmin[j])
    if(diff==0) 
      diff<-1
    X[,j] <- 2*(X[,j] - (Xmax[j]+Xmin[j])/2)/diff
  } 
  
  return(list("X"=X,"Y"=Y,"d"=d))
 }
  