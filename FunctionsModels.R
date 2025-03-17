
# Predict the states for HMM
viterbi_decoding <- function(Xtrain, Xtest, mod) {
  ytrain_pred <- RcppHMM::viterbi(mod, t(as.matrix(Xtrain)))
  ytest_pred <- RcppHMM::viterbi(mod, t(as.matrix(Xtest)))
  
  return(list("ytrain_pred" = as.numeric(ytrain_pred),
              "ytest_pred" = as.numeric(ytest_pred)))
}

# For HMM and predict states on the test and training set
hmm_mod <- function(xtrain, xtest, ytrain, transmat = NULL) {
  states <- sort(unique(ytrain$class))
  N <- length(states)
  M <- length(colnames(xtrain))
  
  # Initial distribution
  Pi <- rep(1/N, N)
  # Transition matrix
  if(is.null(transmat))
  {
    xChar <- as.character(ytrain$class)
    mcX <- markovchainFit(xChar)
    A <- matrix(attr(mcX$estimate, "transitionMatrix"), ncol = N, byrow = FALSE)
  }
  else
  {
    A <- scale_transmat(transmat, states)
  }
  
  # Emission probability
  if(M == 1)
  {
    mu <- c()
    sigma <- c()
    
    for(i in states)
    {
      j <- which(states == i)
      mu[j] <- mean(xtrain[ytrain$class == i, ])
      sigma[j] <- sd(xtrain[ytrain$class == i, ])
    }
    Mu <- matrix(mu, ncol = N)
    Sigma <- array(sigma, dim = c(1, 1, N))
    
    
    
    model <- verifyModel(list( "Model" = "GHMM",
                               "StateNames" = states,
                               "A" = A,
                               "Mu" = Mu,
                               "Sigma" = Sigma,
                               "Pi" = Pi))
    
  }
  if(M > 1)
  {
    Mu <- matrix(NA, nrow = M, ncol = N)
    Sigma <- array(0, dim = c(M, M, N))
    for(i in states)
    {
      j <- which(states == i)
      Mu[, j] <- mlest(xtrain[ytrain$class == i, ])$muhat
      Sigma[, , j] <- mlest(xtrain[ytrain$class == i, ])$sigmahat
    }
    
    model <- verifyModel(list( "Model" = "GHMM",
                               "StateNames" = states,
                               "A" = A,
                               "Mu" = Mu,
                               "Sigma" = Sigma,
                               "Pi" = Pi))
  }
  
  ytrain_pred <- viterbi_decoding(xtrain, xtest, model)$ytrain_pred
  ytest_pred <- viterbi_decoding(xtrain, xtest, model)$ytest_pred
  
  return(list("train" = ytrain_pred,
              "test" = ytest_pred,
              "model" = model))
}

# For GNB and predict states on the test and training set
gnb_mod <- function(xtrain, xtest, ytrain) {
  Xtrain <- as.matrix(xtrain) ; Xtest <- as.matrix(xtest)
  ytrain <- as.factor(as.matrix(ytrain)) 
  model <- gaussian_naive_bayes(Xtrain, ytrain)
  
  ytrain_pred <- as.numeric(as.character(predict(model, newdata = Xtrain, type = "class")))
  ytest_pred <- as.numeric(as.character(predict(model, newdata = Xtest, type = "class")))
  
  return(list("train" = ytrain_pred,
              "test" = ytest_pred))
}

# For RF and predict states on the test and training set
rf_mod <- function(xtrain, xtest, ytrain) {
  ytrain$class <- as.factor(ytrain$class)
  train <- cbind("class" = ytrain, xtrain)
  
  model <- randomForest(class ~ . , data = train)
  
  ytrain_pred <- as.numeric(as.character(predict(model, xtrain)))
  ytest_pred <- as.numeric(as.character(predict(model, xtest)))
  
  return(list("train" = ytrain_pred,
              "test" = ytest_pred))
}

