# Conduct a whole analysis for a given patient id
patient_analysis <- function(data, patient_id, range_train, vec_train, transmat, nMRMR, percVE, tp) {
  
  dfPatient <- patient_data(data, patient_id)
  nObs <- dim(dfPatient)[1]
  y_all <- dfPatient$class
  
  if(is.null(vec_train)) {
    min_train = range_train[1]
    max_train = range_train[2]
    train_s <- floor(nObs*min_train) : floor(nObs*max_train)
  }
  
  if(is.null(range_train)) {
    train_s <- floor(nObs*vec_train)
  }
  
  if(!is.null(range_train) & !is.null(vec_train)) {
    print("Choose one option!")
  }
  
  print(paste0("Patient ", patient_id, ", number of obs: ", nObs))
  print(paste0("Min train size: ", min_train*100, "%, number of obs: ", train_s[1]))
  print(paste0("Max train size: ", max_train*100, "%, number of obs: ", train_s[length(train_s)]))
  
  tns_perc <- round(train_s / nObs * 100, 3)
  tts_perc <- round((nObs - train_s) / nObs * 100, 3)
  lstPredsTrain <- vector("list", length(tns_perc)) |> stats::setNames(tns_perc)
  lstPredsTest <- vector("list", length(tts_perc)) |> stats::setNames(tts_perc)
  
  HMM_ytr <- replicate(3, lstPredsTrain, simplify = F) |> stats::setNames(c("CCA", "PCA", "mRMR"))
  HMMe_ytr <- replicate(3, lstPredsTrain, simplify = F) |> stats::setNames(c("CCA", "PCA", "mRMR"))
  GNB_ytr <- replicate(3, lstPredsTrain, simplify = F) |> stats::setNames(c("CCA", "PCA", "mRMR"))
  RF_ytr <- replicate(3, lstPredsTrain, simplify = F) |> stats::setNames(c("CCA", "PCA", "mRMR"))
  HMM_yte <- replicate(3, lstPredsTest, simplify = F) |> stats::setNames(c("CCA", "PCA", "mRMR"))
  HMMe_yte <- replicate(3, lstPredsTest, simplify = F) |> stats::setNames(c("CCA", "PCA", "mRMR"))
  GNB_yte <- replicate(3, lstPredsTest, simplify = F) |> stats::setNames(c("CCA", "PCA", "mRMR"))
  RF_yte <- replicate(3, lstPredsTest, simplify = F) |> stats::setNames(c("CCA", "PCA", "mRMR"))
  
  train_true_y <- vector("list", length(tns_perc)) |> stats::setNames(tns_perc)
  test_true_y <- vector("list", length(tts_perc)) |> stats::setNames(tts_perc)
  
  for(j in train_s)
  {
    i <- which(j == train_s)
    print(paste0("Train set: ", i, " / ", length(train_s)))
    dfTemp <- data_split(dfPatient, j)
    
    Xtrain <- dfTemp$xtrain ; Xtest <- dfTemp$xtest
    ytrain <- dfTemp$ytrain ; ytest <- dfTemp$ytest
    
    xtr_mrmr <- data_mrmr(Xtrain, ytrain, Xtest, nMRMR)$xtrain
    xte_mrmr <- data_mrmr(Xtrain, ytrain, Xtest, nMRMR)$xtest
    xtr_cca <- data_cca(Xtrain, ytrain, Xtest)$xtrain
    xte_cca <- data_cca(Xtrain, ytrain, Xtest)$xtest
    xtr_pca <- data_pca(Xtrain, Xtest, percVE, NULL)$xtrain
    xte_pca <- data_pca(Xtrain, Xtest, percVE, NULL)$xtest
    
    states <- sort(unique(ytrain$class))

    # Models and predictions
    hmm_cca <- hmm_mod(xtr_cca, xte_cca, ytrain, NULL)
    hmm_pca <- hmm_mod(xtr_pca, xte_pca, ytrain, NULL)
    hmm_mrmr <- hmm_mod(xtr_mrmr, xte_mrmr, ytrain, NULL)
    
    hmme_cca <- hmm_mod(xtr_cca, xte_cca, ytrain, transmat)
    hmme_pca <- hmm_mod(xtr_pca, xte_pca, ytrain, transmat)
    hmme_mrmr <- hmm_mod(xtr_mrmr, xte_mrmr, ytrain, transmat)
    
    gnb_cca <- gnb_mod(xtr_cca, xte_cca, ytrain)
    gnb_pca <- gnb_mod(xtr_pca, xte_pca, ytrain)
    gnb_mrmr <- gnb_mod(xtr_mrmr, xte_mrmr, ytrain)
    
    rf_cca <- rf_mod(xtr_cca, xte_cca, ytrain)
    rf_pca <- rf_mod(xtr_pca, xte_pca, ytrain)
    rf_mrmr <- rf_mod(xtr_mrmr, xte_mrmr, ytrain)
    
    ####
    HMM_ytr[["CCA"]][[i]] <- hmm_cca$train
    HMM_ytr[["PCA"]][[i]] <- hmm_pca$train
    HMM_ytr[["mRMR"]][[i]] <- hmm_mrmr$train
    
    HMM_yte[["CCA"]][[i]] <- hmm_cca$test
    HMM_yte[["PCA"]][[i]] <- hmm_pca$test
    HMM_yte[["mRMR"]][[i]] <- hmm_mrmr$test
    
    HMMe_ytr[["CCA"]][[i]] <- hmme_cca$train
    HMMe_ytr[["PCA"]][[i]] <- hmme_pca$train
    HMMe_ytr[["mRMR"]][[i]] <- hmme_mrmr$train
    
    HMMe_yte[["CCA"]][[i]] <- hmme_cca$test
    HMMe_yte[["PCA"]][[i]] <- hmme_pca$test
    HMMe_yte[["mRMR"]][[i]] <- hmme_mrmr$test
    
    GNB_ytr[["CCA"]][[i]] <- gnb_cca$train
    GNB_ytr[["PCA"]][[i]] <- gnb_pca$train
    GNB_ytr[["mRMR"]][[i]] <- gnb_mrmr$train
    
    GNB_yte[["CCA"]][[i]] <- gnb_cca$test
    GNB_yte[["PCA"]][[i]] <- gnb_pca$test
    GNB_yte[["mRMR"]][[i]] <- gnb_mrmr$test
    
    RF_ytr[["CCA"]][[i]] <- rf_cca$train
    RF_ytr[["PCA"]][[i]] <- rf_pca$train
    RF_ytr[["mRMR"]][[i]] <- rf_mrmr$train
    
    RF_yte[["CCA"]][[i]] <- rf_cca$test
    RF_yte[["PCA"]][[i]] <- rf_pca$test
    RF_yte[["mRMR"]][[i]] <- rf_mrmr$test
    
    train_true_y[[i]] <- ytrain$class
    test_true_y[[i]] <- ytest$class
    
  }
  return(list("HMM_ytrain" = HMM_ytr,
              "HMM_ytest" = HMM_yte,
              "HMMe_ytrain" = HMMe_ytr,
              "HMMe_ytest" = HMMe_yte,
              "GNB_ytrain" = GNB_ytr,
              "GNB_ytest" = GNB_yte,
              "RF_ytrain" = RF_ytr,
              "RF_ytest" = RF_yte,
              "ytrain" = train_true_y,
              "ytest" = test_true_y))
}

# Create a matrix of predictions for plot_states function
pred_data <- function(ytrue, ypred, model, dimred) {
  ypred_ <- ypred[[paste0(model, "_ytest")]][[dimred]]
  testSizes <- as.numeric(names(ypred_))
  patientSize <- length(ytrue)
  
  df_preds <- data.frame(matrix(NA, nrow = patientSize, ncol = length(testSizes) + 1))
  
  for(i in 1:length(testSizes))
  {
    df_preds[, i] <- c(rep(9, patientSize - length(unlist(ypred_[i]))), unlist(ypred_[i]))
  }
  Actual <- ytrue
  df_preds[, tail(colnames(df_preds), 1)] <- Actual
  trainSizes <- 1 - testSizes/100
  colnames(df_preds) <- c(trainSizes, "Actual")
  return(df_preds)
}

# Summarise predictions in terms of EMR and CM for different data splits, models and 
# dimensionality reduction techniques
summarise_performance <- function(pdat) {
  d <- data.frame(matrix(nrow = length(pdat$ytest), ncol = 3*4*2 + 1))
  d[,1] <- as.numeric(names(pdat$ytest))
  colnames(d)[1] <- "TestSize"
  
  mods <- c("HMM", "HMMe", "GNB", "RF")
  dimreds <- c("CCA", "PCA", "mRMR")
  vals <- c("EMR", "CM")
  i <- 2
  for(dr in dimreds)
  {
    for(mod in mods)
    {
      {
        for(vs in vals)
        {
          d[,i] <- as.numeric(lapply(1:length(pdat$ytest), 
                                     function(x) compute_measure(pdat[[paste0("y", "test")]][[x]],
                                                                pdat[[paste0(mod, "_y", "test")]][[dr]][[x]], vs)))
          colnames(d)[i] <- paste0(dr, "_", mod, "_", vs)
          i <- i + 1
        }
      }
    }
  }
  return(d)
}

# Select a subset of data
filter_data <- function(sdat, mods, dimreds, measures) {
  sdatm <- melt(sdat, id.vars = 'TestSize', variable.name = 'Classifier')
  sdatm$Variable <- sdatm$Classifier
  fdata <- separate_wider_delim(sdatm, cols = Classifier, delim = "_", names = c("Method", "Model", "Measure"))
  
  fdata <- fdata[fdata$Method %in% dimreds, ]
  fdata <- fdata[fdata$Model %in% mods, ]
  fdata <- fdata[fdata$Measure %in% measures, ]
  
  return(fdata %>% data.frame())
}

# Create a dataframe of Performance Indexes for various values of K
data_pi_of_k <- function(dat, mods, dimreds, rangeK) {
  df <- data.frame(matrix(ncol = length(mods) * length(dimreds), nrow = length(rangeK)))
  dfl <- data.frame(matrix(ncol = length(mods) * length(dimreds), nrow = length(rangeK)))
  dfu <- data.frame(matrix(ncol = length(mods) * length(dimreds), nrow = length(rangeK)))
  
  i <- 1
  for(m in mods)
  {
    for(dr in dimreds)
    {
      v <- lapply(rangeK, function(x) perf_index(dat, m, dr, "test", x, TRUE, FALSE))
      df[,i] <- as.numeric(lapply(v, function(x) as.numeric(x$pi)))
      dfl[,i] <- as.numeric(lapply(v, function(x) as.numeric(x$lower95)))
      dfu[,i] <- as.numeric(lapply(v, function(x) as.numeric(x$upper95)))
      
      colnames(df)[i] <- paste0(m, "_", dr)
      colnames(dfl)[i] <- paste0(m, "_", dr)
      colnames(dfu)[i] <- paste0(m, "_", dr)
      i <- i + 1
    }
  }
  return(list("mean" = df,
              "lower" = dfl,
              "upper" = dfu,
              "K" = rangeK))
}
