# Reduce dimensionality by applying CCA: resulting data is one dimensional
data_cca <- function(xtrain, ytrain, xtest) {
  xtrain <- as.matrix(xtrain) ; xtest <- as.matrix(xtest)
  ytrain <- as.matrix(ytrain)

  coefs <- data.frame(cc(xtrain, ytrain)$xcoef)
  coefs <- cbind(variable = rownames(coefs), coefs)
  rownames(coefs) <- 1:nrow(coefs)
  colnames(coefs) <- c("Variable", "Coefficient")
  cca_train <- data.frame(rowSums(t(apply(xtrain, 1, function(x) coefs[,2]*x)), na.rm = T))
  cca_test <- data.frame(rowSums(t(apply(xtest, 1, function(x) coefs[,2]*x)), na.rm = T))
  
  colnames(cca_train) <- "value"
  colnames(cca_test) <- "value"
  
  return(list(xtrain = cca_train, xtest = cca_test))
}

# Reduce dimensionality by applying mRMR: the number of features in the resulting
# data is controlled by "size" argument
data_mrmr <- function(xtrain, ytrain, xtest, size) {
  dd <- mRMR.data(data = as.data.frame(cbind("class" = as.numeric(ytrain$class), xtrain)))
  vars_ <- mRMR.classic(data = dd, 
                        target_indices = 1,
                        feature_count = size,
                        continuous_estimator = "kendall")
  vars <- apply(solutions(vars_)[[1]], 2, function(x, y) { return(y[x]) }, y = featureNames(dd))
  return(list(xtrain = xtrain[, vars], 
              xtest = xtest[, vars]))
}

# Reduce dimensionality by applying PCA: the number of features (principal components)
# in the resulting data is controlled by "varExpl" and "nFeatures" arguments:
# nFeatures (int) = number of resulting PCs
# varExpl (float) = min. percentage of explained variance (i.e. find k PCs such
# that they explain at least varExpl% variance)
data_pca <- function(Xtrain, Xtest, varExpl, nFeatures) {
  pc <- prcomp(Xtrain, scale = TRUE)
  cve <- cumsum((pc$sdev)^2 / sum((pc$sdev)^2)) * 100
  
  if(is.null(nFeatures)) { pcaSize <- length(cve[cve < varExpl]) + 1 }
  if(is.null(varExpl)) { pcaSize <- nFeatures }
  if(is.null(nFeatures) & is.null(varExpl)) { print("!!!") }

  PCAtrain <- pc$x[, 1:pcaSize]
  PCAtest <- predict(pc, newdata = scale(Xtest))[, 1:pcaSize]
  
  #print(paste0("Number of PCs: ", pcaSize, 
  #             ", variance explained: ", round(cve[pcaSize], 2), "%"))
  
  return(list(xtrain = PCAtrain, xtest = PCAtest))
}

# Choose functionals of the aggregated data; mean and sd are by default
choose_functionals <- function(data, funs = c("avg", "sd")) {
  cols <- grep(paste(funs, collapse = "|"), colnames(data), value = TRUE)
  return(data[, c("class", cols)])
}

# Split data into train and test sets. Size of the train set (and the corresponding 
# test set) is controlled by "nTrain" and "perc" arguments.
# nTrain (int or float) = a number / fraction of observations in the test set
# perc (bool) = if FALSE, then nTrain is a number of observations in the train set
#               if TRUE, then nTrain is a fraction of observations in the train set
data_split <- function(data, nTrain, perc = FALSE) {
  
  n <- dim(data)[1]
  if(perc) { nTrain <- floor( nTrain * n ) }
  nTest <- n - nTrain
  
  dataTrain <- data[1:nTrain, ]
  dataTest <- data[(nTrain + 1) : n, ]
  
  stopifnot(dim(dataTest)[1] == nTest)
  
  xtrain <- dataTrain[, !colnames(dataTrain) %in% "class"] |> `rownames<-`(NULL)
  xtest <- dataTest[, !colnames(dataTest) %in% "class"] |> `rownames<-`(NULL)
  ytrain <- dataTrain[, "class", drop = FALSE] |> `rownames<-`(NULL)
  ytest <- dataTest[, "class", drop = FALSE] |> `rownames<-`(NULL)

  s_len <- get_states_length(ytrain)$length < 15
  if(sum(s_len > 0)) { 
    print(paste("WARNING! Too few observations from state(s)", which(s_len))) }
  
  return(list("xtrain" = xtrain, "xtest" = xtest,
              "ytrain" = ytrain, "ytest" = ytest))
}

# Investigate what the temporal structure of states looks like 
get_states_length <- function(patient_data) {
  v <- as.numeric(as.character((patient_data$class)))
  states_length <- data.frame("state" = v[v != dplyr::lag(v, default = !v[1])], 
                              "length" = rle(v)$lengths)
  return(states_length)
}

# Subset of the original data for a given id of a patient
patient_data <- function(data, id) {
  pdata <- data[data$patient_id == id, !colnames(data) %in% "patient_id"]
  rownames(pdata) <- NULL
  return(pdata)
}

# Make a transtition matrix a stochastic matrix
scale_transmat <- function(transmat, states) {
  states <- sort(states)
  newtransmat <- transmat[states + 1, states + 1]
  transmat <- newtransmat/rowSums(newtransmat)
  return(transmat)
}


