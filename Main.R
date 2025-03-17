library(mRMRe)
library(praznik)
library(CCA)
library(naivebayes)
library(randomForest)
library(markovchain)
library(mvnmle)
library(reshape2)
library(RcppHMM)
library(dplyr)
library(tidyr)
library(ggplot2)

source("FunctionsAnalyse.R")
source("FunctionsMetrics.R")
source("FunctionsModels.R")
source("FunctionsPreprocessing.R")
source("FunctionsVisualize.R")

### Read data
patients_data <- read.csv("dfp.csv")

### Expert transition matrix based on [12]
expert_transmat <- matrix(c(25, 17, 2, 6,
                   14, 39, 4, 4,
                   4, 2, 1, 2,
                   9, 5, 5, 6), ncol = 4, nrow = 4, 
                 dimnames = list(c("Euthymia", "Depression", "Mixed", "Mania"),
                                 c("Euthymia", "Depression", "Mixed", "Mania")),
                 byrow = TRUE)
expert_transmat <- expert_transmat / rowSums(expert_transmat)

### Available patients
unique(patients_data$patient_id)

### Choose example patient
p_id <- 5736

### Select a subbset of the original data that corresponds to a given patient
pdat <- patient_data(patients_data, p_id)

### How does the structure of hidden states look like?
get_states_length(pdat) # 292 observations from state 1 (depression) followed by 
                        # 28 observations from state 2 (mania)

### Choose a subset of functionals: average and standard deviation
pdat <- choose_functionals(pdat, c("avg", "sd"))

### Split data intro train and test part for C = 0.80
dtemp <- data_split(pdat, 0.80, TRUE)

x_train <- dtemp$xtrain
x_test <- dtemp$xtest
y_train <- dtemp$ytrain
y_test <- dtemp$ytest

### Reduce dimensionality using mRMR algorithm: preserve 3 "best" features
x_mrmr <- data_mrmr(x_train, y_train, x_test, 3)
xtr_mrmr <- x_mrmr$xtrain
xte_mrmr <- x_mrmr$xtest

### Reduce dimensionality using CCA
x_cca <- data_cca(x_train, y_train, x_test)
xtr_cca <- x_cca$xtrain
xte_cca <- x_cca$xtest

### Reduce dimensionality using PCA: choose such number of principal components 
# that allows to explain at least 75% of variance in the data
x_pca <- data_pca(x_train, x_test, 75, NULL)
xtr_pca <- x_pca$xtrain
xte_pca <- x_pca$xtest

### Fit models and make predictions
gnb_mrmr <- gnb_mod(xtr_mrmr, xte_mrmr, y_train)
gnb_cca <- gnb_mod(xtr_cca, xte_cca, y_train)
gnb_pca <- gnb_mod(xtr_pca, xte_pca, y_train)
rf_mrmr <- rf_mod(xtr_mrmr, xte_mrmr, y_train)
rf_cca <- rf_mod(xtr_cca, xte_cca, y_train)
rf_pca <- rf_mod(xtr_pca, xte_pca, y_train)
hmm_mrmr <- hmm_mod(xtr_mrmr, xte_mrmr, y_train, NULL)
hmme_mrmr <- hmm_mod(xtr_mrmr, xte_mrmr, y_train, expert_transmat)
hmm_cca <- hmm_mod(xtr_cca, xte_cca, y_train, NULL)
hmme_cca <- hmm_mod(xtr_cca, xte_cca, y_train, expert_transmat)
hmm_pca <- hmm_mod(xtr_pca, xte_pca, y_train, NULL)
hmme_pca <- hmm_mod(xtr_pca, xte_pca, y_train, expert_transmat)

### Optionally, make prediction based on the fitted model
y_hmm <- viterbi_decoding(xtr_cca, xte_cca, hmm_cca$model)

### Calculate Performance Measure
perf_measure(y_test$class, y_hmm$ytest_pred)

### Conduct an overall analysis
p_analysis <- patient_analysis(patients_data, p_id, c(0.85, 0.99), NULL, NULL, 5, 75, expert_transmat)

### Visualize predictions: here for HMM and CCA
mod <- "HMM"
dimred <- "CCA"
plot_states(pred_data(pdat$class, p_analysis, mod, dimred), paste0("Patient ", p_id, " (", mod, ", ", dimred, ")"))

### Summarise predictions: calculate EMR and CM for each model and each dimred technique
p_analysis_sum <- summarise_performance(p_analysis)

### Choose a subset of models, dimensionality reduction techniques and measures
mod_subset <- c("HMM", "RF")
dimred_subset <- c("CCA")
measure_subset <- c("CM", "EMR")
f_data <- filter_data(p_analysis_sum, mod_subset, dimred_subset, measure_subset)

### Visualize the above
visualize_performance(f_data)

### Calculate PI(K) for a range of K values; add 95%CI
k_range <- seq(-5, 5, by = 0.5)
dat_pik <- data_pi_of_k(p_analysis, mod_subset, dimred_subset, k_range)

### Visualize the above
visualize_pi_of_k(dat_pik)



