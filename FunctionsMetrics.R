# Detect change points locations for a given vector "d"
change_point <- function(d) {
  p <- cumsum(rle(d)$lengths) + 1
  p[-length(p)] }

# Split vector "x" with respect to the given change points locations "pos"
split_vec <- function(x, pos) {
  unname(split(x, cumsum(seq_along(x) %in% pos))) }

# Calculate Jaccard index for two sets "x" and "y"
jaccard <- function(x,y) { length(intersect(x, y)) / length(union(x, y)) }

# Calculate covering metric of true sequence by predicted sequence
covering_metric <- function(true_part, pred_part) {
  assertthat::are_equal(length(true_part), length(pred_part))
  t <- length(true_part)
  
  ### true and pred change points
  true_cp <- change_point(true_part)
  pred_cp <- change_point(pred_part)
  ### split true and pred vectors
  split_true <- split_vec(1:t, true_cp)
  split_pred <- split_vec(1:t, pred_cp)
  
  sumjac <- 0
  for(j in split_true)
  {
    sumjac <- sumjac + length(j) * max(sapply(split_pred, jaccard, y = j))
  }
  return(1/t * sumjac)
}

# Calculate EMR
emr <- function(y_true, y_pred) {
  return( mean(y_true == y_pred) )
}

# Calculate Performance Measure (PM)
perf_measure <- function(y_true, y_pred) {
  return( sqrt(emr(y_true, y_pred) * covering_metric(y_true, y_pred)) )
}

# Calculate Performance Index (PI) for a given K; optionally (bstrap = TRUE), calculate 
# 95% bootstrap confidence interval (based on 1000 samples) around PI
perf_index <- function(pdata, model, dimred, type, K, bstrap, plotDensity) {
  testLen <- length(pdata[[paste0("y", type)]])
  perfMeasure <- as.numeric(lapply(1:testLen, 
                                   function(x) perf_measure(pdata[[paste0("y", type)]][[x]],
                                                            pdata[[paste0(model, "_y", type)]][[dimred]][[x]])))
  C <- as.numeric(names(pdata[[paste0("y", type)]]))
  
  if(bstrap) {
  bsamples <- c()
  for(i in 1:1000)
  {
    idx <- sample(1:length(C), length(C), replace = TRUE)
    bsamples <- c(bsamples, sum(C[idx]^K * perfMeasure[idx]) / sum(C[idx]^K))
  }
  
  if(plotDensity) { plot(density(bsamples), 
                         main = "Distribution of means of bootstrap samples of weighted PM") 
    abline(v = c(quantile(bsamples, 0.025),
                 quantile(bsamples, 0.975)), col = "red", lty = 2)} 
  
  return( list("pi" = as.numeric(C^K %*% perfMeasure / sum(C^K)),
               "lower95" = quantile(bsamples, 0.025),
               "upper95" = quantile(bsamples, 0.925)) )
  }
  else { return(as.numeric(C^K %*% perfMeasure / sum(C^K))) }

}

# Wrapper function for calculating EMR or CM
compute_measure <- function(ytrue, ypred, m) {
  if(m == "EMR") { return(emr(ytrue, ypred)) }
  if(m == "CM") { return(covering_metric(ytrue, ypred)) }
}

