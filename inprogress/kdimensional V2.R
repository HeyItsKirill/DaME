kdimensional = function(data, k){
  library(Rcpp)
  sourceCpp('helpFunction.cpp')
  sourceCpp("helpFunction_1.cpp")
  uniqueTimesList <- lapply(1:k, function(i) {
    ti <- data[[i]]
    deltai <- data[[ncol(data) - k + i]]
    c(0, uniqueTime(ti, deltai)[[1]])
  })
  args <- setNames(uniqueTimesList, paste0("t", 1:k))
  df <- expand.grid(args)
  df <- cbind(df, DabrowskaEstimator(df, data, k))
  return(df)
}