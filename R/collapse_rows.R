collapse_rows <- function(X, row_groups){

	X <- as.data.frame(X)
  X_split_by_gr <- split(X, row_groups)
  max_avg <- lapply(X_split_by_gr, function(x) apply(x, 1, mean)) #rowAverage
  max_avg <- lapply(max_avg, function(x) names(x[which.max(x)])) #names of the highest average
  max_avg <- data.frame(old=unlist(max_avg), new=names(max_avg), stringsAsFactors = F) #map
  ans <- merge(max_avg, X, by.x=1, by.y=0) #merge selected names in relation to novel

  rownames(ans) <- ans$new
  ans <- ans[, -c(1:2)]

  return(ans)

}
