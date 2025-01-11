summary_count <- function(InputMatrix,Col_vec,Col_filter,AddToSum = 0){
  
  if(!(Col_filter %in% names(InputMatrix))){
    stop("Column specified for filter is not present in input matrix")
  }
  
  filter_vec <- InputMatrix[,Col_filter]
  
  if(any(is.na(filter_vec))){
    
    
    warning(paste(sum(is.na(filter_vec)),"rows containing NAs in the filter column has been removed"))
    InputMatrix <- InputMatrix[!is.na(filter_vec),]
    filter_vec <- filter_vec[!is.na(filter_vec)]
  }
  
  Row_vec <- unique(filter_vec)
  
  n <- length(Row_vec)
  m <- length(Col_vec)
  
  df.count <- data.frame(matrix(ncol = m, nrow = n))
  names(df.count) <- Col_vec
  row.names(df.count) <- Row_vec
  
  for(i in 1:n){
    
    filter1 <- filter_vec == Row_vec[i]
    df.temp <- InputMatrix[filter1,]
    
    for(j in 1:m){
      
      df.count[i,j] <- sum(df.temp[,Col_vec[j]])+AddToSum
      
    }
  }
  
  return(df.count)
  
}
