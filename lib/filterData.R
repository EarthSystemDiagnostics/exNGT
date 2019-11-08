#Function to filter annual mean isotope record with different window size
#function for filter
#by Thom Laepple, following Mann et al.,2003
#' @title  Apply a filter to a timeseries
#' @description Apply a filter to a timeseries the timestep provided by ts is
#'   not used!!!! Thus for timeseries with a different spacing than 1, the
#'   filter has to be adapted Using endpoint constrains as describen in  Mann et
#'   al., GRL 2003 no constraint (loss at both ends) (method=0) minimum norm
#'   constraint (method=1) minimum slope constraint (method=2) minimum roughness
#'   constraint (method=3) circular filtering (method=4)
#' @param data Input timeseries (ts object)
#' @param filter vector of filter weights
#' @param method constraint method choice 0-4
#' @return filtered timeseries (ts object)
#' @author Thomas Laepple
#' @export

ApplyFilter <- function(data, filter, method = 0) {
  N <- floor(length(filter)/2)
  if (method == 0) {
    result <- stats::filter(c(data), filter, circular = FALSE)
    return(ts(result, frequency = frequency(data)))
  }
  
  if (method == 4) {
    result <- stats::filter(c(data), filter, circular = TRUE)
    return(ts(result, frequency = frequency(data)))
  }
  
  {
    if (method == 1) # Minimum Norm
    {
      before <- rep(mean(data), N)
      after <- rep(mean(data), N)
    }
    if (method == 2) {
      before <- c(data)[N:1]
      after <- c(data)[length(data):(length(data) - N + 1)]
    }
    if (method == 3) {
      before <- c(data)[N:1]
      after <- c(data)[length(data):(length(data) - N + 1)]
      
      before <- c(data)[1] - (before - mean(before))
      
      after <- c(data)[length(data)] - (after - mean(after))
    }
    result <- stats::filter(c(before, data, after), filter, circular = F)[(N + 1):(N + length(data))]
    return(ts(result, frequency = frequency(data)))
  }
}

#---------------------------------------------

# input vector or dataframe (rows are time!)
# input window length of filter
# output filtered vector or dataframe
# using Thoms Apply Filter for end member treatment
filterData <-function(data, window = 5, method = 2, hasAgeColumn = TRUE) {
  
  if (hasAgeColumn) {
    Year <- data[, 1]
    data <- data[, -1]
  }
  
  if (is.null(dim(data))) {
    
    Datafilt<- ApplyFilter(data, filter = (rep(1 / window, window)), method = method)
    
  } else {
    
    Datafilt<- apply(data, 2, ApplyFilter, filter = (rep(1 / window, window)), method = method)
    
  }
  
  if (hasAgeColumn) Datafilt <- as.data.frame(cbind(Year, Datafilt))
    
  return(Datafilt)
  
}