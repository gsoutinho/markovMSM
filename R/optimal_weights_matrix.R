optimal_weights_matrix <- function(data, id, grid, transition, min_time = 0, 
                                   other_weights = NULL)
{
  # Convert data to etm data
  trans <- attr(data, "trans")
  etmdata <- msdata2etm(data, id)
  trans2 <- to.trans2(trans)
  from <- trans2$from[trans2$transno == transition]
  to <- trans2$to[trans2$transno == transition]
  
  numbers <- sapply(grid, function(x)
    table(factor(etmdata$from)[(etmdata$entry <= x & etmdata$exit > x)]))
  subevent <- sapply(grid, function(x)
    sum(etmdata$from == from & etmdata$to == to & etmdata$exit > x))
  tnumbers <- apply(numbers, 2, sum)
  weights <- sapply(1:dim(numbers)[1], function(x)
    sqrt(subevent * numbers[x, ] * (tnumbers - numbers[x, ]))/tnumbers)
  weights[is.nan(weights)] <- 0
  fn_list <- list()
  for (i in 1:dim(numbers)[1]) {
    # Take into account the distance between grids
    val <- weights[, i] * diff(c(min_time, grid))
    fn_list[[i]] <- list(fn = function(x)
      weighted.mean(abs(x), w = val, na.rm = TRUE))
    if (!is.null(other_weights)) {
      nother <- length(other_weights)
      fn_list[[i]][2:(nother + 1)] <- other_weights
    }
  }
  # Store the weights as an attribute
  attr(fn_list, "weights") <- weights
  fn_list
}
