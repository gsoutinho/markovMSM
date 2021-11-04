optimal_weights_multiple <- function(data, id, grid, transition, min_time = 0)
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
    subevent * numbers[x, ] * (tnumbers - numbers[x, ])/tnumbers^2)
  weights[is.nan(weights)] <- 0
  weight <- apply(weights, 1, max)
  weight * diff(c(min_time, grid))
}