NAfix2 <- function(x, subst=-Inf) {
  ### Written by Christian Hoffmann; propagate last known non-NA value
  ### Input:
  ###     x: numeric vector
  ###     subst: scalar inidicating which value should replace NA
  ###         if x starts with a series of NA's
  ### Output:
  ###     (numeric) vector, with NA's replaced by last known non-NA value,
  ###         or 'subst'
  spec <- max(x[!is.na(x)])+1
  x <- c(spec,x)
  while (any(is.na(x))) x[is.na(x)] <- x[(1:length(x))[is.na(x)]-1]
  x[x==spec] <- subst
  x <- x[-1]
  x
}