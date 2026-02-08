rnormbm <- function(n) {
  # this function is written with normal generation in a for loop so that it is
  # more comparable to our implementation of the Ziggurat. 
  
  # we only need to generate ceiling(n/2) pairs of uniform variables.
  num_pairs <- ceiling(n / 2)
  
  # Generate independent uniform random variables U1 and U2
  u1 <- numeric(num_pairs)
  for (i in 1:n) {
    u1[i] <- runif(1)
  }
  u2 <- numeric(num_pairs)
  for (i in 1:n) {
    u2[i] <- runif(1)
  }
  
  z <- numeric(num_pairs*2)
  # Apply the Box-Muller Transform
  for (i in 1:num_pairs) {
    R <- sqrt(-2 * log(u1[i]))
    theta <- 2 * pi * u2[i]
    
    z1 <- R * cos(theta)
    z2 <- R * sin(theta)
    
    # Combine z1 and z2 into a single vector
    z[i*2-1] <- z1
    z[i*2] <- z2
  }

  # Return exactly n elements (in case n was odd)
  return(z[1:n])
}