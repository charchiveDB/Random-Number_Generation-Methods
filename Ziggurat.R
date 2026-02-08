zigtable <- function(f, f_inv, n_layers) {
  # function to return the breakpoints (x1,x2,...,xn) for a pdf and its inverse function
  n <- n_layers - 1
  
  tail_integral <- function(r) {
    result = integrate(f, lower = r, upper = Inf)
    return(result$value)
  }

  calc_ziggurat_error <- function(r, n) {
    f0 <- f(0)
    # Calculate common area v using the base strip
    v <- r * f(r) + tail_integral(r)
    # Iterate backwards to find x values
    x_current <- r
    for (i in n:2) {
      f_next <- (v / x_current) + f(x_current)
      if (f_next > f0) {
        # Return negative penalty
        return( -1 * ((f_next - f0) + 0.5 * i) )
      }
      x_current <- f_inv(f_next)
    }
    area_top <- x_current * (f0 - f(x_current))
    # Return the discrepancy
    return(area_top - v)
  }
  
  solution <- uniroot(calc_ziggurat_error, interval = c(0.1,10), tol = 1e-12, n =n)
  r_optimal <- solution$root

  v <- r_optimal * f(r_optimal) + tail_integral(r_optimal)
  x <- numeric(n + 1)
  x[n + 1] <- r_optimal # x_255 in paper notation, index shifted for R
  
  for (i in n:2) {
    f_next <- (v / x[i + 1]) + f(x[i + 1])
    x[i] <- f_inv(f_next)
  }
  x[1] <- 0 # x_0
  return(list(x = x, v = v))
}

rnormzig <- function(n_samples, n_layers=256) {
  # set up vector to store samples
  results <- numeric(n_samples)
  
  f <- function(x) {1/sqrt(2*pi) * exp(-x**2/2)}
  f_inv <- function(y) {sqrt(-2*log(y*sqrt(2*pi)))}
  
  # Pre-calculate values that we use multiple times
  pyramid <- zigtable(f, f_inv, n_layers)
  # index 1 is 0, index 256 is the end of layer n-1 and the split between base
  # rectangle and tail
  x_vals <- pyramid$x
  
  fx_vals <- f(x_vals)
  
  v <- pyramid$v # single rectangle area
  r <- x_vals[n_layers]
  base_prop <- r * fx_vals[n_layers] / v # proportion of base rectangle in base region
  
  for (k in 1:n_samples) {
    repeat {
      # Pick a random rectangle index 'i' (usually 2 thru 257)
      # sample is really slow for this
      # in C, this can be extracted from the last 9 bits of a random value
      i <- sample(2:(n_layers + 1), 1)
      
      # If we are in the base layer
      if (i == (n_layers + 1)) {
        if (runif(1) < base_prop) { # generate proportionally from base rectangle and the tail
          x_candidate <- runif(1, 0, x_vals[n_layers]) # if it in the base rectangle, randomly generate a point from it
          if (runif(1) < 0.5) x_candidate <- -x_candidate # make distribution symmetric, other distributions may not have this line
          results[k] <- x_candidate
          break
        }
        # otherwise, sample from tail distribution (this part is specific to the normal)
        repeat {
          x_tail <- -log(runif(1)) / r
          y_tail <- -log(runif(1))
          if (2 * y_tail > x_tail**2) {
            val <- r + x_tail
            if (runif(1) < 0.5) val <- -val
            results[k] <- val
            break
          }
        }
        break # not sure if this is necessary
      }
      # 2. Generate random X coordinate within the current rectangle width
      # x_vals[i] is the width of the current box
      x_candidate <- runif(1) * x_vals[i]
      
      # Fast Path, based on unoptimized 1984 method
      # If x is within the width of the rectangle ABOVE (x_vals[i-1]),
      # it is guaranteed to be under the curve.
      if (x_candidate < x_vals[i-1]) {
        # Apply random sign and move to next iteration
        if (runif(1) < 0.5) x_candidate <- -x_candidate
        results[k] <- x_candidate
        break
      }
      
      # Slow path, adapted from 2000 method
      # should execute less than 5% of the time
      
      # Sampled x is in the "fringe" between x_vals[i-1] and x_vals[i]
      
      # Generate a random Y height relative to the strip
      y_height <- runif(1) * (fx_vals[i-1] - fx_vals[i])
      
      # 5. Exact Density Check
      # Does the random point (x, y) fall below the actual curve?
      if ((fx_vals[i] + y_height) < f(x_candidate)) {
        if (runif(1) < 0.5) x_candidate <- -x_candidate
        results[k] <- x_candidate
        break
      }
    }
  }
  return(results)
}