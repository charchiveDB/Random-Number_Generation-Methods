zigtable <- function(f, f_inv, n) {
  # function to return the breakpoints (x1,x2,...,xn) for a pdf and its inverse function
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
    for (i in (n-1):1) {
      f_next <- (v / x_current) + f(x_current)
      if (f_next > f0) {
        # Return negative penalty
        return( -1 * ((f_next - f0) + 0.5 * i) )
      }
      x_current <- f_inv(f_next)
    }
    area_top <- x_current * (1 - f(x_current))
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

rexpzig <- function(n) {
  u <- runif(n)
  x <- u # add code to transform to pdf here
  return (x)
}

rnormzig <- function(n, levels) {
  results <- numeric(n_samples)
  fx_vals <- f(x_vals)
  
  j <- rdunif(n,1,256)
  u <- runif(n)
  table <- zigtable[[1]]
}