chi_squared_test_norm <- function(x, k){
  t <- length(x) # number of observations
  p <- 1 / 1000000000 # guarantee a very small probability of multiple samples ending up in each bin
  bins <- qnorm(seq(p, 1-p, length.out = k+1), mean = 0, sd = 1) #200 bins with equal prob
  
  Yi <- hist(x, breaks = bins, plot = FALSE)$counts # observed counts
  # expected probabilites (p_i) of standard normal
  p_i <- numeric(k)
  for (i in 1:k) {
    p_i[i] <- pnorm(bins[i+1]) - pnorm(bins[i])
  }
  # using computing chi square using equation from the paper
  chi_sq <- sum((Yi - t * p_i)^2 / (t * p_i))
  df <- k - 1
  return(data.frame(chi_square = chi_sq,
                    df = df))
}

chi_squared_test_exp <- function(x, k){
  t <- length(x) # number of observations
  p <- 1 / 1000000000 # guarantee a very small probability of multiple samples ending up in each bin
  bins <- qexp(seq(p, 1-p, length.out = k+1)) #200 bins with equal prob
  
  Yi <- hist(x, breaks = bins, plot = FALSE)$counts # observed counts
  # expected probabilites (p_i) of standard normal
  p_i <- numeric(k)
  for (i in 1:k) {
    p_i[i] <- pexp(bins[i+1]) - pexp(bins[i])
  }
  # using computing chi square using equation from the paper
  chi_sq <- sum((Yi - t * p_i)^2 / (t * p_i))
  df <- k - 1
  return(data.frame(chi_square = chi_sq,
                    df = df))
}