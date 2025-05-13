####### File to illustrate experiments of the HMC algorithm

target <- function(x){
  # Simple bi-modal target
  return(
    0.6*dnorm(x, mean = 0, sd = 1) + 0.4*dnorm(x, mean = 5, sd = 1)
  )
}

log_target <- function(x){
  # Log of the target
  return(log(target(x)))
}

grad_log_target <- function(x) {
  # Components
  w1 <- 0.6
  w2 <- 0.4
  mu1 <- 0
  mu2 <- 5
  sd1 <- 1
  sd2 <- 1
  
  # Densities
  d1 <- dnorm(x, mean = mu1, sd = sd1)
  d2 <- dnorm(x, mean = mu2, sd = sd2)
  
  # Derivatives of each component
  d1_prime <- - (x - mu1) / (sd1^2) * d1
  d2_prime <- - (x - mu2) / (sd2^2) * d2
  
  # Mixture density
  pi_x <- w1 * d1 + w2 * d2
  pi_prime_x <- w1 * d1_prime + w2 * d2_prime
  
  # Gradient of log-density
  return(pi_prime_x / pi_x)
}


## Visualising target
x <- seq(-3, 8, length.out = 200)
y <- target(x)
plot(x, y, type = "l", lwd = 2, col = "blue", 
     main = "Target Distribution", xlab = "x", ylab = "Density")

HMC_Sampler <- function(L = 10, eps = .1, n = 1e3){
  qt <- numeric(n)
  qt[1] <- 0 # Some starting value
  accept <- 1
  
  # vectorizing this outside to save time
  momentums <- rnorm(n)
  for(k in 2:n)
  {
    # new momentum
    p <- momentums[k]
    q_prop <- qt[k-1]
    
    # one half-Euler step
    p_prop <- p + eps/2 * grad_log_target(q_prop)
    for(l in 1:L) 
    {
      q_prop <- q_prop + eps * p_prop
      if(l != L) p_prop <- p_prop + eps*grad_log_target(q_prop)
    }
    # one last half-Euler step
    p_prop <- p_prop + eps/2 * grad_log_target(q_prop)
    
    # Accept-reject
    energy_prop <- -log_target(q_prop) + (p_prop**2)/2
    energy_pre <- -log_target(qt[k-1]) + (momentums[k]**2)/2
    log.ratio <- -(energy_prop  - energy_pre)
    
    if(log(runif(1)) < log.ratio)
    {
      qt[k] <- q_prop
      accept <- accept + 1
    } else{
      qt[k] <- qt[k-1]
    }
  }
  print(paste("Acceptance = ", accept/n))
  return(qt)
}

chain <- HMC_Sampler(L = 20, eps = .1, n = 1e5)

par(mfrow = c(1, 1))
hist(chain, breaks = 50, probability = TRUE,
     main = "Samples for L = 20, eps = 0.1", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2)
lines(density(chain), col = "red", lwd = 3, lt = 2) # Overlay density estimate
x <- seq(-3, 8, length.out = 200)
y <- target(x)
lines(x, y, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("HMC Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, lty = c(2,1), bty = "n",
       cex = 1.5)

acf(chain, lag.max = 50, main = "ACF of HMC Chain", cex.main = 2,
    col = "lightblue3", lwd = 2)

###############################################################################

## Changing the means to 0 and 8
target <- function(x){
  # Simple bi-modal target
  return(
    0.6*dnorm(x, mean = 0, sd = 1) + 0.4*dnorm(x, mean = 8, sd = 1)
  )
}

grad_log_target <- function(x) {
  # Components
  w1 <- 0.6
  w2 <- 0.4
  mu1 <- 0
  mu2 <- 8
  sd1 <- 1
  sd2 <- 1
  
  # Densities
  d1 <- dnorm(x, mean = mu1, sd = sd1)
  d2 <- dnorm(x, mean = mu2, sd = sd2)
  
  # Derivatives of each component
  d1_prime <- - (x - mu1) / (sd1^2) * d1
  d2_prime <- - (x - mu2) / (sd2^2) * d2
  
  # Mixture density
  pi_x <- w1 * d1 + w2 * d2
  pi_prime_x <- w1 * d1_prime + w2 * d2_prime
  
  # Gradient of log-density
  return(pi_prime_x / pi_x)
}

chain2 <- HMC_Sampler(L = 20, eps = .1, n = 1e5)

par(mfrow = c(1, 1))
hist(chain2, breaks = 50, probability = TRUE,
     main = "Samples for L = 20, eps = 0.1", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2)
lines(density(chain2), col = "red", lwd = 3, lt = 2) # Overlay density estimate
x2 <- seq(-3, 11, length.out = 200)
y2 <- target(x2)
lines(x2, y2, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("HMC Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, lty = c(2,1), bty = "n",
       cex = 1.5)

acf(chain, lag.max = 50, main = "ACF of HMC Chain", cex.main = 2,
    col = "lightblue3", lwd = 2)

###############################################################################

## Changing the means to 0 and 10
target <- function(x){
  # Simple bi-modal target
  return(
    0.6*dnorm(x, mean = 0, sd = 1) + 0.4*dnorm(x, mean = 10, sd = 1)
  )
}

grad_log_target <- function(x) {
  # Components
  w1 <- 0.6
  w2 <- 0.4
  mu1 <- 0
  mu2 <- 10
  sd1 <- 1
  sd2 <- 1
  
  # Densities
  d1 <- dnorm(x, mean = mu1, sd = sd1)
  d2 <- dnorm(x, mean = mu2, sd = sd2)
  
  # Derivatives of each component
  d1_prime <- - (x - mu1) / (sd1^2) * d1
  d2_prime <- - (x - mu2) / (sd2^2) * d2
  
  # Mixture density
  pi_x <- w1 * d1 + w2 * d2
  pi_prime_x <- w1 * d1_prime + w2 * d2_prime
  
  # Gradient of log-density
  return(pi_prime_x / pi_x)
}

chain3 <- HMC_Sampler(L = 20, eps = .1, n = 1e5)

par(mfrow = c(1, 1))
hist(chain3, breaks = 50, probability = TRUE,
     main = "Samples for L = 20, eps = 0.1", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2,
     xlim = c(-4,13))
lines(density(chain3), col = "red", lwd = 3, lt = 2) # Overlay density estimate
x3 <- seq(-3, 13, length.out = 200)
y3 <- target(x3)
lines(x3, y3, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("HMC Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, lty = c(2,1), bty = "n",
       cex = 1.5)

acf(chain, lag.max = 50, main = "ACF of HMC Chain", cex.main = 2,
    col = "lightblue3", lwd = 2)

###############################################################################

## Comparing the 3 plots
par(mfrow = c(3, 1))
hist(chain, breaks = 50, probability = TRUE,
     main = "Samples for mu1 = 0, mu2 = 5", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2,
     xlim = c(-4,13))
lines(density(chain), col = "red", lwd = 3, lt = 2) # Overlay density estimate
lines(x, y, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("HMC Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, lty = c(2,1), bty = "n",
       cex = 1.5)

hist(chain2, breaks = 50, probability = TRUE,
     main = "Samples for mu1 = 0, mu2 = 8", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2,
     xlim = c(-4,13))
lines(density(chain2), col = "red", lwd = 3, lt = 2) # Overlay density estimate
lines(x2, y2, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("HMC Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, lty = c(2,1), bty = "n",
       cex = 1.5)

hist(chain3, breaks = 50, probability = TRUE,
     main = "Samples for mu1 = 0, mu2 = 10", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2,
     xlim = c(-4,13))
lines(density(chain3), col = "red", lwd = 3, lt = 2) # Overlay density estimate
lines(x3, y3, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("HMC Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, lty = c(2,1), bty = "n",
       cex = 1.5)

###############################################################################

### Visualising the engery barrier
x <- seq(-3, 13, length.out = 200)
y <- log_target(x)
par(mfrow = c(1, 1))
plot(x, -y, type = "l", lwd = 2, col = "blue", 
     main = "Potential Energy", xlab = "x", ylab = "Log-Density",
     cex.main = 2)
abline(v = 5, col = "red", lty = 2, lwd = 2)
points(0, -log_target(0), col = "gold2", pch = 19, cex = 3)
legend("topright", legend = c("Energy Barrier", "Particle"), 
       col = c("red", "gold2"), lwd = c(2,NA), lty = c(2, NA), pch = c(NA, 19),
       bty = "n", cex = 1.5)
