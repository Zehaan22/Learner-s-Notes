####### File to illustrate a toy example of the MH algorithm

target <- function(x){
  # Simple bi-modal target
  return(
    0.6*dnorm(x, mean = 0, sd = 1) + 0.4*dnorm(x, mean = 5, sd = 1)
  )
}

## Visualising target
x <- seq(-3, 8, length.out = 200)
y <- target(x)
plot(x, y, type = "l", lwd = 2, col = "blue", 
     main = "Target Distribution", xlab = "x", ylab = "Density")


## MH Sampler Function
MHSampler <- function(n, h = 1){
  chain <- numeric(n)
  curr <- 0 # Starting Value
  
  accept <- 0
  for(i in 1:n){
    proposal <- rnorm(1, mean = curr, sd = h) # Propose a new value
    
    alpha <- min(1, target(proposal) / target(curr)) # Acceptance ratio
    
    if(runif(1) < alpha){ # Accept or reject
      curr <- proposal
      accept <- accept + 1 # Increment acceptance count
    }
    
    chain[i] <- curr # Store the current value
  }
  
  print(paste("Acceptance Rate:", accept / n)) 
  return(chain)
}

chain <- MHSampler(1e5, h = 1)

## Visulising the chain
par(mfrow = c(2, 1))
plot(chain, type = "l", lwd = 1, col = "red", 
     main = "MH Chain trace plot", xlab = "Iteration", ylab = "Value",
     cex.main = 2)
hist(chain, breaks = 50, probability = TRUE,
     main = "MH Chain Samples", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2)
lines(density(chain), col = "red", lwd = 2) # Overlay density estimate
lines(x, y, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("MH Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, bty = "n")

## Visualsing the actual samples
par(mfrow = c(2, 1))
plot(x, y, type = "l", lwd = 2, col = "blue", 
     main = "Samples Visualised", xlab = "x", ylab = "Density",
     cex.main = 2)
points(sample(chain,1e3), (1:1e3)/2e4, pch = 19, col = "red", cex = 0.3)
legend("topright", legend = c("Target Distribution", "MH Samples"), 
       col = c("blue", "red"), pch = 19, bty = "n")

acf(chain, main = "ACF of MH Chain", lag.max = 50, 
    cex.main = 20, col = "lightblue3", lwd = 2)
