##############################################################################
## Code to generate theoretical Physics visuals

## Generating the contour 
t <- seq(0,2*pi, length = 100)
start <- 0
end <- 1.5
t.prime <- seq(start, end, length = 100)
x <- sin(t)
p <- cos(t)

par(mfrow = c(1, 1))
plot(x, p, type = "l", lwd = 2, col = "blue", 
     main = "Position-Momentum Plot", xlab = "x", ylab = "p",
     cex.main = 2,
     asp = 1)
### Plotting a sample trajectory
lines(sin(t.prime), cos(t.prime), col = "green4", lwd = 4, lty = 2)
points(sin(start), cos(start), pch = 19, col = "red", cex = 3)
points(sin(end), cos(end), pch = 19, col = "gold3", cex = 3)
legend("topright", legend = c("Trajectory", "Start", "End"), 
       col = c("green4", "red", "gold3"), pch = 19, bty = "n",
       lty = c(2, NA, NA), lwd = c(4, NA, NA))

### Generating the potential plot
x <- seq(-3, 3, length.out = 100)
plot(x, 0.5*x^2, type = "l", lwd = 2, col = "blue", 
     main = "Potential Energy Plot", xlab = "x", ylab = "U(x)",
     cex.main = 2)
## Plotting the same trajectory
lines(sin(t.prime), 0.5*sin(t.prime)^2, col = "green4", lwd = 4, lty = 2)
points(sin(start), 0.5*sin(start)^2, pch = 19, col = "red", cex = 3)
points(sin(end), 0.5*sin(end)^2, pch = 19, col = "gold3", cex = 3)
legend("top", legend = c("Trajectory", "Start", "End"), 
       col = c("green4", "red", "gold3"), pch = 19, bty = "n",
       lty = c(2, NA, NA), lwd = c(4, NA, NA))

### Plotting the change of state
t <- seq(0,2*pi, length = 100)
start <- 0
end <- 1.5
t.prime <- seq(start, end, length = 100)
x <- sin(t)
p <- cos(t)
par(mfrow = c(1, 1))
plot(1.4*x, 1.4*p, type = "l", lwd = 2, col = "purple4", 
     main = "Momentum Jump Plot", xlab = "x", ylab = "p",
     cex.main = 2,
     asp = 1)
lines(x, p, type = "l", lwd = 2, col = "blue")
lines(1.4*sin(t.prime), 1.4*cos(t.prime), col = "green4", lwd = 4, lty = 2)
lines(c(sin(start), 1.4*sin(start)), c(cos(start), 1.4*cos(start)), col = "coral2", lwd = 4, lty = 2)
points(sin(start), cos(start), pch = 19, col = "red", cex = 3)
points(1.4*sin(start), 1.4*cos(start), pch = 19, col = "red3", cex = 3)
points(1.4*sin(end), 1.4*cos(end), pch = 19, col = "gold3", cex = 3)
legend("bottomleft", 
       legend = c("Trajectory", "Start", "New Momentum", "End"), 
       col = c("green4", "red", "red3", "gold3"), pch = 19, bty = "n",
       lty = c(2, NA, NA), lwd = c(4, NA, NA))


##############################################################################
## Leap Frog Visuals
## Generating the contour
t <- seq(0,2*pi, length = 100)
x <- sin(t)
p <- cos(t)
start <- 0
par(mfrow = c(1, 1))
plot(x, p, type = "l", lwd = 2, col = "purple4", 
     main = "Leap Frog Approximations", xlab = "x", ylab = "p",
     cex.main = 2,
     asp = 1)
## Starting point
points(sin(start), cos(start), pch = 19, col = "red", cex = 3)

## Leap Frog Steps
x.seq <- c()
p.seq <- c()
x.st <- sin(start)
p.st <- cos(start)
eps <- 0.1
L <- 10
for(i in 1:L){
  p.st <- p.st - eps*0.5*x.st
  x.st <- x.st + eps*p.st
  p.st <- p.st - eps*0.5*x.st
  x.seq[i] <- x.st
  p.seq[i] <- p.st
}
## Plotting the leap frog steps
points(x.seq, p.seq, pch = 19, col = "green4", cex = 1)

x.st <- sin(start)
p.st <- cos(start)
eps <- 0.5
L <- 10
for(i in 1:L){
  p.st <- p.st - eps*0.5*x.st
  x.st <- x.st + eps*p.st
  p.st <- p.st - eps*0.5*x.st
  x.seq[i] <- x.st
  p.seq[i] <- p.st
}
## Plotting the leap frog steps
points(x.seq, p.seq, pch = 19, col = "gold2", cex = 1)

x.st <- sin(start)
p.st <- cos(start)
eps <- 0.05
L <- 10
for(i in 1:L){
  p.st <- p.st - eps*0.5*x.st
  x.st <- x.st + eps*p.st
  p.st <- p.st - eps*0.5*x.st
  x.seq[i] <- x.st
  p.seq[i] <- p.st
}
## Plotting the leap frog steps
points(x.seq, p.seq, pch = 19, col = "coral2", cex = 1)

legend("topright", 
       legend = c("Eps = 0.1", "Eps = 0.5", "Eps = 0.05", "Start"), 
       col = c("green4","gold2","coral2", "red"), 
       pch = 19, bty = "n",
       lty = c(2, 2, 2, NA), lwd = c(4, 4, 4, NA))


##############################################################################

#### HMC CODE

## HMC Sampler Function
HMCsampler <- function(s = 1, n = 1e4)
{
  qt <- numeric(length = n)
  qt[1] <- 0 # starting here
  # don't need a starting value of p
  for(k in 2:n)
  {
    # new momentum
    p <- rnorm(1)
    #initial conditions for new H
    r2 <- qt[k-1]**2 + p**2
    
    # choosing +- with probability 1/2
    r <- sample(c(sqrt(r2), -sqrt(r2)), size = 1)
    a <- acos(qt[k-1]/r)
    
    # simulating Hamiltonian forward s time units
    qt[k] <- r*cos(a + s)
    
  }
  return(qt)
}

chain1 <- HMCsampler(s = 5, n = 1e5)
chain2 <- HMCsampler(s = 1, n = 1e5)
chain3 <- HMCsampler(s = 0.1, n = 1e5)

## Visulising the chain
par(mfrow = c(3, 1))
hist(chain1, breaks = 50, probability = TRUE,
     main = "Samples for s = 5", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2)
lines(density(chain1), col = "red", lwd = 3, lt = 2) # Overlay density estimate
x <- seq(-3, 3, length.out = 200)
y <- dnorm(x)
lines(x, y, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("MH Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, bty = "n",
       cex = 1.5)

hist(chain2, breaks = 50, probability = TRUE,
     main = "Samples for s = 1", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2)
lines(density(chain2), col = "red", lwd = 3, lt = 2) # Overlay density estimate
x <- seq(-3, 3, length.out = 200)
y <- dnorm(x)
lines(x, y, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("MH Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, bty = "n",
       cex = 1.5)

hist(chain3, breaks = 50, probability = TRUE,
     main = "Samples for s = 0.1", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2)
lines(density(chain3), col = "red", lwd = 3, lt = 2) # Overlay density estimate
x <- seq(-3, 3, length.out = 200)
y <- dnorm(x)
lines(x, y, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("MH Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, bty = "n",
       cex = 1.5)


## ACF plots
par(mfrow = c(3, 1))
acf(chain1, main = "ACF of HMC Chain",
    lag.max = 10,
     col = "red", lwd = 2)
acf(chain2, main = "ACF of HMC Chain",
    lag.max = 10,
     col = "red", lwd = 2)
acf(chain3, main = "ACF of HMC Chain",
    lag.max = 10,
     col = "red", lwd = 2)

###############################################################################
### Negative ACF plots

chain1 <- HMCsampler(s = 10, n = 1e5)
par(mfrow = c(1, 1))
acf(chain1, main = "ACF of HMC Chain",
    lag.max = 20,
    col = "red", lwd = 2)

hist(chain1, breaks = 50, probability = TRUE,
     main = "Samples for s = 10", xlab = "Value", ylab = "Density",
     col = "lightblue", border = "black",
     cex.main = 2)
lines(density(chain1), col = "red", lwd = 3, lt = 2) # Overlay density estimate
x <- seq(-3, 3, length.out = 200)
y <- dnorm(x)
lines(x, y, col = "blue", lwd = 2) # Overlay target distribution     
legend("topright", legend = c("MH Chain", "Target Distribution"), 
       col = c("red", "blue"), lwd = 2, bty = "n",
       cex = 1.5)