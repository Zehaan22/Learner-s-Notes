## ===============================
## 1. Generate data
## ===============================

set.seed(42)

n <- 150
t <- seq(0, 5, length.out = n)

theta_true <- c(theta1 = 3, theta2 = 1.2, theta3 = 0.5)
sigma <- 0.2

mu <- theta_true[1] * exp(-theta_true[2] * t) + theta_true[3]
y  <- mu + rnorm(n, 0, sigma)

## ===============================
## 2. Model and residuals
## ===============================

model_mean <- function(theta, t) {
  theta[1] * exp(-theta[2] * t) + theta[3]
}

residuals_fn <- function(theta, y, t) {
  y - model_mean(theta, t)
}

## ===============================
## 3. Jacobian (key ingredient)
## ===============================

jacobian <- function(theta, t) {
  theta1 <- theta[1]
  theta2 <- theta[2]
  
  J <- cbind(
    exp(-theta2 * t),                # d mu / d theta1
    -theta1 * t * exp(-theta2 * t),  # d mu / d theta2
    rep(1, length(t))                # d mu / d theta3
  )
  
  J
}

## ===============================
## 4. Gauss–Newton optimizer
## ===============================

gauss_newton <- function(theta_init, y, t, tol = 1e-6, max_iter = 50) {
  
  theta <- theta_init
  history <- matrix(NA, nrow = max_iter, ncol = length(theta))
  
  for (k in 1:max_iter) {
    r <- residuals_fn(theta, y, t)
    J <- jacobian(theta, t)
    
    ## Normal equations: (J^T J) delta = J^T r
    step <- solve(t(J) %*% J, t(J) %*% r)
    
    theta_new <- theta + step
    history[k, ] <- theta_new
    
    if (max(abs(theta_new - theta)) < tol) {
      history <- history[1:k, , drop = FALSE]
      return(list(theta = theta_new,
                  history = history,
                  iter = k))
    }
    
    theta <- theta_new
  }
  
  stop("Gauss–Newton did not converge")
}

## ===============================
## 5. Run Gauss–Newton
## ===============================

theta_init <- c(1, 0.5, 0)   # deliberately rough start
fit <- gauss_newton(theta_init, y, t)

fit$theta
# Estimated parameters

## ===============================
## 6. Visualization
## ===============================

pdf("Images/gauss_newton_diagnostics.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

## (a) Fit vs data
t_grid <- seq(min(t), max(t), length.out = 400)

plot(t, y,
     pch = 16, col = rgb(0, 0, 0, 0.4),
     xlab = "t", ylab = "y",
     main = "Gauss–Newton fit")

lines(t_grid,
      model_mean(theta_true, t_grid),
      col = "steelblue", lwd = 2)

lines(t_grid,
      model_mean(fit$theta, t_grid),
      col = "firebrick", lwd = 2)

legend("topright",
       legend = c("True mean", "Gauss–Newton fit"),
       col = c("steelblue", "firebrick"),
       lwd = 2, bty = "n")

## (b) Parameter trajectories
matplot(fit$history, type = "o", pch = 16,
        col = c("steelblue", "firebrick", "darkgreen"),
        xlab = "Iteration", ylab = "Parameter value",
        main = "Gauss–Newton convergence")

abline(h = theta_true, col = c("steelblue", "firebrick", "darkgreen"),
       lty = 2)

legend("right",
       legend = expression(theta[1], theta[2], theta[3]),
       col = c("steelblue", "firebrick", "darkgreen"),
       pch = 16, bty = "n")
dev.off()
