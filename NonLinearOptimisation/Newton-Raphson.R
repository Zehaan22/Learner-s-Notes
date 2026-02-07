## ===============================
## 1. Generate data
## ===============================

set.seed(123)

n <- 200
t <- seq(0.1, 10, length.out = n)

alpha_true <- 2.0
beta_true  <- 1.5
sigma      <- 0.3

mu <- alpha_true * sin(t) + 1 / (beta_true * t)
y  <- mu + rnorm(n, 0, sigma)

## ===============================
## 2. Define model components
## ===============================

model_mean <- function(theta, t) {
  alpha <- theta[1]
  beta  <- theta[2]
  alpha * sin(t) + 1 / (beta * t)
}

## Least squares objective
objective <- function(theta, y, t) {
  r <- y - model_mean(theta, t)
  sum(r^2)
}

## Gradient
gradient <- function(theta, y, t) {
  alpha <- theta[1]
  beta  <- theta[2]
  
  r <- y - model_mean(theta, t)
  
  d_mu_d_alpha <- sin(t)
  d_mu_d_beta  <- -1 / (beta^2 * t)
  
  grad_alpha <- -2 * sum(r * d_mu_d_alpha)
  grad_beta  <- -2 * sum(r * d_mu_d_beta)
  
  c(grad_alpha, grad_beta)
}

## Hessian
hessian <- function(theta, y, t) {
  alpha <- theta[1]
  beta  <- theta[2]
  
  r <- y - model_mean(theta, t)
  
  d_mu_d_alpha <- sin(t)
  d_mu_d_beta  <- -1 / (beta^2 * t)
  d2_mu_d_beta2 <- 2 / (beta^3 * t)
  
  H_aa <- 2 * sum(d_mu_d_alpha^2)
  H_ab <- 2 * sum(d_mu_d_alpha * d_mu_d_beta)
  H_bb <- 2 * sum(d_mu_d_beta^2 - r * d2_mu_d_beta2)
  
  matrix(c(H_aa, H_ab,
           H_ab, H_bb),
         nrow = 2, byrow = TRUE)
}

## ===============================
## 3. Newton–Raphson optimizer
## ===============================

newton_raphson <- function(theta_init, y, t, tol = 1e-6, max_iter = 50) {
  theta <- theta_init
  history <- matrix(NA, nrow = max_iter, ncol = 2)
  
  for (k in 1:max_iter) {
    g <- gradient(theta, y, t)
    H <- hessian(theta, y, t)
    
    step <- solve(H, g)
    theta_new <- theta - step
    
    history[k, ] <- theta_new
    
    if (max(abs(theta_new - theta)) < tol) {
      history <- history[1:k, , drop = FALSE]
      return(list(theta = theta_new,
                  history = history,
                  iter = k))
    }
    
    theta <- theta_new
  }
  
  stop("Newton–Raphson did not converge")
}

## ===============================
## 4. Run optimization
## ===============================

theta_init <- c(alpha = 1, beta = 1)

fit <- newton_raphson(theta_init, y, t)

fit$theta
#> Estimated parameters

## ===============================
## 5. Diagnostic plots
## ===============================

pdf("Images/Newton_Raphson_Diagnostics.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

## (a) Data + fitted curve
t_grid <- seq(min(t), max(t), length.out = 500)

plot(t, y,
     pch = 16, col = rgb(0, 0, 0, 0.4),
     xlab = "t", ylab = "y",
     main = "Non-linear fit via Newton–Raphson")

lines(t_grid,
      model_mean(c(alpha_true, beta_true), t_grid),
      col = "steelblue", lwd = 5, lty = 2)

lines(t_grid,
      model_mean(fit$theta, t_grid),
      col = "firebrick", lwd = 2)

legend("topright",
       legend = c("True mean", "Newton–Raphson fit"),
       col = c("steelblue", "firebrick"),
       lwd = c(5,2), lty = c(2,1), bty = "n")

## (b) Parameter trajectories
plot(fit$history[,1], type = "o",
     pch = 16, col = "steelblue",
     ylim = range(fit$history),
     xlab = "Iteration", ylab = "Parameter value",
     main = "Newton–Raphson convergence")

lines(fit$history[,2], type = "o",
      pch = 16, col = "firebrick")

abline(h = alpha_true, col = "steelblue", lty = 2)
abline(h = beta_true, col = "firebrick", lty = 2)

legend("right",
       legend = c(expression(alpha), expression(beta)),
       col = c("steelblue", "firebrick"),
       pch = 16, bty = "n")

dev.off()