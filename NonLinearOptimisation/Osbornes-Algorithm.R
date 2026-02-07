## ===============================
## 1. Generate exponential data
## ===============================

set.seed(123)

n <- 60
t <- 0:(n - 1)

a_true <- c(a0 = 0.5, a1 = 2.0, a2 = 1.0)
b_true <- c(b1 = -0.06, b2 = -0.25)
sigma  <- 0.05

y_clean <- a_true["a0"] +
  a_true["a1"] * exp(b_true["b1"] * t) +
  a_true["a2"] * exp(b_true["b2"] * t)

y <- y_clean + rnorm(n, 0, sigma)


## ===============================
## 2. OSBORNE's Algorithm
## ===============================
## Design matrix for fixed b
make_X <- function(t, b) {
  cbind(1, exp(b[1] * t), exp(b[2] * t))
}

## Linear LS for amplitudes (QR-stable)
fit_a_given_b <- function(X, y) {
  qr.coef(qr(X), y)
}

jacobian_r_b <- function(t, y, b, eps = 1e-6) {
  X0 <- make_X(t, b)
  a0 <- fit_a_given_b(X0, y)
  r0 <- y - X0 %*% a0
  
  J <- matrix(0, nrow = length(y), ncol = length(b))
  
  for (j in seq_along(b)) {
    bp <- b
    bp[j] <- bp[j] + eps
    Xp <- make_X(t, bp)
    ap <- fit_a_given_b(Xp, y)
    rp <- y - Xp %*% ap
    J[, j] <- (rp - r0) / eps
  }
  
  J
}

osborne <- function(t, y, b_init,
                    max_iter = 100, tol = 1e-8,
                    lambda0 = 1e-6, verbose = TRUE) {
  
  b <- b_init
  lambda <- lambda0
  history <- matrix(NA, max_iter, length(b))
  
  for (k in 1:max_iter) {
    X <- make_X(t, b)
    a <- fit_a_given_b(X, y)
    r <- y - X %*% a
    ss <- sum(r^2)
    
    Jb <- jacobian_r_b(t, y, b)
    
    ## Damped Gauss–Newton step (QR-safe)
    H <- crossprod(Jb) + lambda * diag(length(b))
    g <- crossprod(Jb, r)
    
    delta <- tryCatch(
      solve(H, g),
      error = function(e) rep(0, length(b))
    )
    
    ## Step halving if needed
    alpha <- 1
    repeat {
      b_new <- b + alpha * delta
      Xn <- make_X(t, b_new)
      an <- fit_a_given_b(Xn, y)
      rn <- y - Xn %*% an
      if (sum(rn^2) < ss || alpha < 1e-6) break
      alpha <- alpha / 2
    }
    
    if (sum(rn^2) < ss) {
      lambda <- max(lambda / 10, 1e-12)
      b <- b_new
    } else {
      lambda <- lambda * 10
    }
    
    history[k, ] <- b
    
    if (verbose) {
      cat(sprintf(
        "iter %3d | SSE = %.6f | step = %.2e | lambda = %.1e\n",
        k, sum(rn^2), alpha, lambda
      ))
    }
    
    if (max(abs(alpha * delta)) < tol) break
  }
  
  X_final <- make_X(t, b)
  a_final <- fit_a_given_b(X_final, y)
  
  list(
    a = a_final,
    b = b,
    fitted = as.numeric(X_final %*% a_final),
    history = history[!is.na(history[,1]), , drop = FALSE]
  )
}


## ===============================
## 3. Obsborne fit
## ===============================

fit_osborne <- osborne(
  t, y,
  b_init = c(b1 = -0.1, b2 = -0.3),
  max_iter = 200,
  verbose = TRUE
)


## ===============================
## 4. Plotting
## ===============================

pdf("Images/Osborne_Diagnostics.pdf", width = 10, height = 5)
par(mfrow = c(1, 1))

## (a) Fit
plot(t, y,
     pch = 16, col = rgb(0,0,0,0.4),
     xlab = "t", ylab = "y",
     main = "Osborne's algorithm fit")

lines(t,
      y_clean,
      col = "steelblue", lwd = 2)

lines(t,
      fit_osborne$fitted,
      col = "firebrick", lwd = 2, lty = 2)

legend("topright",
       legend = c("True signal", "Osborne fit"),
       col = c("steelblue", "firebrick"),
       lwd = 2, lty = c(1,2), bty = "n")

dev.off()
