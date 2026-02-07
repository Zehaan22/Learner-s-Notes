## ===============================
## 1. Generate exponential data
## ===============================

set.seed(123)

n <- 200
t <- 0:(n - 1)

theta_true <- 0.18
A_true <- 2.0
B_true <- -1.2
sigma <- 0.5

y_clean <- A_true * cos(theta_true * t) +
  B_true * sin(theta_true * t)

y <- y_clean + rnorm(n, 0, sigma)

## ===============================
## 2. Alt-Opt Algorithm
## ===============================

fit_amplitudes <- function(theta, t, y) {
  X <- cbind(cos(theta * t), sin(theta * t))
  coef <- qr.coef(qr(X), y)
  fitted <- X %*% coef
  rss <- sum((y - fitted)^2)
  list(A = coef[1], B = coef[2], rss = rss, fitted = fitted)
}

theta_grid <- seq(0.01, pi, length.out = 300)

profile_rss <- numeric(length(theta_grid))

for (i in seq_along(theta_grid)) {
  profile_rss[i] <- fit_amplitudes(theta_grid[i], t, y)$rss
}

theta_init <- theta_grid[which.min(profile_rss)]

refine_theta <- function(theta_start, t, y) {
  opt <- optimize(
    f = function(th) fit_amplitudes(th, t, y)$rss,
    interval = c(theta_start - 0.1, theta_start + 0.1)
  )
  opt$minimum
}

theta_hat <- refine_theta(theta_init, t, y)
fit_final <- fit_amplitudes(theta_hat, t, y)


## ===============================
## 3. Plotting
## ===============================

pdf("Images/Alt_opt_trig_models.pdf", width = 12, height = 9)
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

## (a) Signal and fit
plot(t, y,
     pch = 16, col = rgb(0,0,0,0.4),
     xlab = "t", ylab = "y",
     main = "Trigonometric signal fit")

lines(t, y_clean, col = "steelblue", lwd = 2)
lines(t, fit_final$fitted, col = "firebrick", lwd = 2, lty = 2)

legend("topright",
       legend = c("True signal", "Alt-opt fit"),
       col = c("steelblue", "firebrick"),
       lwd = 2, lty = c(1,2), bty = "n")

## (b) Profile objective (periodogram-style)
plot(theta_grid, profile_rss,
     type = "l", lwd = 2,
     xlab = expression(theta),
     ylab = "Profile RSS",
     main = "Profile objective over frequency")

abline(v = theta_true, col = "steelblue", lty = 2)
abline(v = theta_hat, col = "firebrick", lty = 2)

legend("topright",
       legend = c("True frequency", "Estimated frequency"),
       col = c("steelblue", "firebrick"),
       lty = 2, bty = "n")

## (c) Periodogram (same idea, different scaling)
periodogram <- abs(colSums(y * exp(-1i * outer(t, theta_grid))))^2 / n

plot(theta_grid, periodogram,
     type = "l", lwd = 2,
     xlab = expression(theta),
     ylab = "Periodogram",
     main = "Periodogram view")

abline(v = theta_hat, col = "firebrick", lty = 2)
dev.off()
