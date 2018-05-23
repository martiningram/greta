#' @importFrom gpexp predict_points ard_kernel
compute_gp <- function(x, y, x_to_predict, sigma) {

  sigma_inv <- diag(1 / sigma, 1, 1)
  k_fun <- function (x1, x2) ard_kernel(x1, x2, sigma_inv)

  prediction <- predict_points(as.matrix(x), x_to_predict, sqrt(0.1),
                               as.matrix(y), k_fun, mean_centre = FALSE)

  return(prediction)
}

compute_beta <- function(i, d, delta = 0.1) {
  # i is the iteration
  # d is the dimension
  beta_ip1 <- 2 * log((i + 1) ^ (d / 2 + 2) * pi^2 / (3 * delta))
  return(beta_ip1)
}

compute_u <- function(scaled_means, covariances, beta, p_i) {
  # Compute u using the UCB acquisition function
  sigmas <- sqrt(diag(covariances))
  return(scaled_means + p_i * sqrt(beta) * sigmas)
}

initialise_data <- function() {

  initial_data <- list('i' = 1,
                       'r' = NULL,
                       's' = NULL,
                       'gamma' = NULL)

  return(initial_data)

}

opt_step <- function(gamma, r, data, alpha = 4, k = 100, min_L = 1,
                     max_L = 200, kappa = 0.2) {
  # Gamma is the current setting of the hyperparameters; a column vector
  # r is the "reward" of the current step, assumed scalar
  # data is a list containing:
  # field "i", the current iteration (starting at 1)
  # field "gamma", a matrix of past gamma settings (i x D)
  # field "r", the vector of past rewards (i x 1)
  # field "s", the current scaling factor
  to_predict <- as.matrix(seq(min_L, max_L))

  d <- length(gamma)

  if (data$i == 1) {
    # This is, by definition, the scaling
    # TODO: Think about whether this really makes sense
    max_r <- r
    data$s <- alpha / r
    data$gamma <- matrix(gamma, nrow = 1)
    data$r <- c(r)
  } else {
    max_r <- max(data$r)
    if (r > max_r) {
      data$s <- alpha / r
    }
    # Update the data
    data$gamma <- rbind(data$gamma, gamma)
    data$r <- c(data$r, r)
  }

  u <- runif(1)
  p_i <- max(data$i - k + 1, 1)^(-0.5)
  sigma <- (kappa * (max_L - min_L))^2

  if (u < p_i) {
    # Choose a new gamma
    # First, fit the GP to our data
    gp_results <- compute_gp(data$gamma, data$r, to_predict, sigma)

    # Scale the means using s
    scaled_means <- gp_results$mean * data$s

    # Compute beta
    beta <- compute_beta(data$i, d)

    # Compute the acquisition function
    u <- compute_u(scaled_means, gp_results$cov, beta, p_i)

    # Find its maximum
    argmax_u <- which.max(u)

    # Select the new gamma
    new_gamma <- to_predict[argmax_u]
  } else {
    # Keep current setting
    new_gamma <- gamma
  }

  # Increment i
  data$i <- data$i + 1
  return(list('new_gamma' = new_gamma, 'data' = data))
}
