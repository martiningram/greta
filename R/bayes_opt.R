#' @importFrom gpexp predict_points ard_kernel
compute_gp <- function(x, y, x_to_predict, sigma, obs_var = 0.1) {

  sigma_inv <- sigma
  diag(sigma_inv) <- 1 / diag(sigma)
  k_fun <- function (x1, x2) ard_kernel(x1, x2, sigma_inv)

  prediction <- predict_points(as.matrix(x), x_to_predict, sqrt(obs_var),
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

initialise_tuning_data <- function() {

  initial_data <- list('i' = 1,
                       'r' = NULL,
                       's' = NULL,
                       'gamma' = NULL)

  return(initial_data)

}

#' @importFrom gpexp plot_gp
opt_step <- function(gamma, r, data, alpha = 4, k = 100, min_L = 2,
                     max_L = 200, min_eps = 0.0001, max_eps = 0.1, num_eps = 50,
                     kappa = 0.2) {
  # Gamma is the current setting of the hyperparameters; a column vector
  # r is the "reward" of the current step, assumed scalar
  # data is a list containing:
  # field "i", the current iteration (starting at 1)
  # field "gamma", a matrix of past gamma settings (i x D)
  # field "r", the vector of past rewards (i x 1)
  # field "s", the current scaling factor
  to_predict <- 
    as.matrix(expand.grid(seq(min_L, max_L, 10), 
                          seq(min_eps, max_eps, length.out = num_eps)))

  d <- length(gamma)

  if (data$i == 1) {
    # This is, by definition, the scaling
    # TODO: Think about whether this really makes sense
    max_r <- r
    data$gamma <- t(matrix(gamma))
    data$r <- c(r)
  } else {
    max_r <- max(data$r)
    # Update the data
    data$gamma <- rbind(data$gamma, gamma)
    data$r <- c(data$r, r)
  }

  if (r >= max_r) {
    if (r > 0) {
      data$s <- alpha / r
    } else {
      # Avoid dividing by zero
      data$s <- 1
    }
  }

  u <- runif(1)
  p_i <- max(data$i - k + 1, 1)^(-0.5)
  sigma <- diag(c((kappa * (max_L - min_L))^2, 
                  (kappa * (max_eps - min_eps))^2))
  obs_var <- 0.01 # Observation noise TODO: Check!

  if (u < p_i) {
    # Choose a new gamma
    # First, fit the GP to our data
    gp_results <- compute_gp(data$gamma, data$r, to_predict, sigma, 
                             obs_var = obs_var)

    # Scale the means using s
    scaled_means <- gp_results$mean * data$s

    # Out of interest, plot the gp
    # plot_gp(to_predict, scaled_means, gp_results$cov, sqrt(obs_var), 
    #         x_train = data$gamma, y_train = data$r * data$s,
    #         save_to = paste0('/tmp/gp_plot_', data$i, '.png'))

    # Compute beta
    beta <- compute_beta(data$i, d)

    # Compute the acquisition function
    u <- compute_u(scaled_means, gp_results$cov, beta, p_i)

    # Find its maximum
    argmax_u <- which.max(u)

    # Select the new gamma
    new_gamma <- to_predict[argmax_u, ]

    # Save this for debugging
    saveRDS(list('gp_results' = gp_results,
                 'data' = data,
                 'to_predict' = to_predict,
                 's' = data$s,
                 'beta' = beta,
                 'u' = u,
                 'argmax_u' = argmax_u,
                 'new_gamma' = new_gamma),
            paste0('/tmp/debug_bayes_opt_', data$i, '.Rds'))

  } else {
    # Keep current setting
    new_gamma <- gamma
  }

  # Increment i
  data$i <- data$i + 1
  return(list('new_gamma' = new_gamma, 'data' = data))
}
