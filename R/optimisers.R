adagrad <- function (learning_rate = 0.8,
                     initial_accumulator_value = 0.1,
                     use_locking = TRUE) {

  obj <- list(parameters = list(learning_rate = learning_rate,
                                initial_accumulator_value =
                                  initial_accumulator_value, 
                                use_locking = use_locking),
              class = adagrad_optimiser)

  class(obj) <- c('adagrad optimiser', 'optimiser')
  obj

}
