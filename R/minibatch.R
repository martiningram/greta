from_first_dim <- function (x, idx) {
  # This function picks from the first dimension of an array.
  ndim <- length(dim(x))
  args <- list(x = x, i = idx)
  missings <- replicate(ndim - 1, quote(expr = ))
  args <- c(args, missings)
  args$drop <- FALSE
  do.call('[', args)
}

get_minibatch <- function(data_arrays, batch_size = 8) {

  # The arrays provided all have to be greta data arrays.
  are_data <- vapply(data_arrays, function (x) inherits(x$node, 'data_node'),
                     FUN.VALUE = FALSE)

  if (!all(are_data)) {
    stop('All arrays provided must be data arrays!', call. = FALSE)
  }

  # The arrays must also all have the same first dimension
  first_dim <- vapply(data_arrays, function (x) dim(x)[1],
                      FUN.VALUE = 1)

  if(length(unique(first_dim)) != 1) {
    stop('All arrays provided must have the same first dimension!', 
         call. = FALSE)
  }

  # With that out of the way, we can sample some random indices
  total_num <- dim(data_arrays[[1]])[1]

  # Sample without replacement
  picked <- sample.int(total_num, batch_size)

  batch <- lapply(data_arrays, function(x) from_first_dim(x, picked))

  batch
}
