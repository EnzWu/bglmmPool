# internal cache for compiled Stan models
.stan_cache <- new.env(parent = emptyenv())

#' Get cached Stan model for pooling
#' @keywords internal
get_pooling_stan_model <- function() {
  key <- "Pooling.stan"

  if (exists(key, envir = .stan_cache, inherits = FALSE)) {
    return(get(key, envir = .stan_cache, inherits = FALSE))
  }

  stan_path <- system.file("stan", "Pooling.stan", package = "bglmmPool")
  if (!nzchar(stan_path)) stan_path <- file.path("inst", "stan", "Pooling.stan")
  if (!file.exists(stan_path)) stop("Cannot find Stan file: ", stan_path)

  sm <- rstan::stan_model(file = stan_path)
  assign(key, sm, envir = .stan_cache)
  sm
}

