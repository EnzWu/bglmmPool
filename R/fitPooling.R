#' Fit pooling model (rstan)
#'
#' @name fit_pooling
#' @title Fit pooling model (rstan)
#' @description Fits the pooled/blinded BGLMM using the Stan program shipped with the package.
#'
#' @param data Named list for Stan data
#' @param iter Total iterations per chain (warmup + sampling)
#' @param warmup Warmup iterations per chain
#' @param chains Number of MCMC chains
#' @param seed Sampling seed
#' @param show_messages Display rstan::sampling() message or not
#' @return A stanfit object
#' @export
fit.pooling <- function(data,
                        iter = 40000,
                        warmup = 20000,
                        chains = 4,
                        seed = 43,
                        show_messages = FALSE,
                        check = FALSE) {
  stan_path <- system.file("stan", "Pooling.stan", package = "bglmmPool")
  if (!nzchar(stan_path)) {
    # dev fallback when running from source without install
    stan_path <- file.path("inst", "stan", "Pooling.stan")
  }
  if (!file.exists(stan_path)) stop("Can't find Stan file: ", stan_path)

  sm <- get_pooling_stan_model()

  fit <- rstan::sampling(
    object = sm,
    data = data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed,
    show_messages = show_messages,
    check = check
  )
  fit
}
