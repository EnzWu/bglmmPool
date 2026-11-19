library(rstan)
studies <- read.csv(file.choose())
studies <- read.csv("CaseData.csv")

is_blinded <- is.na(studies$NumCtrl) & is.na(studies$NumTrtm)
unblinded <- studies[!is_blinded, ]
blinded <- studies[is_blinded, ]
r_unblinded <- unblinded$AllocationRatio / (1 + unblinded$AllocationRatio)
n_trt <- round(unblinded$SampleSize * r_unblinded)
n_ctrl <- unblinded$SampleSize - n_trt

# Unblinded study data
data_u <- list(
  N_unblind = nrow(unblinded),
  N_trt = n_trt,
  N_ctrl = n_ctrl,
  y_trt = unblinded$NumTrtm,
  y_ctrl = unblinded$NumCtrl,
  fixed_prior_sd = 5,
  re_prior_scale = 1
)

# Combined study data
data_c <- c(data_u, list(
  N_blind = nrow(blinded),
  N_total = blinded$SampleSize,
  y_total = blinded$NumPooled,
  r = blinded$AllocationRatio / (1 + blinded$AllocationRatio)
))

# Stan models
unblinded_model <- "
// Bivariate logistic-normal model for the case study.
data {
  int<lower=1> N_unblind;
  array[N_unblind] int<lower=0> y_trt;
  array[N_unblind] int<lower=0> N_trt;
  array[N_unblind] int<lower=0> y_ctrl;
  array[N_unblind] int<lower=0> N_ctrl;
  real<lower=0> fixed_prior_sd;
  real<lower=0> re_prior_scale;
}
transformed data {
  int N_studies = N_unblind;
}
parameters {
  // Separate intercepts: beta0 for control, beta1 for treatment.
  real beta0;
  real beta1;
  vector<lower=0>[2] sigma;
  real z_rho;
  matrix[N_studies, 2] z;
}
transformed parameters {
  real rho = tanh(z_rho);
  matrix[2,2] Sigma;
  matrix[2,2] L_cov;
  matrix[N_studies,2] b;
  Sigma[1,1] = square(sigma[1]);
  Sigma[2,2] = square(sigma[2]);
  Sigma[1,2] = rho * sigma[1] * sigma[2];
  Sigma[2,1] = Sigma[1,2];
  L_cov = cholesky_decompose(Sigma);
  for (i in 1:N_studies)
    b[i,] = (L_cov * to_vector(z[i,]))';
}
model {
  beta0 ~ normal(0, fixed_prior_sd);
  beta1 ~ normal(0, fixed_prior_sd);
  sigma ~ normal(0, re_prior_scale);
  z_rho ~ normal(0, 1);
  to_vector(z) ~ normal(0, 1);
  for (i in 1:N_unblind) {
    y_ctrl[i] ~ binomial_logit(N_ctrl[i], beta0 + b[i,1]);
    y_trt[i] ~ binomial_logit(N_trt[i], beta1 + b[i,2]);
  }

}
generated quantities {
  real C = 16 * sqrt(3) / (15 * pi());
  real denom_ctrl = sqrt(1 + square(C * sigma[1]));
  real denom_trt = sqrt(1 + square(C * sigma[2]));
  real p_ctrl_overall = inv_logit(beta0 / denom_ctrl);
  real p_trt_overall = inv_logit(beta1 / denom_trt);
  real RD = p_trt_overall - p_ctrl_overall;
  real RR_overall = p_trt_overall / p_ctrl_overall;
  vector[N_studies] p_ctrl_i;
  vector[N_studies] p_trt_i;
  vector[N_studies] RR_i;
  for (i in 1:N_studies) {
    p_ctrl_i[i] = inv_logit(beta0 + b[i,1]);
    p_trt_i[i] = inv_logit(beta1 + b[i,2]);
    RR_i[i] = p_trt_i[i] / p_ctrl_i[i];
  }
}
"

combined_model <- "
// Bivariate logistic-normal model for the case study.
data {
  int<lower=1> N_unblind;
  int<lower=0> N_blind;
  array[N_unblind] int<lower=0> y_trt;
  array[N_unblind] int<lower=0> N_trt;
  array[N_unblind] int<lower=0> y_ctrl;
  array[N_unblind] int<lower=0> N_ctrl;
  array[N_blind] int<lower=0> y_total;
  array[N_blind] int<lower=0> N_total;
  array[N_blind] real<lower=0,upper=1> r;
  real<lower=0> fixed_prior_sd;
  real<lower=0> re_prior_scale;
}
transformed data {
  int N_studies = N_unblind + N_blind;
}
parameters {
  // Separate intercepts: beta0 for control, beta1 for treatment.
  real beta0;
  real beta1;
  vector<lower=0>[2] sigma;
  real z_rho;
  matrix[N_studies, 2] z;
}
transformed parameters {
  real rho = tanh(z_rho);
  matrix[2,2] Sigma;
  matrix[2,2] L_cov;
  matrix[N_studies,2] b;
  Sigma[1,1] = square(sigma[1]);
  Sigma[2,2] = square(sigma[2]);
  Sigma[1,2] = rho * sigma[1] * sigma[2];
  Sigma[2,1] = Sigma[1,2];
  L_cov = cholesky_decompose(Sigma);
  for (i in 1:N_studies)
    b[i,] = (L_cov * to_vector(z[i,]))';
}
model {
  beta0 ~ normal(0, fixed_prior_sd);
  beta1 ~ normal(0, fixed_prior_sd);
  sigma ~ normal(0, re_prior_scale);
  z_rho ~ normal(0, 1);
  to_vector(z) ~ normal(0, 1);
  for (i in 1:N_unblind) {
    y_ctrl[i] ~ binomial_logit(N_ctrl[i], beta0 + b[i,1]);
    y_trt[i] ~ binomial_logit(N_trt[i], beta1 + b[i,2]);
  }
  if (N_blind > 0) {
    for (j in 1:N_blind) {
      int i = N_unblind + j;
      real p_mix = r[j] * inv_logit(beta1 + b[i,2])
                   + (1-r[j]) * inv_logit(beta0 + b[i,1]);
      y_total[j] ~ binomial(N_total[j], p_mix);
    }
  }
}
generated quantities {
  real C = 16 * sqrt(3) / (15 * pi());
  real denom_ctrl = sqrt(1 + square(C * sigma[1]));
  real denom_trt = sqrt(1 + square(C * sigma[2]));
  real p_ctrl_overall = inv_logit(beta0 / denom_ctrl);
  real p_trt_overall = inv_logit(beta1 / denom_trt);
  real RD = p_trt_overall - p_ctrl_overall;
  real RR_overall = p_trt_overall / p_ctrl_overall;
  vector[N_studies] p_ctrl_i;
  vector[N_studies] p_trt_i;
  vector[N_studies] RR_i;
  for (i in 1:N_studies) {
    p_ctrl_i[i] = inv_logit(beta0 + b[i,1]);
    p_trt_i[i] = inv_logit(beta1 + b[i,2]);
    RR_i[i] = p_trt_i[i] / p_ctrl_i[i];
  }
}
"

# Fit both models
fit_u <- stan(model_code = unblinded_model, data = data_u, seed = 40001011,
              chains = 4, iter = 40000, warmup = 20000)
fit_c <- stan(model_code = combined_model, data = data_c, seed = 40001012,
              chains = 4, iter = 40000, warmup = 20000)

# Summarize the posterior draws and sampler diagnostics
fits <- list("Unblinded only" = fit_u, "Combined" = fit_c)
result_rows <- list()

for (model in names(fits)) {
  fit <- fits[[model]]
  draws <- rstan::extract(fit, pars = c("p_ctrl_overall", "p_trt_overall", "RD"))
  rd <- draws$RD
  interval <- quantile(rd, c(0.025, 0.975), names = FALSE)
  prob <- mean(rd > 0)
  
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- 0
  post_warmup <- 0
  
  for (chain in sampler) {
    divergences <- divergences + sum(chain[, "divergent__"])
    post_warmup <- post_warmup + nrow(chain)
  }
  
  result_rows[[model]] <- data.frame(
    model = model,
    control_risk = mean(draws$p_ctrl_overall),
    treatment_risk = mean(draws$p_trt_overall),
    RD_mean = mean(rd),
    RD_median = median(rd),
    RD_lower = interval[1],
    RD_upper = interval[2],
    Pr_RD_gt_0 = prob,
    exceed_080 = prob > 0.80,
    exceed_090 = prob > 0.90,
    exceed_095 = prob > 0.95,
    divergences = divergences,
    post_warmup = post_warmup
  )
}

results <- rbind(result_rows[["Unblinded only"]], result_rows[["Combined"]])
print(results, row.names = FALSE, digits = 4)
