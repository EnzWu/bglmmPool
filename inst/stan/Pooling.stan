data {
  // Dimensions of data
  int<lower=0> N_unblind;
  int<lower=0> N_blind;
  int<lower=0> N_studies;
  int<lower=0> y_trt[N_unblind];
  int<lower=0> N_trt[N_unblind];
  int<lower=0> y_ctrl[N_unblind];
  int<lower=0> N_ctrl[N_unblind];
  int<lower=0> y_total[N_blind];
  int<lower=0> N_total[N_blind];
  real<lower=0,upper=1> r[N_blind];
}
parameters {
  real beta0;
  real beta1;
  vector<lower=0>[2] sigma;
  real z_rho;
  matrix[N_studies, 2] z; // Non-centered random effect
}
transformed parameters {
  real rho = tanh(z_rho);

  matrix[2,2] Sigma;
  matrix[2,2] L_cov;
  matrix[N_studies, 2] b;

  Sigma[1,1] = square(sigma[1]);
  Sigma[2,2] = square(sigma[2]);
  Sigma[1,2] = rho * sigma[1] * sigma[2];
  Sigma[2,1] = Sigma[1,2];
  L_cov = cholesky_decompose(Sigma);

  for (i in 1:N_studies) {
    b[i,] = (L_cov * to_vector(z[i,]))'; // sample random effect using non-centered parameterization
  }
}
model {
  // Priors
  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  sigma ~ normal(0, 1);
  z_rho ~ normal(0, 1);
  to_vector(z) ~ normal(0, 1);

  // Likelihood for unblinded studies
  for (i in 1:N_unblind) {
    real logit_ctrl = beta0 + b[i, 1];
    real logit_trt = beta1 + b[i, 2];
    y_ctrl[i] ~ binomial_logit(N_ctrl[i], logit_ctrl);
    y_trt[i] ~ binomial_logit(N_trt[i],  logit_trt);
  }
  // Likelihood for blinded studies
  for (j in 1:N_blind) {
    int idx = N_unblind + j;
    real p_treat = inv_logit(beta1 + b[idx, 2]);
    real p_ctrl = inv_logit(beta0 + b[idx, 1]);
    real p_mix = r[j] * p_treat + (1 - r[j]) * p_ctrl;
    y_total[j] ~ binomial(N_total[j], p_mix);
  }
}
generated quantities {
  // Constant for logistic-normal mean approximation
  real C = 16*sqrt(3) / (15*pi());

  // Marginal probabilities integrating out RE variance
  real denom_ctrl = sqrt(1 + square(C*sigma[1]));
  real denom_trt= sqrt(1 + square(C*sigma[2]));

  real p_ctrl_overall = inv_logit(beta0 / denom_ctrl);
  real p_trt_overall = inv_logit(beta1 / denom_trt);

  // Overall effect measures
  real RR_overall_all = p_trt_overall / p_ctrl_overall;
  real logOR_overall_all = logit(p_trt_overall) - logit(p_ctrl_overall);
  real RD_overall_all = p_trt_overall - p_ctrl_overall;

  // Study-specific probabilities and RRs with random effects
  vector[N_studies] p_ctrl_i;
  vector[N_studies] p_trt_i;
  vector[N_studies] RD_i;

  for (i in 1:N_studies) {
    p_ctrl_i[i] = inv_logit(beta0 + b[i, 1]);
    p_trt_i[i] = inv_logit(beta1 + b[i, 2]);
    RD_i[i] = p_trt_i[i] - p_ctrl_i[i];
  }
}
