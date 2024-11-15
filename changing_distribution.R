set.seed(528)

library(ggplot2)
library(purrr)
library(tidyr)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)

numCores = detectCores()

# simulation parameters
n_test = 100000
reps = 100
eta = 2
alphas = c(0.1, 0.15, 0.2)
n_alphas = length(alphas)
deltas = c(0.1)
n_deltas = length(deltas)
pred_fun = mean
residual_fun = function(z_true, z_pred) {return(abs(z_true - z_pred))}
prob_content = function(lower, upper) {return(pnorm(upper) - pnorm(lower))}

# log pmf of discretized log normal
h.p_fun = function(t, mu=6, sigma=1) {
  return(plnorm(t + 1, meanlog=mu, sdlog=sigma) - plnorm(t, meanlog=mu, sdlog=sigma))
}

# upper tail of discretized log normal
h.c_fun = function(t, mu=6, sigma=1) {
  return(plnorm(t, lower.tail=FALSE, meanlog=mu, sdlog=sigma))
}

dynamic_skorski = function(t, k, t_k, alpha, delta, h.p, h.c, g.p) {
  curr_size = t-ceiling(eta^k)+1
  beta = (curr_size+1)*(1-alpha-(4*(2*alpha-1)*log(h.p(curr_size)))/(3*(curr_size+3))+
                          sqrt(-2*alpha*(1-alpha)*log(h.p(curr_size))/(curr_size+2))+
                          sqrt((2*alpha*(1-alpha)*log((h.c(t_k) - h.c(floor(eta^(k+2))))/g.p(k)))/
                                 (curr_size + 2))+sqrt((pi*alpha*(1-alpha))/(2*(curr_size+2))))
  return(ifelse(beta >= (curr_size + 1)*(1-alpha), ceiling(beta), Inf))
}

dynamic_skorski_fun = partial(dynamic_skorski, 
                              h.p=partial(h.p_fun, mu=11), 
                              h.c=partial(h.c_fun, mu=11), 
                              g.p=partial(h.p_fun, mu=log(16)))

start_time = Sys.time()

registerDoParallel(numCores)
res_rep = foreach (r=1:reps, 
                   .combine=rbind, 
                   .errorhandling='remove') %dopar% {
                     z_test = c(rnorm(n_test/2), rnorm(n_test/2, mean=1)) # simulate normal data points
                     pred_prev = 0
                     pred_curr = 0
                     residuals_prev = c()
                     residuals_curr = c()
                     t_0s_prev = rep(1, n_alphas * n_deltas)
                     t_0s_curr = rep(1, n_alphas * n_deltas)
                     k = 0
                     t_k_record = c()
                     betas_prime = c()
                     q_hats_prime = c()
                     preds_prime = c()
                     ks = c()
                     for(t in 1:n_test) {
                       if (t >= eta^(k+1)) {
                         k = k + 1
                         pred_prev = pred_curr
                         pred_curr = pred_fun(z_test[1:(t-1)])
                         t_0s_prev = t_0s_curr
                         t_0s_curr = rep(t, n_alphas * n_deltas)
                         residuals_prev = residuals_curr
                         residuals_curr = c()
                       }
                       t_k_record = c(t_k_record, t_0s_prev, t_0s_curr)
                       z = z_test[t]
                       residual_prev = residual_fun(z, pred_prev)
                       residuals_prev = sort(c(residuals_prev, residual_prev))
                       beta_prev = mapply(function(a, d) {
                         return(dynamic_skorski_fun(t=t, 
                                                    k=k-1,
                                                    t_k=t_0s_prev[n_deltas * ((a-1) %% 3) + d],
                                                    alpha=alphas[a], 
                                                    delta=deltas[d]))}, 
                         rep(1:n_alphas, each=n_deltas), 
                         rep(1:n_deltas, n_alphas))
                       q_hat_prev = ifelse(!is.na(beta_prev) & (beta_prev <= (t - ceiling(eta^(k-1)) + 1)),
                                           residuals_prev[beta_prev], Inf)
                       t_0s_prev = ifelse(!is.na(beta_prev) & (beta_prev <= (t - ceiling(eta^(k-1)) + 1)),
                                          t_0s_prev, t+1) 
                       residual_curr = residual_fun(z, pred_curr)
                       residuals_curr = sort(c(residuals_curr, residual_curr))
                       beta_curr = mapply(function(a, d) {
                         return(dynamic_skorski_fun(t=t, 
                                                    k=k,
                                                    t_k=t_0s_curr[n_deltas * ((a-1) %% 3) + d],
                                                    alpha=alphas[a], 
                                                    delta=deltas[d]))}, 
                         rep(1:n_alphas, each=n_deltas), 
                         rep(1:n_deltas, n_alphas))
                       q_hat_curr = ifelse(!is.na(beta_curr) & (beta_curr <= (t - ceiling(eta^k) + 1)),
                                           residuals_curr[beta_curr], Inf)
                       t_0s_curr = ifelse(!is.na(beta_curr) & (beta_curr <= (t - ceiling(eta^k) + 1)),
                                          t_0s_curr, t+1)
                       betas_prime = c(betas_prime, beta_prev, beta_curr)
                       q_hats_prime = c(q_hats_prime, q_hat_prev, q_hat_curr)
                       preds_prime = c(preds_prime, rep(pred_prev, n_alphas * n_deltas), rep(pred_curr, n_alphas * n_deltas))
                       ks = c(ks, k)
                     }
                     return(data.frame(r=rep(r, n_test*2*n_alphas*n_deltas), 
                                       t=rep(1:n_test, each=2*n_alphas*n_deltas),
                                       k=rep(ks, each=2*n_alphas*n_deltas),
                                       t_k=t_k_record,
                                       alpha=rep(rep(alphas, each=n_deltas), 2*n_test),
                                       delta=rep(deltas, n_alphas*2*n_test),
                                       pred_prime=preds_prime,
                                       beta_prime=betas_prime,
                                       q_hat_prime=q_hats_prime,
                                       beta_name_prime=rep(rep(c("Previous", "Current"), each=n_alphas*n_deltas), n_test)))
                   }
stopImplicitCluster()

end_time = Sys.time()
print(end_time - start_time)

res_rep = res_rep %>%
  left_join(res_rep %>%
              group_by(beta_name_prime, k, alpha, delta) %>%
              filter(is.infinite(q_hat_prime)) %>%
              summarise(t_max_prime=max(t))) %>%
  mutate(q_hat_prime_correct=ifelse(!is.na(t_max_prime) & (t <= t_max_prime), Inf, q_hat_prime))

save(res_rep, file="numerical_illustration.Rda")
