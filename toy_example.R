set.seed(528)

library(parallel)
library(foreach)
library(doParallel)
library(spatstat.utils)

numCores = detectCores()

# simulation parameters
n_train = 100 
n_test = 100000
n_rep = 100
alphas = c(0.1, 0.15, 0.2)
n_alphas = length(alphas)
delta = 0.1
prob_content = function(lower, upper) {return(pnorm(upper) - pnorm(lower))}

# log pmf of discretized log normal
h.p = function(t, mu=6, sigma=1) {
  return(plnorm(t + 1, meanlog=mu, sdlog=sigma) - plnorm(t, meanlog=mu, sdlog=sigma))
}

# upper tail of discretized log normal
h.c = function(t, mu=6, sigma=1) {
  return(plnorm(t, lower.tail=FALSE, meanlog=mu, sdlog=sigma))
}

# entropy of Bernoulli random variable
psi = function(x, p) {
  return(p*log(p/x) + (1-p)*log((1-p)/(1-x)))
}

# computes probability content of split conformal interval at time t and coverage level 1-alpha
split_conformal = function(residuals, t, pred, alpha) {
  beta = (t+1)*(1-alpha)
  q_hat = ifelse(beta <= t,
                 orderstats(residuals, ceiling(beta)),
                 Inf)
  lower = pred - q_hat
  upper = pred + q_hat
  return(prob_content(lower, upper))
}

# computes order statistic for constructing split tuc interval at time t and coverage level 1-alpha
split_tuc = function(t, t_0, alpha) {
  beta = (t+1)*(1-alpha-4*(2*alpha-1)*log(h.p(t))/(3*(t+3))+
                  sqrt(-2*alpha*(1-alpha)*log(h.p(t))/(t+2))+
                  sqrt(pi*alpha*(1-alpha)/(2*(t+2)))*h.c(t_0+1))
  beta = ifelse(beta >= (t+1)*(1-alpha), ceiling(beta), Inf)
  return(beta)
}

# computes probability content of (alpha, delta)-split sequential inductive prediction interval at time t
split_si = function(residuals, t, pred, alpha, delta) {
  g_t = 0.85 * sqrt((log(log(exp(1) * t)) + 0.8 * log(1612/delta))/t)
  low_prob = alpha/2 - g_t
  high_prob = 1 - alpha/2 + g_t
  q_hat_low = ifelse(low_prob >= 0,
                     quantile(residuals, probs=low_prob, type=1),
                     -Inf)
  q_hat_high = ifelse(high_prob <= 1,
                      quantile(residuals, probs=high_prob, type=1),
                      Inf)
  lower = q_hat_low + pred
  upper = q_hat_high + pred
  return(prob_content(lower, upper))
}

split_tupac = function(t, t_0, alpha, delta, x_grid_size=10) {
  s = max(t_0, ceiling((t+1)*(1-alpha)))
  if (t < s) {
    return(Inf)
  } else {
    for (beta in s:t) {
      p = beta/(t+1)
      if (any(psi(seq(1-alpha, p, length.out=x_grid_size), p) >= 
              (log(2*h.c(t_0+1)/delta) - log(h.p(t)))/(t+1))) {
        return(beta)
      }
    }
  } 
  return(Inf)
}

start_time = Sys.time()

registerDoParallel(numCores)
res_rep = foreach (r=1:n_rep, 
                   .combine=rbind, 
                   .errorhandling='remove') %dopar% {
                     print(r)
                     z_seq = rnorm(n_train + n_test) # simulate standard normal data points
                     z_train = z_seq[1:n_train] # construct training dataset
                     z_test = z_seq[(n_train+1):(n_train+n_test)] # construct test dataset
                     pred = mean(z_train) # train model
                     tuc_t_0s = rep(1, n_alphas) # initialize t_0 for each alpha level for split TUC
                     tupac_t_0s = rep(1, n_alphas) # initialize t_0 for each alpha level for split TUPAC
                     residuals = c() 
                     abs_residuals = c()
                     split_conformal_probs = c()
                     split_si_probs = c()
                     split_tuc_probs = c()
                     split_tupac_probs = c()
                     for (t in 1:n_test) {
                       z = z_test[t]
                       residuals = c(residuals, z - pred) # update residuals
                       abs_residuals = c(abs_residuals, abs(z-pred)) # update absolute residuals
                       
                       # compute probability content of split conformal prediction intervals
                       split_conformal_prob = sapply(1:n_alphas,
                                                     function(a) {
                                                       return(split_conformal(residuals=abs_residuals,
                                                                              t=t,
                                                                              pred=pred,
                                                                              alpha=alphas[a]))})
                       split_conformal_probs = c(split_conformal_probs, split_conformal_prob)
                       
                       # compute probability content of split sequential inductive prediction intervals
                       split_si_prob = sapply(1:n_alphas,
                                              function(a) {
                                                return(split_si(residuals=residuals,
                                                                t=t,
                                                                pred=pred,
                                                                alpha=alphas[a],
                                                                delta=delta))
                                              })
                       split_si_probs = c(split_si_probs, split_si_prob)
                       
                       # compute probability content of split TUC prediction intervals
                       split_tuc_beta = sapply(1:n_alphas,
                                               function(a) {
                                                 return(split_tuc(t=t,
                                                                  t_0=tuc_t_0s[a],
                                                                  alpha=alphas[a]))
                                               })
                       split_tuc_q_hat = sapply(split_tuc_beta,
                                                function(beta) {return(ifelse(!is.na(beta) & 
                                                                                (beta <= t),
                                                                              orderstats(abs_residuals, beta),
                                                                              Inf))})
                       tuc_t_0s = ifelse(!is.na(split_tuc_beta) &
                                           (split_tuc_beta <= t),
                                         tuc_t_0s,
                                         t+1)
                       if (t > 1) { # correct probability content if t_0 is updated
                         for (a in 1:n_alphas) {
                           if (is.infinite(split_tuc_q_hat[a])) {
                             split_tuc_probs[n_alphas * c(0:(t-2)) + a ] = 1
                           }
                         }
                       }
                       split_tuc_prob = prob_content(pred - split_tuc_q_hat, pred + split_tuc_q_hat)
                       split_tuc_probs = c(split_tuc_probs, split_tuc_prob)
                       
                       # compute probability content of split TUPAC prediction intervals
                       split_tupac_beta = sapply(1:n_alphas,
                                                 function(a) {
                                                   return(split_tupac(t=t,
                                                                      t_0=tupac_t_0s[a],
                                                                      alpha=alphas[a],
                                                                      delta=delta))
                                               })
                       split_tupac_q_hat = sapply(split_tupac_beta,
                                                  function(beta) {return(ifelse(!is.na(beta) &
                                                                                  (beta <= t),
                                                                                orderstats(abs_residuals, beta),
                                                                                Inf))})
                       tupac_t_0s = ifelse(!is.na(split_tupac_beta) &
                                             (split_tupac_beta <= t),
                                           tupac_t_0s,
                                           t+1)
                       if (t > 1) { # correct probability content if t_0 is updated
                         for (a in 1:n_alphas) {
                           if (is.infinite(split_tupac_q_hat[a])) {
                             split_tupac_probs[n_alphas * c(0:(t-2)) + a ] = 1
                           }
                         }
                       }
                       split_tupac_prob = prob_content(pred - split_tupac_q_hat, pred + split_tupac_q_hat)
                       split_tupac_probs = c(split_tupac_probs, split_tupac_prob)
                     }
                     return(data.frame(r=rep(r, times=n_test*n_alphas),
                                       alpha=rep(alphas, n_test),
                                       t=rep(1:n_test, each=n_alphas),
                                       split_conformal_probs=split_conformal_probs,
                                       split_si_probs=split_si_probs,
                                       split_tuc_probs=split_tuc_probs,
                                       split_tupac_probs=split_tupac_probs))
                   }

stopImplicitCluster()

end_time = Sys.time()
print(end_time - start_time)

# save(res_rep, file="toy_example.Rda") # save data