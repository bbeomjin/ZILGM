zipoisson = function(x, y, par, n) {
  pi = par[1]
  bvec = par[-1]
  mu = exp(x %*% bvec)
  obj = -1 / n * sum(log(pi * (y == 0) + (1 - pi) * (exp(-mu) * mu^y) / gamma(y + 1)))
  return(obj)
}

zipoisson_grad = function(x, y, par, n) {
  pi = par[1]
  bvec = par[-1]
  mu = exp(x %*% bvec)
  grad1 = -1 / n * (sum(((y == 0) - (exp(-mu) * mu^y) / gamma(y + 1)) / (pi * (y == 0) + (1 - pi) * (exp(-mu) * mu^y) / gamma(y + 1))))
  grad2 = -1 / n * (t(x) %*% ((1 - pi) * exp(-mu) * (y * mu^(y - 1) - mu^y) * mu / gamma(y + 1) / (pi * (y == 0) + (1 - pi) * (exp(-mu) * mu^y) / gamma(y + 1))))
  grad = c(grad1, grad2)
  return(grad)
}


zinegbin = function(x, y, par, theta, n) {
  lambda = par[1]
  bvec = par[-1]
  mu = exp(x %*% bvec)
  pi = exp(lambda) / (1 + exp(lambda))
  obj = -1 / n * sum(log(pi * (y == 0) + (1 - pi) * dnbinom(x = y, size = theta, mu = mu, log = FALSE)))
  return(obj)
}

zinegbin_grad = function(x, y, par, theta, n) {
  lambda = par[1]
  bvec = par[-1]
  mu = exp(x %*% bvec)
  pi = exp(lambda) / (1 + exp(lambda))
  grad1 = 1 / n * crossprod(x[y == 0, , drop = FALSE], ((mu * (1 + mu / theta)^{-1 -theta}) / (pi / (1 - pi) + (1 + mu / theta)^{-theta}))[y == 0])
  grad2 = -1 / n * crossprod(x[y != 0, , drop = FALSE], ((y - mu) / (1 + mu / theta))[y != 0])
  grad_bvec = grad1 + grad2
  grad_pi = -1 / n * (sum(((y == 0) - dnbinom(x = y, size = theta, mu = mu)) / (pi * (y == 0) + (1 - pi) * dnbinom(x = y, size = theta, mu = mu))))
  grad = c(grad_pi, grad_bvec)
  return(grad)
}




zigm_coord_poisson2 = function(y, x, lambda, alpha = 1, weight = NULL, update_type = c("IRLS", "MM", "lbfgs"),
                               thres = 1e-6, conv_eps = 1e-6, MaxIter = 300, inner_eps = 1e-7, inner_iter = 1e+5,
                               fixed_theta = FALSE)
{
  x1 = cbind(1, x)
  p = ncol(x1)
  system.time((poisson_fit = lbfgs(zipoisson, zipoisson_grad, x = x1, y = y, n = n, vars = rep(0, p + 1), invisible = 1,
                                   orthantwise_c = lambda, linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING",
                                   orthantwise_start = 2, orthantwise_end = p + 1)))
  poisson_fit
}


zigm_coord_negbin_inner = function(y, x, lambda, theta0, weight = NULL, bvec0, prob0) {
  x1 = cbind(1, x)
  p = ncol(x1)
  negbin_fit = lbfgs(zinegbin, zinegbin_grad, x = x1, y = y, n = n, vars = rep(0, p + 1), invisible = 1, theta = theta0,
                     orthantwise_c = lambda, linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING",
                     orthantwise_start = 2, orthantwise_end = p + 1)
  lambda = negbin_fit$par[1]
  prob = exp(lambda) / (1 + exp(lambda))
  bvec = negbin_fit$par[-1]
  mu = exp(x1 %*% bvec)
  eta = log(mu)
  
  return(list(bvec = bvec, prob = prob, mu = mu, eta = eta))
}
