
zilgm_negbin2 = function(y, x, lambda, weights = NULL, update_type = c("IRLS", "MM"), penalty.factor = NULL,
                        fixed_sigma = FALSE, tol = 1e-6, EM_tol = 1e-6, EM_iter = 300, thresh = 1e-6, maxit = 1e+3)
{
  update_type = match.arg(update_type)
  fun_call = match.call()
  
  if ((p == 1) & (update_type == "MM")) {update_type = "onecol_MM"}
  if ((p == 1) & (update_type == "IRLS")) {update_type = "onecol_IRLS"}
  
  out = list()
  update_fun = switch(update_type,
                      onecol_MM = wlasso_p,
                      onecol_irls = glm_p,
                      MM = pglm_p_mm,
                      IRLS = pglm_p_irls)
  
  n = NROW(x)
  p = NCOL(x)
  pos_zero = (y == 0)
  pos_nzero = !pos_zero
  z = rep(1e-6, n)
  
  if (is.null(penalty.factor)) {
    penalty.factor = rep(1, p)
  }
  
  if (is.null(weights)) {
    weights = rep(1, n)
  }
  
  if (length(unique(y)) == 1) {
    param = list(bvec = rep(0, p + 1), sigma = 0, prob = 0, pos_zero = which(pos_zero), iter = 0)
    return(param)
  }
  
  weights = weights / sum(weights)
  
  mu0 = rep(mean(y[y > 0]), n)
  eta0 = log(mu0)
  bvec0 = c(eta0[1], rep(0, p))
  # sigma0 = sigma_ml(y = y, mu = mu0)
  sigma0 = 1e-4
  prob0 = (sum(pos_zero) - sum(dNBII(0, mu = mu0, sigma = sigma0, log = FALSE))) / n
  prob0 = ifelse(prob0 < 1e-10, 1e-10, ifelse(prob0 > 1, 1, prob0))
  
  erisk_prev = 1e+150
  
  if (sum(pos_zero) == 0) {
    sol_bvec = update_fun(y = y, x = x, weights = weights, penalty.factor = penalty.factor,
                               bvec0 = bvec0, eta0 = eta0, mu0 = mu0, lambda = lambda,
                               thresh = tol, maxit = maxit, n = n, p = p)
    bvec = sol_bvec$bvec
    eta = sol_bvec$eta
    mu = sol_bvec$mu
    
    prob = prob0
    iter = 0
    erisk = erisk_prev
    sigma = sigma0
  } else {
    
    for (iter in 1:EM_iter) {
      
      # E-step
      tmp_z = prob0 / (prob0 + (1 - prob0) * dNBII(0, sigma = sigma0, mu = mu0, log = FALSE))
      tmp_z[is.nan(tmp_z)] = 1
      tmp_z = ifelse(tmp_z >= (1 - 1e-6), 1 - 1e-6, tmp_z)
      z[pos_zero] = tmp_z[pos_zero]
      
      prob = sum(z) / n
      prob = ifelse(prob < 1e-10, 1e-10, ifelse(prob > 1, 1, prob))
      
      # M-step
      sol_bvec = update_fun(y = y, x = x, weights = weights * (1 - z), penalty.factor = penalty.factor,
                            bvec0 = bvec0, eta0 = eta0, mu0 = mu0, lambda = lambda, thresh = tol, 
                            maxit = maxit, n = n, p = p)
      
      bvec = sol_bvec$bvec
      eta = sol_bvec$eta
      mu = sol_bvec$mu
      
      if (fixed_sigma) {
        sigma = sigma0
      } else {
        sigma = sigma_ml(y, mu = mu, weights = weights * (1 - z))
      }
      
      erisk = nb2_objective(y = y, prob = prob, bvec = bvec, mu = mu, lambda = lambda, 
                            weights = weights, penalty.factor = penalty.factor, sigma = sigma, posz = pos_zero)
      
      if ((abs((erisk_prev - erisk) / (erisk_prev + 1)) < EM_tol)) {
        bvec = bvec
        sigma = sigma
        prob = prob
        z = z
        break
        # } else if (erisk > erisk_prev + 1e-10) {
        #   bvec = bvec0
        #   sigma = sigma
        #   prob = prob0
        #   break
      } else {
        erisk_prev = erisk
        bvec0 = bvec
        eta0 = eta
        mu0 = mu
        sigma0 = sigma
        prob0 = prob
      }
    }
  }
  flag = abs(bvec) < thresh
  bvec[flag] = 0
  
  out$bvec = bvec
  out$sigma = sigma
  out$prob = prob
  out$pos_zero = which(pos_zero)
  out$iterations = iter
  out$loglik = erisk
  out$call = fun_call
  
  return(out)
}
