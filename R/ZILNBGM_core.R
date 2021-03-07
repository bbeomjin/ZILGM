# NB regression with l1 regularization for x with 1 columns
# glmreg_fit = mpath:::glmreg_fit
wlasso_nb = function(y, x, weights, penalty.factor = NULL, eta0 = NULL, mu0 = NULL, theta0 = NULL,
                     lambda, thresh = 1e-6, maxit = 100, n = NROW(x), p = NCOL(x))
{
  fun_call = match.call()
  obj_prev = 1e+150

  if (is.null(eta0)) {
    eta0 = rep(0, n)
  }

  if (is.null(mu0)) {
    mu0 = rep(1, n)
  }

  for (i in 1:maxit) {

    sig = max(n * weights * ((1 + 1 / theta0 * y) * mu0)/(1 + 1 / theta0 * mu0)^2)

    bobj = wlasso(X = x, y = eta0 + n * weights * ((y - mu0) / (1 + 1 / theta0 * mu0)) / sig,
                  eta0 = lambda / sig, wID = rep(1, n), weight = penalty.factor,
                  maxStep = maxit, eps = thresh, stand.scale = FALSE, trace = FALSE)
    bvec = bobj$coefficients
    eta = drop(bvec[1] + x %*% bvec[-1])
    mu = exp(eta)

    obj = nb_bvec_obj(y = y, weights = weights, bvec = bvec, mu = mu, lambda = lambda,
                      penalty.factor = penalty.factor, theta = theta0)
    if (is.nan(obj) | is.infinite(obj)) {obj = obj_prev}
    if (abs((obj_prev - obj) / obj_prev) < thresh) {
      bvec = bvec
      mu = mu
      eta = eta
      break
    } else if (obj > obj_prev + 1e-10) {
      bvec = bvec0
      mu = mu0
      eta = eta0
      break
    } else {
      obj_prev = obj
      bvec0 = bvec
      mu0 = mu
      eta0 = eta
    }
  }
  return(list(bvec = bvec, mu = mu, eta = eta, theta = theta0, iter = i))
}

glm_nb = function(y, x, weights, penalty.factor = NULL, eta0 = NULL, mu0 = NULL, theta0 = NULL,
                  lambda, thresh = 1e-6, maxit = 100, n = NROW(x), p = NCOL(x))
{
  bobj = glm.fit(x = cbind(1, x), y = y, family = negative.binomial(theta = theta0), intercept = TRUE, weights = n * weights,
                 control = list(epsilon = thresh, maxit = maxit), etastart = eta0, mustart = mu0)
  bvec = bobj$coefficients
  mu = bobj$fitted.values
  eta = log(mu)
  return(list(bvec = bvec, mu = mu, eta = eta, theta = theta0, iter = 0))
}


irls_nb = function(y, x, weights, penalty.factor = NULL, eta0 = NULL, mu0 = NULL, theta0 = NULL,
                   lambda, thresh = 1e-6, maxit = 100, n = NROW(x), p = NCOL(x))
{
  fun_call = match.call()
  obj_prev = 1e+150

  if (is.null(eta0)) {
    eta0 = rep(0, n)
  }

  if (is.null(mu0)) {
    mu0 = rep(1, n)
  }

  for (i in 1:maxit) {
    w = mu0 / (1 + mu0 / theta0)
    z = eta0 + (y - mu0) / mu0
     
    bobj = try((glmnet(x = x, y = z, family = "gaussian", weights = w * weights, lambda = lambda / sum(w * weights),
                  standardize = FALSE, alpha = 1, thresh = thresh, maxit = 10 * maxit, nlambda = 1, penalty.factor = penalty.factor)), silent = TRUE)
    
	if (inherits(bobj, "try-error")) {
	  bvec = rep(0, ncol(x) + 1)
	  mu = rep(1e-8, length(y))
      eta = log(mu)
	} else {
	  bvec = drop(coefficients(bobj))
      eta = drop(bvec[1] + x %*% bvec[-1])
      eta = ifelse(eta > log(1e+4), log(1e+4), eta)
      mu = exp(eta)
	}
    mu = ifelse(mu < 1e-6, 1e-6, mu)
    obj = nb_bvec_obj(y = y, weights = weights, bvec = bvec, mu = mu, lambda = lambda, penalty.factor = penalty.factor, theta = theta0)
    if (is.nan(obj) | is.infinite(obj)) {obj = obj_prev}
    if (abs((obj_prev - obj) / obj_prev) < thresh) {
      bvec = bvec
      mu = mu
      eta = eta
      break
      # } else if (obj > obj_prev + 1e-10) {
      #   bvec = bvec0
      #   mu = mu0
      #   eta = eta0
      #   break
    } else {
      obj_prev = obj
      bvec0 = bvec
      mu0 = mu
      eta0 = eta
    }
  }
  return(list(bvec = bvec, mu = mu, eta = eta, theta = theta0, iter = i))
}

# NB regression with l1 regularization using MM algorithm
pglm_nb_mm = function(y, x, weights, penalty.factor = NULL, bvec0 = NULL, eta0 = NULL, mu0 = NULL, theta0 = NULL,
                 lambda, thresh = 1e-6, maxit = 100, n = NROW(x), p = NCOL(x))
{
  fun_call = match.call()
  obj_prev = 1e+150

  if (is.null(eta0)) {
    eta0 = rep(0, n)
  }

  if (is.null(mu0)) {
    mu0 = rep(1, n)
  }

  for (i in 1:maxit) {

    sig = max(n * weights * ((1 + 1 / theta0 * y) * mu0)/(1 + 1 / theta0 * mu0)^2)
    bobj = glmnet(x = x, y = eta0 + n * weights * ((y - mu0) / (1 + 1 / theta0 * mu0)) / sig, family = "gaussian", alpha = 1,
                  lambda = lambda / sig, penalty.factor = penalty.factor, maxit = 10 * maxit, thresh = thresh, standardize = FALSE)
    bvec = drop(coefficients(bobj, s = lambda / sig))
    eta = drop(bvec[1] + x %*% bvec[-1])
    mu = exp(eta)

    obj = nb_bvec_obj(y = y, weights = weights, bvec = bvec, mu = mu, lambda = lambda,
                      penalty.factor = penalty.factor, theta = theta0)
    if (is.nan(obj) | is.infinite(obj)) {obj = obj_prev}
	
    if (abs((obj_prev - obj) / obj_prev) < thresh) {
      bvec = bvec
      mu = mu
      eta = eta
      break
    } else if (obj > obj_prev + 1e-10) {
      bvec = bvec0
      mu = mu0
      eta = eta0
      break
    } else {
      obj_prev = obj
      bvec0 = bvec
      mu0 = mu
      eta0 = eta
    }
  }
  return(list(bvec = bvec, mu = mu, eta = eta, theta = theta0, iter = i))
}


# NB regression with l1 regularization using IRLS algorithm
pglm_nb_irls = function(y, x, weights, theta0 = NULL, bvec0 = NULL, eta0 = NULL, mu0 = NULL,
                   lambda, penalty.factor = rep(1, NCOL(x)), thresh = 1e-6, maxit = 1e+3, n = NROW(x), p = NCOL(x))
{
  fun_call = match.call()
  negbin_fit = try((glmreg(y = y, x = x, weights = weights, lambda = lambda, alpha = 1, theta = theta0,
                           family = "negbin", thresh = thresh, maxit = maxit, penalty.factor = penalty.factor,
                           start = bvec0, mustart = mu0, etastart = eta0, standardize = FALSE, penalty = "enet",
                           x.keep = FALSE, y.keep = FALSE, trace = FALSE)), silent = TRUE)
  if (inherits(negbin_fit, "try-error")) {
    negbin_fit = try((irls_nb(y = y, x = x, weights = weights, lambda = lambda, theta0 = theta0,
                          thresh = thresh, maxit = maxit, penalty.factor = penalty.factor, eta0 = eta0, mu0 = mu0)), silent = TRUE)
    if (inherits(negbin_fit, "try-error")) {
      bvec = rep(0, ncol(x) + 1)
      mu = rep(1e-8, length(y))
      eta = log(mu)
    } else {
      bvec = negbin_fit$bvec
      eta = negbin_fit$eta
      mu = negbin_fit$mu
    }
  } else {
    bvec = drop(c(negbin_fit$b0, negbin_fit$beta))
    mu = negbin_fit$fitted.values
    eta = log(mu)
  }

  return(list(bvec = bvec, mu = mu, eta = eta, theta = theta0))
}



zilgm_negbin = function(y, x, lambda, weights = NULL, update_type = c("IRLS", "MM"), penalty.factor = NULL,
                        thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2, init_theta = NULL)
{
  update_type = match.arg(update_type)
  fun_call = match.call()
  out = list()

  n = NROW(x)
  p = NCOL(x)

  if ((p == 1) & (update_type == "MM")) {update_type = "onecol_MM"}
  if ((p == 1) & (update_type == "IRLS")) {update_type = "onecol_IRLS"}

  update_fun = switch(update_type,
                      onecol_MM = wlasso_nb,
                      onecol_irls = glm_nb,
                      MM = pglm_nb_mm,
                      IRLS = pglm_nb_irls)

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
    param = list(bvec = rep(0, p + 1), theta = 1e+8, prob = 0, pos_zero = which(pos_zero), iter = 0)
    return(param)
  }

  weights = weights / sum(weights)

  mu0 = rep(mean(y[y > 0]), n)
  eta0 = log(mu0)
  bvec0 = c(eta0[1], rep(0, p))

  if (is.null(init_theta)) {
    theta0 = 1e+8
  } else {
    theta0 = init_theta
  }
  
  theta0 = 1e+8
  prob0 = (sum(pos_zero) - sum(dNBI(0, mu = mu0, theta = theta0, log = FALSE))) / n
  prob0 = ifelse(prob0 < 1e-10, 1e-10, ifelse(prob0 > 1, 1, prob0))

  erisk_prev = 1e+150

  if (sum(pos_zero) == 0) {

    for (iter in 1:EM_iter) {

      sol_bvec = update_fun(y = y, x = x, weights = weights, penalty.factor = penalty.factor,
                            bvec0 = bvec0, eta0 = eta0, mu0 = mu0, lambda = lambda, theta0 = theta0,
                            thresh = tol, maxit = maxit, n = n, p = p)
      bvec = sol_bvec$bvec
      eta = sol_bvec$eta
      mu = sol_bvec$mu
      theta = sol_bvec$theta

      if (!is.null(init_theta)) {
        theta = init_theta
      } else {
        theta = theta_ml(y = y, mu = mu, weights = weights)
      }

      erisk = nb_objective(y = y, prob = prob, bvec = bvec, mu = mu, lambda = lambda,
                           weights = weights, penalty.factor = penalty.factor, theta = theta, posz = pos_zero)
      if (is.infinite(erisk) | is.nan(erisk)) {erisk = erisk_prev}
      if ((abs((erisk_prev - erisk) / (erisk_prev + 1)) < EM_tol)) {
        bvec = bvec
        theta = theta
        prob = 0
        break
      } else if (erisk > erisk_prev + 1e-10) {
        bvec = bvec0
        theta = theta
        prob = 0
        break
      } else {
        erisk_prev = erisk
        bvec0 = bvec
        eta0 = eta
        mu0 = mu
        theta0 = theta
        prob = 0
      }
    }
  } else {

    for (iter in 1:EM_iter) {

      # E-step
      tmp_z = prob0 / (prob0 + (1 - prob0) * dNBI(0, theta = theta0, mu = mu0, log = FALSE))
      tmp_z[is.nan(tmp_z)] = 1
      tmp_z = ifelse(tmp_z >= (1 - 1e-6), 1 - 1e-6, tmp_z)
      z[pos_zero] = tmp_z[pos_zero]

      prob = sum(z) / n
      prob = ifelse(prob < 1e-10, 1e-10, ifelse(prob > 1, 1, prob))

      # M-step
      sol_bvec = update_fun(y = y, x = x, weights = weights * (1 - z), penalty.factor = penalty.factor,
                            bvec0 = bvec0, eta0 = eta0, mu0 = mu0, lambda = lambda, theta0 = theta0,
                            thresh = tol, maxit = maxit, n = n, p = p)

      bvec = sol_bvec$bvec
      eta = sol_bvec$eta
      mu = sol_bvec$mu

      if (!is.null(init_theta)) {
        theta = init_theta
      } else {
        theta = theta_ml(y = y, mu = mu0, weights = weights * (1 - z))
      }

      erisk = nb_objective(y = y, prob = prob, bvec = bvec, mu = mu, lambda = lambda,
                           weights = weights, penalty.factor = penalty.factor, theta = theta, posz = pos_zero)
	  if (is.infinite(erisk) | is.nan(erisk)) {erisk = erisk_prev}
	  
      if ((abs((erisk_prev - erisk) / (erisk_prev + 1)) < EM_tol)) {
        bvec = bvec
        theta = theta
        prob = prob
        break
      # } else if (erisk > erisk_prev + 1e-10) {
      #   bvec = bvec0
      #   theta = theta
      #   prob = prob0
      #   break
      } else {
        erisk_prev = erisk
        bvec0 = bvec
        eta0 = eta
        mu0 = mu
        theta0 = theta
        prob0 = prob
      }
    }
  }
  flag = abs(bvec) < thresh
  bvec[flag] = 0

  out$bvec = bvec
  out$theta = theta
  out$prob = prob
  out$pos_zero = which(pos_zero)
  out$iterations = iter
  out$loglik = erisk
  out$call = fun_call
  class(out) = "zilgm"
  return(out)
}
