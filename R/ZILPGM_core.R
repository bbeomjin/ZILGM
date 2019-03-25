# Poisson regression with l1 regularization for x with 1 columns using MM algorithm
wlasso_p = function(y, x, weights, penalty.factor = NULL, eta0 = NULL, mu0 = NULL,
                    lambda, thresh = 1e-7, maxit = 100, n = NROW(x), p = NCOL(x))
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

    sig = max(n * weights * mu0)

    bobj = wlasso(X = x, y = eta0 + n * weights * (y - mu0) / sig,
                  eta0 = lambda / sig, wID = rep(1, n), weight = penalty.factor,
                  maxStep = maxit, eps = thresh, stand.scale = FALSE, trace = FALSE)
    bvec = bobj$coefficients
    eta = drop(bvec[1] + x %*% bvec[-1])
    mu = exp(eta)

    obj = p_bvec_obj(y = y, weights = weights, bvec = bvec, mu = mu, lambda = lambda, penalty.factor = penalty.factor)

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
  return(list(bvec = bvec, mu = mu, eta = eta, iter = i))
}


# Poisson regression with l1 regularization for x with 1 columns using IRLS algorithm
glm_p = function(y, x, weights, penalty.factor = NULL, eta0 = NULL, mu0 = NULL,
                 lambda, thresh = 1e-7, maxit = 100, n = NROW(x), p = NCOL(x))
{
  bobj = glm.fit(x = cbind(1, x), y = y, family = "poisson", intercept = TRUE, weights = n * weights,
                 control = list(epsilon = thresh, maxit = maxit), etastart = eta0, mustart = mu0)
  bvec = bobj$coefficients
  mu = bobj$fitted.values
  eta = log(mu)
  return(list(bvec = bvec, mu = mu, eta = eta, iter = 0))
}


irls_p = function(y, x, weights, penalty.factor = NULL, eta0 = NULL, mu0 = NULL,
                  lambda, thresh = 1e-7, maxit = 100, n = NROW(x), p = NCOL(x))
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
    w = mu0
    z = eta0 + (y - mu0) / mu0

    bobj = glmnet(x = x, y = z, family = "gaussian", weights = w * weights, lambda = lambda / sum(w * weights),
                  standardize = FALSE, alpha = 1, thresh = thresh, maxit = 10 * maxit, nlambda = 1)

    bvec = drop(coefficients(bobj))
    eta = drop(bvec[1] + x %*% bvec[-1])
    eta = ifelse(eta > log(1e+4), log(1e+4), eta)
    mu = exp(eta)

    obj = p_bvec_obj(y = y, weights = weights, bvec = bvec, mu = mu, lambda = lambda, penalty.factor = penalty.factor)

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
  return(list(bvec = bvec, mu = mu, eta = eta, iter = i))
}


pglm_p_mm = function(y, x, weights, penalty.factor = NULL, bvec0 = NULL, eta0 = NULL, mu0 = NULL,
                      lambda, thresh = 1e-7, maxit = 100, n = NROW(x), p = NCOL(x))
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

    sig = max(n * weights * mu0)
    bobj = glmnet(x = x, y = eta0 + n * weights * (y - mu0) / sig, family = "gaussian", alpha = 1,
                  lambda = lambda / sig, penalty.factor = penalty.factor, maxit = 10 * maxit, thresh = thresh, standardize = FALSE)
    bvec = drop(coefficients(bobj, s = lambda / sig))
    eta = drop(bvec[1] + x %*% bvec[-1])
    mu = exp(eta)

    obj = p_bvec_obj(y = y, weights = weights, bvec = bvec, mu = mu, lambda = lambda,
                      penalty.factor = penalty.factor)

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
  return(list(bvec = bvec, mu = mu, eta = eta, iter = i))
}

pglm_p_irls = function(y, x, weights, bvec0 = NULL, eta0 = NULL, mu0 = NULL,
                        lambda, penalty.factor = NULL, thresh = 1e-7, maxit = 1e+5, n = NROW(x), p = NCOL(x))
{
  fun_call = match.call()
  poisson_fit = try((penalized_glm(y = y, x = x, weights = weights, lambda = lambda, alpha = 1, family = "poisson",
                                   thresh = thresh, maxit = maxit, penalty.factor = penalty.factor,
                                   start = bvec0, mustart = mu0, etastart = eta0)), silent = FALSE)
  if (inherits(poisson_fit, "try-error")) {
    poisson_fit = irls_p(y = y, x = x, weights = weights, lambda = lambda, thresh = thresh, maxit = maxit,
                          penalty.factor = penalty.factor, eta0 = eta0, mu0 = mu0)
    bvec = poisson_fit$bvec
    eta = poisson_fit$eta
    mu = poisson_fit$mu
  } else {
    bvec = drop(c(poisson_fit$b0, poisson_fit$beta))
    mu = poisson_fit$fitted.values
    eta = log(mu)
  }

  return(list(bvec = bvec, mu = mu, eta = eta))
}

zilgm_poisson = function(y, x, lambda, weights = NULL, update_type = c("IRLS", "MM"), penalty.factor = NULL,
                         thresh = 1e-6, EM_tol = 1e-6, EM_iter = 500, tol = 1e-6, maxit = 1e+3)
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
  z0 = z = rep(1e-6, n)

  if (is.null(penalty.factor)) {
    penalty.factor = rep(1, p)
  }

  if (is.null(weights)) {
    weights = rep(1, n)
  }

  if (length(unique(y)) == 1) {
    param = list(bvec = rep(0, p + 1), prob = 0, pos_zero = which(pos_zero), iter = 0)
    return(param)
  }

  weights = weights / sum(weights)

  # mu0 = rep(mean(y[y > 0]), n)
  mu0 = rep(mean(y), n)
  eta0 = log(mu0)
  bvec0 = c(eta0[1], rep(0, p))

  prob0 = (sum(pos_zero) - sum(exp(-mu0))) / n
  prob0 = ifelse(prob0 < 1e-10, 1e-10, ifelse(prob0 > 1, 1, prob0))

  erisk_prev = 1e+150

  if (sum(pos_zero) == 0) {
    sol_bvec = update_fun(y = y, x = x, weights = weights, bvec0 = bvec0, eta0 = eta0, mu0 = mu0, lambda = lambda,
                          penalty.factor = penalty.factor, thresh = tol, maxit = maxit, n = n, p = p)
    bvec = sol_bvec$bvec
    eta = sol_bvec$eta
    mu = sol_bvec$mu

    prob = prob0
    iter = 0
    erisk = 1e+150

  } else {

    for (iter in 1:EM_iter) {

      # E-step
      tmp_z = prob0 / (prob0 + (1 - prob0) * dP(0, mu0, log = FALSE))
      tmp_z[is.nan(tmp_z)] = 1
      tmp_z = ifelse(tmp_z >= (1 - 1e-6), 1 - 1e-6, tmp_z)
      z[pos_zero] = tmp_z[pos_zero]

      prob = sum(z) / n
      prob = ifelse(prob < 1e-10, 1e-10, ifelse(prob > 1, 1, prob))

      # M-step
      sol_bvec = update_fun(y = y, x = x, weights = weights * (1 - z), bvec0 = bvec0, eta0 = eta0, mu0 = mu0, lambda = lambda,
                            penalty.factor = penalty.factor, thresh = tol, maxit = maxit, n = n, p = p)
      bvec = sol_bvec$bvec
      eta = sol_bvec$eta
      mu = sol_bvec$mu

      erisk = p_objective(y = y, weights = weights, prob = prob, bvec = bvec, mu = mu, lambda = lambda,
                          penalty.factor = penalty.factor, posz = pos_zero)


      if ((abs((erisk_prev - erisk) / (erisk_prev + 1)) < EM_tol)) {
        bvec = bvec
        erisk = erisk
        prob = prob
        z = z
        break
        # } else if (erisk > erisk_prev + 1e-10) {
        #   bvec = bvec0
        #   erisk = erisk_prev
        #   prob = prob0
        #   z = z0
        #   break
      } else {
        erisk_prev = erisk
        bvec0 = bvec
        eta0 = eta
        mu0 = mu
        prob0 = prob
        z0 = z
      }
    }
  }
  flag = abs(bvec) < thresh
  bvec[flag] = 0

  out$bvec = bvec
  out$prob = prob
  out$pos_zero = which(pos_zero)
  out$iterations = iter
  out$loglik = erisk
  out$z = z
  out$call = fun_call
  class(out) = "zilgm"
  return(out)
}
