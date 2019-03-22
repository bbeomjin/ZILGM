library(flux)

dNBI = function(y, mu, theta, log = FALSE)
{
  density = lgamma(y + theta + 1e-10) - lgamma(y + 1) - lgamma(theta + 1e-10) + theta * (log(theta + 1e-10) - log(theta + mu + 1e-10)) + y * (log(mu + 1e-10) - log(theta + mu + 1e-10))
  if (log == FALSE) {density = exp(density)}
  return(density)
}

dP = function(y, mu, log = FALSE)
{
  density = -mu + y * log(mu + 1e-10) - lgamma(y + 1)
  if (log == FALSE) {density = exp(density)}
  return(density)
}

dNBII = function(y, mu, sigma, log = FALSE) 
{
  density = dNBI(y, mu = mu, theta = mu / sigma, log = log)
  return(density)
}


Adjacency = function (mat)
{
  # absolute matrix
  diag(mat) = 0
  idx = which(abs(mat) != 0, arr.ind = TRUE)
  mat[idx] = 1
  return(mat)
}

find_lammax = function (X)
{
  tmp = t(X) %*% X
  lammax = 1/nrow(X) * max(abs(tmp[upper.tri(tmp)]))
  return(lammax)
}

cal_ebic_inflation = function (results, X, gamval = 1)
{
  p = ncol(X) # p by length of lambda
  n = nrow(X)  # n by length of lambda
  nlam = ncol(results[[1]]$Bmat)
  risk = rep(0, nlam); df = rep(0, nlam)
  
  for (k in 1:nlam) {
    ll = matrix(0, n, p)
    dfj = 0
    
    for (j in 1:p) {
      eta = results[[j]]$b0[k] + drop(X %*% results[[j]]$Bmat[, k])
      mu = exp(eta)
      prob0 = results[[j]]$prob0[k]
      dfj = dfj + sum(results[[j]]$Bmat[,k] != 0)
      
      flag0 = X[,j] == 0
      ll[,j] = (1 - prob0)*exp(-mu)*mu^X[, j]/factorial(X[, j])
      ll[flag0, j] = ll[flag0, j] + prob0
    }
    risk[k] = sum(log(ll))
    df[k] = dfj
  }
  risk = -2 * risk
  npair = p * (p-1)/2
  
  #npair <- 1
  bic = risk + log(n)*(df)
  ebic = risk + (log(n)+2*gamval*log(npair))*df
  
  return(list(bic = bic, ebic = ebic, risk = risk, df = df))
}

cal_ebic = function (results, X, gamval = 1)
{
  p = ncol(X) # p by length of lambda
  n = nrow(X)  # n by length of lambda
  nlam = ncol(results[[1]]$Bmat)
  risk = rep(0, nlam); df = rep(0, nlam)
  for (k in 1:nlam) {
    ll = matrix(0, n, p)
    dfj = 0
    
    for (j in 1:p) {
      eta = results[[j]]$b0[k] + drop(X %*% results[[j]]$Bmat[, k])
      mu = exp(eta)
      prob0 = results[[j]]$prob0[k]
      dfj = dfj + sum(results[[j]]$Bmat[, k] != 0)
      ll[,j] = exp(-mu) * mu^X[, j]/factorial(X[, j])
    }
    risk[k] = sum(log(ll))
    df[k] = dfj
  }
  risk = -2 * risk
  npair = p * (p - 1)/2
  #npair <- 1
  bic = risk + log(n) * (df)
  ebic = risk + (log(n) + 2 * gamval * log(npair)) * df
  
  return(list(bic = bic, ebic = ebic, risk = risk, df = df))
}

get_A = function (results, pos)
{
  p = length(results)
  A = matrix(0, p, p)
  
  for (j in 1:p) {
    flag = results[[j]]$Bmat[,pos] != 0
    A[flag, j] = 1 
  }
  return(A)
}  

to_upper = function (tX) {tX[lower.tri(tX, diag = F)]}

vec2tri = function (k, p)
{
  i = ceiling(0.5 * (2 * p - 1 - sqrt((2 * p - 1)^2 - 8 * k)))
  j = k - p * (i - 1) + i * (i - 1)/2 + i
  return(as.matrix(cbind(i, j)))
}

tri2vec = function(i, j, p){
  return(p * (i - 1) - i * (i - 1)/2 + j - i)
}


p_bvec_obj = function(y, weights, bvec, mu, lambda, penalty.factor)
{
  penalty = lambda * sum(abs(penalty.factor * bvec[-1]))
  pnl = -sum(weights * dP(y = y, mu = mu, log = FALSE) + 1e-10)
  return(pnl + penalty)
}

nb_bvec_obj = function(y, weights, bvec, mu, theta = NULL, lambda, penalty.factor)
{
  penalty = lambda * sum(abs(penalty.factor * bvec[-1]))
  pnl = -sum(weights * log(dNBI(y = y, theta = theta, mu = mu, log = FALSE) + 1e-10))
  return(pnl + penalty)
}


p_objective = function(y, weights, prob, bvec, mu, lambda, penalty.factor, posz)
{
  penalty = lambda * sum(abs(penalty.factor * bvec[-1]))
  pnl = -sum(weights * log(prob * posz + (1 - prob) * dP(y = y, mu = mu, log = FALSE) + 1e-10))
  return(pnl + penalty)
}  

nb_objective = function(y, weights, prob, bvec, mu, theta = NULL, lambda, penalty.factor, posz)
{
  penalty = lambda * sum(abs(penalty.factor * bvec[-1]))
  pnl = -sum(weights * log(prob * posz + (1 - prob) * dNBI(y = y, theta = theta, mu = mu, log = FALSE) + 1e-10))
  return(pnl + penalty)
}

nb2_objective = function(y, weights, prob, bvec, mu, sigma = NULL, lambda, penalty.factor, posz)
{
  penalty = lambda * sum(abs(penalty.factor * bvec[-1]))
  pnl = -sum(weights * log(prob * posz + (1 - prob) * dNBII(y = y, sigma = sigma, mu = mu, log = FALSE) + 1e-10))
  return(pnl + penalty)
}

# thresholding matrix
thresholding_mat = function(mat, thres = 0.1)
{
  #Bmat <- matrix(0, nrow(mat), ncol(mat))
  Bmat = Matrix(0, nrow(mat), ncol(mat), sparse = TRUE)
  for (i in 1:nrow(mat)) {
    flag = abs(mat[i, ]) > thres
    Bmat[i, flag] = mat[i, flag]
  }
  return(Bmat)
}

f_K_fold = function(Nobs, K = 5)
{
  rs = runif(Nobs)
  id = seq(Nobs)[order(rs)]
  k = as.integer(Nobs * seq(1, K - 1)/K)
  k = matrix(c(0, rep(k, each = 2), Nobs), ncol = 2, byrow = TRUE)
  k[, 1] = k[, 1] + 1
  l = lapply(seq.int(K), function(x, k, d) {list(train = d[!(seq(d) %in% seq(k[x, 1],k[x, 2]))],
                                                 test = d[seq(k[x, 1], k[x, 2])])},
             k = k, d = id)
  return(l)
}
#################### modifying
split_to_subset = function(n, folds = 5)
{
  size = floor(n/folds)
  tsize = 0
  pset = vector(mode = "list", length = folds)
  for (k in 1:folds) {
    tsize = tsize + size
    pset[[k]] = (size * (k - 1) + 1):(size * k)
  }
  if (tsize < n) {
    pset[[folds]] = ((size * (folds - 1)) + 1):n
  }
  return(pset)
}


hat_net = function(coef_mat, thresh = 1e-6, type = c("AND", "OR"))
{
  type = match.arg(type)
  
  tmp_mat = abs(coef_mat) > thresh
  
  if (type == "AND") {
    res_mat = tmp_mat * t(tmp_mat)
  }
  
  if (type == "OR") {
    res_mat = (tmp_mat + t(tmp_mat) > 0) * 1
  }
  return(res_mat)
}


theta_ml = function(y, mu, weights = NULL) {
  n = length(y)
  if (is.null(weights)) {weights = rep(1, n)}
  nb_theta = function(theta, mu, y, weights) {
    return(sum(weights * dNBI(y = y, theta = theta, mu = mu, log = TRUE)))
  }
  fit = optimize(nb_theta, y = y, mu = mu, weights = weights, interval = c(1e-4, 5e+3), maximum = TRUE)
  theta = ifelse(fit$maximum > 1e+3, 1e+8, fit$maximum)
  return(theta)
}


theta_ml2 = function(y, mu, posz, prob0, bvec, lambda, weights, penalty.factor) {
  n = length(y)
  fit = optimize(nb_objective, y = y, posz = posz, bvec = bvec, lambda = lambda, penalty.factor = penalty.factor,
                 prob = prob0, mu = mu, weights = weights, interval = c(1e-4, 1e+4), maximum = FALSE)
  theta = fit$minimum
  return(theta)
}


theta_mm = function(y, mu, weights, nparam = NULL) {
  w = sum(weights)
  if (nparam >= w) {theta = Inf} else {
    theta = theta.mm(y = y, mu = mu, weights = weights, dfr = length(y) - nparam)
    theta = ifelse(theta > 1e+4, Inf, ifelse(theta <= 0, Inf, theta))
  }
  return(theta)
}

theta_md = function(y, mu, weights, nparam = NULL) {
  w = sum(weights)
  if (nparam >= w) {theta = Inf} else {
    theta = theta.md(y = y, mu = mu, weights = weights, dfr = sum(weights) - nparam)
    theta = ifelse(theta > 1e+4, Inf, ifelse(theta <= 0, Inf, theta))
  }
  return(theta)
}

pNBII = function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(mu <= 0)) 
    stop(paste("mu must be greater than 0 ", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(q < 0)) 
    stop(paste("y must be >=0", "\n", ""))
  if (length(sigma) > 1) 
    cdf = ifelse(sigma > 1e-04, pnbinom(q, size = mu / sigma, 
                                        mu = mu, lower.tail = lower.tail, log.p = log.p), 
                 ppois(q, lambda = mu, lower.tail = lower.tail, log.p = log.p))
  else cdf = if (sigma < 1e-04) 
    ppois(q, lambda = mu, lower.tail = lower.tail, log.p = log.p)
  else pnbinom(q, size = mu / sigma, mu = mu, lower.tail = lower.tail, 
               log.p = log.p)
  return(cdf)
}

qNBII = function(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(mu <= 0)) 
    stop(paste("mu must be greater than 0 ", "\n", ""))
  if (any(sigma < 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  if (length(sigma) > 1) 
    q = ifelse(sigma > 1e-04, qnbinom(p, size = mu / sigma, 
                                      mu = mu, lower.tail = lower.tail, log.p = log.p), 
               qpois(p, lambda = mu, lower.tail = lower.tail, log.p = log.p))
  else q = if (sigma < 1e-04) 
    qpois(p, lambda = mu, lower.tail = lower.tail, log.p = log.p)
  else qnbinom(p, size = mu/sigma, mu = mu, lower.tail = lower.tail, 
               log.p = log.p)
  return(q)
}


rNBII = function(n, mu = 1, sigma = 1) 
{
  if (any(mu <= 0)) 
    stop(paste("mu must be greater than 0 ", "\n", ""))
  if (any(sigma <= 0)) 
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  n = ceiling(n)
  p = runif(n)
  r = qNBII(p, mu = mu, sigma = sigma)
  return(r)
}


sigma_ml = function(y, mu, weights = NULL)
{
  n = length(y)
  if (is.null(weights)) {weights = rep(1 / n, n)}
  NB2_theta = function(sigma, mu, y, weights) {
    return(sum(n * weights * dNBII(y = y, sigma = sigma, mu = mu, log = TRUE)))
  }
  # start = c(0.01)
  fit = optimize(NB2_theta, y = y, mu = mu, weights = weights, interval = c(1e-6, 1000), maximum = TRUE)
  sigma = fit$maximum
  # sigma = ifelse(sigma <= 5e-5, 0, sigma)
  return(sigma)
}
