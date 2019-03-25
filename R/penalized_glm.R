penalized_glm = function(x, y, weights, start = NULL, etastart = NULL, mustart = NULL, lambda = NULL,
                         alpha = 1, gamma = 3, rescale = TRUE, penalty.factor = rep(1, p), thresh = 1e-6, 
                         eps.bino = 1e-5, maxit = 1000, eps = .Machine$double.eps, theta, family = c("gaussian", "binomial", "poisson", "negbin"),
                         trace = FALSE)
{
  family = match.arg(family)
  if(!family %in% c("gaussian", "binomial", "poisson", "negbin")){
    print(family)
    stop("'family' not recognizied\n")
  }
  if(family == "gaussian") rescale = FALSE
  
  if (alpha < 0 || alpha > 1){
    cat("alpha=", alpha)
    stop("alpha must be greater than 0 and less than or equal to 1")
  }
  
  if (alpha == 0 && is.null(lambda)){
    stop("not designed for alpha=0 and lambda=NULL\n")
  }
  
  if(!is.null(etastart) && !is.null(mustart)){
    if((length(etastart) != length(mustart)) || length(etastart) != length(y))
      stop("length of etastart and mustart should be the same as y\n") 
  }
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  np = dim(x)
  n = np[1]
  p = np[2]
  
  if(missing(weights)) weights = rep(1, n)
  weights = as.vector(weights)
  
  if(!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  
  ## check weights and offset
  if( !is.null(weights) && any(weights < 0) ){
    stop("negative weights not allowed")
  }
  if(family=="binomial"){
    if(is.factor(y))
      y = as.integer(y) - 1
    if(any(y > 1))
      stop("response should be between 0 and 1 for family=binomial")
  }
  
  if(family == "negbin"){
    if(missing(theta))
      stop("theta has to be provided for family='negbin'\n")
    else if(length(theta) > 1)
      stop("theta has to be one scale parameter for family='negbin'\n")
    else if(theta <= 0){
      cat("theta=", theta, "\n")
      stop("theta has to be positive for family='negbin'\n")
    }
  }
  if(family != "negbin") theta = 1 ### this theta is not useful but required as input in Fortran glmlink subroutine
  tracetype = 0
  if(trace) tracetype = 1
  famtype = switch(family,
                    "gaussian"=1,
                    "binomial"=2,
                    "poisson"=3,
                    "negbin"=4)
  pentype = 1
  stantype = as.numeric(FALSE)
  #  warning("alp = 1/theta in glm negative.binomial function\n")
  if(length(penalty.factor) != p) stop("number of penalty.factor should be the same as nvars")
  im = inactive = seq(p)
  ### compute the pathwise coordinate descent, cf: section 2.5 in Friedman et al., JSS 2010, 33(1)
  if(is.null(weights)) weights = rep(1, n)
  wt = weights / sum(weights)
  if(is.null(mustart) || is.null(etastart)){
    mu = rep(1, n)
    eta = rep(0, n)
  }
  else{
    mu = mustart
    eta = etastart
  }
  if(is.null(lambda)){
    stop("Penalty parameter must be given \n")
  }
  else nlambda = length(lambda)
  ### ref: McCullagh and Nelder, 2nd edition, 1989, page 121
  if(length(mu) != n)
    mu = rep(mu, n)
  nulldev = .Fortran("deveval",
                     n = as.integer(n),
                     y = as.double(y),
                     mu = as.double(mu),
                     theta = as.double(theta),
                     weights = as.double(weights),
                     family = as.integer(famtype),
                     dev = as.double(0), PACKAGE = "ZILGM")$dev
  beta = matrix(0, ncol = 1, nrow = p)
  if(is.null(start))
    startv = 0
  if(!is.null(start)){
    if(length(start) != (p + 1))
      stop("length of start doesn't match x dimension\n")
    else{
      startv = 1
      beta = start[-1]
      b0 = start[1]
    }
  }
  else{
    b0 = 0
    start = rep(0, p + 1)
  }
  resdev = rep(0, 1)
  yhat = matrix(0, n, 1)
  penfac = penalty.factor / sum(penalty.factor) * p
  lam = penfac %o% lambda
  
  if(family == "gaussian") {
    innermaxit = maxit
    maxit = 1
  }
  else innermaxit = 1
  wtnew = weights / sum(weights)
  meanx = drop(wtnew %*% x)
  
  tmp = .Fortran("outloop",
                x = as.double(x),
                y = as.double(y),
                weights = as.double(weights),
                wt = as.double(wt),
                n = as.integer(n),
                m = as.integer(p),
                penalty = as.integer(pentype),
                nlambda = as.integer(nlambda),
                lam = as.double(lam),
                alpha = as.double(alpha),
                gam = as.double(gamma),
                theta = as.double(theta),
                rescale = as.integer(rescale),
                mu = as.double(mu),
                eta = as.double(eta),
                family = as.integer(famtype),
                standardize = as.integer(stantype),
                nulldev = as.double(nulldev),
                thresh = as.double(thresh),
                maxit = as.integer(maxit),
                innermaxit = as.integer(innermaxit),
                eps = as.double(eps),
                trace = as.integer(tracetype),
                start = as.double(start),
                startv = as.integer(startv),
                b = as.double(beta),
                bz = as.double(b0),
                resdev = as.double(resdev),
                ypre = as.double(yhat),
                convout = as.integer(rep(0, 1)), 
                satu = as.integer(0),
                good = as.integer(nlambda),
                ep = as.double(eps.bino),
                outpll = as.double(matrix(0, maxit, 1)), PACKAGE = "ZILGM") 
  if(nlambda > 1 || tmp$satu == 0)
    good = 1:tmp$good # only those non-saturated models
  else if(tmp$satu == 1)
    return(RET = list(satu = 1))
  ### Names
  if(is.null(colnames(x))) varnames = paste("X", 1:ncol(x), sep = "")
  else varnames = colnames(x)
  beta = matrix(tmp$b, ncol = nlambda)[,good]
  beta = as.matrix(beta)
  b0 = tmp$bz[good]
  ### note: pll was from standardized beta values if standardize=TRUE
  if(trace)
    pll = matrix(tmp$outpll, ncol = nlambda)
  else pll = NULL
  normx = NULL
  resdev = tmp$resdev[good]
  yhat = matrix(tmp$ypre, ncol = 1)[,good] 
  theta = theta
  if(is.null(dim(beta)))
    names(beta) = colnames(x)
  else{
    rownames(beta) = colnames(x)
    colnames(beta) = round(lambda, digits = 4)
  }
  RET = list(family = family, satu = tmp$satu, lambda = lambda[good], beta = beta, b0 = b0,
             theta = theta[good], pll = pll, fitted.values = yhat, converged = tmp$convout[good],
             alpha = alpha)
  class(RET) = c("glmreg", "zilgm")
  ###penalized log-likelihood function value for rescaled beta
  
  penval = .Fortran("penGLM", 
                    start=as.double(beta[, 1]),
                    m = as.integer(p),
                    lambda = as.double(RET$lambda * penalty.factor),
                    alpha = as.double(alpha),
                    gam = as.double(gamma),
                    penalty = as.integer(pentype),
                    pen = as.double(0.0), PACKAGE = "ZILGM")$pen
  RET$penval = penval
  RET$df = sum(abs(beta) > 0) ### df= number of nonzero coefficients (intercept excluded)
  return(RET)
}


