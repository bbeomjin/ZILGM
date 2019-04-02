require(network)

generate_random_network = function(node, prob) {
  # Generate matrix values, sampling 0 or 1 with given probabilities
  matvals = sample(c(0, 1), node * (node - 1)/2, 
                   replace = TRUE, prob = c(1 - prob, prob))
  
  # From the values above, generate a symmetric matrix
  networkmat = matrix(rep(0, node * node), ncol = node)
  mv = 1
  for (i in 1:node) {
    for (j in 1:node) {
      if (i > j) {
        networkmat[i, j] = networkmat[j, i] = matvals[mv]
        mv = mv + 1
      }
    }
  }
  return(networkmat)
}

zilgm_sim = function(prob, n, p, zlvs, family = c("poisson", "negbin"),
                     signal, theta = NULL, noise)
{
  family = match.arg(family)
  is.symm = TRUE
  
  npair = p * (p - 1) / 2
  Y = E = matrix(0, nrow = n, ncol = p)
  Ypair = matrix(0, n, npair)
  
  if (family == "poisson") {
    for (j in 1:p)
    {
      Y[, j] = rpois(n = n, lambda = signal) # rZeroPoisson(p+pair, w=pi0, mu=mu_t)
      E[, j] = rpois(n = n, lambda = noise) #rZeroPoisson(m, w=pi0, mu=mu_n) # rpois(p, mu_n)
    }
    
    for(j in 1:npair)
    {
      Ypair[, j] = rpois(n = n, lambda = signal) # rZeroPoisson(p + pair, w = pi0, mu = mu_t)
    }
  }
  
  if (family == "negbin") {
    if (is.null(theta)) {stop("Need a theta")}
    
    for (j in 1:p)
    {
      Y[, j] = rnegbin(n = n, mu = signal, theta = theta) # rZeroPoisson(p+pair, w=pi0, mu=mu_t)
      E[, j] = rpois(n = n, lambda = noise) #rZeroPoisson(m, w=pi0, mu=mu_n) # rpois(p, mu_n)
    }
    
    for (j in 1:npair)
    {
      Ypair[, j] = rnegbin(n = n, mu = signal, theta = theta) # rZeroPoisson(p + pair, w = pi0, mu = mu_t)
    }
  }
  
  Y = cbind(Y, Ypair)
  
  # B is real signal
  tri_A = t(A)[lower.tri(t(A), diag = FALSE)]
  B = matrix(0, nrow = npair, ncol = p)
  ix = vec2tri(1:npair, p)
  adj = cbind(ix, tri_A)
  
  nadj = (1:npair)[adj[, 3] != 0]
  for(i in nadj)
  {
    B[i, adj[i, 2]] = 1
    
    if(is.symm) {B[i, adj[i, 1]] = 1}
  }
  B = rbind(diag(1, p), B)
  X = Y %*% B + E
  
  for(j in 1:p)
  {
    X[runif(n) <= zlvs, j] = 0
  }
  return(list(A = A, X = X))
}
