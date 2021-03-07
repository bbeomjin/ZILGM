# require(MASS)
# require(glmnet)
# require(mpath)
# require(pscl)

zilgm = function(X, lambda = NULL, nlambda = 50, family = c("Poisson", "NBI", "NBII"), update_type = c("IRLS", "MM"),
                sym = c("AND", "OR"), theta = NULL, thresh = 1e-6, weights_mat = NULL, penalty_mat = NULL,
                do_boot = FALSE, boot_num = 10, beta = 0.05, lambda_min_ratio = 1e-4,
                init_select = FALSE, nCores = 1, ...)
{

  family = match.arg(family)
  update_type = match.arg(update_type)
  sym = match.arg(sym)
  fun_call = match.call()


  if (!any(class(X) %in% "matrix")) {
    X = as.matrix(X)
  }

  if (!any(class(X) %in% "matrix")) {
    stop("X must be a matrix")
  }

  if (any(lambda < 0)) {
    stop("lambda must be non-negative values")
  }

  n = NROW(X)
  p = NCOL(X)

  if (p < 2) {stop("X must be a matrix with 2 or more columns")}

  penalty = "LASSO"

  cat("learning for ", family, " graphical model \n",
      "nlambda : ", nlambda, "\n",
      "penalty function : ", penalty, "\n",
      "update type : ", update_type, "\n", sep = "")

  if (is.null(lambda)) {
    cat("\n Searching lambda \n")

    rho_max = find_lammax(X)
    rho_min = lambda_min_ratio * rho_max
    tmp_lams = c(exp(seq(log(rho_max),log(rho_min), length = 15)))

    tmp_net = zigm_network(X = X, lambda = tmp_lams, family = family, update_type = update_type, sym = sym, theta = theta,
                           thresh = thresh, weights_mat = weights_mat, penalty_mat = penalty_mat,
                           init_select = init_select, nCores = nCores, n = n, p = p, ...)

    nOfEdge = unlist(lapply(tmp_net$hat_net, function(x) sum(x != 0)))
    s_lam = tmp_lams[which.max(nOfEdge > 1)]
    e_lam = tmp_lams[which.max(nOfEdge)]
    lambda = seq(s_lam, e_lam, length = nlambda)
    rm(tmp_net); gc();
    cat("Complete \n")
  } else {
    nlambda = length(lambda)
  }

  out = list()

  if (do_boot) {
    if (n < 250) {m = round(0.632 * n)} else {m = round(10 * sqrt(n))}

    boot_tmp = vector(mode = "list", length = nlambda)
    for (i in 1:nlambda) {boot_tmp[[i]] = Matrix(0, p, p)}

    for (b in 1:boot_num) {
      cat(paste("Conducting sampling in progress : ", floor(100 * (b/boot_num)), "%", collapse = ""), "\r")
      flush.console()

      # bad_sample = 1
      # while (bad_sample) {
      #   sub_ind = sample(1:n, m, replace = FALSE)
      #   if (sum(apply(X[sub_ind, , drop = FALSE], 2, function(x) (sum(x != 0) < 2))) == 0) {
      #     bad_sample = 0
      #   }
      # }
      sub_ind = sample(1:n, m, replace = FALSE)

      boot_net = zigm_network(X = X[sub_ind, , drop = FALSE], lambda = lambda, family = family, update_type = update_type,
                              sym = sym, theta = theta, thresh = thresh, weights_mat = weights_mat, penalty_mat = penalty_mat,
                              init_select = init_select, nCores = nCores, n = m, p = p, ...)

      for (l in 1:nlambda) {
        boot_tmp[[l]] = boot_tmp[[l]] + boot_net$hat_net[[l]]
      }
    }

    v = rep(0, nlambda)
    for (l in 1:nlambda) {
      gv = as.matrix(boot_tmp[[l]] / boot_num)
      gv_tmp = 2 * gv * (1 - gv)
      v[l] = mean(gv_tmp[upper.tri(gv_tmp)])
    }
    rm(boot_tmp); gc();

    opt_index = max(which.max(v >= beta)[1] - 1, 1)
    opt_lambda = lambda[opt_index]

    out$v = v
    out$opt_index = opt_index
    out$opt_lambda = opt_lambda
  }

  net = zigm_network(X = X, lambda = lambda, family = family, update_type = update_type,
                     sym = sym, theta = theta, thresh = thresh, weights_mat = weights_mat, penalty_mat = penalty_mat,
                     init_select = init_select, nCores = nCores, n = n, p = p, ...)

  out$network = net$hat_net
  out$coef_network = net$coef_net
  out$lambda = lambda
  out$call = fun_call
  return(out)
}



zigm_network = function(X, lambda = NULL, family = c("Poisson", "NBI", "NBII"), update_type = c("IRLS", "MM"), sym = c("AND", "OR"), theta = NULL,
                        thresh = 1e-6, weights_mat = NULL, penalty_mat = NULL, init_select = FALSE, nCores = 1, n, p, ...)
{
  family = match.arg(family)
  update_type = match.arg(update_type)
  sym = match.arg(sym)

  coord_fun = switch(family,
                     Poisson = zilgm_poisson,
                     NBI = zilgm_negbin,
                     NBII = zilgm_negbin2)

  nlambda = length(lambda)
  coef_mat = array(dim = c(p, p, nlambda))

  if (is.null(weights_mat)) {
    weights_mat = matrix(1, n, p)
  }

  if (any(weights_mat < 0)) {"The elements in weights_mat must have non-negative values"}
  if ((NROW(weights_mat) != n) | (NCOL(weights_mat) != p)) {"The number of elements in weights_mat not equal to the number of rows and columns on X"}

  coef_tmp = mclapply(1:p, FUN = function(j) {zigm_wrapper(jth = j, X = X, lambda = lambda, family = family, update_type = update_type, theta = theta,
                                                           thresh = thresh, weights = weights_mat[, j], penalty.factor = penalty_mat[, j],
                                                           init_select = init_select, fun = coord_fun, n = n, p = p, nlambda = nlambda, ...)},
                      mc.cores = nCores, mc.preschedule = FALSE)

  for (j in 1:p) {
    coef_mat[, j, ] = as.matrix(coef_tmp[[j]]$Bmat)
  }

  ghat = lapply(1:nlambda, FUN = function(l) hat_net(coef_mat[, , l], thresh = thresh, type = sym))
  gs = lapply(1:nlambda, FUN = function(l) Matrix(ghat[[l]]))

  return(list(hat_net = gs, coef_net = coef_mat))
}


zigm_wrapper = function(jth, X, lambda, family, update_type, theta, weights, penalty.factor, init_select, fun,
                        n, p, nlambda, thresh, ...)
{
  seqP = 1:p
  Bmat = Matrix(0, p, nlambda, sparse = TRUE)
  b0 = rep(0, nlambda)

  if (init_select) {

    fit0 = glmnet(x = X[, -jth, drop = FALSE], y = X[, jth], standardize = FALSE,
                  family = "poisson", nlambda = 100, dfmax = p)
    bic = (1 - fit0$dev.ratio) * fit0$nulldev + 2 * fit0$df
    p0.b = which.min(bic[-1])
    # lam_ind = ifelse(p0.b >= length(fit0$lambda), p0.b, p0.b + 1)
    lam_ind = p0.b
    coeff = drop(predict.glmnet(fit0, s = fit0$lambda[lam_ind], type = "coefficients"))
    nset = seqP[-jth][which(abs(coeff[-1]) > (thresh / 100))]

    wthres = thresh / 100
    for (init_iter in 1:100) {
      if (length(nset) == 0) {
        wthres = wthres/10
        nset = seqP[-jth][which(abs(coeff[-1]) > wthres)]
        # bvec_init = c(coeff[1],coeff[-1][abs(coeff[-1]) > wthres])/2
      } else {
        # bvec_init = c(coeff[1],coeff[-1][abs(coeff[-1]) > wthres])/2
        break
      }
    }
  } else {
    nset = seqP[-jth]
  }

  # elastic net; if alpha = 1, LASSO penalty, if alpha = 0, ridge penalty
  if (length(nset) == 0) {
    Bmat = Bmat
    b0 = b0
  } else {
    for (iter in 1:nlambda) {
      cat("lambda = ", lambda[iter], ", ", jth, "/", p, "th node learning \n", sep = "")
      coef_res = fun(x = X[, nset, drop = FALSE], y = X[, jth], lambda = lambda[iter], init_theta = theta, weights = weights,
                     update_type = update_type, penalty.factor = penalty.factor, thresh = thresh, ...)

      Bmat[nset, iter] = coef_res$bvec[-1]
      b0[iter] = coef_res$bvec[1]
    }
  }
  return(list(b0 = b0, Bmat = Bmat))
}
