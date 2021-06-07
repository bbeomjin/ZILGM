# devtools::install_github("bbeomjin/ZILGM")
require(ZILGM)

cal_adj_pfms = ZILGM:::cal_adj_pfms

pars_set = expand.grid(n = c(50, 100, 200), theta = c(1e+8, 1, 0.5, 0.25), zlvs = c(0.1))

nlam = 50
signal = 1.5; noise = 0.0

nb_p_result = list()
nb_res_mat = nb2_res_mat = p_res_mat = lpgm_res_mat = list()
nb_roc = nb2_roc = p_roc = lpgm_roc = vector(mode = "list", length = nrow(pars_set))
nb_time_list = nb2_time_list = p_time_list = lpgm_time_list = vector("list", length = nrow(pars_set))
nb_aic_list = nb2_aic_list = p_aic_list = vector("list", length = nrow(pars_set))
nb_bic_list = nb2_bic_list = p_bic_list = vector("list", length = nrow(pars_set))
nb_ll_list = nb2_ll_list = p_ll_list = vector("list", length = nrow(pars_set))
nb_coef_list = nb2_coef_list = p_coef_list = vector("list", length = nrow(pars_set))

i = 1
j = 10

for (j in 1:nrow(pars_set)) {
  cat("start ==========", j, "================" , "\n")
  kk = j
  n = pars_set[kk, 1]; 
  zlvs = pars_set[kk, 3];
  theta = pars_set[kk, 2]
  
  for (i in 1:100) {
    set.seed(i)
    A = generate_network_structure2(node = p, prob = prob, NofHub = 3, type = "scale-free")
    mdat = zilgm_sim(A = A, n = n, p = p, zlvs = zlvs, signal = signal, noise = noise,
                     theta = theta, family = "negbin")
    
    lam_max = find_lammax(mdat$X)
    lams = exp(seq(log(lam_max), log(1e-4 * lam_max), length.out = 50))
    
    nb_time = system.time((
      simul_nb_result = zilgm(X = mdat$X, lambda = lams, family = "NBI", thresh = 1e-6, EM_tol = 1e-4,
                              EM_iter = 1e+2, nCores = numWorkers, update_type = "IRLS", weights_mat = NULL, do_boot = TRUE,
                              boot_num = 30, init_select = FALSE, beta = 0.05, sym = "OR", maxit = 50)
    ))[3]
    
    nb2_time = system.time((
      simul_nb2_result = zilgm(X = mdat$X, lambda = lams, family = "NBII", thresh = 1e-6, EM_tol = 1e-4,
                               EM_iter = 1e+2, nCores = numWorkers, update_type = "IRLS", weights_mat = NULL, do_boot = TRUE,
                               boot_num = 30, init_select = FALSE, beta = 0.05, sym = "OR", maxit = 50)
    ))[3]
    
    p_time = system.time((
      simul_p_result = zilgm(X = mdat$X, lambda = lams, family = "Poisson", thresh = 1e-6, EM_tol = 1e-4,
                             EM_iter = 1e+2, nCores = numWorkers, update_type = "IRLS", weights_mat = NULL, do_boot = TRUE,
                             boot_num = 30, init_select = FALSE, beta = 0.05, sym = "OR", maxit = 50)
    ))[3]
    
    # lpgm_time = system.time((
    #   simul_lpgm_result = LPGM1(t(mdat$X), method = "LPGM", N = 30, beta = 0.05, stability = "star", lambda.path = lams,
    #                             th = 1e-6)
    # ))[3]
    
    nb_net = simul_nb_result$network[[simul_nb_result$opt_index]]
    nb_coef_net = simul_nb_result$coef_network[, , simul_nb_result$opt_index]
    nb_coef_list[[j]][[i]] = nb_coef_net
    # nb_ll_list[[j]][[i]] = rowMeans(simul_nb_result$loglik)[simul_nb_result$opt_index]
    # nb_aic_list[[j]][[i]] = rowMeans(simul_nb_result$aic)[simul_nb_result$opt_index]
    # nb_bic_list[[j]][[i]] = rowMeans(simul_nb_result$bic)[simul_nb_result$opt_index]
    
    nb2_net = simul_nb2_result$network[[simul_nb2_result$opt_index]]
    nb2_coef_net = simul_nb2_result$coef_network[, , simul_nb2_result$opt_index]
    nb2_coef_list[[j]][[i]] = nb2_coef_net
    # nb2_ll_list[[j]][[i]] = rowMeans(simul_nb2_result$loglik)[simul_nb2_result$opt_index]
    # nb2_aic_list[[j]][[i]] = rowMeans(simul_nb2_result$aic)[simul_nb2_result$opt_index]
    # nb2_bic_list[[j]][[i]] = rowMeans(simul_nb2_result$bic)[simul_nb2_result$opt_index]
    
    p_net = simul_p_result$network[[simul_p_result$opt_index]]
    p_coef_net = simul_p_result$coef_network[, , simul_p_result$opt_index]
    p_coef_list[[j]][[i]] = p_coef_net
    # p_ll_list[[j]][[i]] = rowMeans(simul_p_result$loglik)[simul_p_result$opt_index]
    # p_aic_list[[j]][[i]] = rowMeans(simul_p_result$aic)[simul_p_result$opt_index]
    # p_bic_list[[j]][[i]] = rowMeans(simul_p_result$bic)[simul_p_result$opt_index]
    
    # lpgm_net = simul_lpgm_result$network[[simul_lpgm_result$opt.index]]
    
    nb_roc[[j]][[i]] = cal_adj_pfms(A, nb_net)
    nb2_roc[[j]][[i]] = cal_adj_pfms(A, nb2_net)
    p_roc[[j]][[i]] = cal_adj_pfms(A, p_net)
    # lpgm_roc[[j]][[i]] = cal_adj_pfms(A, lpgm_net)
    nb_time_list[[j]][[i]] = nb_time
    nb2_time_list[[j]][[i]] = nb2_time
    p_time_list[[j]][[i]] = p_time
    # lpgm_time_list[[j]][[i]] = lpgm_time
  }
  i = 1
  
  nb_roc_dat = lapply(nb_roc[[j]], ROC_fun)
  nb2_roc_dat = lapply(nb2_roc[[j]], ROC_fun)
  p_roc_dat = lapply(p_roc[[j]], ROC_fun)
  # lpgm_roc_dat = lapply(lpgm_roc[[j]], ROC_fun)
  
  nb_res_mat[[j]] = do.call(rbind.data.frame, nb_roc_dat)
  nb2_res_mat[[j]] = do.call(rbind.data.frame, nb2_roc_dat)
  p_res_mat[[j]] = do.call(rbind.data.frame, p_roc_dat)
  # lpgm_res_mat[[j]] = do.call(rbind.data.frame, lpgm_roc_dat)
  # nb_p_result[[j]] = t(data.frame(
  #   NB = colMeans(nb_res_mat[[j]]),
  #   NB2 = colMeans(nb2_res_mat[[j]]),
  #   P = colMeans(p_res_mat[[j]]),
  #   LPGM = colMeans(lpgm_res_mat[[j]])))
}