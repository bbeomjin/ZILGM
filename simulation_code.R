rm(list = ls())
gc()

# If the ZILGM package does not installed, please run below line
# devtools::install_github("bbeomjin/ZILGM")
require(ZILGM)

# Define function to compute true positive rate (TPR) and precision
ROC_fun = function(x) {
  TPR = x[6] / x[4]
  FNR = x[5] / x[4]
  TNR = x[2] / x[1]
  FPR = x[3] / x[1]
  precision = x[6] / (x[3] + x[6])
  return(list(TPR = TPR, FNR = FNR, TNR = TNR, FPR = FPR, Precision = precision))
}

# Define the combination of experimental settings
# n is the number of samples, theta is the over-dispersion parameter and zlvs is the probability of zero
pars_set = expand.grid(n = c(50, 100, 200), theta = c(1e+8, 1, 0.5, 0.25), zlvs = c(0.1))

# The number of nodes
p = 30

# The number of cores. If the OS is windows, only the value 1 is supported
numWorkers = 1

# Grid for regularization parameter lambda
nlam = 50

# The probability of randomly connected to each nodes
prob = 2 / p

# True mu and noise mu
signal = 1.5; noise = 0.0

nb_p_result = list()
nb_res_mat = nb2_res_mat = p_res_mat = list()
nb_roc = nb2_roc = p_roc = vector(mode = "list", length = nrow(pars_set))
nb_time_list = nb2_time_list = p_time_list = vector("list", length = nrow(pars_set))

for (j in 1:nrow(pars_set)) {
  cat("start ==========", j, "================" , "\n")
  kk = j
  n = pars_set[kk, 1]; 
  zlvs = pars_set[kk, 3];
  theta = pars_set[kk, 2]
  
  for (i in 1:100) {
    set.seed(i)
    
    # Generate network structure. If type == "hub", the hub network is generated and if type == "random", the random network is generated.
    A = generate_network(node = p, prob = prob, type = "scale-free")
    
    # Generate simulated data for the graph structure
    mdat = zilgm_sim(A = A, n = n, p = p, zlvs = zlvs, signal = signal, noise = noise,
                     theta = theta, family = "negbin")
    
    # Find maximum of regularization parameter
    lam_max = find_lammax(mdat$X)
    
    # Define the sequence of regularization parameters
    lams = exp(seq(log(lam_max), log(1e-4 * lam_max), length.out = nlam))
    
    # Fitting the ZILNBGM-I
    nb_time = system.time((
      simul_nb_result = zilgm(X = mdat$X, lambda = lams, family = "NBI", update_type = "IRLS", 
                              do_boot = TRUE, boot_num = 30, beta = 0.05, sym = "OR", nCores = numWorkers)
    ))[3]
    
    # Fitting the ZILNBGM-II
    nb2_time = system.time((
      simul_nb2_result = zilgm(X = mdat$X, lambda = lams, family = "NBII", update_type = "IRLS", 
                               do_boot = TRUE, boot_num = 30, beta = 0.05, sym = "OR", nCores = numWorkers)
    ))[3]
    
    # Fitting the ZILPGM
    p_time = system.time((
      simul_p_result = zilgm(X = mdat$X, lambda = lams, family = "Poisson", update_type = "IRLS",
                             do_boot = TRUE, boot_num = 30, beta = 0.05, sym = "OR", nCores = numWorkers)
    ))[3]
    
    # Select the optimal estimated network
    nb_net = simul_nb_result$network[[simul_nb_result$opt_index]]
    nb2_net = simul_nb2_result$network[[simul_nb2_result$opt_index]]
    p_net = simul_p_result$network[[simul_p_result$opt_index]]
    
    # Compute true positive (TP) and false positive (FP) to compute the ROC
    nb_roc[[j]][[i]] = ZILGM:::cal_adj_pfms(A, nb_net)
    nb2_roc[[j]][[i]] = ZILGM:::cal_adj_pfms(A, nb2_net)
    p_roc[[j]][[i]] = ZILGM:::cal_adj_pfms(A, p_net)
    
    nb_time_list[[j]][[i]] = nb_time
    nb2_time_list[[j]][[i]] = nb2_time
    p_time_list[[j]][[i]] = p_time
  }
  nb_roc_dat = lapply(nb_roc[[j]], ROC_fun)
  nb2_roc_dat = lapply(nb2_roc[[j]], ROC_fun)
  p_roc_dat = lapply(p_roc[[j]], ROC_fun)
  
  nb_res_mat[[j]] = do.call(rbind.data.frame, nb_roc_dat)
  nb2_res_mat[[j]] = do.call(rbind.data.frame, nb2_roc_dat)
  p_res_mat[[j]] = do.call(rbind.data.frame, p_roc_dat)
}