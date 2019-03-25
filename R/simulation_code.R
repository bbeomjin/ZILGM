setwd("~/beomfile/simulation_result/")
rm(list = ls())
gc()

source("~/beomfile/ZINBGM2/network_gen.R")
source("~/beomfile/ZINBGM2/pfm_fun.R")
source("~/beomfile/ZINBGM2/zigm_funcs.R")
# source("~/beomfile/ZINBGM2/ZIGM_core.R")
source("~/beomfile/ZINBGM2/ZIPGM_core.R")
source("~/beomfile/ZINBGM2/ZINBGM_core.R")
source("~/beomfile/ZINBGM2/ZINB2GM_core.R")
source("~/beomfile/ZINBGM2/ZIGM.R")
source("~/beomfile/ZINBGM2/wlasso.R")
source("~/beomfile/ZINBGM2/LPGM1.R")
source("~/beomfile/ZINBGM2/penalized_glm.R")

require(parallel)
require(MASS)
require(XMRF)
require(igraph)

# load("./STAR_select_simulation_p30_4_180808.RData")

numWorkers = detectCores()

p = 30
prob = 0.05
# prob = 2 / p

set.seed(111)
A = generate_network(p, prob)
# A = barabasi.game(n = p, power = 0.01, zero.appeal = p, directed = FALSE) %>%
#   get.adjacency(., type = "both") %>%
#   as.matrix()
pars_set = expand.grid(n = c(50, 100, 200), theta = c(1e+8, 1, 0.5, 0.25), zlvs = c(0.1))

nlam = 100
signal = 1.5; noise = 0.0

nb_p_result = list()
nb_res_mat = nb2_res_mat = p_res_mat = lpgm_res_mat = list()
nb_roc = nb2_roc = p_roc = lpgm_roc = vector(mode = "list", length = nrow(pars_set))
# p = 30, j = 8, i = 48~50 돌리기
for (j in 1:nrow(pars_set)) {
  cat("start ==========", j, "================" , "\n")
  kk = j
  n = pars_set[kk, 1]; 
  zlvs = pars_set[kk, 3];
  theta = pars_set[kk, 2]
  
  i = 1
  for (i in i:50) {
    set.seed(i)
    mdat = ex_rgraph(A = A, n = n, p = p, zlvs = zlvs, signal = signal, noise = noise,
                     theta = theta, family = "negbin")
    
    # simul_nb_result = zigm(X = mdat$X, nlambda = 30, family = "NBI", thres = 1e-6, EM_thresh = 1e-6,
    #                        EM_iter = 100, nCores = numWorkers, update_type = "IRLS", weights_mat = NULL, do_boot = TRUE,
    #                        boot_num = 30, init_select = FALSE, beta = 0.05, edge_type = "OR",
    #                        thresh = 1e-6, maxit = 1e+3)
    
    # system.time((simul_snb_result = zigm(X = mdat$X, nlambda = 50, family = "NBI", thres = 1e-6, EM_thresh = 1e-6,
    #                         EM_iter = 300, nCores = numWorkers, update_type = "IRLS", weights_mat = NULL, do_boot = TRUE,
    #                         boot_num = 30, init_select = FALSE, beta = 0.05, fixed_theta = T, edge_type = "OR",
    #                         thresh = 1e-6, maxit = 1e+3, rescaled = TRUE)))
    
    simul_nb2_result = zigm(X = mdat$X, nlambda = 30, family = "NBII", thres = 1e-6, EM_thresh = 1e-6,
                            EM_iter = 1000, nCores = numWorkers, update_type = "IRLS", weights_mat = NULL, do_boot = TRUE,
                            boot_num = 30, init_select = FALSE, beta = 0.05, edge_type = "OR", thresh = 1e-6, maxit = 1e+3)
    
    # simul_p_result = zigm(X = mdat$X, nlambda = 30, family = "Poisson", thres = 1e-6, EM_thresh = 1e-6,
    #                       EM_iter = 500, nCores = numWorkers, update_type = "IRLS", weights_mat = NULL, do_boot = TRUE,
    #                       boot_num = 30, init_select = FALSE, beta = 0.05, thresh = 1e-6, maxit = 1e+3, edge_type = "OR")
    
    # lpgm_lams = c(exp(seq(log(rho.max + 2),log(rho.min), length = nlam)))
    simul_lpgm_result = LPGM1(t(mdat$X), method = "LPGM", N = 30, beta = 0.05, stability = "star", nlams = 30,
                              th = 1e-6)
    
    # simul_pgm_result = XMRF:::PGM(t(mdat$X), method = "LPGM", N = 30, beta = 0.05, stability = "star", nlams = 50)
    
    # nb_net = simul_nb_result$network[[simul_nb_result$opt_index]]
    # snb_net = simul_snb_result$network[[simul_snb_result$opt_index]]
    nb2_net = simul_nb2_result$network[[simul_nb2_result$opt_index]]
    # p_net = simul_p_result$network[[simul_p_result$opt_index]]
    lpgm_net = simul_lpgm_result$network[[simul_lpgm_result$opt.index]]
    # pgm_net = simul_pgm_result$network[[simul_pgm_result$opt.index]]
    
    # nb_roc[[j]][[i]] = cal_adj_pfms(A, nb_net)
    # snb_roc[[j]][[i]] = cal_adj_pfms(A, snb_net)
    nb2_roc[[j]][[i]] = cal_adj_pfms(A, nb2_net)
    # p_roc[[j]][[i]] = cal_adj_pfms(A, p_net)
    lpgm_roc[[j]][[i]] = cal_adj_pfms(A, lpgm_net)
    # cal_adj_pfms(A, pgm_net)
    save.image("STAR_select_simulation_p30_1910225(OR,subsamp632).RData")
  }
  
  nb_roc_dat = lapply(nb_roc[[j]], function(x) {
    TPR = x[6] / x[4]
    FNR = x[5] / x[4]
    TNR = x[2] / x[1]
    FPR = x[3] / x[1]
    precision = x[6] / (x[3] + x[6])
    return(list(TPR = TPR, FNR = FNR, TNR = TNR, FPR = FPR, Precision = precision))
  })
  
  # snb_roc_dat = lapply(snb_roc[[j]], function(x) {
  #   TPR = x[6] / x[4]
  #   FNR = x[5] / x[4]
  #   TNR = x[2] / x[1]
  #   FPR = x[3] / x[1]
  #   precision = x[6] / (x[3] + x[6])
  #   return(list(TPR = TPR, FNR = FNR, TNR = TNR, FPR = FPR, Precision = precision))
  # })
  
  nb2_roc_dat = lapply(nb2_roc[[j]], function(x) {
    TPR = x[6] / x[4]
    FNR = x[5] / x[4]
    TNR = x[2] / x[1]
    FPR = x[3] / x[1]
    precision = x[6] / (x[3] + x[6])
    return(list(TPR = TPR, FNR = FNR, TNR = TNR, FPR = FPR, Precision = precision))
  })
  
  p_roc_dat = lapply(p_roc[[j]], function(x) {
    TPR = x[6] / x[4]
    FNR = x[5] / x[4]
    TNR = x[2] / x[1]
    FPR = x[3] / x[1]
    precision = x[6] / (x[3] + x[6])
    return(list(TPR = TPR, FNR = FNR, TNR = TNR, FPR = FPR, Precision = precision))
  })
  
  
  lpgm_roc_dat = lapply(lpgm_roc[[j]], function(x) {
    TPR = x[6] / x[4]
    FNR = x[5] / x[4]
    TNR = x[2] / x[1]
    FPR = x[3] / x[1]
    precision = x[6] / (x[3] + x[6])
    return(list(TPR = TPR, FNR = FNR, TNR = TNR, FPR = FPR, Precision = precision))
  })
  
  nb_res_mat[[j]] = do.call(rbind.data.frame, nb_roc_dat)
  # snb_res_mat[[j]] = do.call(rbind.data.frame, snb_roc_dat)
  nb2_res_mat[[j]] = do.call(rbind.data.frame, nb2_roc_dat)
  p_res_mat[[j]] = do.call(rbind.data.frame, p_roc_dat)
  lpgm_res_mat[[j]] = do.call(rbind.data.frame, lpgm_roc_dat)
  nb_p_result[[j]] = t(data.frame(
    NB = colMeans(nb_res_mat[[j]]),
    # SNB = colMeans(snb_res_mat[[j]]),
    NB2 = colMeans(nb2_res_mat[[j]]),
    P = colMeans(p_res_mat[[j]]),
    LPGM = colMeans(lpgm_res_mat[[j]])))
}


nb_p_result[[8]]
t(cbind(cbind(apply(nb2_res_mat[[8]], 2, sd) / sqrt(50), apply(p_res_mat[[8]], 2, sd) / sqrt(50)),
        apply(lpgm_res_mat[[8]], 2, sd) / sqrt(50)))

sd(nb2_res_mat[[11]]$TPR - p_res_mat[[11]]$TPR) / sqrt(50)

wilcox.test(snb_res_mat1[[13]]$TPR, snb_res_mat2[[13]]$TPR, paired = T)
wilcox.test(p_res_mat[[13]]$Precision, snb_res_mat1[[13]]$Precision, paired = T)

sum(snb_res_mat[[13]]$Precision > p_res_mat[[13]]$Precision)
sum(snb_res_mat[[13]]$TPR < p_res_mat[[13]]$TPR)


colMeans(nb2_res_mat[[10]])


colMeans(p_res_mat[[8]])
apply(p_res_mat[[8]], 2, sd) / sqrt(200)

colMeans(nb2_res_mat[[8]])
apply(nb2_res_mat[[8]], 2, sd) / sqrt(100)

colMeans(p_res_mat[[4]])
apply(p_res_mat[[4]], 2, sd) / sqrt(50)

f1_score = function(data)
{
  return(2 * data$TPR * data$TNR / (data$TPR + data$TNR))
}

f1_nb_3 = f1_score(nb_res_mat[[3]])
mean(f1_nb_3)
sd(f1_nb_3) / sqrt(length(f1_nb_3))

f1_p_3 = f1_score(p_res_mat[[3]])
mean(f1_p_3)
sd(f1_p_3) / sqrt(length(f1_p_3))

f1_nb_4 = f1_score(nb_res_mat[[4]])
mean(f1_nb_4)
sd(f1_nb_4) / sqrt(length(f1_nb_4))

f1_p_4 = f1_score(p_res_mat[[4]])
mean(f1_p_4)
sd(f1_p_4) / sqrt(length(f1_p_4))


source("~/beomfile/ZINBGM2/ZIPGM_core.R")
source("~/beomfile/ZINBGM2/ZINBGM_core.R")
source("~/beomfile/ZINBGM2/ZINB2GM_core.R")
source("~/beomfile/ZINBGM2/zigm_funcs.R")
# source("~/beomfile/ZINBGM2/ZIGM_core_tmp2.R")
m = round(0.632 * n)
set.seed(i)
inx = 10
sub_ind = sample(1:n, m, replace = FALSE)

y = mdat$X[sub_ind, inx]
x = mdat$X[sub_ind, -inx]

y = mdat$X[, inx]
x = mdat$X[, -inx]


lambda = simul_snb_result1$opt_lambda
lambda = 1
update_type = "IRLS"
fixed_theta = TRUE
rescaled = TRUE
tmp_type = "zip"
EM_iter = 100
EM_thresh = 1e-6
thres = 1e-2
thresh = 1e-6
maxit = 1e+3

zeroinfl(y ~ x|1, dist = "negbin")

system.time((fit1 = zigm_coord_negbin(y = y, x = x, lambda = 0.1, update_type = "IRLS", fixed_theta = FALSE, rescaled = FALSE,
                                      EM_iter = 100, EM_thresh = 1e-4, thres = 1e-6, thresh = 1e-5, maxit = 3e+2)))


system.time((fit2 = zigm_coord_negbin2(y = y, x = x, lambda = 0.1, update_type = "IRLS", 
                                       EM_iter = 2000, EM_thresh = 1e-4, thres = 1e-6, thresh = 1e-5, maxit = 3e+2)))


system.time((fit3 = zigm_coord_negbin(y = y, x = x, lambda = 0.1, update_type = "IRLS", fixed_theta = FALSE, rescaled = FALSE,
                                      EM_iter = 100, EM_thresh = 1e-6, thres = 1e-6, thresh = 1e-6, maxit = 1e+5)))


system.time((fit4 = zigm_coord_poisson(y = y, x = x, lambda = 0.1, update_type = "IRLS", 
                                       EM_iter = 100, EM_thresh = 1e-4, thres = 1e-6, thresh = 1e-5, maxit = 3e+2)))
