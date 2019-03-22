cal_adj_pfms = function (net, hat_net)
{
  hat_net = abs(hat_net != 0)
  p = ncol(net)
  pfms = matrix(0, p, 6)
  colnames(pfms) = c("ZE","TN", "FP", "NZ", "FN", "TP")
  rownames(pfms) = 1:p
  
  for (j in 1:p) {
    # netative cases
    flag = net[, j] == 0
    flag[j] = FALSE
    pfms[j, 1] = sum(flag)
    pfms[j, 2] = sum(net[flag, j] == hat_net[flag, j]) # TN
    pfms[j, 3] = sum(net[flag, j] != hat_net[flag, j]) # FP
    
    # positive cases
    flag = net[, j] != 0
    pfms[j, 4] = sum(flag)
    pfms[j, 5] = sum(net[flag, j] != hat_net[flag, j]) # FN
    pfms[j, 6] = sum(net[flag, j] == hat_net[flag, j]) # TP 
  }
  return(colSums(pfms))
}

roc = function(hat_net, net)
{
  lams_len = length(hat_net)
  rocm = matrix(nrow = lams_len, ncol = 2)
  
  for (i in 1:lams_len) {
    vec = cal_adj_pfms(net, hat_net[[i]])
    rocm[i, ] = c(vec[3]/vec[1], vec[6]/vec[4])
  }
  return(rocm)
}

# roc
ftpr_to_roc = function (rocm, upper_fpr = 0.1)
{
  if (max(rocm[, 1]) < upper_fpr) {
    upper_fpr = max(rocm[, 1])
  }  
  
  flag = rocm[, 1] <= upper_fpr
  rocm = rocm[flag, , drop = FALSE]
  n_pos = nrow(rocm)  
  
  if (rocm[n_pos, 1] < upper_fpr) {
    rocm = rbind(rocm, c(upper_fpr, rocm[n_pos, 2]))
  }
  
  if (rocm[1, 1] > 0 & rocm[1, 2] > 0) {
    rocm = rbind(0, rocm)
  }
  
  #val <- auc(rocm[,1],rocm[,2])
  #if(val<0)
  val = auc2(rocm)
  return(val)
}
auc2 = function (rocm)
{
  n = nrow(rocm)
  tot = 0
  
  for (i in 1:(n - 1)) {
    tot = tot + (rocm[i + 1, 1] - rocm[i, 1]) * (rocm[i + 1, 2] + rocm[i, 2])/2
  }  
  return(tot)
}

ftpr_to_roc2 = function (rocm, upper_fpr = 0.1)
{
  
  adj_sum = function (vec)
  {
    n = length(vec)
    return(vec[1:(n - 1)]+c(0, vec[1:(n - 2)]))
  }
  
  if (max(rocm[, 1]) < upper_fpr) {upper_fpr = max(rocm[, 1])}
  
  flag = rocm[, 1] <= upper_fpr
  
  if (sum(flag) == 0) {flag[c(1, 2, 3)] = TRUE}
  
  rocm = rocm[flag, ,drop = FALSE]
  n_pos = nrow(rocm)
  rocm = rbind(rocm, c(upper_fpr, rocm[n_pos, 2]))
  
  if (rocm[1, 1] > 0 & rocm[1, 2] > 0) {rocm = rbind(0, rocm)}
  
  l = diff(rocm[, 1], differences = 1)
  h = adj_sum(rocm[, 2])
  ps = l * h/2
  return(sum(ps))
}

