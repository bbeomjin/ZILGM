wlasso <- function(X, y, eta0=0, wID=rep(1,nrow(X)), weight=rep(1,ncol(X)), 
                   maxStep=1e3, eps=1e-10, stand.scale=FALSE, trace=FALSE){
  #trace <- TRUE
  n = length(y)
  p = ncol(X)    # number of components
  my <- mean(y)
  y <- y-my
  mx <- apply(X,2,mean)  
  
  X <- X - matrix(1,nrow=n,ncol=1) %*% mx 
  
  if(stand.scale){
    sdx <- apply(X,2,sd)
    X <- sweep(X, 2, sdx, "/")  
  }  
  tX <- t(X)
  
  seqN = 1:n
  seqP = 1:p
  
  BetaMatr     = matrix(0,maxStep, p)
  ObjValTrace  = rep(0,maxStep)
  LamTrace     = rep(0,maxStep)
  
  allPredIncluded = FALSE    # if all the predictors are included, then it is set to TRUE
  
  ####### Initial Solution ######
  conv = FALSE
  Step = 1
  wmat <- diag(wID)
  weight[abs(weight)< 1e-10] <- 1e-10
  
  wgcorr = -drop(tX %*% wmat %*% y)/weight
  Beta  = rep(0,p)
  eta   = max(abs(wgcorr))
  
  if(eta<=eta0){
    conv = TRUE
    BetaMatr[1,] <- 0
    LamTrace[1] <- eta0 
    ObjValueTrace = 1
  }
  
  V = which(eta==abs(wgcorr))
  Sign = rep(0,p)
  Sign[V] = -sign(wgcorr[V])
  
  ########## Start Loop ###########
  while(!conv){
    #cat("eta   >> ", eta, "\t wgcorr >> ", max(abs(wgcorr)), "\n")
    
    fx = drop(X[,V,drop=F] %*% Beta[V])
    Residual = y-fx
    wgcorr = -drop(t(Residual) %*% wmat %*% X)/weight
    
    ObjValTrace[Step] = getObjective(Residual, wID, Beta, eta, weight)
    BetaMatr[Step,]   = Beta
    LamTrace[Step]    = eta
    if(trace){
      current_status(Step, V, Beta, eta, Sign, wgcorr, Residual, ObjValTrace)
    }
    if(eta<=eta0){
      cat(">>>>eta", eta, "\n")
      conv = TRUE
      break
    }else if(Step > maxStep){   
      #cat("Check 3 : Limited iteration!!!\n")
      conv = FALSE
      break
    }
    #cat("############################Step ",Step, "th EVENT SEARCH########################################################\n")
    ################################################################################
    # Step 1 : Get the right derivatives : dbeta0/ds, dbeta/ds, deta/ds, dbeta^c/ds
    ################################################################################
    ##### Solve the linear system for derivatives : A%*%sol = b
    nV = length(V)
    d_beta  = rderiv(XV=X[,V,drop=FALSE], wmat, nV, weight[V]*Sign[V], eps=eps)
    if(trace){
      cat(d_beta, "\n")
    }
    if(sum(abs(d_beta))<eps){
      cat("coefficients don't move with s\n")
      break       
    }
    # Step 1 : Compute the residual for every points
    d_wgcorr = 1/weight*drop(tX %*% wmat %*% X[,V,drop=FALSE] %*% d_beta)
    d_wgcorr[abs(d_wgcorr)<eps] <- 0
    ################################################################################
    # Step 2 : Compute how much increase of s is needed to get to each type of event
    ################################################################################
    # Find which event happens first, 
    Events = Find.Event(p, V, Beta, eta, eta0, wgcorr, d_wgcorr, Sign, d_beta, trace=trace, eps=eps)
    EPO = Events$EPO
    ds  = Events$ds
    
    d_S1_var = Events$d_S1_var
    d_S2_var = Events$d_S2_var
    
    d_lamT1 = Events$d_lamT1
    d_lamT2 = Events$d_lamT2
    
    # Update the current solutions and step size
    Beta[V] = Beta[V] - ds*d_beta
    eta     = eta - ds
    #cat("BETA",Beta,eta,"\n")        
    ########################################################
    # Step 4 : Update V.
    if(EPO==1){
      # a predictor reduces to 0
      removedVar = (seqP[V])[ds==d_S1_var]
      if(Sign[removedVar] > 0){
        if(trace){
          cat("################## Event 1: variable", removedVar, "is removed from V+\n")
        }
        removedSign = 1
      }else if(Sign[removedVar]<0){
        if(trace){
          cat("################## Event 1: variable", removedVar, "is removed from V-\n")
        }
        removedSign = -1
      }else{
        if(trace){
          cat("error: removedVar is not in V\n")
        }
      }
      V = setdiff(V, removedVar)
      Sign[removedVar] = 0
    }else if(EPO==2){                
      temp = min(c(1:(p-nV))[abs(ds-d_S2_var)<eps])
      newVar = seqP[-V][temp]
      if(abs(d_lamT1[temp] - ds)<eps){
        Sign[newVar] = 1
        if(trace){
          cat("################# Event 2:  New variable", newVar, "is added into V+\n")
        }
      }else if(abs(d_lamT2[temp]-ds)<eps){
        Sign[newVar] = -1
        if(trace){
          cat("################# Event 2: New variable", newVar, "is added into V-\n")
        }
      }
      V = sort(union(V, newVar))
    }else if(EPO==3){ 
      if(trace){
        cat("reach at ", eta0, "\n")
      }
      Step <- Step + 1
      fx = drop(X[,V,drop=F] %*% Beta[V])
      Residual = y-fx   
      ObjValTrace[Step] = getObjective(Residual, wID, Beta, eta, weight)
      BetaMatr[Step,]   = Beta     
      LamTrace[Step]    = eta         
      # break
    }else{
      cat("Error: EPO is not correct value.!!!\n")
    }        
    Step = Step + 1
  }#while
  
  ##################### Loop Ends  #################################
  BetaMatr <- BetaMatr[1:Step,]
  if(!stand.scale){
    sdx = rep(1, p)
  }
  BetaMatr = BetaMatr %*% diag(1/sdx, p, p)
  
  Beta0 = rep(0, Step)
  for(j in 1:Step){
    Beta0[j] <- my-sum(BetaMatr[j,]*mx)
    #print(sum(y+my-Beta0[j]-oX%*%BetaMatr[j,]))
  }
  
  ###
  return(list(Beta = BetaMatr,
              Beta0 = Beta0,  
              coefficients = c(Beta0[Step],BetaMatr[Step,]),
              fitted.values = drop(X %*% BetaMatr[Step,]),
              LamTrace = LamTrace[1:Step], 
              ObjValueTrace = ObjValTrace[1:Step],
              Step = Step,
              conv = conv))
}
current_status <- function(Step, V, Beta, eta, Sign, wgcorr, Residual, ObjValTrace){
  cat("#######################################################################################################\n")
  cat("############################ Step ",Step, "th Status########################################################\n")
  cat("#######################################################################################################\n")
  cat("V >>> ", V, "Sign[",V,"]", Sign[V],"\n")
  
  cat("Beta[",V,"] >>", Beta[V], "\t L1 norm >>", sum(Sign[V]*Beta[V]), "\n")            
  cat("eta   >> ", eta, "\t wgcorr >> ", wgcorr, "\n")
  cat("resid >>> ", round(Residual,3), "\n")
  cat("ObjValTrace >>> ", round(ObjValTrace[Step],4), "\n")
}
getObjective <- function(res, wID, Beta, eta, weight){
  loss = sum(res^2*wID)/2+eta*sum(weight*abs(Beta))
  return(loss)
}
# right derivatives of d_beta, d_eta
rderiv <- function(XV, wmat, nV, wSignV, eps=1e-8){
  A = t(XV) %*% wmat %*% XV
  b = -wSignV
  tmpqr = qr(A)
  
  if(tmpqr$rank < nV){
    sol = rep(0,nV)
    sol[nV] = 1e5 # return delta which has a big postive value.
    # break
  }
  else{
    sol = array(qr.solve(tmpqr,b)) #gamma_Active
  }        
  return(sol)
}
Find.Event <- function(p, V, Beta, eta, eta0, wgcorr, d_wgcorr, Sign, d_beta, trace=FALSE, eps = 1e-10){
  nV = length(V)
  # 1: Active variable becomes inactive.
  index1 = (Sign[V]>0 & d_beta > eps) | (Sign[V]<0 & d_beta < -eps) 
  d_S1_var = rep(Inf, nV)
  if(sum(index1)>0){
    d_S1_var[index1] = (Beta[V])[index1]/d_beta[index1]
    d_S1_var[d_S1_var<=0] = Inf
  }
  d_S1 = min(d_S1_var)
  
  # 2: Inactive variable joins the active set.
  index2 = rep(FALSE, p)
  index2[-V] = TRUE
  d_lamT1 = (eta+wgcorr[index2])/(1+d_wgcorr[index2])
  d_lamT2 = (eta-wgcorr[index2])/(1-d_wgcorr[index2])
  
  d_lamT1[d_lamT1 <= eps] = Inf
  d_lamT2[d_lamT2 <= eps] = Inf        
  
  d_S2_var = apply(cbind(d_lamT1, d_lamT2), 1, min)
  d_S2_which = apply(cbind(d_lamT1, d_lamT2), 1, which.min) 
  
  if(sum(d_S2_var == -Inf)==p-nV | sum(d_S2_var < -eps) ==p-nV){
    d_S2_var = rep(Inf,p-nV)
  }
  if(sum(d_S2_var==Inf)==p-nV){
    d_S2 = Inf
  }else{   
    d_S2 = min(d_S2_var)
  }
  # 3: The generalized correlation of active variables reduces to zero
  d_S3 = ifelse(eta-eta0<0, Inf, eta-eta0)
  
  ds = min(c(d_S1, d_S2, d_S3))
  
  if(ds==Inf){
    #cat("ds is Inf !!! therefore no further update!!!\n")
  }
  EPO = which.min(c(d_S1, d_S2, d_S3))
  
  if(trace){        
    #if(1){
    cat("d_S1 >>> ", d_S1, "\t")
    cat("d_S2 >>> ", d_S2, "\t")
    cat("d_S3 >>> ", d_S3, "\t")
    cat("Event ", EPO, "occurs!!!\n")            
  }
  return( list(EPO = EPO, 
               ds = ds, 
               d_S1_var = d_S1_var,
               d_S2_var = d_S2_var,
               d_lamT1 = d_lamT1,
               d_lamT2 = d_lamT2
  )
  )
}
