C     
C     #coordinate descent algorithm
C     #ref: Regularization paths for generalized linear models via coordinate descent, Friedman et al., JSS, 2010, 33(1)
C     output: updated
C     beta
C     b0
C     yhat
C     jj
C     input: others including following
C     beta
C     b0
C     mu
C     w
C     z
C     X
C     y

      subroutine loop_glm(x,y,z,n,m,w,mu, penalty,thresh, eps, 
     +     standardize, family, beta, b0, lambda, alpha, gam,
     +     wtold,trace,jj,rescale,converged,theta,pll,activeset, jk)
      implicit none
      integer trace, n, m, standardize,family, jj, i,j,penalty,
     +     converged, rescale, activeset(m), jk, ii
      double precision x(n, m), y(n), thresh, eps,beta(m),beta_old(m) 
      double precision lambda(m), alpha, gam, wtold(n)
      double precision z(n), b0, b0_old, xd(m), yhat(n), w(n),
     +     r(n),xwx,xwr, wtnew(n), wsum,mu(n), theta, pll_old, pll

      wsum = 0
      do i=1, n
         wsum = wsum + wtold(i)
      enddo
      do i=1, n
         if(standardize .EQ. 1)then
            wtnew(i) = wtold(i)/wsum
         else 
            wtnew(i) = wtold(i)
         endif
      enddo
C     current penalized log-likehood value
      if(trace .EQ. 1) then
         call evalpll(y,x,n,m,beta,b0,family, theta, wtnew, alpha,
     +        gam,lambda, penalty, pll)
      else 
         pll=0
      endif     
      pll_old = pll
      jj = 1
      do 100 j = 1, m
         beta_old(j) = beta(j)
 100  continue
      b0_old = b0
      do 10 i = 1, n
         yhat(i) = b0
         do 20 ii = 1, jk
            j = activeset(ii)
            yhat(i) = yhat(i) + x(i,j) * beta(j)
 20      continue         
 10   continue
C     for binomial and poisson distributions
C     family=1 is not really used here
      if(family .EQ. 1)then
         do i = 1, n
            r(i) = y(i) - mu(i)
         enddo
      else 
         if(family .EQ. 2 .OR. family .EQ. 3)then
            do i = 1, n
               r(i) = (y(i) - mu(i))/w(i)
            enddo
         else 
            if(family .EQ. 4)then
C     for negative binomila distribution
               do i = 1, n
                  r(i) = (y(i) - mu(i))/mu(i)
               enddo
            endif
         endif
      endif 
C     intercept 
      xwr =0 
      xwx =0
      do 200 i=1, n
         xwr = xwr + wtnew(i) * w(i) * r(i)
         xwx = xwx + wtnew(i) * w(i)
 200  continue
      b0 = xwr/xwx + b0_old
      do 210 i=1, n
         r(i) = r(i) - (b0 - b0_old)
 210  continue
C     calculate response z
C     cycle through in the active set (begin with if statement)
      do 40 ii = 1, jk
         j = activeset(ii)
         z = 0.0D0
         xwr = 0
         xwx = 0
         do 50 i = 1, n
            xwr = xwr + wtnew(i) * x(i,j) * w(i) * r(i)
            xwx = xwx + wtnew(i) * x(i,j) * w(i) * x(i,j)
 50      continue 
         z = xwr+ xwx * beta_old(j)
         xd(j) = xwx   
         if(penalty.EQ.1) then
            call enet(z,xd(j), lambda(j)*alpha, lambda(j)*(1-alpha),
     +           beta(j))
         else if(penalty.EQ.2) then
            call mcp(z, xd(j), lambda(j)*alpha, lambda(j)*(1-alpha), 
     +           gam, rescale, beta(j))
         else if(penalty.EQ.3) then
            call scad(z, xd(j), lambda(j)*alpha, 
     +           lambda(j)*(1-alpha), gam, rescale, beta(j))
         endif

C     check if beta(j) strictly increases penalized loglikehood function 
         if(trace .EQ. 1)then
            call evalpll(y,x,n,m,beta,b0,family, theta, wtnew, alpha,
     +           gam,lambda, penalty, pll) 
         else 
            pll=0
         endif    
         if(dabs(beta(j)-beta_old(j)) .GT. eps)then
            do 65 i = 1, n
               r(i) = r(i) - x(i, j) * (beta(j) - beta_old(j))
 65         continue
         endif
C     70: cycle through in the active set (end with endif statement) 
 40   continue
      call checkConvergence(m, beta, beta_old, eps, thresh, converged,
     +     activeset, jk)

      return
      end

C     compute penalized log-likelihood value, this is only practical useful if rescale=FALSE
      subroutine evalpll(y,x,n,m,beta,b0,family, theta, weights, alpha,
     +     gam,lambda, penalty, pll)     
      implicit none
      integer penalty, n,m,i,j, family 
      double precision x(n, m), y(n), weights(n), lambda(m),eta(n),
     +     yhat(n), alpha, gam, b0, beta(m), mu(n), theta, pll, ll

      do 220 i = 1, n
         yhat(i) = b0
         do 230 j = 1, m
            yhat(i) = yhat(i) + x(i,j) * beta(j)
 230     continue
 220  continue
      call DCOPY(n, yhat, 1, eta, 1)
      call linkinv(n, eta, family, mu)
      call loglikFor(n, y, mu, theta, weights, family, ll)
C     compute penalty value pll
      call penGLM(beta, m, lambda, alpha, gam,
     +     penalty, pll)
      pll = ll - pll

      return
      end
	  