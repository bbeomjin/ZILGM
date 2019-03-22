CCC   ONLY used for family=1; compare with midloopGLM for other family
C     Middle loop: Update the quadratic approximation likelihood
C     input:
C     n
C     m
C     x
C     y
C     mu
C     family
C     nlambda
C     lamk
C     alpha
C     gam
C     maxit
C     del
C     thresh
C     wt              used only for family=1
C     trace
C     innermaxit: not used except for family=1      
C     maxit: not used (and maxit=1) for family=1      
C     output:
C     beta
C     b0
C     yhat
C     dev

      subroutine midloop(n,m,x,y,xold,yold,weights, mu, eta, family,
     +     penalty,lamk, 
     +     alpha, gam, theta, rescale, standardize,eps,innermaxit,
     +     maxit, thresh, nulldev, wt, beta, b0,yhat,dev,trace,convmid, 
     +     satu, ep, pll, normx, xd, avg)
      
      implicit none
      integer standardize, trace, penalty, maxit, i, j, jj, nmid, n, 
     +     family, innermaxit, m,converged,convmid, satu, rescale
      double precision x(n,m),y(n), mu(n), z(n), eta(n), wt(n), w(n), 
     +     del,olddev,weights(n),xold(n,m), yold(n),normx(m),xd(m), 
     +     thresh, nulldev, dev, theta, wtw(n),lamk(m),alpha, 
     +     gam, eps, beta(m), betaold(m), b0, b0old, yhat(n),avg, ep, 
     +     pll(maxit)

      innermaxit = maxit
      maxit = 1
      
C      do 5 jj=1, maxit
C         pll(jj)=0
C 5    enddo
      dev = nulldev
      call glmlink(n,mu,family,theta,w, ep)
      call zeval(n, y, eta, mu, w, family, z)

      do 10 i=1, n
         wtw(i)=wt(i) * w(i)
 10   continue
      call preprocess(x, z, n, m, wtw, family, standardize,
     +     normx, xd, avg)
      call lmnetGaus(x, z, n, m, wtw, lamk, alpha, gam, thresh, 
     +     innermaxit, eps, standardize, penalty, xd, 
     +     beta, b0, avg, nmid,rescale, converged)
      do 220 i = 1, n
         yhat(i) = b0
         do 230 j = 1, m
            if(family .EQ. 1)then
               yhat(i) = yhat(i) + xold(i,j) * beta(j)
            else
               yhat(i) = yhat(i) + x(i,j) * beta(j)
            endif
 230     continue
 220  continue
      call DCOPY(n, yhat, 1, eta, 1)
      call linkinv(n, eta, family, mu)
      olddev = dev
C     compute deviance dev
      call deveval(n, yold, mu, theta, weights, family, dev)
      del = dabs(dev - olddev)
      convmid = converged
      if(trace.EQ.1)then
         call dblepr("deviance difference at the end of middle loop "
     +        , -1, del, 1)
      endif
      
      return
      end
	  