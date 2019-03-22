C     outer loop: decrement lambda
C     input:
C     mu: a vector
C     output
C     b: beta
C     bz: b0
C     resdev: residual deviance
C     ypre: yhat

      subroutine outloop(x,y,weights, wt, n,m, penalty, nlambda, lam, 
     +     alpha,gam,theta,rescale,mu,eta, family,standardize, nulldev, 
     +     thresh, maxit, innermaxit, eps, trace, start, startv,b, bz,
     +     resdev,ypre, convout, satu, good, ep, outpll)

      implicit none
      integer n,m,i,j,k,kk,penalty, nlambda, family, standardize, maxit,
     +     innermaxit, trace,convmid,convout(nlambda), startv, rescale, 
     +     satu,good,convact,fullset(m),activeset(m),activesetold(m),jk
      double precision x(n,m), y(n), wt(n), lam(m, nlambda),alpha,
     +     gam, theta,mu(n),eta(n), nulldev,thresh, eps, b(m,nlambda),
     +     bz(nlambda),xold(n,m), yold(n), start(m+1), resdev(nlambda), 
     +     v(n), ypre(n,nlambda), lamk(m), beta(m), b0,dev,
     +     weights(n),yhat(n),ep, pll(maxit), outpll(maxit, nlambda),
     +     normx(m),xd(m),avg 

      if(family .NE. 1)then
         call preprocess(x, y, n, m, weights, family, standardize,
     +        normx, xd, avg)
      endif
C     keep a copy of x and y in case x and y will be changed in subroutine lmnet if standardize = 1
      do 100 j=1, m
         do 110 i=1, n
            xold(i,j)=x(i,j)
 110     continue
 100  continue
      call DCOPY(n, y, 1, yold, 1)
      if(startv .EQ. 0)then
         b0 = eta(1)
         do j=1, m
            beta(j) = 0
         enddo
      else 
         b0 = start(1)
         do j=1, m
            beta(j) = start(j+1)
         enddo
      endif

      k = 1
      satu = 0

 1000 if(k .LE. nlambda .AND. satu .EQ. 0)then
         if(trace.EQ.1)then
            call dblepr("", -1, 1,0)
            call dblepr("Outer loop: Decrement lambda", -1, 1,0)
            call intpr("  lambda iteration", -1, k, 1)
            call dblepr("  lambda value", -1, lam(1,k), 1)
         endif
         do 10 j=1,m
            lamk(j) = lam(j,k)
 10      continue
C     active set here only applies for family!=1. For family=1, active set
C     is implemented in lmnetGaus
C     70: begin if the block 
         if(family .EQ. 1)then
            call midloop(n,m,x,y, xold,yold,weights,mu, eta, family, 
     +           penalty,lamk,alpha,gam,theta,rescale,standardize, eps,
     +           innermaxit, maxit, thresh, nulldev, wt, beta, b0, yhat,
     +           dev, trace, convmid,satu,ep,pll,normx,xd,avg)
         else 
            do 401 j=1, m
               activeset(j)=j
               fullset(j)=j
 401        continue
            convact=0
C     some loop, if no change to the active set, stop
            kk = 1
            do 2000 while (kk .LT. 100 .AND. convact .EQ. 0)
            do 501 j=1, m
                  activesetold(j)=activeset(j)
 501           continue
C     set maxit=1, and have a complete cycle through all the variables
               call midloopGLM(n,m,x,y,xold,yold,weights,mu,eta,family, 
     +              penalty,lamk,alpha,gam,theta,rescale,standardize,
     +              eps,innermaxit,1,thresh,nulldev,wt,beta,b0,yhat, 
     +              dev,trace,convmid,satu,ep,pll,normx,xd,avg,fullset,
     +              m)
C     determine the active set with only non-zero coefficients 
C     jk: number of variables in active set
               jk = 0
               do 601 j=1, m
                  if(dabs(beta(j)) .GT. eps)then
                     jk=jk+1
                     activeset(jk)=j
                  endif
 601           continue
C     it is possible, jk=0 if beta=0, like intercept-only model for
C     large lambda value
               if(jk .EQ. 0)then
                  convact=1
                  exit
               endif
C     check if the active set was changed--begin
               if(kk .GT. 1)then
                  do 901 j=1, m
                     if(activesetold(j) .NE. activeset(j))then
                        exit
                     endif
                       if(j .EQ. m)then
                               convact = 1
                       endif
 901              continue
                  if(convact .EQ. 1)then
                     exit
                  endif
               endif
C     check if the active set was changed--end
C     now cycle through only the active set
               call midloopGLM(n,m,x,y,xold,yold,weights,mu,eta,family, 
     +              penalty,lamk,alpha,gam,theta,rescale,standardize,
     +              eps,innermaxit,maxit,thresh,nulldev,wt,beta,b0,yhat,
     +              dev,trace,convmid,satu,ep,pll,normx,xd,avg,
     +              activeset, jk)
               kk=kk+1
 2000       continue
         endif
C     70: end if the block
         if(satu .EQ. 1)then
            good = k - 1
         endif
         convout(k)=convmid
         if(family .NE. 1)then
            do 15 i=1, maxit
               outpll(i,k) = pll(i)
 15         continue
         endif
         
         do 20 j=1, m
            b(j, k) = beta(j)
 20      continue
         bz(k) = b0
         resdev(k) = dev
         call linkinv(n, yhat, family, v)
         do 30 i=1, n
            ypre(i,k) = v(i)
 30      continue
         k = k + 1
         if(k .LE. nlambda .AND. satu .EQ. 0)then
            do 40 j=1,m
               b(j, k) = b(j, k-1)
 40         continue
         endif
         goto 1000
      endif
      
      return
      end  