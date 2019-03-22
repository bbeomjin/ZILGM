      subroutine preprocess(x,y,n,m,weights,family,standardize,normx,xd,
     +     mu)
      implicit none
      integer i,j,n,m, standardize,family
      double precision x(n,m),y(n),weights(n),normx(m), xd(m), 
     +     wtnew(n),meanx(m),mu,wsum,ddot,xtmp(n,m)

      do i=1, n
         do j=1,m
            xtmp(i,j) = x(i,j)
         enddo
      enddo
      wsum=0      
C     compute weighted means sum(weights_i*y_i)
      mu = ddot(n, y, 1, weights, 1)
      do 20 i=1, n
         wsum = wsum + weights(i)
 20   continue
      mu = mu/wsum
      
C###  Center x and y, and scale x, and save the means and sds 
C     for now, except for gaussian family, ignore weights
      do 30 i=1, n
         wtnew(i) = weights(i)/wsum
 30   continue
C     compute weighted column averages meanx = x^(transpose) * wtnew

      call DGEMV('T',n, m, 1.0D0, x, n, wtnew, 1, 0.0D0, meanx, 1)

      if(standardize.EQ.1)then
         do 70 j=1, m
            do 60 i=1, n
               x(i,j) = x(i,j) - meanx(j)
C     deal the weights later
               x(i,j) = dsqrt(wtnew(i)) * x(i,j)
 60         continue
 70      continue
         if(family .EQ. 1)then
            do 80 i=1, n
               y(i) = y(i) - mu
               y(i)= dsqrt(wtnew(i)) * y(i)
 80         continue
         endif
         do 110 j=1, m
            normx(j)=0
            xd(j) =1
            do 100 i=1, n
               normx(j) = normx(j) + x(i,j) * x(i,j)
 100        continue
            normx(j) = dsqrt(normx(j))
 110     continue
      else
         do 120 j=1, m
            xd(j) = 0
            normx(j) = 1
            do 125 i=1, n
               xd(j) = xd(j) + weights(i) * x(i,j) * x(i,j)
 125        continue
 120     continue
      endif
      if(standardize.EQ.1 .AND. family .EQ. 1)then
         do 140 j=1, m
            do 130 i=1, n
               x(i,j) = x(i,j)/normx(j)
 130        continue
 140     continue
      endif
      if(standardize .EQ. 1 .AND. family .NE. 1)then
         do 150 j=1, m
            do 160 i=1, n
               x(i,j) = (xtmp(i,j)-meanx(j))/normx(j)
 160        continue
 150     continue
      endif
      
      return
      end

C     #coordinate descent algorithm
C     #ref: Regularization paths for generalized linear models via coordinate descent, Friedman et al., JSS, 2010, 33(1)
C     output:
C     beta
C     b0
C     mu
C     yhat
C     jj
C     input: others except the above mentioned output
      subroutine lmnetGaus(x, y, n, m, weights, lambda, alpha, 
     +     gam,thresh, maxit, eps, standardize, penalty, xd,beta, b0, 
     +     avg, jj, rescale, converged)

      implicit none
      integer maxit, standardize, penalty, n,m,i,j, jj,converged,
     +     rescale,fullset(m),k,activeset(m),activesetold(m),convact,jk 
      double precision x(n, m), y(n), weights(n), lambda(m),
     +     meanx(m), xd(m), resid(n),
     +     alpha, gam, thresh, eps, avg, wsum, mu(n),
     +     wtnew(n), b0, beta(m), ddot

      do i=1, n
         resid(i) = y(i)
      enddo
      if(standardize .EQ. 1)then
         b0 = 0
      else 
         b0 = avg
      endif
      wsum=0
C     compute weighted means sum(weights_i*y_i)
      mu = ddot(n, y, 1, weights, 1)
      do 20 i=1, n
         wsum = wsum + weights(i)
 20   continue
      do 30 i=1, n
         wtnew(i) = weights(i)/wsum
 30   continue
C     compute weighted column averages meanx = x^(transpose) * wtnew
      call DGEMV('T',n, m, 1.0D0, x, n, wtnew, 1, 0.0D0, meanx, 1)
      do 40 j=1, m
C     activeset(j)=1
C     fullset(j)=1
         activeset(j)=j
         fullset(j)=j
 40   continue
      convact=0
C     some loop, if no change to the active set, stop
      k = 1
      do 2000 while (k .LT. 100 .AND. convact .EQ. 0)
         do 50 j=1, m
            activesetold(j)=activeset(j)
 50      continue
C     set maxit=1, and have a complete cycle through all the variables
         call loop_gaussian(x,y, n,m,penalty,thresh,eps,1,standardize,
     +        beta,b0,resid,xd,lambda,alpha,gam,weights,avg,meanx, 
     +        jj,rescale, converged, fullset, m)
C     determine the active set with only non-zero coefficients 
C     jk: number of variables in active set
         jk = 0
C     it is possible, jk=0 if beta=0, like intercept-only model
         do 60 j=1, m
            if(dabs(beta(j)) .GT. eps)then
               jk=jk+1
               activeset(jk)=j
C     activeset(j)=1
C     else 
C     activeset(j)=0
            endif
 60      continue
C     it is possible, jk=0 if beta=0, like intercept-only model
         if(jk .EQ. 0)then
            convact=1
            exit
         endif
C     check if the active set was changed--begin
         if (k .GT. 1)then
            do 90 j=1, m
               if(activesetold(j) .NE. activeset(j))then
                  convact=0
                  exit
               else 
                  if(j .EQ. m)then
                     convact = 1
                  endif
               endif 
 90         continue
            if(convact .EQ. 1)then
               exit
            endif
         endif
C     check if the active set was changed--end
C     now cycle through only the active set
         call loop_gaussian(x,y, n,m,penalty,thresh,eps,maxit,
     +        standardize,beta,b0,resid,xd,lambda,alpha,gam,weights,
     +        avg,meanx, jj,rescale, converged, activeset, jk)
         k = k + 1
 2000 continue
C     used to match outside loop in subroutine midloop
      jj = jj - 1
C     if(standardize.EQ.1)then
C     do 200 j=1, m
C     beta(j) = beta(j)/normx(j)
C     if(dabs(beta(j)) .LT. 1e-10)then
C     beta(j) = 0
C     endif
C     200 continue
C     endif
CCCC  changed 4/22/2015
C     b0 = 0
C     do 210 j=1, m
C     b0 = b0 + meanx(j) * beta(j)
C     210 continue
C     b0 = avg - b0

      return
      end

!     coordinate descent algorithm

      subroutine enet(z, t, lone, ltwo, res)
      implicit none
      double precision z, t, lone, ltwo, res

      call soth(z, lone, res)
      res= res/(t + ltwo)
      return
      end
      
      subroutine mcp(z, t, lone, ltwo, gam, rescale, res)
      implicit none
      integer rescale
      double precision z, t, v, lone, ltwo, gam, res

      if(rescale .EQ. 1) then
         v = 1
      else
         v = t
      endif
      if(dabs(z) .LE. v*gam * lone * (1+ltwo)) then
         call soth(z, lone, res)                
         if(rescale .EQ. 1) then
            res =res/(t*(1 + ltwo - 1/gam))
         else 
            res =res/(t + ltwo - 1/gam)
         endif
      else 
         if(rescale .EQ. 1) then
            res = z/(t*(1 + ltwo))
         else 
            res = z/(t + ltwo)
         endif
      endif
      return
      end

      subroutine scad(z, t, lone, ltwo, gam, rescale, res)
      implicit none
      integer rescale
      double precision z, t, lone, ltwo, gam, res
      
      if(rescale .EQ. 1) then
         if(dabs(z) .LE. lone   + lone * (1+ltwo)) then
            call soth(z, lone, res)                
C     if(rescale .EQ. 1) then
            res =res/(t*(1 + ltwo))
C     else 
C     res =res/(t + ltwo)
C     endif
         else if(dabs(z) .LE. gam * lone * (1+ltwo)) then
            call soth(z, gam*lone/(gam-1), res)
            res=res/(t*(1-1/(gam-1)+ltwo))
         else 
            res=z/(t*(1+ltwo))
         endif     
      else
         if(dabs(z) .LE. lone   + lone * (t+ltwo)) then
            call soth(z, lone, res)                
            res =res/(t + ltwo)
         else if(dabs(z) .LE. gam * lone * (t+ltwo)) then
            call soth(z, gam*lone/(gam-1), res)
            res=res/(t-1/(gam-1)+ltwo)
         else 
            res=z/(t+ltwo)
         endif
      endif     
      return
      end   

!     input variables
!     x
!     y
!     n
!     m
!     eps
!     maxit
!     threshold
!     standardize
!     beta
!     b0
!     xd
!     lambda
!     alpha
!     gam
!     
!     output variables
!     beta
!     jj       


      subroutine loop_gaussian(x, y, n, m, penalty, thresh, eps, maxit, 
     +     standardize, beta, b0, resid,xd, lambda, alpha, gam,wtold,
     +     avg, meanx, jj, rescale, converged, activeset, jk)
      implicit none
      integer n, m, maxit, standardize, jj, i,j,penalty,converged
      integer rescale, activeset(m), jk, ii
      double precision x(n, m), y(n), thresh, eps,beta(m),beta_old(m) 
      double precision lambda(m), alpha, gam, wtold(n), avg, meanx(m)
      double precision z, b0, xd(m), yhat(n), resid(n),b0_old
      
      jj = 1
      converged = 0
 1000 if(jj .LE. maxit .AND. converged .EQ.0)then
         do 100 ii = 1, jk
            j = activeset(ii)
            beta_old(j) = beta(j)
 100     continue
         b0_old = b0
         do 10 i = 1, n
            yhat(i) = b0
            do 20 ii = 1, jk
               j = activeset(ii)
               yhat(i) = yhat(i) + x(i,j) * beta(j)
 20         continue         
 10      continue
         do 30 i = 1, n
            resid(i) = y(i) - yhat(i)
 30      continue
!     When cycling through the variables, we could restrict ourselves to the current 
!     active set and visit all the variables e.g. every 10th iteration to update 
!     the active set.
         do 40 ii = 1, jk
            j = activeset(ii)
            z = 0.0D0
            if(standardize.EQ.1)then
               do 50 i = 1, n
                  z = z + x(i,j) * resid(i)
 50            continue 
               z = z + beta(j)    
            else 
               do 60 i = 1, n
                  z= z + wtold(i)*x(i,j) * (resid(i) + x(i,j)*beta(j)) 
 60            continue
            endif
            if(penalty.EQ.1) then
               call enet(z, xd(j), lambda(j)*alpha,lambda(j)*(1-alpha),
     +              beta(j))
            else if(penalty.EQ.2) then
               call mcp(z, xd(j), lambda(j)*alpha, lambda(j)*(1-alpha),
     +              gam, rescale, beta(j))
            else if(penalty.EQ.3) then
               call scad(z, xd(j), lambda(j)*alpha, 
     +              lambda(j)*(1-alpha), gam, rescale, beta(j))
            endif
            if(dabs(beta(j) - beta_old(j)) .GT. eps)then
               do 70 i = 1, n
                  resid(i) = resid(i) - x(i, j) * (beta(j)-beta_old(j))
 70            continue
            endif
 40      continue
         if(standardize.EQ.0)then
            b0 = 0.0D0
            do 80 ii = 1, jk
               j = activeset(ii)
               b0 = b0 + meanx(j) * beta(j)
 80         continue
            b0 = avg - b0
         endif
         call checkConvergence(m, beta, beta_old, eps, thresh, 
     +        converged, activeset, jk)
         jj = jj + 1
         goto 1000
      endif

      return
      end

!     soft-threshold=g operator
      subroutine soth(z, g, res)
      implicit none
      double precision res, z, g
      
      if(z .GT. g)then
         res=z-g
      else if (dabs(z) .LE. g) then
         res=0
      else if (z .LT. -g) then
         res=z+g
      endif
      return
      end

      subroutine checkConvergence(m, beta, betaold,eps,thresh,converged,
     +     activeset, jk)
      implicit none
      integer m, j, jk, ii, converged, activeset(m)
      double precision beta(m), betaold(m), eps, thresh

      converged=1
      do 10 ii=1, jk
         j = activeset(ii)
C     20: only check active set      
         if(dabs(beta(j)) .GT. eps .AND. dabs(betaold(j)) .GT. eps) then
            if(dabs((beta(j)-betaold(j))/betaold(j)) .GT. thresh) then
               converged=0
               exit
            endif
         else if(dabs(beta(j)) .LE. eps .AND. dabs(betaold(j)) .GT. eps)
     +           then
            converged=0
            exit
         else if(dabs(beta(j)) .GT. eps .AND. dabs(betaold(j)) .LE. eps)
     +           then
            converged=0
            exit
         endif
C     20: only check active set      
 10   continue

      return
      end
	  