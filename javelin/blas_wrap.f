      subroutine dcopy_wrap(x,y,nx)
cf2py intent(inplace) y
cf2py intent(hide) nx
cf2py threadsafe
      double precision x(nx), y(nx)
      external dcopy
!       dcopy(n,dx,incx,dy,incy)

      CALL dcopy(nx,x,1,y,1)
      RETURN
      END



      SUBROUTINE dtrmm_wrap(M,N,A,B,SIDE,TRANSA,UPLO,ALPHA)
cf2py intent(inplace)::B
cf2py integer intent(hide)::M
cf2py integer intent(hide)::N
cf2py threadsafe

      
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(M,M),B(M,N)
      
      EXTERNAL DTRMM
! DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*     X = alpha*op( A )*B,   or   X = alpha*B*op( A ),
      
      DIAG = 'N'
      LDA = M
      LDB = M

      
      CALL DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      RETURN
      END      


      SUBROUTINE dtrsm_wrap(M,N,A,B,SIDE,TRANSA,UPLO,ALPHA)
cf2py intent(inplace)::B
cf2py integer intent(hide)::M
cf2py integer intent(hide)::N
cf2py threadsafe
      
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(M,M),B(M,N)
      
      EXTERNAL DTRSM
! DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
      
      DIAG = 'N'
      LDA = M
      LDB = M
      
      CALL DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      RETURN
      END


      SUBROUTINE dpotrf_wrap(A,n,info)

cf2py double precision dimension(n,n), intent(inplace) :: A
cf2py integer depend(A),intent(hide)::n=shape(A,0)
cf2py integer intent(out)::info
cf2py threadsafe

      DOUBLE PRECISION A(n,n)
      INTEGER n, info, i, j

      EXTERNAL DPOTRF
! DPOTRF( UPLO, N, A, LDA, INFO ) Cholesky factorization
      
!     C <- cholesky(C)      
      call DPOTRF( 'U', n, A, n, info )
      do i=2,n
        do j=1,i-1
          A(i,j)=0.0D0
        enddo
      enddo
      
      return
      END

      
      SUBROUTINE dpotrf2_wrap(A,n,info)

cf2py double precision dimension(n,n), intent(inplace) :: A
cf2py integer depend(A),intent(hide)::n=shape(A,0)
cf2py integer intent(out)::info
cf2py threadsafe

      DOUBLE PRECISION A(n,n)
      INTEGER n, info, i, j

      EXTERNAL DPOTRF
! DPOTRF( UPLO, N, A, LDA, INFO ) Cholesky factorization
      
!     C <- cholesky(C)      
      call DPOTRF( 'L', n, A, n, info )
      do i=1,n
        !do j=1,i-1
        do j=i+1,n
          A(i,j)=0.0D0
        enddo
      enddo
      
      return
      END


      subroutine dchdc_wrap(a,p,work,piv,info)

cf2py intent(inplace) a
cf2py intent(hide) p
cf2py intent(hide) work       
cf2py intent(out) piv
cf2py intent(out) info
cf2py threadsafe

      integer lda,p,piv(p),job,info
      double precision a(p,p),work(p)      

c
c     dchdc computes the cholesky decomposition of a positive definite
c     matrix.  a pivoting option allows the user to estimate the
c     condition of a positive definite matrix or determine the rank
c     of a positive semidefinite matrix.
c
c     on entry
c
c         a      double precision(lda,p).
c                a contains the matrix whose decomposition is to
c                be computed.  onlt the upper half of a need be stored.
c                the lower part of the array a is not referenced.
c
c         lda    integer.
c                lda is the leading dimension of the array a.
c
c         Since f2py will make sure that a is column-major contiguous,
c         lda is set to p in this wrapper. 
c
c         p      integer.
c                p is the order of the matrix.
c
c         work   double precision.
c                work is a work array.
c
c         piv   integer(p).
c                piv contains integers that control the selection
c                of the pivot elements, if pivoting has been requested.
c                each diagonal element a(k,k)
c                is placed in one of three classes according to the
c                value of piv(k).
c
c                   if piv(k) .gt. 0, then x(k) is an initial
c                                      element.
c
c                   if piv(k) .eq. 0, then x(k) is a free element.
c
c                   if piv(k) .lt. 0, then x(k) is a final element.
c
c                before the decomposition is computed, initial elements
c                are moved by symmetric row and column interchanges to
c                the beginning of the array a and final
c                elements to the end.  both initial and final elements
c                are frozen in place during the computation and only
c                free elements are moved.  at the k-th stage of the
c                reduction, if a(k,k) is occupied by a free element
c                it is interchanged with the largest free element
c                a(l,l) with l .ge. k.  piv is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c        In this wrapper, job=1.
c
c     on return
c
c         a      a contains in its upper half the cholesky factor
c                of the matrix a as it has been permuted by pivoting.
c
c         piv   piv(j) contains the index of the diagonal element
c                of a that was moved into the j-th position,
c                provided pivoting was requested.
c
c         info   contains the index of the last positive diagonal
c                element of the cholesky factor.
c
c     for positive definite matrices info = p is the normal return.
c     for pivoting with positive semidefinite matrices info will
c     in general be less than p.  however, info may be greater than
c     the rank of a, since rounding error can cause an otherwise zero
c     element to be positive. indefinite systems will always cause
c     info to be less than p.
c
c     linpack. this version dated 08/14/78 .
c     j.j. dongarra and g.w. stewart, argonne national laboratory and
c     university of maryland.
c
c
c     blas daxpy,dswap
c     fortran dsqrt
c
c     internal variables
c
      integer pu,pl,plp1,i,j,jp,jt,k,kb,km1,kp1,l,maxl
      double precision temp
      double precision maxdia
      logical swapk,negk
c
      pl = 1
      pu = 0
      info = p
      job=1
      lda=p
      if (job .eq. 0) go to 160
c
c        pivoting has been requested. rearrange the
c        the elements according to piv.
c
         do 70 k = 1, p
            swapk = piv(k) .gt. 0
            negk = piv(k) .lt. 0
            piv(k) = k
            if (negk) piv(k) = -piv(k)
            if (.not.swapk) go to 60
               if (k .eq. pl) go to 50
                  call dswap(pl-1,a(1,k),1,a(1,pl),1)
                  temp = a(k,k)
                  a(k,k) = a(pl,pl)
                  a(pl,pl) = temp
                  plp1 = pl + 1
                  if (p .lt. plp1) go to 40
                  do 30 j = plp1, p
                     if (j .ge. k) go to 10
                        temp = a(pl,j)
                        a(pl,j) = a(j,k)
                        a(j,k) = temp
                     go to 20
   10                continue
                     if (j .eq. k) go to 20
                        temp = a(k,j)
                        a(k,j) = a(pl,j)
                        a(pl,j) = temp
   20                continue
   30             continue
   40             continue
                  piv(k) = piv(pl)
                  piv(pl) = k
   50          continue
               pl = pl + 1
   60       continue
   70    continue
         pu = p
         if (p .lt. pl) go to 150
         do 140 kb = pl, p
            k = p - kb + pl
            if (piv(k) .ge. 0) go to 130
               piv(k) = -piv(k)
               if (pu .eq. k) go to 120
                  call dswap(k-1,a(1,k),1,a(1,pu),1)
                  temp = a(k,k)
                  a(k,k) = a(pu,pu)
                  a(pu,pu) = temp
                  kp1 = k + 1
                  if (p .lt. kp1) go to 110
                  do 100 j = kp1, p
                     if (j .ge. pu) go to 80
                        temp = a(k,j)
                        a(k,j) = a(j,pu)
                        a(j,pu) = temp
                     go to 90
   80                continue
                     if (j .eq. pu) go to 90
                        temp = a(k,j)
                        a(k,j) = a(pu,j)
                        a(pu,j) = temp
   90                continue
  100             continue
  110             continue
                  jt = piv(k)
                  piv(k) = piv(pu)
                  piv(pu) = jt
  120          continue
               pu = pu - 1
  130       continue
  140    continue
  150    continue
  160 continue
      do 270 k = 1, p
c
c        reduction loop.
c
         maxdia = a(k,k)
         kp1 = k + 1
         maxl = k
c
c        determine the pivot element.
c
         if (k .lt. pl .or. k .ge. pu) go to 190
            do 180 l = kp1, pu
               if (a(l,l) .le. maxdia) go to 170
                  maxdia = a(l,l)
                  maxl = l
  170          continue
  180       continue
  190    continue
c
c        quit if the pivot element is not positive.
c
         if (maxdia .gt. 0.0d0) go to 200
            info = k - 1
c     ......exit
            go to 280
  200    continue
         if (k .eq. maxl) go to 210
c
c           start the pivoting and update piv.
c
            km1 = k - 1
            call dswap(km1,a(1,k),1,a(1,maxl),1)
            a(maxl,maxl) = a(k,k)
            a(k,k) = maxdia
            jp = piv(maxl)
            piv(maxl) = piv(k)
            piv(k) = jp
  210    continue
c
c        reduction step. pivoting is contained across the rows.
c
         work(k) = dsqrt(a(k,k))
         a(k,k) = work(k)
         if (p .lt. kp1) go to 260
         do 250 j = kp1, p
            if (k .eq. maxl) go to 240
               if (j .ge. maxl) go to 220
                  temp = a(k,j)
                  a(k,j) = a(j,maxl)
                  a(j,maxl) = temp
               go to 230
  220          continue
               if (j .eq. maxl) go to 230
                  temp = a(k,j)
                  a(k,j) = a(maxl,j)
                  a(maxl,j) = temp
  230          continue
  240       continue
            a(k,j) = a(k,j)/work(k)
            work(j) = a(k,j)
            temp = -a(k,j)
            call daxpy(j-k,temp,work(kp1),1,a(kp1,j),1)
  250    continue
  260    continue
  270 continue
  280 continue
  
C Zero the lower triangle 
      do i=2,info
        do j=1,i-1
          a(i,j) = 0.0D0
        enddo
      enddo

C decrement the pivots.
      do i=1,p
        piv(i) = piv(i)-1
      enddo
      return
      end


      SUBROUTINE dpotrs_wrap(chol_fac, b, info, n, m, uplo)

cf2py double precision dimension(n,n), intent(in) :: chol_fac
cf2py integer depend(chol_fac), intent(hide)::n=shape(chol_fac,0)
cf2py integer depend(b), intent(hide)::m=shape(b,1)
cf2py optional character intent(in):: uplo='U'
cf2py integer intent(out)::info
cf2py double precision intent(inplace), dimension(n,m)::b
cf2py threadsafe

      DOUBLE PRECISION chol_fac(n,n), b(n,m)
      INTEGER n, info,m
      
      CHARACTER uplo

      EXTERNAL DPOTRS
! DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO ) Solves triangular system

      call DPOTRS(uplo,n,m,chol_fac,n,b,n,info)
      
      return
      END
