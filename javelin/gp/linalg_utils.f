! Copyright (c) Anand Patil, 2007

! TODO: Multithread zeroing the lower triangle, do it straight from Covariance.

      SUBROUTINE remove_duplicates(x,N,Nd,Nr,rf,rt,Nu,xu,ui)
cf2py intent(out) rt
cf2py intent(out) rf
cf2py intent(out) Nr
cf2py intent(out) xu
cf2py intent(out) Nu
cf2py intent(out) ui
cf2py intent(hide) N
cf2py intent(hide) Nd
cf2py threadsafe
      INTEGER N, Nd, i, j, k, rf(N), rt(N), Nr, ui(N)
      INTEGER Nu
      DOUBLE PRECISION x(N,Nd), xu(N,Nd)
      LOGICAL match
      
      Nr = 0
      Nu = 1
      do k=1,Nd
        xu(1,k) = x(1,k)
      end do
      ui(1)=0
      match=.FALSE.
      do i=2,N
        do j=1,i-1
          match=.TRUE.
          do k=1,Nd
            if(x(i,k).NE.x(j,k)) then
              match=.FALSE.
              go to 10
            end if
          end do
   10   if (match) then
          Nr=Nr+1
          rt(Nr)=i-1
          rf(Nr)=j-1
          go to 20
        end if
        end do
   20   if (.NOT.match) then
          Nu=Nu+1
          ui(Nu)=i-1
          do k=1,Nd
            xu(Nu,k)=x(i,k)
          end do
        end if
      end do
          
      RETURN
      END


      SUBROUTINE check_repeats(x, x_sofar, f_sofar, N, N_dim, N_sofar, 
     +f, new_indices, N_new_indices)
cf2py double precision dimension(N,N_dim), intent(in) :: x
cf2py double precision dimension(N_sofar, N_dim), intent(in) :: x_sofar
cf2py double precision dimension(N_sofar), intent(in) :: f_sofar
cf2py integer intent(hide), depend(x):: N = shape(x,0)
cf2py integer intent(hide), depend(x):: N_dim = shape(x,1)
cf2py integer intent(hide), depend(x_sofar):: N_sofar = shape(x_sofar,0)
cf2py double precision intent(out), dimension(N) :: f
cf2py integer intent(out), dimension(N):: new_indices
cf2py integer intent(out):: N_new_indices
cf2py threadsafe
      INTEGER N, N_dim, N_sofar,
     +N_new_indices, new_indices(N), i, j, k
      DOUBLE PRECISION x(N,N_dim), x_sofar(N_sofar, N_dim),
     +f(N), f_sofar(N_sofar)
      LOGICAL match
      
      N_new_indices = 0
      match=.FALSE.
      do i=1,N
!         N_new_indices = N_new_indices + 1
!         new_indices(N_new_indices) = i-1
        do j=1,N_sofar

          match=.TRUE.

          do k=1,N_dim
            if (x(i,k) .NE. x_sofar(j,k)) then
              match=.FALSE.
              GO TO 10
            endif
          enddo

   10     continue  
        
          if (match) then
            GO TO 20
          endif        
        enddo

   20   continue
   
        if (match) then
          f(i) = f_sofar(j)
        else
          N_new_indices = N_new_indices+1
          new_indices(N_new_indices) = i-1
        endif

      enddo
          
      RETURN
      END


      subroutine diag_call(x,n,ndim,V,cov_fun)
cf2py intent(hide) n
cf2py intent(hide) ndim
cf2py intent(out) V
        integer i, n, ndim, j
        double precision x(n,ndim), xe(1,ndim), V(n)
        external cov_fun
cf2py double precision q
cf2py q = cov_fun(xe,ndim)
        
        do i=1,n
             do j=1,ndim
                 xe(1,j) = x(i,j)
             enddo

            V(i) = cov_fun(xe,ndim)
!             print *,xe,cov_fun(xe,ndim)
        enddo
        
        return
        END 

      SUBROUTINE basis_diag_call(basis_x, V, n, nbas)
cf2py intent(hide) n
cf2py intent(hide) nbas
cf2py intent(out) V
        integer n, i, j
        double precision V(n), basis_x(nbas,n)
        
        do i=1,n
            V(i) = 0
            do j=1,nbas
                V(i) = V(i) + basis_x(j,i) ** 2
            enddo
        enddo
        
        return
        END
        
      SUBROUTINE gp_array_logp(x, mu, sig, n, like, info)

cf2py intent(copy) x, mu
cf2py intent(in) sig
cf2py intent(out) like
cf2py intent(hide) info, n
cf2py threadsafe

      DOUBLE PRECISION sig(n,n), x(n), mu(n), like
      INTEGER n, info, i
      DOUBLE PRECISION twopi_N, log_detC, gd
      DOUBLE PRECISION infinity
      PARAMETER (infinity = 1.7976931348623157d308)      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793238462643d0)

      EXTERNAL DTRSV
! DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
!       EXTERNAL DPOTRS
! DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO ) Solves triangular system
      EXTERNAL DAXPY
! DAXPY(N,DA,DX,INCX,DY,INCY) Adding vectors
      EXTERNAL DCOPY
! DCOPY(N,DX,INCX,DY,INCY) copies x to y
! NB DDOT from ATLAS, compiled with gfortran 4.2 on Ubuntu Gutsy,
! was producing bad output- hence the manual dot product.
      
!     x <- (x-mu)      
      call DAXPY(n, -1.0D0, mu, 1, x, 1)

!       mu <- x
!       call DCOPY(n,x,1,mu,1)
      
!     x <- sig ^-1 * x
!       call DPOTRS('L',n,1,sig,n,x,n,info)
      call DTRSV('U','T','N',n,sig,n,x,1)

      gd=0.0D0
      do i=1,n
          gd=gd+x(i)*x(i)
      end do
      
!     like <- .5 dot(x,mu) (.5 (x-mu) C^{-1} (x-mu)^T)
      like = -0.5D0 * gd
!       print *, like
      
      twopi_N = 0.5D0 * N * dlog(2.0D0*PI)
!       print *, twopi_N
      
      log_detC = 0.0D0
      do i=1,n
        log_detC = log_detC + log(sig(i,i))
      enddo
!       print *, log_detC
      
      like = like - twopi_N - log_detC
      
      return
      END


c
      SUBROUTINE asqs(C,S,nx,ny,cmin,cmax)

cf2py intent(in) C
cf2py intent(inplace) S
cf2py integer intent(in), optional :: cmin = 0
cf2py integer intent(in), optional :: cmax = -1
cf2py intent(hide) nx,ny
cf2py threadsafe

      DOUBLE PRECISION C(nx,ny), cn, S(ny)
      INTEGER nx, ny, i, j, cmin, cmax

      EXTERNAL DSCAL

      if (cmax.EQ.-1) then
          cmax = ny
      end if


        do j=cmin+1,cmax
            S(j) = 0.0D0
            do i=1,nx
                cn = C(i,j)
                S(j) = S(j) + cn * cn
            end do
 !          CALL DSCAL(nx,a,C(1,j),1)
        enddo


      RETURN
      END