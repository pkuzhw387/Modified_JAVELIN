! Last-modified: 25 Feb 2013 01:46:47 AM
! modified by ZHW

MODULE spear_covfunc
implicit none

contains
SUBROUTINE covmat_bit(mat,jd1,jd2,id1,id2,sigma,tau,slagarr,swidarr,scalearr,nx,ny,ncurve,cmin,cmax,symm)
implicit none
!f2py intent(inplace) mat
!f2py intent(in) jd1,jd2,id1,id2
!f2py intent(hide) nx,ny,ncurve
!f2py logical intent(in), optional :: symm=0
!f2py integer intent(in), optional :: cmin=0
!f2py integer intent(in), optional :: cmax=-1
!f2py intent(in) 
!f2py threadsafe
INTEGER(kind=4)  :: nx,ny,ncurve,cmin,cmax
REAL(kind=8), DIMENSION(nx,ny) :: mat
REAL(kind=8), DIMENSION(nx) :: jd1
REAL(kind=8), DIMENSION(ny) :: jd2
INTEGER(kind=4), DIMENSION(nx) :: id1
INTEGER(kind=4), DIMENSION(ny) :: id2
!REAL(kind=8) :: A,gama
REAL(kind=8) :: sigma,tau
REAL(kind=8), DIMENSION(ncurve) :: slagarr,swidarr,scalearr
LOGICAL :: symm
INTEGER(kind=4)  :: i,j
REAL(kind=8) :: slag1,swid1,scale1,slag2,swid2,scale2

if (cmax .eq. -1) then
    cmax = ny
endif

if (symm) then
    do j = cmin+1,cmax
        slag2 = slagarr(id2(j))
        swid2 = swidarr(id2(j))
        scale2=scalearr(id2(j))
        do i=1,j
            slag1 = slagarr(id1(i))
            swid1 = swidarr(id1(i))
            scale1=scalearr(id1(i))
            call covmatij(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),sigma,tau,slag1,swid1,scale1,slag2,swid2,scale2)
        enddo
    enddo
else
    do j = cmin+1,cmax
        slag2 = slagarr(id2(j))
        swid2 = swidarr(id2(j))
        scale2=scalearr(id2(j))
        do i=1,nx
            slag1 = slagarr(id1(i))
            swid1 = swidarr(id1(i))
            scale1=scalearr(id1(i))
            call covmatij(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),sigma,tau,slag1,swid1,scale1,slag2,swid2,scale2)
        enddo
    enddo
endif
return
END SUBROUTINE covmat_bit


SUBROUTINE covmatpmap_bit(mat,jd1,jd2,id1,id2,sigma,tau,slagarr,swidarr,scalearr,&
scale_hidden,model,nx,ny,ncurve,nline,cmin,cmax,symm)
implicit none
!f2py intent(inplace) mat
!f2py intent(in) jd1,jd2,id1,id2
!f2py intent(hide) nx,ny,ncurve,nline
!f2py logical intent(in), optional :: symm=0
!f2py integer intent(in), optional :: cmin=0
!f2py integer intent(in), optional :: cmax=-1
!f2py intent(in) 
!f2py threadsafe
INTEGER(kind=4)  :: nx,ny,ncurve,nline,cmin,cmax,model
REAL(kind=8), DIMENSION(nx,ny) :: mat
REAL(kind=8), DIMENSION(nx) :: jd1
REAL(kind=8), DIMENSION(ny) :: jd2
INTEGER(kind=4), DIMENSION(nx) :: id1
INTEGER(kind=4), DIMENSION(ny) :: id2
REAL(kind=8) :: sigma,tau
! here ncurve is not the actual number --> we count the line band flux as two.
REAL(kind=8), DIMENSION(ncurve) :: slagarr,swidarr,scalearr
REAL(kind=8), DIMENSION(nline) :: scale_hidden
LOGICAL :: symm
INTEGER(kind=4)  :: i,j,nlines
REAL(kind=8) :: slag1,swid1,scale1,slag2,swid2,scale2

if (cmax .eq. -1) then
    cmax = ny
endif

! scale_hidden is the scale of the hidden continua under line bands.

nlines = size(scale_hidden)

! print*, scale_hidden, model

if (symm) then
    do j = cmin+1,cmax
        ! in the single line version, the id can only be 1 or 2, so the lag, wid, scale of line+cont and cont bands can both be accessed using id(i), but in the multi-line version, id == 1 and other cases must be handled differently.
        ! e.g. scalearr = [cont_scale, line1_scale, line_cont1_scale, line2_scale, line_cont2_scale,...].
        ! In this case, we scalearr(1) <-> id == 1; scalearr(2) <-> id == 2; but!!! scalearr(3) != id == 3!!!!
        if (id2(j) .eq. 1) then
            slag2 = slagarr(1)
            swid2 = swidarr(1)
            scale2=scalearr(1)
        else
            slag2 = slagarr(2 * id2(j) - 2)
            swid2 = swidarr(2 * id2(j) - 2)
            scale2=scalearr(2 * id2(j) - 2)
        endif
        do i=1,j
            if (id1(i) .eq. 1) then
                slag1 = slagarr(1)
                swid1 = swidarr(1)
                scale1=scalearr(1)
            else
                slag1 = slagarr(2 * id1(i) - 2)
                swid1 = swidarr(2 * id1(i) - 2)
                scale1=scalearr(2 * id1(i) - 2)
            endif
            ! print*, model, (model .eq. 1)
            if (model .eq. 1) then
                call covmatpmapij(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),&
                    sigma,tau,slag1,swid1,scale1,&
                    &slag2,swid2,scale2,scale_hidden,nlines)
            elseif (model .eq. 2) then
                call covmatpmapij2(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),&
                    sigma,tau,slag1,swid1,scale1,&
                    &slag2,swid2,scale2,scale_hidden,nlines)
            endif
        enddo
    enddo
else
    do j = cmin+1,cmax
        if (id2(j) .eq. 1) then
            slag2 = slagarr(1)
            swid2 = swidarr(1)
            scale2=scalearr(1)
        else
            slag2 = slagarr(2 * id2(j) - 2)
            swid2 = swidarr(2 * id2(j) - 2)
            scale2=scalearr(2 * id2(j) - 2)

        endif
        do i=1,nx
            if (id1(i) .eq. 1) then
                slag1 = slagarr(1)
                swid1 = swidarr(1)
                scale1=scalearr(1)
            else
                slag1 = slagarr(2 * id1(i) - 2)
                swid1 = swidarr(2 * id1(i) - 2)
                scale1=scalearr(2 * id1(i) - 2)
            endif
            if (model .eq. 1) then
                call covmatpmapij(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),&
                    sigma,tau,slag1,swid1,scale1,&
                    &slag2,swid2,scale2,scale_hidden,nlines)
            elseif (model .eq. 2) then
                call covmatpmapij2(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),&
                    sigma,tau,slag1,swid1,scale1,&
                    &slag2,swid2,scale2,scale_hidden,nlines)
            endif
        enddo
    enddo
endif
return
END SUBROUTINE covmatpmap_bit


SUBROUTINE covmatpmap_bit_baldwin(mat,jd1,jd2,id1,id2,sigma,tau,slagarr,swidarr,scalearr,&
scale_hidden,model,nx,ny,ncurve,nline,nepoch,cmin,cmax,symm)
implicit none
!f2py intent(inplace) mat
!f2py intent(in) jd1,jd2,id1,id2
!f2py intent(hide) nx,ny,ncurve,nline,nepoch
!f2py logical intent(in), optional :: symm=0
!f2py integer intent(in), optional :: cmin=0
!f2py integer intent(in), optional :: cmax=-1
!f2py intent(in) 
!f2py threadsafe
INTEGER(kind=4)  :: nx,ny,ncurve,nline,nepoch,cmin,cmax,model
REAL(kind=8), DIMENSION(nx,ny) :: mat
REAL(kind=8), DIMENSION(nx) :: jd1
REAL(kind=8), DIMENSION(ny) :: jd2
INTEGER(kind=4), DIMENSION(nx) :: id1
INTEGER(kind=4), DIMENSION(ny) :: id2
REAL(kind=8) :: sigma,tau
! here ncurve is not the actual number --> we count the line band flux as two.
REAL(kind=8), DIMENSION(ncurve) :: slagarr,swidarr
REAL(kind=8), DIMENSION(nline + 1) :: scale_hidden
REAL(kind=8), DIMENSION(nline, nepoch) :: scalearr
LOGICAL :: symm
INTEGER(kind=4)  :: i,j,nlines,jj
REAL(kind=8) :: slag1,swid1,scale1,slag2,swid2,scale2

if (cmax .eq. -1) then
    cmax = ny
endif

! scale_hidden is the scale of the hidden continua under line bands.

nlines = size(scale_hidden)

! print*, scale_hidden, model

if (symm) then
    do j = cmin+1,cmax
        ! in the single line version, the id can only be 1 or 2, so the lag, wid, scale of line+cont and cont bands can both be accessed using id(i), but in the multi-line version, id == 1 and other cases must be handled differently.
        ! e.g. scalearr = [cont_scale, line1_scale, line_cont1_scale, line2_scale, line_cont2_scale,...].
        ! In this case, we scalearr(1) <-> id == 1; scalearr(2) <-> id == 2; but!!! scalearr(3) != id == 3!!!!
        jj = j / (nlines + 1)
        if (id2(j) .eq. 1) then
            slag2 = slagarr(1)
            swid2 = swidarr(1)
            scale2=scale_hidden(1)
        else
            slag2 = slagarr(2 * id2(j) - 2)
            swid2 = swidarr(2 * id2(j) - 2)
            scale2=scalearr(jj, id2(j) - 1)
        endif
        do i=1,j
            if (id1(i) .eq. 1) then
                slag1 = slagarr(1)
                swid1 = swidarr(1)
                scale1=scale_hidden(1)
            else
                slag1 = slagarr(2 * id1(i) - 2)
                swid1 = swidarr(2 * id1(i) - 2)
                scale1=scalearr(jj, id1(i) - 1)
            endif
            ! print*, model, (model .eq. 1)
            if (model .eq. 1) then
                call covmatpmapij(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),&
                    sigma,tau,slag1,swid1,scale1,&
                    &slag2,swid2,scale2,scale_hidden,nlines)
            elseif (model .eq. 2) then
                call covmatpmapij2(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),&
                    sigma,tau,slag1,swid1,scale1,&
                    &slag2,swid2,scale2,scale_hidden,nlines)
            endif
        enddo
    enddo
else
    do j = cmin+1,cmax
        if (id2(j) .eq. 1) then
            slag2 = slagarr(1)
            swid2 = swidarr(1)
            scale2=scale_hidden(1)
        else
            slag2 = slagarr(2 * id2(j) - 2)
            swid2 = swidarr(2 * id2(j) - 2)
            scale2=scalearr(jj, id2(j) - 1)

        endif
        do i=1,nx
            if (id1(i) .eq. 1) then
                slag1 = slagarr(1)
                swid1 = swidarr(1)
                scale1=scale_hidden(1)
            else
                slag1 = slagarr(2 * id1(i) - 2)
                swid1 = swidarr(2 * id1(i) - 2)
                scale1=scalearr(jj, id1(i) - 1)
            endif
            if (model .eq. 1) then
                call covmatpmapij(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),&
                    sigma,tau,slag1,swid1,scale1,&
                    &slag2,swid2,scale2,scale_hidden,nlines)
            elseif (model .eq. 2) then
                call covmatpmapij2(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),&
                    sigma,tau,slag1,swid1,scale1,&
                    &slag2,swid2,scale2,scale_hidden,nlines)
            endif
        enddo
    enddo
endif
return
END SUBROUTINE covmatpmap_bit_baldwin





SUBROUTINE covmatpmap_bit_single(mat,jd1,jd2,id1,id2,sigma,tau,slagarr,swidarr,scalearr,&
nx,ny,ncurve,cmin,cmax,symm)
implicit none
!f2py intent(inplace) mat
!f2py intent(in) jd1,jd2,id1,id2
!f2py intent(hide) nx,ny,ncurve
!f2py logical intent(in), optional :: symm=0
!f2py integer intent(in), optional :: cmin=0
!f2py integer intent(in), optional :: cmax=-1
!f2py intent(in) 
!f2py threadsafe
INTEGER(kind=4)  :: nx,ny,ncurve,cmin,cmax
REAL(kind=8), DIMENSION(nx,ny) :: mat
REAL(kind=8), DIMENSION(nx) :: jd1
REAL(kind=8), DIMENSION(ny) :: jd2
INTEGER(kind=4), DIMENSION(nx) :: id1
INTEGER(kind=4), DIMENSION(ny) :: id2
REAL(kind=8) :: sigma,tau
! here ncurve is not the actual number --> we count the line band flux as two.
REAL(kind=8), DIMENSION(ncurve) :: slagarr,swidarr,scalearr
LOGICAL :: symm
INTEGER(kind=4)  :: i,j
REAL(kind=8) :: slag1,swid1,scale1,slag2,swid2,scale2,scale_hidden

if (cmax .eq. -1) then
    cmax = ny
endif

! scale_hidden is the scale of the hidden continuum under line band.
scale_hidden =scalearr(3)

if (symm) then
    do j = cmin+1,cmax
        slag2 = slagarr(id2(j))
        swid2 = swidarr(id2(j))
        scale2=scalearr(id2(j))
        do i=1,j
            slag1 = slagarr(id1(i))
            swid1 = swidarr(id1(i))
            scale1=scalearr(id1(i))
            call covmatpmapij_single(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),sigma,tau,&
                slag1,swid1,scale1,slag2,swid2,scale2,scale_hidden)
        enddo
    enddo
else
    do j = cmin+1,cmax
        slag2 = slagarr(id2(j))
        swid2 = swidarr(id2(j))
        scale2=scalearr(id2(j))
        do i=1,nx
            slag1 = slagarr(id1(i))
            swid1 = swidarr(id1(i))
            scale1=scalearr(id1(i))
            call covmatpmapij_single(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),sigma,tau,&
                slag1,swid1,scale1,slag2,swid2,scale2,scale_hidden)
        enddo
    enddo
endif
return
END SUBROUTINE covmatpmap_bit_single

SUBROUTINE covmatpmapij_single(covij,id1,id2,jd1,jd2,sigma,tau,slag1,swid1,scale1,&
slag2,swid2,scale2,scale_hidden)
implicit none
REAL(kind=8),intent(out) :: covij
INTEGER(kind=4),intent(in) :: id1,id2
REAL(kind=8),intent(in)  :: jd1,jd2
REAL(kind=8),intent(in)  :: sigma,tau
REAL(kind=8),intent(in)  :: slag1,swid1,scale1,slag2,swid2,scale2,scale_hidden
REAL(kind=8) :: twidth,twidth1,twidth2
REAL(kind=8) :: tgap,tgap1,tgap2
INTEGER(kind=4) :: imax,imin

imax = max(id1,id2)
imin = min(id1,id2)

! here imin(imax) is either one or two.
if (imin .le. 0) then
    print*,"ids can not be smaller than 1"
    covij = -1.0D0
    return
endif

if (imin .eq. imax) then
    ! between two epochs of the same light curve
    if (imin .eq. 1) then
        ! id1 = id2 = 1
        ! continuum band auto: cov(c0i,c0j)
        covij = exp(-abs(jd1 - jd2) / tau)
    else
        ! id1 = id2 = 2
        ! line band auto: cov(c1i,c1j) + cov(c1i, lj) + cov(li, c1j) + cov(li, lj)
        ! Both (slag1, swid1, scale1) and (slag2, swid2, scale2) are
        ! referring to the line properties.
        covij = scale_hidden * scale_hidden * exp(-abs(jd1 - jd2) / tau) +&
                scale_hidden * scale2 * exp(-abs(jd1 - jd2 + slag2) / tau) +&
                scale_hidden * scale1 * exp(-abs(jd1 - jd2 - slag1) / tau) +&
                scale1 * scale2 * exp(-abs(jd1 - jd2 - slag1 + slag2) / tau)
    endif
else
    ! between two epochs of different light curves 
    !XXX now just the continuum and line bands
    ! continuum band and line band cross cov(c0i, c1j) + cov(c0i, lj)
    ! cov(c0i, c1j)
    covij = scale_hidden * exp(-abs(jd1 - jd2) / tau)
    ! cov(c0i, lj)
    if ((id1.eq.1).and.(id2.ge.2)) then
        covij = covij + scale2 * exp(-abs(jd1 - jd2 + slag2) / tau)
    elseif ((id2.eq.1).and.(id1.ge.2)) then
        covij = covij + scale1 * exp(-abs(jd2 - jd1 + slag1) / tau)
    endif
endif

covij = 0.5D0*sigma*sigma*covij
return
END SUBROUTINE covmatpmapij_single



SUBROUTINE covmat_bit2(mat,jd1,jd2,id1,id2,A,gama,slagarr,swidarr,scalearr,nx,ny,ncurve,cmin,cmax,symm)
implicit none
!f2py intent(inplace) mat
!f2py intent(in) jd1,jd2,id1,id2
!f2py intent(hide) nx,ny,ncurve
!f2py logical intent(in), optional :: symm=0
!f2py integer intent(in), optional :: cmin=0
!f2py integer intent(in), optional :: cmax=-1
!f2py intent(in) 
!f2py threadsafe
INTEGER(kind=4)  :: nx,ny,ncurve,cmin,cmax
REAL(kind=8), DIMENSION(nx,ny) :: mat
REAL(kind=8), DIMENSION(nx) :: jd1
REAL(kind=8), DIMENSION(ny) :: jd2
INTEGER(kind=4), DIMENSION(nx) :: id1
INTEGER(kind=4), DIMENSION(ny) :: id2
REAL(kind=8) :: A,gama
REAL(kind=8), DIMENSION(ncurve) :: slagarr,swidarr,scalearr
LOGICAL :: symm
INTEGER(kind=4)  :: i,j
REAL(kind=8) :: slag1,swid1,scale1,slag2,swid2,scale2

if (cmax .eq. -1) then
    cmax = ny
endif

if (symm) then
    do j = cmin+1,cmax
        slag2 = slagarr(id2(j))
        swid2 = swidarr(id2(j))
        scale2=scalearr(id2(j))
        do i=1,j
            slag1 = slagarr(id1(i))
            swid1 = swidarr(id1(i))
            scale1=scalearr(id1(i))
            call covmatij(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),A,gama,slag1,swid1,scale1,slag2,swid2,scale2)
        enddo
    enddo
else
    do j = cmin+1,cmax
        slag2 = slagarr(id2(j))
        swid2 = swidarr(id2(j))
        scale2=scalearr(id2(j))
        do i=1,nx
            slag1 = slagarr(id1(i))
            swid1 = swidarr(id1(i))
            scale1=scalearr(id1(i))
            call covmatij(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),A,gama,slag1,swid1,scale1,slag2,swid2,scale2)
        enddo
    enddo
endif
return
END SUBROUTINE covmat_bit2




!XXX this is deprecated.
SUBROUTINE covmat(mat,npt,ncurve,idarr,jdarr,sigma,tau,slagarr,swidarr,scalearr)
implicit none
REAL(kind=8), DIMENSION(npt,npt),intent(out)  :: mat
INTEGER(kind=4),intent(in)  :: npt,ncurve
INTEGER(kind=4), DIMENSION(npt),intent(in)  :: idarr
REAL(kind=8), DIMENSION(npt),intent(in)  :: jdarr
REAL(kind=8), DIMENSION(ncurve),intent(in)  :: slagarr,swidarr,scalearr
REAL(kind=8),intent(in)  :: sigma,tau
INTEGER(kind=4) :: m,n,id1,id2
REAL(kind=8) :: jd1,jd2,slag1,swid1,scale1,slag2,swid2,scale2

do m=1,npt
    do n=m,npt
        id1=idarr(m)
        id2=idarr(n)
        jd1=jdarr(m)
        jd2=jdarr(n)
        slag1 =slagarr(id1)
        swid1 =swidarr(id1)
        scale1=scalearr(id1)
        slag2 =slagarr(id2)
        swid2 =swidarr(id2)
        scale2=scalearr(id2)
        call covmatij(mat(m,n), id1,id2,jd1,jd2,sigma,tau,slag1,swid1,scale1,slag2,swid2,scale2)
        if (m .ne. n) then
            mat(n,m) = mat(m,n)
        endif
    enddo
enddo

END SUBROUTINE covmat

SUBROUTINE covmatij(covij,id1,id2,jd1,jd2,sigma,tau,slag1,swid1,scale1,slag2,swid2,scale2)
implicit none
REAL(kind=8),intent(out) :: covij
INTEGER(kind=4),intent(in) :: id1,id2
REAL(kind=8),intent(in)  :: jd1,jd2
REAL(kind=8),intent(in)  :: sigma,tau
REAL(kind=8),intent(in)  :: slag1,swid1,scale1,slag2,swid2,scale2
REAL(kind=8) :: twidth,twidth1,twidth2
REAL(kind=8) :: tgap,tgap1,tgap2
INTEGER(kind=4) :: imax,imin

imax = max(id1,id2)
imin = min(id1,id2)

if (imin .le. 0) then
    print*,"ids can not be smaller than 1"
    covij = -1.0D0
    return
endif

if (imin .eq. imax) then
    ! between two epochs of the same light curve
    if (imin .eq. 1) then
        ! continuum auto
        covij = getcmat_delta(id1,id2,jd1,jd2,tau,slag1,scale1,slag2,scale2)
    else
        ! line auto
        ! Why does the value of swid1 matter when we do this?
        if(swid1 .le. 0.01D0) then
            covij = getcmat_delta(id1,id2,jd1,jd2,tau,slag1,scale1,slag2,scale2)
        else
            covij = getcmat_lauto(id1,jd1,jd2,tau,slag1,swid1,scale1)
        endif
    endif
else
    ! between two epochs of different light curves
    if (imin .eq. 1) then
        ! continuum and line cross
        ! assume swid of the continuum is 0.0
        twidth = max(swid1, swid2)
        if (twidth .le. 0.01D0) then
            covij = getcmat_delta(id1,id2,jd1,jd2,tau,slag1,scale1,slag2,scale2)
        else
            covij = getcmat_lc(id1,id2,jd1,jd2,tau,slag1,swid1,scale1,slag2,swid2,scale2)
        endif
    else 
        ! line1 and line2 cross
        twidth1 = swid1
        twidth2 = swid2
        if((twidth1.le.0.01D0).and.(twidth2.le.0.01D0)) then
            covij = getcmat_delta(id1,id2,jd1,jd2,tau,slag1,scale1,slag2,scale2)
        else if((twidth1 .le. 0.01D0).or.(twidth2 .le. 0.01D0)) then
            covij = getcmat_lc(id1,id2,jd1,jd2,tau,slag1,swid1,scale1,slag2,swid2,scale2)
        else
            covij = getcmat_lcross(id1,id2,jd1,jd2,tau,slag1,swid1,scale1,slag2,swid2,scale2)
        endif
    endif
endif
covij = sigma*sigma*covij
return
END SUBROUTINE covmatij



! pmap version of covmatij, now allows more than two light curves.
!covmatpmapij(mat(i,j), id1(i),id2(j),jd1(i),jd2(j),sigma,tau,slag1,swid1,scale1,slag2,swid2,scale2,scale_hidden)
SUBROUTINE covmatpmapij(covij,id1,id2,jd1,jd2,sigma,tau,slag1,swid1,scale1,&
slag2,swid2,scale2,scale_hidden,nlines)
implicit none
REAL(kind=8),intent(out) :: covij
INTEGER(kind=4),intent(in) :: id1,id2
REAL(kind=8),intent(in)  :: jd1,jd2
REAL(kind=8),intent(in)  :: sigma,tau
REAL(kind=8),intent(in)  :: slag1,swid1,scale1,slag2,swid2,scale2
REAL(kind=8) :: twidth,twidth1,twidth2
REAL(kind=8) :: tgap,tgap1,tgap2
INTEGER(kind=4) :: imax,imin,nlines
REAL(kind=8), DIMENSION(nlines) :: scale_hidden

imax = max(id1,id2)
imin = min(id1,id2)

! here imin(imax) is either one, two,...
if (imin .le. 0) then
    print*,"ids can not be smaller than 1"
    covij = -1.0D0
    return
endif

! the following is the delta transfer function version.

if (imin .eq. 1) then
    if (imax .eq. 1) then
        ! id1 = id2 = 1
        ! continuum auto correlation
        covij = exp(-abs(jd1 - jd2) / tau)
        ! print*,covij
    else
        covij = scale_hidden(imax) * exp(-abs(jd1 - jd2) / tau)
        ! between two epochs of continuum and one of the lines
        if ((id1 .eq. 1) .and. (id2 .ge. 2)) then
            ! print*,scale_hidden(id2 - 1), id1, id2, scale1, scale2
            covij = covij + scale2 * exp(-abs(jd1 - jd2 + slag2) / tau)
            ! print*,covij
        elseif ((id2 .eq. 1) .and. (id1 .ge. 2)) then
!             print*,scale_hidden(id1 - 1), id1, id2, scale1, scale2
            covij = covij + scale1 * exp(-abs(jd2 - jd1 + slag1) / tau)
            ! print*,covij
        endif
    endif

else
!     print*, slag1, slag2
    ! id1, id2 >= 2
    ! line band cross correlation: cov(cmi,cnj) + cov(cmi,lnj) + cov(lmi,cnj) + cov(lmi,lnj)
    covij = scale_hidden(id1) * scale_hidden(id2) * exp(-abs(jd1 - jd2) / tau) +&
            scale_hidden(id1) * scale2 * exp(-abs(jd1 - jd2 + slag2) / tau) +&
            scale_hidden(id2) * scale1 * exp(-abs(jd1 - jd2 - slag1) / tau) +&
            scale1 * scale2 * exp(-abs(jd1 - jd2 - slag1 + slag2) / tau)
endif
! print*, covij
covij = 0.5D0*sigma*sigma*covij
return
END SUBROUTINE covmatpmapij




SUBROUTINE covmatpmapij2(covij,id1,id2,jd1,jd2,A,gama,slag1,swid1,scale1,&
slag2,swid2,scale2,scale_hidden,nlines)
implicit none
REAL(kind=8),intent(out) :: covij
INTEGER(kind=4),intent(in) :: id1,id2
REAL(kind=8),intent(in)  :: jd1,jd2
REAL(kind=8),intent(in)  :: A,gama
REAL(kind=8),intent(in)  :: slag1,swid1,scale1,slag2,swid2,scale2
! index 1 and 2 correspond to id1 and id2.
REAL(kind=8) :: twidth,twidth1,twidth2
REAL(kind=8) :: tgap,tgap1,tgap2
INTEGER(kind=4) :: imax,imin,nlines
REAL(kind=8), DIMENSION(nlines) :: scale_hidden


imax = max(id1,id2)
imin = min(id1,id2)

! here imin(imax) is either one or two.
if (imin .le. 0) then
    print*,"ids can not be smaller than 1"
    covij = -1.0D0
    return
endif


! the following is the delta transfer function version.

if (imin .eq. 1) then
    if (imax .eq. 1) then
        ! id1 = id2 = 1
        ! continuum auto correlation
        covij = (10.0D0)**gama - 0.5D0 * (abs(jd1 - jd2) / 365.25D0)**gama

    else
        ! between two epochs of continuum and one of the lines
        covij = scale_hidden(imax) * ((10.0D0)**gama - 0.5D0 * (abs(jd1 - jd2) / 365.25D0)**gama) 
        if ((id1 .eq. 1) .and. (id2 .ge. 2)) then
            covij = covij + scale2 * ((10.0D0)**gama-(abs(jd1 - jd2 + slag2) / 365.25D0)**gama)
        elseif ((id2 .eq. 1) .and. (id1 .ge. 2)) then
            covij = covij + scale1 * ((10.0D0)**gama-(abs(jd2 - jd1 + slag1) / 365.25D0)**gama)
        endif
    endif
else
    ! id1, id2 >= 2
    ! line band cross correlation: cov(cmi,cnj) + cov(cmi,lnj) + cov(lmi,cnj) + cov(lmi,lnj)
    covij = scale_hidden(id1) * scale_hidden(id2) * ((10.0D0)**gama - 0.5D0 * (abs(jd1 - jd2) / 365.25D0)**gama) +&
            scale_hidden(id1) * scale2 * ((10.0D0)**gama-(abs(jd1 - jd2 + slag2) / 365.25D0)**gama) +&
            scale_hidden(id2) * scale1 * ((10.0D0)**gama-(abs(jd2 - jd1 + slag1) / 365.25D0)**gama) +&
            scale1 * scale2 * ((10.0D0)**gama-(abs(jd2 - jd1 + slag1 - slag2) / 365.25D0)**gama)

endif

covij = A*A*covij
!print*, A, gama, covij
return
END SUBROUTINE covmatpmapij2



FUNCTION expcov(djd,tau)
implicit none
REAL(kind=8) :: expcov,djd,tau
expcov = exp(-dabs(djd)/tau)
return
END FUNCTION expcov

FUNCTION getcmat_delta2(id1,id2,jd1,jd2,gama,tspike1,scale1,tspike2,scale2)
implicit none
REAL(kind=8) :: getcmat_delta2
INTEGER(kind=4) ::  id1,id2
REAL(kind=8) ::  jd1,jd2,tspike1,tspike2,gama
REAL(kind=8) ::  scale1,scale2
getcmat_delta2 = 10D0**(2.0D0*gama)-abs((jd1-jd2-tspike1+tspike2)/365.25D0)**(2.0D0*gama)
getcmat_delta2 = abs(scale1*scale2)*getcmat_delta2
return
END FUNCTION getcmat_delta2



FUNCTION getcmat_delta(id1,id2,jd1,jd2,tau,tspike1,scale1,tspike2,scale2)
implicit none
REAL(kind=8) :: getcmat_delta
INTEGER(kind=4) ::  id1,id2
REAL(kind=8) ::  jd1,jd2,tspike1,tspike2,tau
REAL(kind=8) ::  scale1,scale2
getcmat_delta = exp(-abs(jd1-jd2-tspike1+tspike2)/tau)
getcmat_delta = abs(scale1*scale2)*getcmat_delta
return
END FUNCTION getcmat_delta





FUNCTION getcmat_lc(id1,id2,jd1,jd2,tau,slag1,swid1,scale1,slag2,swid2,scale2)
implicit none
REAL(kind=8) :: getcmat_lc
INTEGER(kind=4) ::  id1,id2
REAL(kind=8) ::  jd1,jd2,tau
REAL(kind=8) ::  twidth,emscale,tlow,thig
REAL(kind=8) ::  slag1,swid1,scale1,slag2,swid2,scale2

!1 is the continuum, 2 is line or line-continuum
!t2 = slag2 - 0.5swid2 t1 = slag2 + 0.5swid2
if ((id1.eq.1).and.(id2.ge.2)) then
    tlow = jd2-jd1-slag2-0.5D0*swid2
    thig = jd2-jd1-slag2+0.5D0*swid2
    emscale = dabs(scale2)
    twidth  = swid2
else if ((id2.eq.1).and.(id1.ge.2)) then
    tlow = jd1-jd2-slag1-0.5D0*swid1
    thig = jd1-jd2-slag1+0.5D0*swid1
    emscale = dabs(scale1)
    twidth  = swid1
    ! XXX the following code is inherited from the old spear code where the
    ! DOUBLE_HAT mode used to call getcmat_lc from the cross-correlation
    ! between two lines with one of their widths zero (DEPRECATED, but no harm
    ! if kept).
else if((id1.ge.2).and.(id2.ge.2)) then
    if (swid1.le.0.01D0) then
        tlow = jd2-(jd1-slag1)-slag2-0.5D0*swid2
        thig = jd2-(jd1-slag1)-slag2+0.5D0*swid2
        emscale = dabs(scale2*scale1)
        twidth  = swid2
    else if(swid2.le.0.01D0) then
        tlow = jd1-(jd2-slag2)-slag1-0.5D0*swid1
        thig = jd1-(jd2-slag2)-slag1+0.5D0*swid1
        emscale = dabs(scale2*scale1)
        twidth  = swid1
    endif
endif
if (thig.le.0.0D0) then
    getcmat_lc = exp( thig/tau)-exp( tlow/tau)
else if (tlow.ge.0.0D0) then
    getcmat_lc = exp(-tlow/tau)-exp(-thig/tau)
else 
    getcmat_lc = 2.0D0-exp(tlow/tau)-exp(-thig/tau)
endif
getcmat_lc = tau*(emscale/twidth)*getcmat_lc
return
END FUNCTION getcmat_lc




FUNCTION getcmat_lc2(id1,id2,jd1,jd2,gama,slag1,swid1,scale1,slag2,swid2,scale2)
implicit none
REAL(kind=8) :: getcmat_lc2
INTEGER(kind=4) ::  id1,id2
REAL(kind=8) ::  jd1,jd2,gama
REAL(kind=8) ::  twidth,emscale,tlow,thig
REAL(kind=8) ::  slag1,swid1,scale1,slag2,swid2,scale2

!1 is the continuum, 2 is line or line-continuum
!t2 = slag2 - 0.5swid2 t1 = slag2 + 0.5swid2
if ((id1.eq.1).and.(id2.ge.2)) then
    tlow = jd2-jd1-slag2-0.5D0*swid2
    thig = jd2-jd1-slag2+0.5D0*swid2
    emscale = dabs(scale2)
    twidth  = swid2
else if ((id2.eq.1).and.(id1.ge.2)) then
    tlow = jd1-jd2-slag1-0.5D0*swid1
    thig = jd1-jd2-slag1+0.5D0*swid1
    emscale = dabs(scale1)
    twidth  = swid1
    ! XXX the following code is inherited from the old spear code where the
    ! DOUBLE_HAT mode used to call getcmat_lc from the cross-correlation
    ! between two lines with one of their widths zero (DEPRECATED, but no harm
    ! if kept).
else if((id1.ge.2).and.(id2.ge.2)) then
    if (swid1.le.0.01D0) then
        tlow = jd2-(jd1-slag1)-slag2-0.5D0*swid2
        thig = jd2-(jd1-slag1)-slag2+0.5D0*swid2
        emscale = dabs(scale2*scale1)
        twidth  = swid2
    else if(swid2.le.0.01D0) then
        tlow = jd1-(jd2-slag2)-slag1-0.5D0*swid1
        thig = jd1-(jd2-slag2)-slag1+0.5D0*swid1
        emscale = dabs(scale2*scale1)
        twidth  = swid1
    endif
endif
if (thig.le.0.0D0) then
    getcmat_lc2 = 10D0**(2.0D0*gama)*twidth + (1/(1+2.0D0*gama))*((-thig)**(1+2.0D0*gama) -&
 (-tlow)**(1+2.0D0*gama))/(365.25D0)**(2.0D0*gama)
else if (tlow.ge.0.0D0) then
    getcmat_lc2 = 10D0**(2.0D0*gama)*twidth + (1/(1+2.0D0*gama))*((tlow)**(1+2.0D0*gama) - &
(thig)**(1+2.0D0*gama))/(365.25D0)**(2.0D0*gama)
else 
    getcmat_lc2 = 10D0**(2.0D0*gama)*twidth + (1/(1+2.0D0*gama))*(-(-tlow)**(1+2.0D0*gama) - &
(thig)**(1+2.0D0*gama))/(365.25D0)**(2.0D0*gama)
endif
getcmat_lc2 = (emscale/twidth)*getcmat_lc2
return
END FUNCTION getcmat_lc2





FUNCTION getcmat_lauto(id,jd1,jd2,tau,slag,swid,scale)
implicit none
REAL(kind=8) :: getcmat_lauto
INTEGER(kind=4) :: id
REAL(kind=8) :: jd1,jd2,tau
REAL(kind=8) :: slag,swid,scale
REAL(kind=8) :: twidth,emscale,tlow,tmid,thig

twidth  = swid
emscale = scale

tlow = jd1-jd2-twidth
tmid = jd1-jd2
thig = jd1-jd2+twidth

if((thig.le.0.0D0).or.(tlow.ge.0.0D0))then
    getcmat_lauto = exp(-abs(tmid)/tau)*(exp(0.5D0*twidth/tau)-exp(-0.5D0*twidth/tau))**2
else
    getcmat_lauto = -2.0D0*exp(-abs(tmid)/tau)+exp(-twidth/tau)*(exp(-tmid/tau)+exp(tmid/tau))
    if(tmid.ge.0.0D0)then
        getcmat_lauto = getcmat_lauto+2.0D0*(twidth-tmid)/tau
    else
        getcmat_lauto = getcmat_lauto+2.0D0*(twidth+tmid)/tau
    endif
endif
getcmat_lauto = ((emscale*tau/twidth)**2)*getcmat_lauto
return
END FUNCTION getcmat_lauto



FUNCTION getcmat_lauto2(id,jd1,jd2,gama,slag,swid,scale)
implicit none
REAL(kind=8) :: getcmat_lauto2
INTEGER(kind=4) :: id
REAL(kind=8) :: jd1,jd2,gama
REAL(kind=8) :: slag,swid,scale
REAL(kind=8) :: twidth,emscale,tlow,tmid,thig

twidth  = swid
emscale = scale

tlow = jd1-jd2-twidth
tmid = jd1-jd2
thig = jd1-jd2+twidth

if((thig.le.0.0D0).or.(tlow.ge.0.0D0))then
    getcmat_lauto2 = 10D0**(2.0D0*gama)*twidth**2.0D0 + 1/((1.0D0+2.0D0*gama)*(2.0D0+2.0D0*gama)*365.25D0**(2.0D0*gama))*&
(2.0D0*abs(tmid)**(2.0D0*gama+2.0D0)-abs(tlow)**(2.0D0*gama+2.0D0)-abs(thig)**(2.0D0*gama+2.0D0))
else
    if(tmid.ge.0.0D0)then
        getcmat_lauto2 = 10D0**(2.0D0*gama)*twidth**2.0D0 - 1/((1.0D0+2.0D0*gama)*(2.0D0+2.0D0*gama)*365.25D0**(2.0D0*gama))*&
((thig)**(2.0D0*gama+2.0D0) - 2.0D0*(twidth)**(2.0D0*gama+2.0D0) + (-tlow)**(2.0D0*gama+2.0D0))
    else
        getcmat_lauto2 = 10D0**(2.0D0*gama)*twidth**2.0D0 - 1/((1.0D0+2.0D0*gama)*(2.0D0+2.0D0*gama)*365.25D0**(2.0D0*gama))*&
((thig)**(2.0D0*gama+2.0D0) + (-tlow)**(2.0D0*gama+2.0D0))
    endif
endif
getcmat_lauto2 = ((emscale/twidth)**2)*getcmat_lauto2
return
END FUNCTION getcmat_lauto2

FUNCTION getcmat_lcross(id1,id2,jd1,jd2,tau,slag1,swid1,scale1,slag2,swid2,scale2)
implicit none
REAL(kind=8) :: getcmat_lcross
INTEGER(kind=4) ::  id1,id2
REAL(kind=8) ::  jd1,jd2,tau
REAL(kind=8) ::  slag1,swid1,scale1,slag2,swid2,scale2
REAL(kind=8) ::  twidth1,twidth2,bottleneck
REAL(kind=8) :: t1,t2,t3,t4,ti,tj,tlow,tmid1,tmid2,thig

twidth1 = swid1
twidth2 = swid2

if(twidth1.ge.twidth2) then
    t1 = slag1-0.5D0*twidth1
    t2 = slag1+0.5D0*twidth1
    t3 = slag2-0.5D0*twidth2
    t4 = slag2+0.5D0*twidth2
    bottleneck = twidth2
    ti = jd1
    tj = jd2
else
    t1 = slag2-0.5D0*twidth2
    t2 = slag2+0.5D0*twidth2
    t3 = slag1-0.5D0*twidth1
    t4 = slag1+0.5D0*twidth1
    bottleneck = twidth1
    ti = jd2
    tj = jd1
endif

tlow  = (ti-tj)-(t2-t3)
tmid1 = (ti-tj)-(t2-t4)
tmid2 = (ti-tj)-(t1-t3)
thig  = (ti-tj)-(t1-t4)

if((thig.le.0.0D0).or.(tlow.ge.0.0D0)) then
    getcmat_lcross = dexp(-dabs(tlow)/tau) +dexp(-dabs(thig)/tau)&
                    -dexp(-dabs(tmid1)/tau)-dexp(-dabs(tmid2)/tau)
else 
    getcmat_lcross = dexp(tlow/tau)+dexp(-thig/tau)&
                    -dexp(-dabs(tmid1)/tau)-dexp(-dabs(tmid2)/tau)
    if(tmid2.le.0.0D0) then
        getcmat_lcross = getcmat_lcross+2.0D0*thig/tau
    else if(tmid1.le.0.0D0) then
        getcmat_lcross = getcmat_lcross+2.0D0*bottleneck/tau
    else if(tlow .lt.0.0D0) then
        getcmat_lcross = getcmat_lcross-2.0D0*tlow/tau
    endif
endif

getcmat_lcross = (tau*tau*scale1*scale2/(twidth1*twidth2))&
                 *getcmat_lcross
RETURN
END FUNCTION getcmat_lcross

END MODULE spear_covfunc