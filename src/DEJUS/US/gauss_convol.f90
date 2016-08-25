        module gauss_convol

        contains
        function gauss_convol_func(e,spec,ns,sige,ierror,nsigma) result(spec1)

! FUNCTIONAL DESCRIPTION:
! Routine to make a Gaussian convolution of array "spec," which contains ns elements.
! Width of Gaussian is given by "sigma" sige in the same unit as energy array "e."
! Equidistant step size in energy is required.
! One of two convolution methods may be selected by changing the code -- see method 1 and 2 below.
! Method 1 is currently selected which sets the results array to zero at the ends.
! RJD ASD/ANL, June 30, 2014.
! Updated to use same definitions of "conv" and abscissa "x" as in econ.pro and gaus_convol.pro, respectively.
! Differs slightly from the definitions used in the convolution routine used in the code tc.f but either
! result has enough accuracy.
! RJD ASD/ANL, August 11, 2014.
! Removed trailing white spaces and expanded tabs into spaces.
! RJD ASD/ANL, October 19, 2014.
! Updated to use "precision_standard." New name with same declarations to avoid potential
! conflict with Fortran 90 intrinsic function with the name "precision."
! RJD ASD/ANL, November 9, 2014.
! Changed variable name NPPSIGMA to NPPSIGMA_MIN.
! RJD ASD/ANL, November 10, 2014.

        use precision_standard
        implicit none

        integer(kind=i4b),intent(in) :: ns
        real(kind=dp),dimension(ns),intent(in) :: e,spec
        real(kind=dp),intent(in) :: sige
        integer(kind=i4b),intent(out) :: ierror
        integer(kind=i4b),optional,intent(in) :: nsigma
        real(kind=dp),dimension(ns) :: spec1

!  Declarations of scalars:
        integer(kind=i4b) :: ip,ie,nsigma2,np,ne1,ne2
        real(kind=dp) :: conv,sigx,sum_gauss

!  Declarations of arrays:
        real(kind=dp),dimension(:),allocatable :: x,gs,spec2

!  Labeled constants:
        real(kind=dp),parameter :: ZERO=0.0D0
        integer(kind=i4b),parameter :: NSIGMA_DEFAULT=3,NPPSIGMA_MIN=6

        ierror = 0
        if (present(nsigma) .and. nsigma > 0 .and. nsigma < 7) then
          nsigma2 = 2*NSIGMA
        else
          nsigma2 = 2*NSIGMA_DEFAULT
        end if

!  generate Gaussian with correct sigma in units of x-axis (energy axis), equidistant step-size required
!       conv   = e(2) -e(1)                             ! step-size in units of x-axis; result is sensitive to this value
       conv = (maxval(e) -minval(e))/(size(e) -1)      ! use same as in econ.pro

        sigx = sige/conv                        ! sigma in units of x-axis
        if (sigx < NPPSIGMA_MIN -1) goto 900

        np   = dble(nsigma2)*sigx +1.0d0
        if (np/2*2 == np) np = np +1            ! make odd

! allocate storage
        allocate (x(np),gs(np))                 ! Gaussian
        allocate (spec2(ns+np-1))               ! results of convolution; works for both methods

! old style save for future reference; replaced by array constructor
!       do ip=1,np
!         x(ip) = ((ip-1)/dble(np-1) -0.5d0)*dble(np-1)
!       end do

! make Gaussian
!        x = ((/(ip, ip=0,np-1)/)/dble(np-1) -0.5d0)*dble(np-1)
        x = ((/(ip, ip=0,np-1)/)/dble(np-1) -0.5d0)*dble(np)    ! use same as in gaus_convol.pro
        gs = exp(-x**2/2.0d0/sigx**2)           ! whole array
        sum_gauss = sum(gs)
        gs = gs/sum_gauss                       ! normalize

! old style, make convolution using explicit do loops, safe for future reference
!       ne1 = np/2 +1
!       ne2 = ns +1 -ne1
!       out: do ie=ne1,ne2
!         spec2(ie) = ZERO
!         in: do ip=1,np
!           spec2(ie) = spec2(ie) +spec(ie+ne1-ip)*gs(ip)
!         end do in
!       end do out

! new style, make convolution (using forall ... sometimes gives poor performance but we will use it)
        spec1 = ZERO    ! whole array operation
        spec2 = ZERO    ! whole array operation
        ne1 = np/2 +1
        ne2 = ns +1 -ne1

! method 1 (zeros at the end of array spec for half-width of Gaussian)
        forall (ie=ne1:ne2, ip=1:np) spec2(ie) = spec2(ie) +spec(ie+ne1-ip)*gs(ip)      ! compiler warning ok
        forall (ie=ne1:ne2) spec1(ie) = spec2(ie)                                       ! array spec1 has zero at the ends

! method 2 (convolution with the end of array spec non-zero but not complete for half-width of Gaussian)
!        forall (ie=1:ns, ip=1:np) spec2(ie+ip-1) = spec2(ie+ip-1) +spec(ie)*gs(ip)     ! no compiler warning
!        forall (ie=1:ns) spec1(ie) = spec2(ie+ne1-1)                                   ! shift array spec2 and save in array spec1

! debug
!        write (13,fmt='(" ",f12.6,1pe14.6,0p)') (e(ie),spec1(ie), ie=1,ns)
!        write (14,fmt='(" ",f12.6,1p1e18.10,0p)') (x(ip),gs(ip), ip=1,np)

! deallocate storage
        deallocate (x,gs,spec2)
        return

200     format(' ',8a)
213     format(' ',3(a,f12.7))

!+ Error returns
!  -----------------------------------------------------------------------------
900     continue
        print 200,'&gauss_convol_func-F-SIGERR, Sigma error'
        print 200,'- too few data points for Gaussian convolution, check input data'
        print 200,'- rerun with a finer x-grid or make "sige" larger'
        print 213,'- sige:',sige, ', x-grid:',conv, ', number of points/sigma:',sigx
        ierror = -1
        goto 999

999     continue
        return
        end function gauss_convol_func

        end module gauss_convol
