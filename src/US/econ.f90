        module econ

! External routines may be placed here instead of inside the procedure.
! Routines are made private using the private statement (also private if inside the procedure).
!        use precision_standard
!        use gauss_convol
!        use linint

	private        		! Everything is private unless explicitly made public 
	public :: econ_func 	! Declare public routines to be called from outside this module 

        contains

!        subroutine econ_sub(e,spec,ns,sige,ierror,econ_func)
        function econ_func(e,spec,ns,sige,ierror)

! FUNCTIONAL DESCRIPTION:
! Routine for applying the effect of beam energy spread to a synchrotron radiation spectrum  
! (linear dependency on energy).
! The routine follows the "recipe" given in routine econ.pro.
! RJD ASD/ANL, August 11, 2014.
! Checked subroutine procedure instead of function procedure. Either one works fine and we
! will keep the function procedure.
! RJD ASD/ANL, October 19, 2014.
! Updated to use "precision_standard." New name with same declarations to avoid potential
! conflict with Fortran 90 intrinsic function with the name "precision."
! RJD ASD/ANL, November 9, 2014.
! Checked location of use statements. We want all procedures to be private and this can be accomplished
! by placing them in the specification part of the module or inside the procedure. Here we keep them
! inside the procedure.
! RJD ASD/ANL, November 22, 2014.

        use precision_standard
        use gauss_convol
        use linint

        implicit none

        integer(kind=i4b),intent(in) :: ns
        real(kind=dp),dimension(ns),intent(in) :: e,spec
        real(kind=dp),intent(in) :: sige
        real(kind=dp),dimension(ns) :: econ_func
        integer(kind=i4b),intent(out) :: ierror

!  Declarations of scalars:
        integer(kind=i4b) :: i,c1,c2,ng
        real(kind=dp) :: t,q,sigf

!  Declarations of arrays:
        real(kind=dp),dimension(ns) :: r,cf
        real(kind=dp),dimension(:),allocatable :: cg,specg,specc

!  Labeled constants:
        real(kind=dp),parameter :: ZERO=0.0D0
        integer(kind=i4b),parameter :: NPPSIGMA=6

        ierror = 0

        c1 = 1
        c2 = ns
        t  = 1.0d0/log(e(ns)/e(1))
        r  = log(e/e(1))*t
        q  = (c2 -c1)*t                 ! dim. less
        cf = (c2 -c1)*r                 ! new distribution "cf" with unequal step size
        sigf = 2.0d0*sige*q             ! sige = DE/E = constant;  CONSTANT rms width of Gaussian on new distribution (channels)

        if (sigf < 0.01d0) then
          econ_func = spec              ! return original array if small value of sigf
          return
        endif

        ng = ceiling(NPPSIGMA*ns/sigf)  ! number of channels for new distribution with equidistant step size
        if (ng < ns) ng = ns

        allocate (cg(ng),specg(ng),specc(ng))           ! allocate storage

        cg = (/(i, i=0,ng-1)/)*dble(ns -1)/dble(ng-1)   ! generate channel distribution with equidistant step size

        specg = linint_func(cf,spec,ns,cg,ng)   ! new function interpolated onto new distribution with equidistant step size
        specc = gauss_convol_func(cg,specg,ng,sigf,ierror,nsigma=3)     !  convolve with Gaussian with constant width
        econ_func = linint_func(cg,specc,ng,cf,ns)      ! revert to "cf" distribution (variable step size)

! debug
! write (2,fmt='(" ",f15.6,1p1e18.10,0p)') (e(i),econ_func(i), i=1,ns)

        deallocate (cg,specg,specc)                     ! deallocate storage
        end function econ_func
!        end subroutine econ_sub

        end module econ
