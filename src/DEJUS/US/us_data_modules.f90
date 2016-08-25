! US Data Modules
! Latest update: Sun Nov 23 12:17:10 CST 2014 

  module us_size
    use precision_standard
    implicit none

!  Array sizes
!  E_SZ: max size of energy related arrays (E_SZ+1 max size used internally for method = 4)
!  PHI_SZ: max size of arrays for angle phi (for 1st quadrant; full range of 2 pi uses size 4*PHI_SZ)
!  ALPHA_SZ: max size of arrays for angle alpha (= gamma*theta)
!  Brightness array of index phi uses max PHI_SZ+1 elements and of index alpha uses max ALPHA_SZ elements
!  P_SZ-1 max size of pinhole (aperture) arrays (same in both x and y)
!  PW_SZ: max size of power density function arrays (same in both x and y) and will be reset to this value
!  if the internally generated number of points exceeds this value. Not necessary to give warnings to the users.
!  Array sizes are checked inside the code us.f90 and warnings are given if arrays are dynamically allocated
!  beyond those limits. There are no memory limitations imposed by the code.

    integer(kind=i4b),parameter :: E_SZ=100001,PHI_SZ=200,ALPHA_SZ=200,P_SZ=201,PW_SZ=1000

  end module us_size

  module prt
    use precision_standard
    implicit none
    save

    !  Declarations of scalars:
    logical(kind=4) :: LANG
    real(kind=dp) :: D,EMINU,EMAXU,XPC,YPC,XPS,YPS
  end module prt

  module calc
    use precision_standard
    implicit none
    save

    !  Declarations of scalars:
    integer(kind=i4b) :: MODE,METHOD,IHARM,NPER,NSIGMA
    real(kind=dp) :: K,KX,KY,K3,NPI,PE,GAMMA_US,ER,EW,DEW,LDEV
  end module calc

  module factor
    use precision_standard
    implicit none
    save

    !  Declarations of scalars:
    real(kind=dp) :: FAC,C1,C2,C3,C4,C5
    real(kind=dp),parameter :: BW=1.0D-3 ! 0.1%
  end module factor

  module beam
    use precision_standard
    implicit none
    save

    !  Declarations of scalars:
    real(kind=dp) :: AP2MIN,AP2MAX,AP2CNT,ARGMAX,FU,FV,SIGX2,SIGY2,SIGE
    real(kind=dp) :: SIGU2,SIGV2,SIGU,SIGV
  end module beam

  module pinhole
    use precision_standard
    use us_size
    implicit none
    save

    !  Declarations of scalars:
    integer(kind=i4b) :: NXP,NYP
    real(kind=dp) :: DXP,DYP

    ! Declarations of arrays:
    real(kind=dp),dimension(:),allocatable :: XP,YP,CX,CY		! max size P_SZ
  end module pinhole

  module angle_phi
    use precision_standard
    use us_size
    implicit none
    save

    !  Declarations of scalars:
    integer(kind=i4b) :: NPHI4,NPHI,NPHI_BRIGHT
    real(kind=dp) :: DPHI

    ! Declarations of arrays:
    integer(kind=i4b),dimension(:),allocatable :: INDEX_PHI		! max size 4*PHI_SZ
    real(kind=dp),dimension(:),allocatable :: COSPHI,SINPHI,S2SIGN	! max size 4*PHI_SZ
  end module angle_phi

  module angle_alp
    use precision_standard
    use us_size
    implicit none
    save

    !  Declarations of scalars:
    integer(kind=i4b) :: NALPHA
    real(kind=dp) :: CALPHA2
  end module angle_alp

  module energy_mod
    use precision_standard
    use us_size
    implicit none
    save

    !  Declarations of scalars:
    integer(kind=i4b) :: NE,NE1,NE2,NEU
    real(kind=dp) :: DE

    ! Declarations of arrays:
    integer(kind=i4b),dimension(:),allocatable :: I1,I2				! max size E_SZ
    real(kind=dp),dimension(:),allocatable :: E,EU,SPEC0,SPEC1,SPEC2,SPEC3	! max size E_SZ
  end module energy_mod

  module step ! method = 4
    use precision_standard
    use us_size
    implicit none
    save

    ! Declarations of arrays:
    real(kind=dp),dimension(:),allocatable :: HE,HA			! max size E_SZ
  end module step

  module line_shape ! method = 14
    use precision_standard
    use us_size
    implicit none
    save

    !  Declarations of scalars:
    integer(kind=i4b) :: NW
    real(kind=dp) :: DW

    ! Declarations of arrays:
    real(kind=dp),dimension(:),allocatable :: HW			! max size E_SZ
  end module line_shape

  module spectra
    use precision_standard
    use us_size
    implicit none
    save

    ! Declarations of arrays:
    real(kind=dp),dimension(:,:),allocatable :: BR0,BR1,BR2,BR3		! max size PHI_SZ+1,ALPHA_SZ
    real(kind=dp),dimension(:,:),allocatable :: RA0,RA1,RA2,RA3		! max size P_SZ,P_SZ
  end module spectra

  module power
    use precision_standard
    use us_size
    implicit none
    save

    ! Declarations of arrays:
    real(kind=dp),dimension(:,:),allocatable :: PW0,PW1,PW2,PW3		! max size PW_SZ,PW_SZ
  end module power

  module physical_constants
    ! Fundamental physical constants; Physics Today Aug. 1990:
    use precision_standard
    implicit none

    real(kind=dp),parameter :: C    =2.99792458D8	! Speed of light [m/s]
    real(kind=dp),parameter :: ME   =9.1093897D-31	! Electron rest mass [kg]
    real(kind=dp),parameter :: MEE  =0.51099906D0	! Electron rest mass [MeV]
    real(kind=dp),parameter :: EC   =1.60217733D-19	! Elementary charge [C]
    real(kind=dp),parameter :: H    =6.6260755D-34	! Planck's constant [Js]
    real(kind=dp),parameter :: HBAR =1.05457266D-34	! Planck's constant/2Pi [Js]
    real(kind=dp),parameter :: MUZ  =1.2566370614D-6	! Permeability of vacuum [NA-2]
    real(kind=dp),parameter :: EPSZ =8.854187817D-12	! Permittivity of vacuum [Fm-1]
  end module physical_constants

  module pi_constants
    use precision_standard
    implicit none

    real(kind=dp),parameter :: PI    =3.141592653589793238462643D0
    real(kind=dp),parameter :: PIHALF=1.570796326794896619231322D0
    real(kind=dp),parameter :: TWOPI= 6.283185307179586476925287D0
  end module pi_constants
