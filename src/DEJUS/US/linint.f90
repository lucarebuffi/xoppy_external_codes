        module linint

        contains
        pure function linint_func(za,ya,na,z,nz)

! FUNCTIONAL DESCRIPTION:
! Function for linear interpolations.
! Given arrays za, ya with na elements and given array z with nz elements for which interpolated values
! are sought the routine returns the array y with requested values.
! Abscissas za and z need to be monotonically increasing (or decreasing).
! Requested y values for z beyond the range of za are set to 0.
! RJD, ASD/ANL, August 11, 2014.
! Updated to use "precision_standard." New name with same declarations to avoid potential
! conflict with Fortran 90 intrinsic function with the name "precision."
! RJD ASD/ANL, November 9, 2014.

        use precision_standard
        implicit none

        integer(kind=i4b),intent(in) :: na,nz
        real(kind=dp),dimension(na),intent(in) :: za,ya
        real(kind=dp),dimension(nz),intent(in) :: z

        real(kind=dp),dimension(nz) :: linint_func      ! automatic array

!  Declarations of scalars:
        integer(kind=i4b) :: i,j,io,ierror
        real(kind=dp) :: dz1,dz2,dz

!  Declarations of arrays:

!  Labeled constants:
        real(kind=dp),parameter :: ZERO=0.0D0

        ierror = 0
        forall (j=1:nz) linint_func(j) = ZERO

        io = 1
        out: do j=1,nz
          in: do i=io,na -1
            dz1 = z(j) -za(i)
            dz2 = z(j) -za(i+1)
            dz  = dz1 -dz2
            if (dz1*dz2 <= ZERO) then
              linint_func(j) = (-ya(i)*dz2 +ya(i+1)*dz1)/dz
              io = i
              exit
            endif
          end do in
        end do out
        return

200     format(' ',8a)
210     format(' ',3(a,f12.7))

!+ Error returns
!  -----------------------------------------------------------------------------
900     continue
!        print 200,'&linint_func-F-ERR, Error ... for future use' ! must not use pure function if used
        ierror = -1
        goto 999

999     continue
        return
        end function linint_func

        end module linint
