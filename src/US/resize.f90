    module resize

      contains

      subroutine resize_array(array,n)

! SUBROUTINE DESCRIPTION:
! Routine to dynamically extend the size of an allocatable array (rank 1).
! The new array size is given by the parameter n. It must be larger than the original size.
! Note: Allocatable arrays in a procedure (subroutine or function) need an explicit interface
! or to be embedded in a module procedure.
! RJD ASD/ANL, November 23, 2014.

        use precision_standard

        implicit none

        integer(kind=i4b),intent(in) :: n
        integer(kind=i4b) :: na,lb
        real(kind=dp),dimension(:),allocatable :: array
        real(kind=dp),dimension(:),allocatable :: tmp_array

	na = size(array) ; lb = lbound(array,1)
	! check new size
	if (n .lt. na) then
	  print '(" ",8a)','&RESIZE_ARRAY-F-SZERR, Resize error'
	  print '(" ",a,i8,a,i8)','- new array size ',n,' needs to be greater than ',na
	  call exit (0)
	end if

        allocate(tmp_array(n))
        tmp_array(1:na) = array
        deallocate(array)
        allocate(array(lb:lb+n-1))	! set same lower bound
        array = tmp_array		! copy whole array

      end subroutine resize_array

    end module resize
