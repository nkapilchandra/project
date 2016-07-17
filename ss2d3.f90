! =====================================================================
! STEADY STATE SOLUTION WITH 3 TEMPERATURES AND 1 INSULTED BOUNDARY
! Size of the matrix is taken from -1 to accomodate insulation criterion
! It is solved for 0 - n in y-direction whereas 1 - n in x-direction
! =====================================================================

program ss2d3
  implicit none
  real, dimension(0:4,-1:4) :: C,D ! D is a second matrix to declared to contain immediate previous values of C while iterating
  integer :: i,j
  real :: e,error

! Boundary conditions
  C(0,1:3) = 75
  C(4,1:3) = 50
  C(1:3,4) = 100
  C(0,0) = 75
  C(4,0) = 50

  e = 10
  error = 0.00001

! NOT NECESSARY --- Initialising all the elements to zero
  do i = 1,3
    do j = 0,3
      C(i,j) = 0
    end do
  end do

! calculating the solution until error falls below a point
  do while (e>error)
! Taking old values of C into D to calculate error later
    do i=0,4
      do j=0,4
        D(i,j) = C(i,j)
      end do
    end do
! Condition which is obtained when an insulated boundary condition is applied
    do i = 1,3
      C(i,-1) = C(i,1)
    end do
! Solving it using iteration on the generic equation obtained
    do i=1,3
      do j=0,3
        C(i,j) = (C(i+1,j) + C(i,j+1) + C(i-1,j) + C(i,j-1))/4.0
      end do
    end do

! Calculation of error
    e = ABS(C(2,2)-D(2,2))*100/C(2,2)

  end do

  do i = 1,3
    do j= 0,3
      print *,C(i,j)
    end do

  end do

end program
