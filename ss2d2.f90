!========================================================================================
! STEADY STATE SOLUTION IN 2-D WITH 2 TEMPERATURE AND 2 INSULATED BOUNDARY CONDITIONS
! It solves for 10 elements each in x- and y-direction
!=========================================================================================

program ss2d2
  implicit none

  real, dimension(0:11,0:11) :: C,D ! D is a matrix declared to store the immediate previous values of C to calculate error
  integer :: i,j
  real :: e,error
! Boundary conditions
  C(0,1:10) = 75
  C(11,1:10) = 0

  e = 10
  error = 0.0001

! NOT NECESSARY -- Initialising all the elements to 0
  do i = 1,10
    do j = 1,10
      C(i,j) = 0
    end do
  end do


  do i = 1,10
    C(i,0) = C(i,2)
    C(i,11) = C(i,9)
  end do

  do while (e>error)

! Storing old values of C in D
    do i=0,11
      do j=0,11
        D(i,j) = C(i,j)
      end do
    end do

! Solving the generic equation for all the elements
    do i=1,10
      do j=1,10
        C(i,j) = (C(i+1,j) + C(i,j+1) + C(i-1,j) + C(i,j-1))/4.0
      end do
    end do

! Conditions obtained for insulated boundary conditions on the two walls
    do i = 1,10
      C(i,0) = C(i,2)
      C(i,11) = C(i,9)
    end do

    e = ABS(C(5,5)-D(5,5))*100.0/C(5,5)
    !print *,e
  end do

  !print *,C(5,5)


end program
