!======================================================================
! STEADY STATE SOLUTION IN 2-D WITH 4 TEMPERATURE BOUNDARY CONDITIONS
! First and last elements are used for boundary conditions
!======================================================================

program ss2d1
  implicit none
  real, dimension(0:4,0:4) :: C,D ! D is a matrix declared to store the immediate previous values of C to calculate error
  integer :: i,j,counter
  real :: e,error

! Boundary conditions
  C(0,1:3) = 75
  C(1:3,0) = 0
  C(1:3,4) = 100
  C(4,1:3) = 50

  counter=0 ! To check the number of steps taken

  e=10
  error = 0.00001

  do i=0,4
    do j=0,4
      D(i,j) = C(i,j)
    end do
  end do

! NOT NECESSARY -- Initialising all the elements to 0
  do i = 1,3
    do j=1,3
      C(i,j) = 0
    end do
  end do

do while (e>error)
counter = counter +1
! Storing old values of C in D
  do i=0,4
    do j=0,4
      D(i,j) = C(i,j)
    end do
  end do

! Solving the generic equation
  do i = 1, 3
    do j = 1,3
      C(i,j) = (C(i+1,j)+C(i,j+1)+C(i-1,j)+C(i,j-1))/4.0
    end do
  end do

  e = ABS(C(2,2)-D(2,2))*100/C(2,2)

end do

!print *,e
!print *,C(2,2)
!print *,counter

end program
