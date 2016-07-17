!==================================================================
! EXPLICIT METHOD --- 1-D TRANSIENT SOLUTION
! First element and last element in the rod are used for boundary conditions
! The solution is only for element 1 to size-1.
!==================================================================

program ts1dex
implicit none

integer :: i,j,counter,size,size2 ! size--LENGTH; size2--TIME
real :: del_x,del_t,length,lambda,time
real, allocatable, dimension(:,:) :: C

!Inputs from the user
print *,"enter delta x:"
read *,del_x ! For delta x

print *,"enter delta t:"
read *,del_t ! For delta t

print *,"enter length of the rod:"
read *,length

print *,"enter the time period:"
read *,time

print *,"enter lambda value:"
read *,lambda

size = int(length/del_x) ! Number of elements along the length of rod
size2 = int(time/del_t)  ! Number of time elements

allocate(C(0:size2,0:size)) ! Allocating size from 0 to accomodate boundary and initial conditions

!Boundary conditions

C(0:size2,0) = 100
C(0:size2,size) = 50

!Initial condition

C(0,1:size-1) = 0

! NOT NECESSARY - Initialising all the values in the matrix to zero
do i = 0,size2
  do j = 1,size-1
    C(i,j) = 0
  end do
end do

! Solving the generic equation, dynamically storing the values in the matrix
do i = 0,size2
  do j = 1,size-1
    C(i+1,j) = C(i,j) + lambda*(C(i,j+1)-2*C(i,j) + C(i,j-1))
  end do
end do


do j = 1,size-1
  print *,C(size2,j)
end do

end program
