!==========================================================================
! IMPLICIT METHOD --- 1-D TRANSIENT SOLUTION
! First and last elements in the rod are used for boundary conditions
! The solution is only for element 1 to size-1
!==========================================================================

program ts1dim
  implicit none

  integer :: i,j,k,counter,size,size2 ! size--LENGTH; size2--TIME
  real :: del_x,del_t,length,lambda,time
  real, allocatable, dimension(:,:) :: C,D
  real, allocatable,dimension(:) :: S,b

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

  size = int(length/del_x)
  size2 = int(time/del_t)

  allocate(C(0:size2,0:size)) ! Number of elements along the length of rod
  allocate(D(size,size)) ! Matrix of contants obtained from implcit method
  allocate(S(size)) ! Matrix to take the solution from the subroutine
  allocate(b(size-1)) ! Matrix of contants to feed into the subroutine

! Creating the matrix
  do i = 1,size
    do j =1,size
      D(i,j) = 0
    end do
  end do
  do i = 1,size
    D(i,i) = 1+2*lambda
  end do
  do i = 2,size
    if (i /= size) then
    D(i,i+1) = -lambda
    end if
    D(i,i-1) = -lambda
  end do

  !Boundary conditions
  C(0:size2,0) = 100
  C(0:size2,size) = 50

  !Initial conditions
  C(0,1:size-1) = 0

! Taking values in the constant matrix from C
  do j =1,size2
    b(1) = C(j-1,1) + lambda*C(j,0)
    b(size-1) = C(j-1,size-1) + lambda*C(j,size) ! Separately taking values for first and last element in b because they are different from other elements
    do k = 2,size-2
    b(k) = C(j-1,k)
    end do

      call gauss(D,b,size,S) ! S stores the solution

! Storing the solution in C
      do k = 1, size-1
      C(j,k) = S(k)
    end do
  end do

  do i =1,size2
    do j= 1,size -1
      print *,C(i,j)
    end do
  end do


end program

!!!!!!!!!
!subroutine for gauss elimination taken from other file
!!!!!!!!!

subroutine gauss(A,B,n,X)
  implicit none

  integer :: n,i,j,k
  real,dimension(n,n) :: A,C
  real,dimension(n,n+1) :: D
  real,dimension(n) :: B,B1
  real,dimension(n) :: X
  real :: f,sum

  do i = 1,n
    do j = 1,n
      D(i,j) = A(i,j)
    end do
  end do

  do i = 1,n
    D(i,n+1) = B(i)
  end do

do k = 1,n-1
  do i = k+1,n
    f = D(i,k)/D(k,k)
    do j = k,n+1
      D(i,j) = D(i,j) - f*D(k,j)
    end do
  end do
end do

do i = 1,n
  do j =1,n
    C(i,j) = D(i,j)
  end do
end do

do i = 1,n
  B1(i) = D(i,n+1)
end do

X(n) = B1(n)/C(n,n)
do i = n-1,1,-1
  sum = 0
  do j = n,i+1,-1
    sum = sum + X(j)*C(i,j)
  end do

X(i) = (B1(i)-sum)/C(i,i)

end do


return
end subroutine
