!==========================================================================
! IMPLICIT METHOD --- 1-D TRANSIENT SOLUTION
! First and last elements in the rod are used for boundary conditions
! The solution is only for element 1 to size-1
!==========================================================================

program trans2
  implicit none

  integer :: i,j,k,l,size,size2,n ! size--LENGTH; size2--TIME
  real :: del_x,del_t,length,lambda,time,kv,f,sum
  real, allocatable, dimension(:,:) :: C,D1,M,A,D
  real, allocatable,dimension(:) :: b,b1,X

  !Inputs from the user
  print *,"enter delta x:"
  read *,del_x ! For delta x

  print *,"enter delta t:"
  read *,del_t ! For delta t

  print *,"enter length of the rod:"
  read *,length

  print *,"enter the time period:"
  read *,time

  print *,"enter k value:"
  read *,kv

  lambda = kv*del_t/(del_x**2)
  print *,lambda

  size = int(length/del_x)
  size2 = int(time/del_t)

  allocate(C(0:size2,0:size)) ! Number of elements along the length of rod
  allocate(D1(size-1,size-1)) ! Matrix of contants obtained from implcit method
  allocate(b(size-1)) ! Matrix of constants to feed into the subroutine
  allocate(b1(size-1))
  allocate(M(size-1,size))
  allocate(A(size-1,size-1))
  allocate(X(size-1))
  allocate(D(size-1,size))

! Creating the matrix
  do i = 1,size-1
    do j =1,size-1
      D1(i,j) = 0
    end do
  end do
  do i = 1,size-1
    D1(i,i) = 1+2*lambda
  end do
  do i = 1,size-2
      D1(i,i+1) = -lambda
  end do
do i = 2,size-1
    D1(i,i-1) = -lambda
  end do
do i =1,size-1
  do j = 1,size-1
  !print *,D1(i,j)
end do
end do
  !Boundary conditions
  C(0:size2,0) = 100
  C(0:size2,size) = 50

  !Initial conditions
  C(0,1:size-1) = 0
  n=size-1

! Taking values in the constant matrix from C
  do l =1,size2
    b(1) = C(l-1,1) + lambda*C(l,0)
    b(size-1) = C(l-1,size-1) + lambda*C(l,size) ! Separately taking values for first and last element in b because they are different from other elements
    do k = 2,size-2
    b(k) = C(l-1,k)
    end do

!Gauss elimination method starts

!Creating an augemented matrix
    do i=1,n
      do j= 1,n
        D(i,j) = D1(i,j)
      end do
    end do
    do i =1,n
      D(i,n+1) = b(i)
    end do

!Applying elementary row transformations
    do i = 1,n-1
      do j = i+1,n
        f = D(j,i)/D(i,i)
        do k = 1,n+1
          D(j,k) = D(j,k) - f*D(i,k)
        end do
      end do
    end do
!Breaking the augmented matrix
    do i =1,n
      do j=1,n
        A(i,j) = D(i,j)
      end do
    end do
    do i =1,n
      b1(i) = D(i,n+1)
    end do

!Using backward substitution to get the solution
    X(n) = b1(n)/A(n,n)
    do i = n-1,1,-1
      sum =0
      do j = i+1,n
        sum = sum + X(j)*A(i,j)
    end do
    X(i) = (b1(i)-sum)/A(i,i)
  end do

print *,X
  do k = 1,n
    C(l,k) = X(k)
  end do
end do

!Writing the data to a file 
open(unit=1,file='data.dat')
do i = 0,n+1
write(1,*),i*del_x,' ',C(size2,i)
end do
close(1)

end program

