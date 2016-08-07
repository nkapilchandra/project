program check
  real, dimension(5,5) :: A
  real,dimension(5) :: B
  integer :: i,j

  print *,'A:'
  do i =1,5
    do j =1,5
      read *,A(i,j)
    end do
  end do

  print *,'B:'
  do i =1,5
    read *,B(i)
  end do

call gauss(A,B,5)

end program

subroutine gauss(A,B,n)
  implicit none
  real, dimension(n,n) :: A,A1
  real, dimension(n) :: B,B1,X
  real,dimension(n,n+1) :: D
  integer ::n,i,j,k
  real ::f,sum


!Making the augmented matrix
  do i=1,n
    do j= 1,n
      D(i,j) = A(i,j)
    end do
  end do
  do i =1,n
    D(i,n+1) = B(i)
  end do


!Elementary row transformations
  do i = 1,n-1
    do j = i+1,n
      f = D(j,i)/D(i,i)
      do k = 1,n+1
        D(j,k) = D(j,k) - f*D(i,k)
      end do
    end do
  end do

  !Breaking the obtained augmented matrix
    do i =1,n
      do j=1,n
        A1(i,j) = D(i,j)
      end do
    end do
    do i =1,n
      B1(i) = D(i,n+1)
    end do

  !Backward substitution
    X(n) = B1(n)/A1(n,n)
    do i = n-1,1,-1
      sum =0
      do j = i+1,n
        sum = sum + X(j)*A1(i,j)
    end do
    X(i) = (B1(i)-sum)/A1(i,i)
  end do

print *,X

end subroutine
