program adi
  implicit none

  real :: tk,lambda,del_x,del_t,length,time
  integer :: i,j,k,l,size1,size2,n
  real, allocatable, dimension(:,:,:) :: T
  real, allocatable, dimension(:,:) :: M,T1
  real, allocatable,dimension(:) :: P,buff

  !Inputs from the user
  print *,"enter delta x:"
  read *,del_x ! For delta x

  !print *,"enter delta t:"
  read *,size2 ! For delta t

  print *,"enter length of the rod:"
  read *,length

  !print *,"enter the time period:"
  !read *,time

  print *,"enter thermal diffusivity:"
  read *,lambda

  !lambda = tk*del_t/(del_x**2)
  size1 = int(length/del_x)+1
  !size2 = 2*int(time/del_t)

  allocate(T(0:size2,0:size1,0:size1))
  allocate(T1(0:size1,0:size1))
  allocate(M(size1-1,size1-1))
  allocate(P(size1-1))
  allocate(buff(size1-1))

!Making the constant matrix
  do i = 1,size1-1
    do j =1,size1-1
      M(i,j) = 0
    end do
  end do

  do i = 1,size1-1
    M(i,i) = 2*(1+lambda)
  end do
  do i = 1,size1-2
      M(i,i+1) = -lambda
  end do
do i = 2,size1-1
    M(i,i-1) = -lambda
  end do

  !Initial conditions
  T(0,1:size1-1,1:size1-1) = 0

  !Boundary conditions
  T(0:size2,0,1:size1-1) = 75 !Left surface
  T(0:size2,size1,1:size1-1) = 50 !Right surface
  T(0:size2,1:size1-1,size1) = 100 !Top surface
  T(0:size2,1:size1-1,0) = 0 !Bottom surface

  do i = 1,size2
    do j = 1,size1-1
      do k = 1,size1-1
        T(i,j,k) = 0
      end do
    end do
  end do

  do i = 1,size1-1
    do k = 1,size1-1
      T1(i,k) = T(0,i,k)
    end do
  end do

  T1(0,1:size1-1) = 75 !Left surface
  T1(size1,1:size1-1) = 50 !Right surface
  T1(1:size1-1,size1) = 100 !Top surface
  T1(1:size1-1,0) = 0 !Bottom surface


  do l = 1, size2
    do i = 0,size1
      do j = 0,size1
        T(l,i,j) = T(l-1,i,j)
      end do
    end do
    if (MOD(l,2) .ne. 0) then
    do i = 1,size1-1
    P(1) = lambda*T1(i-1,1) + 2*(1-lambda)*T1(i,1) + lambda*T1(i+1,1) + lambda*T1(i,0)
    P(size1-1) = lambda*T1(i-1,size1-1) + 2*(1-lambda)*T1(i,size1-1) + lambda*T1(i+1,size1-1) + lambda*T1(i,size1)
    do k = 2,size1-2
      P(k) = lambda*T1(i-1,k) + 2*(1-lambda)*T1(i,k) + lambda*T1(i+1,k)
    end do
!print *,P
    call gauss(M,P,size1-1,buff)
    do k = 1,size1-1
      T(l,i,k) = buff(k)
    end do
!print *,buff
do j = 0,size1
  do k = 0,size1
    T1(j,k) = T(l,j,k)
  end do
end do

  end do

else
  do i = 0,size1
    do j = 0,size1
      T(l,i,j) = T(l-1,i,j)
    end do
  end do
  do i = 1,size1-1
  P(1) = lambda*T1(1,i-1) + 2*(1-lambda)*T1(1,i) + lambda*T1(1,i+1) + lambda*T1(0,i)
  P(size1-1) = lambda*T1(size1-1,i-1) + 2*(1-lambda)*T1(size1-1,i) + lambda*T1(size1-1,i+1) + lambda*T1(size1,i)
  do k = 2,size1-2
    P(k) = lambda*T1(k,i-1) + 2*(1-lambda)*T1(k,i) + lambda*T1(k,i+1)
  end do
!print *,P
  call gauss(M,P,size1-1,buff)
  do k = 1,size1-1
    T(l,k,i) = buff(k)
  end do
!print *,buff
do j = 0,size1
do k = 0,size1
  T1(j,k) = T(l,j,k)
end do
end do

if (i .eq. 1) then
do j = 1,size1-1
do k = 1,size1-1
  !print *,T1(j,k)
end do
end do
end if


end do
end if

if (l .eq. 20) then
do j = 1,size1-1
  do i = 1,size1-1
    print *,T(l,i,j)
  end do
end do
end if

end do

end program

subroutine gauss(A,B,n,X)
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

!print *,X

end subroutine
