!======================================================================
! STEADY STATE SOLUTION IN 2-D WITH 4 TEMPERATURE BOUNDARY CONDITIONS
! First and last elements are used for boundary conditions
!======================================================================

program solve
  implicit none
  real, allocatable, dimension(:,:) :: C,D,Fx,Fy,F,A ! D is a matrix declared to store the immediate previous values of C to calculate error
  integer :: i,j,counter,x,y
  real :: e,error,lx,ly,del_x,del_y,tk

print *,'length:'
read *,lx
print *,'breadth'
read *,ly
print *,'thermal conductivity'
read *,tk
print *,'number of points in x-direction:'
read *,x
print *,'number of points in y-direction'
read *,y

del_x = lx/(1+x)
del_y = ly/(1+y)

allocate(C(0:(1+x),0:(1+y)))
allocate(D(0:(1+x),0:(1+y)))
allocate(F(0:(1+x),0:(1+y)))
allocate(Fx(0:(1+x),0:(1+y)))
allocate(Fy(0:(1+x),0:(1+y)))
allocate(A(0:(1+x),0:(1+y)))

! Boundary conditions
do i = 1,y
  C(0,i) = 75
end do
do i = 1,x
  C(i,0) = 0
end do
do i = 1,x
  C(i,1+y) = 100
end do
do i = 1,y
  C(1+x,i) = 50
end do


  counter=0 ! To check the number of steps taken

  e=10
  error = 0.000001

  do i=0,1+x
    do j=0,1+y
      D(i,j) = C(i,j)
    end do
  end do


! NOT NECESSARY -- Initialising all the elements to 0
  do i = 1,x
    do j=1,y
      C(i,j) = 0
      Fx(i,j) = 0
      Fy(i,j) = 0
      F(i,j) = 0
    end do
  end do

do while (e>error)
counter = counter +1
! Storing old values of C in D
  do i=0,1+x
    do j=0,1+y
      D(i,j) = C(i,j)
    end do
  end do

! Solving the generic equation
  do i = 1, x
    do j = 1,y
      C(i,j) = (C(i+1,j)+C(i,j+1)+C(i-1,j)+C(i,j-1))/4.0
    end do
  end do

  e = ABS(C(x-1,x-1)-D(x-1,x-1))*100/C(x-1,x-1)

end do


do i = 1,x
  do j =1,y
    Fx(i,j) = (-tk)*(C(i+1,j)-C(i-1,j))/(2*del_x)
    Fy(i,j) = (-tk)*(C(i,j+1)-C(i,j-1))/(2*del_y)
  end do
end do

do i = 1,x
  do j = 1,y
    F(i,j) = sqrt((Fx(i,j)**2) + (Fy(i,j)**2))
  end do
end do

do i = 1,x
  do j =1,y
    A(i,j) = atan(Fy(i,j)/Fx(i,j))*180/(4*atan(1.0))
  end do
end do

do i =1,x
  do j = 1,y
print *,A(i,j)
end do
end do

end program
