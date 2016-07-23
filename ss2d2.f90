!========================================================================================
! STEADY STATE SOLUTION IN 2-D WITH 2 TEMPERATURE AND 2 INSULATED BOUNDARY CONDITIONS
! It solves for 10 elements each in x- and y-direction
!=========================================================================================

program solve1
  implicit none

  real,allocatable, dimension(:,:) :: C,D,F,Fx,Fy,A ! D is a matrix declared to store the immediate previous values of C to calculate error
  integer :: i,j,x,y
  real :: e,error,tk,del_x,del_y,lx,ly

  print *,'length:'
  read *,lx
  print *,'breadth'
  read *,ly
  print *,'thermal diffusivity'
  read *,tk
  print *,'number of points in x-direction'
  read *,x
  print *,'number of points in y-direction'
  read *,y

allocate(C(0:(1+x),0:(1+y)))
allocate(D(0:(1+x),0:(1+y)))
allocate(F(0:(1+x),0:(1+y)))
allocate(Fx(0:(1+x),0:(1+y)))
allocate(Fy(0:(1+x),0:(1+y)))
allocate(A(0:(1+x),0:(1+y)))

del_x = lx/(1+x)
del_y = ly/(1+y)

! Boundary conditions
  C(0,1:y) = 100
  C(1+x,1:y) = 0

  e = 10
  error = 0.0000001

! NOT NECESSARY -- Initialising all the elements to 0
  do i = 1,x
    do j = 1,y
      C(i,j) = 0
    end do
  end do


  do i = 1,x
    C(i,0) = C(i,2)
    C(i,1+y) = C(i,y-1)
  end do

  do while (e>error)

! Storing old values of C in D
    do i=0,1+x
      do j=0,1+y
        D(i,j) = C(i,j)
      end do
    end do

! Solving the generic equation for all the elements
    do i=1,x
      do j=1,y
        C(i,j) = (C(i+1,j) + C(i,j+1) + C(i-1,j) + C(i,j-1))/4.0
      end do
    end do

! Conditions obtained for insulated boundary conditions on the two walls
    do i = 1,x
      C(i,0) = C(i,2)
      C(i,1+y) = C(i,y-1)
    end do

    e = maxval(ABS(C(1:x,1:y)-D(1:x,1:y)))
  end do

do i =1,x
  do j =2,y-1
    Fx(i,j) = (-tk)*(C(i+1,j)-C(i-1,j))/(2*del_x)
  end do
end do


do i =1,x
  do j =2,y-1
    Fy(i,j) = (-tk)*(C(i,j+1)-C(i,j-1))/(2*del_y)
  end do
end do

do i =1,x
  Fy(i,1) = 0
  Fy(i,y) = 0
end do

do i = 1,x
  do j =1,y
    F(i,j) = sqrt((Fx(i,j)**2) + (Fy(i,j)**2))
  end do
end do

do i =1,x
  do j =1,y
    if (Fx(i,j).ne.0) then
    A(i,j) = atan(Fy(i,j)/Fx(i,j))*180/(4*atan(1.0))
  else
    A(i,j) = 2*atan(1.0)
  end if
  end do
end do

do i =1,x
  do j= 1,y
    print *,Fy(i,j)
  end do
end do

open(unit=1,file='data.xls')
do i = 1,x
  do j =1,y
    write(1,*), C(i,j)
  end do
end do
close(1)

open(unit=2,file='data1.xls')
do i = 1,x
  do j =1,y
    write(2,*), Fx(i,j)
  end do
end do
close(2)

open(unit=3,file='data2.xls')
do i = 1,x
  do j =1,y
    write(3,*), Fy(i,j)
  end do
end do
close(3)


end program
