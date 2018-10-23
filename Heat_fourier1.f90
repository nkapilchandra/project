program fft
implicit none

integer :: nnn,i
real,dimension(:),allocatable :: x,k,k2,uo,ut,q
real :: h,l,dx,alpha,dk,pi,z,t,left,right

pi = 3.1416

print*, 'Enter the half length of the rod'
read*, l

print*, 'Enter the no. of points'
read*, nnn

dx = 2*l/nnn

allocate(x(2*nnn))

do i = 1,2*nnn

x(i) = 0

end do

z = 1-nnn/2

do i = 1,nnn

x(2*i-1) = z*dx

z = z + 1

end do

print*, 'Enter the thermal diffusivity'
read*, alpha

dk = pi/l

allocate(k(2*nnn))

do i = 1,2*nnn

k(i) = 0

end do

z = 0

do i = 1,nnn/2+1

k(2*i-1) = z*dk

z = z + 1

end do

z = 1-nnn/2

do i = nnn/2+2,nnn

k(2*i-1) = z*dk

z = z + 1

end do

allocate(k2(2*nnn))

do i = 1,2*nnn

k2(i) = k(i)**2

end do

allocate(uo(2*nnn))

do i = 1,2*nnn

uo(i) = 0

end do

print*, 'Enter the temperature at initial point'
read*, left

print*, 'Enter the temperature at final point'
read*, right

uo(1) = left
uo(2*nnn-1) = right

allocate(q(2*nnn))

allocate(ut(2*nnn))

call four1(uo,nnn,-1)

do i = 1,2*nnn

q(i) = 0

end do

do t = 0.1,1,0.1

do i = 1,nnn

q(2*i-1) = exp(-alpha*k2(2*i-1)*t)

end do

do i = 1,2*nnn

ut(i) = q(i)*uo(i)

end do

call four1(ut,nnn,1)

open(unit = 20,file = 'Temp.txt')

do i = 1,nnn

write(20,*) ut(2*i-1)/nnn

end do

end do

end program

SUBROUTINE four1(data,nn,isign)
INTEGER isign,nn
REAL data(2*nn)
INTEGER i,istep,j,m,mmax,n
REAL tempi,tempr,theta,wi,wpi,wpr,wr,wtemp

n=2*nn
j=1
do i=1,n,2
if(j.gt.i)then
tempr=data(j)
tempi=data(j+1)
data(j)=data(i)
data(j+1)=data(i+1)
data(i)=tempr
data(i+1)=tempi
endif
m=n/2
1 if ((m.ge.2).and.(j.gt.m)) then
j=j-m
m=m/2
goto 1
endif
j=j+m
end do

mmax=2
2 if (n.gt.mmax) then
istep=2*mmax
theta=6.2831/(isign*mmax)
wpr=-2*sin(0.5*theta)**2
wpi=sin(theta)
wr=1
wi=0
do m=1,mmax,2

do i=m,n,istep
j=i+mmax
tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
data(j)=data(i)-tempr
data(j+1)=data(i+1)-tempi
data(i)=data(i)+tempr
data(i+1)=data(i+1)+tempi
end do
wtemp=wr
wr=wr*wpr-wi*wpi+wr
wi=wi*wpr+wtemp*wpi+wi
end do
mmax=istep
goto 2
endif

return
end subroutine four1
