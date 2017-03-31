program Euler
  implicit none

  real(8) ::Tmax=100
  real(8), parameter :: PI=3.14
  integer, parameter ::Num=100000
  real(8) :: h,r,finish,start
  real(8), dimension(Num+1)::vx,vy,x,y
  integer :: i

  h=Tmax/Num

!Initial Conditions from JPS at Midnight on 3/16/17
  
  x(1)=-0.988251
  y(1)=0.08499779
  vx(1)=(365.25*(-.0068024))
  vy(1)=(365.25*(-0.01719988))

! Euler Method Implementation
  call cpu_time(start)
  do i=1,Num
     x(i+1) = x(i) + h*vx(i)
     y(i+1) = y(i) + h*vy(i)
     vx(i+1) = vx(i)-((h*4*PI*PI*x(i))/(r(x(i),y(i))**3.d0))
     vy(i+1) = vy(i)-((h*4*PI*PI*y(i))/(r(x(i),y(i))**3.d0))
     print*,r(x(i),y(i))

  end do
  call cpu_time(finish)

  print*,'time',finish-start
  open(1,file="Euler.dat")
  do i=1,Num+1
     write(1,*), x(i), y(i)
  end do
  close(1)

  
  return
end program Euler


function r(x,y)
  implicit none
  real(8)::r,x,y

  r=(x*x+y*y)**(1.d0/2.d0)
  
  return
end function r
