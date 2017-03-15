program Euler
  implicit none

  real(8) ::Tmax=10
  real(8), parameter :: PI=3.14
  integer, parameter ::Num=1000
  real(8) :: h,r
  real(8), dimension(Num+1)::vx,vy,x,y
  integer :: i

  h=Tmax/Num

  x(1)=-0.9864
  y(1)=0.10218
  vx(1)=(365.25*(-.0019768))
  vy(1)=(365.25*(-0.017172))

  do i=1,Num
     vx(i+1) = vx(i)-((h*4*PI*PI*x(i))/(r(x(i),y(i))**3.d0))   
     x(i+1) = x(i) + h*vx(i)
     vy(i+1) = vy(i)-((h*4*PI*PI*y(i))/(r(x(i),y(i))**3.d0))
     y(i+1) = y(i) + h*vy(i)

  end do
  
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
