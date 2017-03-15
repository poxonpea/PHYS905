program Velver
  implicit none

  real(8) ::Tmax=1
  real(8), parameter :: PI=3.14
  integer, parameter ::Num=100
  real(8) :: h,r
  real(8), dimension(Num+1)::vx,vy,x,y
  integer :: i

  h=Tmax/Num

  x(1)=-0.9864
  y(1)=0.10218
  vx(1)=(365.25*(-.0019768))
  vy(1)=(365.25*(-0.017172))

  do i=1,Num
     x(i+1)=x(i)+h*vx(i)-(h*h/2)*((4.d0*pi*pi*x(i))/(r(x(i),y(i))**3))
     y(i+1)=y(i)+h*vy(i)-(h*h/2)*((4.d0*pi*pi*y(i))) 
     vx(i+1)=vx(i)-(h/2)*((4.d0*pi*pi*x(i+1))/(r(x(i+1),y(i+1))**3))-(h/2)*((4.d0*pi*pi*x(i))/(r(x(i),y(i))**3))
     vy(i+1)=vy(i)-(h/2)*((4.d0*pi*pi*y(i+1))/(r(x(i+1),y(i+1))**3))-(h/2)*((4.d0*pi*pi*y(i))/(r(x(i),y(i))**3))
  end do
  
  open(1,file="vel_verl.dat")
  do i=1,Num+1
     write(1,*), x(i), y(i)
  end do
  close(1)

  
  return
end program VelVer


function r(x,y)
  implicit none
  real(8)::r,x,y

  r=(x*x+y*y)**(1.d0/2.d0)
  
  return
end function r
