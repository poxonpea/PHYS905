program Euler
  implicit none

  real(8) ::Tmax=100
  real(8), parameter :: PI=3.14
  integer, parameter ::Num=100
  real(8) :: h,r
  real(8), dimension(Num+1)::vx,vy,x,y
  integer :: i

  h=Tmax/Num

  x(1)=1.d0
  y(1)=0.d0
  vx(1)=0.d0
  vy(1)=1.d0

  do i=1,Num
 !    vx(i+1) = vx(i)-((h*4*PI*PI*x(i))/((x(i)*x(i)+y(i)*y(i))**(3.d0/2.d0)))
 !    if (i==5) then
 !       print*, vx(i+1),vx(i),x(i),y(i)       
 !    end if   
 !    x(i+1) = x(i) + h*vx(i)
 !    vy(i+1) = vy(i)-((h*4*PI*PI*y(i))/((x(i)*x(i)+y(i)*y(i))**(3.d0/2.d0)))
 !    y(i+1) = y(i) + h*vy(i)
     x(i+1)=x(i)+h*vx(i)+(h*h/2)*((4*pi*pi*x(i))/(r(x(i),y(i))**3))
     vx(i+1)=vx(i)+(h/2)*((4*pi*pi*x(i+1))/(r(x(i+1),y(i+1))**3))+(h/2)*((4*pi*pi*x(i))/(r(x(i),y(i))**3))
     y(i+1)=y(i)+h*vy(i)+(h*h/2)*((4*pi*pi*y(i))/(r(x(i),y(i))**3))
     vy(i+1)=vy(i)+(h/2)*((4*pi*pi*y(i+1))/(r(x(i+1),y(i+1))**3))+(h/2)*((4*pi*pi*y(i))/(r(x(i),y(i))**3))

  end do
  
  open(1,file="x.dat")
  do i=1,Num+1
     write(1,*), h*i, x(i)
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
