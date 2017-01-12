program precisiondouble
  implicit none
  ! The program calculates the limits on step size for 2 and 3 point formulas
  real(8)::taninv,f2c,f3c,calcepsilon
  real(8)::x, exact
  integer::i
  integer::numhpower=10
  real(8), dimension(:),Allocatable ::h,twopoint,threepoint,epsilon2,epsilon3

  Allocate(h(numhpower))
  Allocate(twopoint(numhpower))
  Allocate(threepoint(numhpower))
  Allocate(epsilon2(numhpower))
  Allocate(epsilon3(numhpower))

  x=sqrt(2.d0)
  exact=1.d0/3.d0
  do i=1,numhpower
     h(i)=1/10.d0**i
     twopoint(i)=f2c(x,h(i))
     threepoint(i)=f3c(x,h(i))
     !print *, h(i),twopoint(i),threepoint(i)
     epsilon2(i) = calcepsilon(twopoint(i),exact)
     epsilon3(i)=calcepsilon(threepoint(i),exact)
  end do

  call printing(h,twopoint,threepoint,epsilon2,epsilon3,numhpower)

end program precisiondouble

function taninv(x)
  ! This function calculates taninv(x)
  implicit none
  real(8)::x
  real(8)::taninv

  taninv=atan(x)
  !print*, "tan^-1(", x, ") is", taninv
  return
end function taninv


function f2c(x,h)
  ! This is the function that calcualtes the derivative in the 2 point approx
  implicit none
  real(8):: taninv
  real(8):: x
  real(8):: h
  real(8):: f2c

  f2c=(taninv(x+h)-taninv(x))/h
  !print *, f2c

  return
end function f2c

function f3c(x,h)
  ! This is the function that calcualtes the derivative in the 2 point approx
  implicit none
  real(8):: taninv
  real(8):: x
  real(8):: h
  real(8):: f3c

  f3c=(taninv(x+h)-taninv(x-h))/(2.d0*h)
  !print *, f3c

  return
end function f3c


function calcepsilon(comp,exact)
  !This function computes the relative error
  implicit none
  real(8):: calcepsilon
  real(8):: comp
  real(8):: exact

  calcepsilon=log10( abs( comp - exact) / exact)

  return
end function calcepsilon

subroutine printing(h,twopoint,threepoint,epsilon2,epsilon3,numhpower)
  implicit none
  ! The program prints the results of the above calculations
  integer::i
  integer::numhpower
  real(8), dimension(numhpower) ::h,twopoint,threepoint,epsilon2,epsilon3
  
  open (unit=1,file='DoubleTwoPoint.txt')
  open (unit=2,file='DoubleThreePoint.txt')
  open (unit=5,file='DoubleTwoPointep.txt')
  open (unit=6,file='DoubleThreePointep.txt')
  do i=1,numhpower
     write(1,"(2f)")h(i),twopoint(i)
     write(2,"(2f)")h(i),threepoint(i)
     write(5,"(2f)")log10(h(i)),epsilon2(i)
     write(6,"(2f)")log10(h(i)),epsilon3(i)
  end do
  close(1)
  close(2)
  close(5)
  close(6)

end subroutine printing

