program precisionsingle
  implicit none
  ! The program calculates the limits on step size for 2 and 3 point formulas
  real(4)::taninv,f2c,f3c,calcepsilon
  real(4)::x, exact
  integer::i
  integer::numhpower=10
  real(4), dimension(:),Allocatable ::h,twopoint,threepoint,epsilon2, epsilon3

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
  
end program precisionsingle

  function taninv(x)
    ! This function calculates taninv(x)
    implicit none
    real(4)::x
    real(4)::taninv

    taninv=atan(x)
    !print*, "tan^-1(", x, ") is", taninv
    return
  end function taninv
  

  function f2c(x,h)
    ! This is the function that calcualtes the derivative in the 2 point approx
    implicit none
    real(4):: taninv
    real(4):: x
    real(4):: h
    real(4):: f2c

    f2c=(taninv(x+h)-taninv(x))/h
    !print *, f2c
    
    return
  end function f2c
  
  function f3c(x,h)
    ! This is the function that calcualtes the derivative in the 2 point approx
    implicit none
    real(4):: taninv
    real(4):: x
    real(4):: h
    real(4):: f3c

    f3c=(taninv(x+h)-taninv(x-h))/(2.d0*h)
    !print *, f3c
    
    return
  end function f3c
  


function calcepsilon(comp,exact)
  !This function computes the relative error
  implicit none
  real(4):: calcepsilon
  real(4):: comp
  real(4):: exact

  calcepsilon=log10( abs( comp - exact) / exact)

  return
end function calcepsilon

subroutine printing(h,twopoint,threepoint,epsilon2,epsilon3,numhpower)
  implicit none
  ! The program prints the results of the above calculations
  integer::i
  integer::numhpower
  real(4), dimension(numhpower) ::h,twopoint,threepoint,epsilon2,epsilon3
  
  open (unit=3,file='SingleTwoPoint.txt')
  open (unit=4,file='SingleThreePoint.txt')
  open (unit=7,file='SingleTwoPointep.txt')
  open (unit=8,file='SingleThreePointep.txt')
  do i=1,numhpower
     write(3,"(2f)")h(i),twopoint(i)
     write(4,"(2f)")h(i),threepoint(i)
     write(7,"(2f)")log10(h(i)),epsilon2(i)
     write(8,"(2f)")log10(h(i)),epsilon3(i)
  end do
  close(3)
  close(4)
  close(7)
  close(8)

end subroutine printing
