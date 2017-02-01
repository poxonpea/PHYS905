program GaussElimination
  implicit none
  ! Declare Variables
  real(8)::SourceValue, ExactValue
  real(8)::StepSize
  integer:: NumGridPoints,Power
  real(8), Dimension(:),Allocatable :: u,d,f,e,error,exact
  integer:: i,NumPower
  real(8):: start, finish
  
  !Ask user for sizes of matrix and allocate Matricies
  print*,"What is the maximum power (base 10) of the matrix dimension you want?"
  read (*,*) Power

  do NumPower=1, Power
    NumGridPoints=10**NumPower

    !print*,NumGridPoints
    
    Allocate(u(NumGridPoints+1))
    Allocate(d(NumGridPoints+1))
    Allocate(f(NumGridPoints+1))
    Allocate(e(NumGridPoints+1))
    Allocate(error(NumGridPoints+1))
    Allocate(exact(NumGridPoints+1))
  
   !Initialize Variables
    StepSize = 1.d0/(NumGridPoints)
    u(1)=0.d0
    u(NumGridPoints+1)=0.d0
    do i=1,NumGridPoints+1
       f(i)=SourceValue((i-1)*StepSize)*StepSize**2.d0
       d(i)=2.d0
       e(i)=-1.d0
    end do

    call cpu_time(start)
    
    !Forward Substitution
    do i=3, NumGridPoints
       d(i)= d(i) - (e(i-1) * e(i-1))/ d(i-1)
      !print*,f(i)
       f(i)=f(i) - (e(i-1) * f(i-1))/ d(i-1)
      !print *, (i-1)*StepSize,d(i), f(i),f(i-1),d(i-1)
    end do

    !Backward Substitution
    u(NumGridPoints) = f(NumGridPoints)/d(NumGridPoints)

    do i=NumGridPoints-1, 2, -1
       u(i)= (f(i) - e(i-1) * u(i+1) )/ d(i)
   !    print*,i
    end do

    call cpu_time(finish)
    print*, power, finish-start
    
    !Calculate relative error
    do i=1, NumGridPoints+1
      exact(i)=ExactValue((i-1)*StepSize)
      error(i)= log10( abs( u(i) - exact(i)) / exact(i))
    end do

    call printing(NumPower,u,stepsize,error,exact,NumGridPoints)


    Deallocate(u)
    Deallocate(d)
    Deallocate(f)
    Deallocate(e)  
    Deallocate(error)
    Deallocate(exact)

    end do

end program GaussElimination

function SourceValue(x)
  implicit none
  real(8)::SourceValue, x

  SourceValue = 100.d0 * Exp(- 10.d0 * x)
  return

end function SourceValue	

function ExactValue(x)
  implicit none
  real(8)::ExactValue, x

  ExactValue = 1.d0 - (1.d0 - exp(-10.d0)) * x - exp(-10.d0 * x)
  return

end function ExactValue

subroutine printing(power, u,stepsize, error, exact,NumGridPoints)
  implicit none
  ! The program prints the results of the above calculations
  character(len=10) :: fmt 
  character(len=10) ::x1
  integer::power, NumGridPoints,i
  real(8) :: stepsize
  real(8), dimension(NumGridPoints+1)::u,error,exact
  character(len=25)::solnfile
  character(len=4)::number

  fmt = '(I1.1)' ! an integer of width 5 with zeros at the left

  write (x1,fmt) power ! converting integer to string using a 'internal file'
  solnfile='solution'//trim(x1)//'.dat' 
  open(power,file=solnfile)
  if (power == 3) then
     do i=1,NumGridPoints+1
        write(power,*) (i-1)*stepsize, exact(i),u(i)
     end do
  else if (power /= 3) then
     do i=1,NumGridPoints+1
        write(power,*) (i-1)*stepsize, u(i)
     end do
  end if

  close(power)
  !print *, solnfile


end subroutine printing
