program GaussElimination
  implicit none
  ! Declare Variables
  real(8)::SourceValue, ExactValue
  real(8)::StepSize
  integer:: NumGridPoints,Power
  real(8), Dimension(:),Allocatable :: u,d,f,e,error
  integer:: i,NumPower

  !Ask user for sizes of matrix and allocate Matricies
  print*,"What is the maximum power (base 10) of the matrix dimension you want?"
  read (*,*) Power

  do NumPower=1, Power
    NumGridPoints=10**NumPower

    print*,NumGridPoints
    
    Allocate(u(NumGridPoints+1))
    Allocate(d(NumGridPoints+1))
    Allocate(f(NumGridPoints+1))
    Allocate(e(NumGridPoints+1))
    Allocate(error(NumGridPoints+1))
  
   !Initialize Variables
    StepSize = 1.d0/(NumGridPoints)
    u(1)=0.d0
    u(NumGridPoints+1)=0.d0
    do i=1,NumGridPoints+1
       f(i)=SourceValue((i-1)*StepSize)*StepSize**2.d0
       d(i)=2.d0
       e(i)=-1.d0
    end do

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

    !Calculate relative error
    do i=1, NumGridPoints+1
      error(i)= log10( abs( u(i) - ExactValue((i-1)*StepSize)) / ExactValue((i-1)*StepSize))
    end do

    call printing(power,u,error,NumGridPoints)


    Deallocate(u)
    Deallocate(d)
    Deallocate(f)
    Deallocate(e)  
    Deallocate(error)

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

subroutine printing(power, u, error,NumGridPoints)
  implicit none
  ! The program prints the results of the above calculations
  integer::i
  integer::power, NumGridPoints
  real(8), dimension(NumGridPoints+1)::u,error
  character(len=8) :: filename
  
do i=1, power
  write (filename, '("myfile" i ".txt")' ) 
  open (file=filename,unit=16)
  write (16,'(5E12.4)') u(i),error(i)
  close (16)
end do

end subroutine printing
