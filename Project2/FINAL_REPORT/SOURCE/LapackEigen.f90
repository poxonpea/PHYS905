program LapackEigenValue
  implicit none
  integer, parameter :: Num=200
  real(8), parameter :: rhoMax= 10.d0
  real(8), dimension(Num-1,Num-1)::a
  real(8)::Stepsize,start,finish,potential
  integer:: i,j,info,lwork
  real(8), dimension(Num-1)::b
  real(8), dimension((num-1)*(2+(num-1)/2)) :: work


  StepSize=rhoMax/Num
  lwork=(num-1)*(2+(num-1)/2)
  
  !Fill Matrix
  do i =1,Num-1
     do j=1,Num-1
        if (i==j) then
           a(j,i)=2.d0/(stepsize**2) + Potential(i*StepSize)
        else if (j==i+1) then
           a(j,i)=-1.d0/(stepsize**2)
        else if (j==i-1) then
           a(j,i)=-1.d0/(stepsize**2)
        else
           a(j,i)=0.d0
        end if
     end do      
!     b(i)=SourceValue(i*StepSize)*StepSize**2.d0
  end do


  !Call LAPACK Solver

  call cpu_time(start)
  
  call DSYEV('V','U',Num-1,A,Num-1,b,work,lwork,info)
  print*,info

  call cpu_time(finish)
  print*, finish-start
  
  open(12,file="Eigen.dat")
  do i=1,Num-1
     write(12,*), b(i)
  end do
  close(12)
  
  end program LapackEigenValue

function Potential(x)
  implicit none
  real(8)::Potential,x

  Potential = x**2
  return

end function Potential
