program Ising
  implicit none

  integer:: NumSpins,NumTemp,MCS,i,j,idum
  integer, dimension(:,:),Allocatable:: SpinMatrix
  real(8)::w(-8:8),average(5)
!  real(8), dimension(5)::average
  real(8):: InTemp,FinTemp,E,M,TempStep,Temp

  ! Ask the user for information
  print*, "How large (N) is the spin matrix? (Form N by N)"
  read (*,*) NumSpins
  print*, "What is the initial and final temp? How many temp steps? (initial final steps)"
  read (*,*) InTemp, FinTemp, NumTemp
  print*, "How many Monte Carlo cycles do you want to run?"
  read (*,*) MCS

  Allocate(SpinMatrix(NumSpins,NumSpins))

  TempStep = (FinTemp-InTemp)/NumTemp
  idum=-1
  Temp=intemp
  SpinMatrix = 1

  open(Unit=1,file="MonteCarloResults.txt")
  do i=1, NumTemp
     Temp = Temp + TempStep
     E=0.d0;
     M=0.d0;

     do j=-8,8
        w(j)=0.d0
     end do
     
     do j=-8,8,4
        w(j)=exp(-j/temp)
     end do
     
     do j=1,5
        average(j)=0.d0
     end do

     call initialize(NumSpins,Temp,SpinMatrix,E,M)

     do j=1, MCS
        call Metropolis(NumSpins, idum, SpinMatrix, E,M,w)
        average(1)=average(1)+E
        average(2)=average(2)+(E*E)
        average(3)=average(3)+M
        average(4)=average(4)+(M*M)
        average(5)=average(5)+abs(M)
     end do
     call output(NumSpins,MCS,temp,average)
  end do
  
  
  Deallocate(SpinMatrix)
  
end program Ising

subroutine  initialize(numspins, temperature, spinmatrix, E, M)
  implicit none
  
  integer :: numspins
  Integer :: spinmatrix(numspins, numspins)
  real(8) :: temperature
  real(8) :: E, M
  INTEGER :: x, y, right, left, up, down

  ! setup initial energy and magnetization
  do y =1, numspins
     DO x= 1, numspins
        right = x+1
        if (x == numspins  ) then
           right = 1
        end if
        left = x-1
        if (x == 1  ) then
           left = numspins 
        end if        
        up = y+1
        if(y == numspins  ) then
           up = 1
        end if
        down = y-1
        if (y == 1  ) then
           down = numspins
        end if
        
        E= E-spinmatrix(x,y)*(spinmatrix(right,y)+&
             spinmatrix(left,y)+spinmatrix(x,up)+ &
             spinmatrix(x,down) )
        M = M + spinmatrix(x,y) 
     end do
  end do
  E = E*0.5

end subroutine  initialize

subroutine Metropolis(NumSpins, idum, SpinMatrix, E, M, w)
  implicit none

  integer::NumSpins,x,y,ix,iy,deltaE,right,left,up,down,idum
  integer, dimension(NumSpins,NumSpins)::SpinMatrix
  real(8)::E,M
  real(8),dimension(17)::w
  do y=1,NumSpins
     do x=1,NumSpins
        ix=int(rand(idum)*NumSpins)+1
        iy=int(rand(idum)*NumSpins)+1
        right = ix+1
        if (ix == numspins  ) then
           right = 1
        end if        
        left = ix-1
        if (ix == 1  )then
           left = numspins
        end if
        up = iy+1
        if (iy == numspins  ) then
           up = 1
        end if
        down = iy-1
        if (iy == 1  ) then
           down = numspins
        end if
        
        deltae = 2*spinmatrix(ix,iy)*(spinmatrix(right,iy)+&
             spinmatrix(left,iy)+spinmatrix(ix,up)+ &
             spinmatrix(ix,down) )

        if ( rand(idum) <= w(deltae) ) then
           spinmatrix(ix,iy) = -spinmatrix(ix,iy)  !flip one spin and accept new spin config
           M = M+2*spinmatrix(ix,iy)
           E = E+deltaE
           
        end if
     end do
  end do
  
end subroutine Metropolis

subroutine output(numspins, mcs, temperature, average)
  implicit none
  integer :: numspins, mcs
  real(8) :: temperature, average(5)
  real(8) :: norm, Eaverage, E2average, Maverage, M2average, Mabsaverage
  real(8) ::Evariance, Mvariance

  norm = 1.0/mcs  ! divided by total number of cycles 
  Eaverage = average(1)*norm
  E2average = average(2)*norm
  Maverage = average(3)*norm
  M2average = average(4)*norm
  Mabsaverage = average(5)*norm
  ! all expectation values are per spin, divide by 1/n_spins/n_spins
  Evariance = (E2average- Eaverage*Eaverage)/numspins/numspins
  Mvariance = (M2average - Mabsaverage*Mabsaverage)/numspins/numspins

  write(1,*)temperature, Eaverage/numspins/numspins, Evariance/temperature/temperature, &
       Maverage/numspins/numspins, Mvariance/temperature, Mabsaverage/numspins/numspins

end subroutine output
