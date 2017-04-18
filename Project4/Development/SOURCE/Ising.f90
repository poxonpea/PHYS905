program Ising
  implicit none

  integer:: NumSpins,NumTemp,MCS,i,j,idum,orientation,accept,acceptance
  integer, dimension(:,:),Allocatable:: SpinMatrix
  real(8)::w(-8:8),average(5)
  real(8):: InTemp,FinTemp,E,M,TempStep,Temp,random

  ! Ask the user for information
  print*, "How large (N) is the spin matrix? (Form N by N)"
  read (*,*) NumSpins
  print*, "What is the initial and final temp? How many temp steps? (initial final steps)"
  read (*,*) InTemp, FinTemp, NumTemp
  print*, "How many Monte Carlo cycles do you want to run?"
  read (*,*) MCS
  print*, "Should I start with an ordered (1) or random (2) orientation?"
  read (*,*) Orientation

  Allocate(SpinMatrix(NumSpins,NumSpins))

  TempStep = (FinTemp-InTemp)/NumTemp
  idum=-1
  Temp=intemp
  accept = 0
  
  call init_random_seed()
  if (orientation ==1) then
     SpinMatrix = 1
  else if (orientation ==2) then
     do i=1,numspins
        do j=1, numspins
           if (rand() .ge. 0.5) then
              SpinMatrix(j,i)=1
           else
              SpinMatrix(j,i)=-1
           end if
        end do
     end do
  end if
  
  
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
        call Metropolis(NumSpins, idum, SpinMatrix, E,M,w,accept)
        average(1)=average(1)+E
        average(2)=average(2)+(E*E)
        average(3)=average(3)+M
        average(4)=average(4)+(M*M)
        average(5)=average(5)+abs(M)
        !        acceptance = acceptance + accept
!       if (j .ge. 200000 .and. mod(j,10)==0) then
!          call output(NumSpins,j,temp,average,acceptance)
!        end if
     end do
       call output(NumSpins,MCS,temp,average,acceptance)
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
     do x= 1, numspins
        right = x+1
        if (x == numspins  ) then
           right = 1
        end if
        left = x-1
        if (x == 1  ) then
           left = numspins 
        end if        
        up = y+1
        if(y == numspins ) then
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

subroutine Metropolis(NumSpins, idum, SpinMatrix, E, M, w,accept)
  implicit none

  integer::NumSpins,x,y,ix,iy,deltaE,right,left,up,down,idum,accept
  integer, dimension(NumSpins,NumSpins)::SpinMatrix
  real(8)::E,M
  real(8)::w(-8:8)
  accept=0
  do y=1,NumSpins
     do x=1,NumSpins
        ix=int(rand()*NumSpins)+1
        iy=int(rand()*NumSpins)+1
        right = ix+1
        if (ix == numspins) then
           right = 1
        end if        
        left = ix-1
        if (ix == 1)then
           left = numspins
        end if
        up = iy+1
        if (iy == numspins) then
           up = 1
        end if
        down = iy-1
        if (iy == 1) then
           down = numspins
        end if    
        deltae = 2*spinmatrix(ix,iy)*(spinmatrix(right,iy)+&
             spinmatrix(left,iy)+spinmatrix(ix,up)+ &
             spinmatrix(ix,down) )
        if ( rand() <= w(deltae) ) then
           accept=accept+1
           spinmatrix(ix,iy) = -spinmatrix(ix,iy)  
           M = M+2*spinmatrix(ix,iy)
           E = E+deltaE           
        end if
     end do
  end do  
end subroutine Metropolis

subroutine output(numspins, mcs, temperature, average,accept)
  implicit none
  integer :: numspins, mcs,accept
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

  write(1,*)temperature, mcs,Eaverage/numspins/numspins, Evariance/temperature/temperature, &
       Maverage/numspins/numspins, Mvariance/temperature, Mabsaverage/numspins/numspins

end subroutine output


 subroutine init_random_seed()
      integer i,n,clock
      integer, dimension(:), allocatable :: seed
 
      call random_seed(size = n)
      allocate(seed(n))
 
      call system_clock(count=clock)
 
      seed = clock + 37 *(/ (i-1, i=1, n) /)
      call random_seed(put=seed)
 
      deallocate(seed)
 
end subroutine
