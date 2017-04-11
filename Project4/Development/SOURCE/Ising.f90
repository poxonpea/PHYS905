program Ising
  implicit none

  integer:: NumSpins,NumTemp,MCS,i,j,idum
  integer, dimension(:,:),Allocatable:: SpinMatrix
  real(8), dimension(17)::w
  real(8), dimension(5)::average
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

  do i=1, NumTemp
     Temp = Temp + TempStep
     E=0.d0;
     M=0.d0;

     do j=1,17
        w(j)=0.d0
     end do
     
     do j=1,17,4
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

     !print*,average(1)
     
  end do
  
  
  Deallocate(SpinMatrix)
  
end program Ising


subroutine initialize(NumSpins, temp, SpinMatrix,E,M)
  implicit none

  Integer::NumSpins,y,x,period,periodicx,periodicy
  Integer,dimension(NumSpins, NumSpins)::SpinMatrix
  real(8)::temp,E,M
  
  do y=1, NumSpins
     do x=1,NumSpins
        if (temp .lt. 1.5) then
           spinmatrix(x,y)=1
           M = M + float(spinmatrix(x,y))
        end if
     end do
  end do

  do y=1,NumSpins
     do x=1,NumSpins
        periodicy= period(y,numspins,-1)
        periodicx=period(x,numspins,-1)
        E = E-float(spinmatrix(x,y))*(spinmatrix(x,periodicy)+spinmatrix(periodicx,y))
     end do
  end do

end subroutine initialize

function period(i,limit,add)
  implicit none

  integer:: i, limit, add
  integer:: period

  period=(i+limit+add)

  return
end function period

subroutine Metropolis(NumSpins, idum, SpinMatrix, E, M, w)
  use ifport
  implicit none

  integer::NumSpins,x,y,ix,iy,deltaE,period,pperiodicx,pperiodicy,nperiodicx,nperiodicy,idum
  integer, dimension(NumSpins,NumSpins)::SpinMatrix
  real(8)::E,M,ixtest,randtest,iytest
  real(8),dimension(17)::w
  do y=1,NumSpins
     do x=1,NumSpins
        ixtest=rand()*NumSpins
        print*,"ixtest is",ixtest
        ix=int(ixtest)+1
        print*,"ixint is", ix
        iytest=rand()*NumSpins
        print*,"iytest is", iytest
        iy=int(iytest)+1
        print*, "iyint is ", iy
        pperiodicx=period(x,NumSpins,1)
        pperiodicy=period(y,NumSpins,1)
        nperiodicx=period(x,NumSpins,-1)
        nperiodicy=period(y,NumSpins,-1)
        deltaE= 2*SpinMatrix(ix,iy)*(SpinMatrix(nperiodicx,iy)+Spinmatrix(ix,nperiodicy)&
             +Spinmatrix(pperiodicx,iy)+spinmatrix(ix,pperiodicy))
        randtest = rand()
        if (randtest .le. w(deltaE+8)) then
           spinMatrix(ix,iy) = SpinMatrix(ix,iy)
           M = M+2*spinmatrix(ix,iy)
           E = E +deltaE
        end if
     end do
  end do
  
end subroutine Metropolis
