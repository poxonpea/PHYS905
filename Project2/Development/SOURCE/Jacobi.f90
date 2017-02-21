program Jacobi

  implicit none
  integer, parameter :: Num=300
  integer, parameter :: testnum=4
  real(8), parameter :: rhoMax= 10.d0
  real(8), dimension(Num-1,Num-1)::a,r
  real(8)::Stepsize,start,finish,potential,testmax,max,tolerance,test
  integer:: i,j,maxi,maxj,count
  real(8), dimension(testnum,testnum)::testmat,testr

  testmat(1,1)=1.d0
  testmat(1,2)=2.d0
  testmat(1,3)=3.d0
  testmat(1,4)=4.d0
  testmat(2,1)=2.d0
  testmat(2,2)=6.d0
  testmat(2,3)=7.d0
  testmat(2,4)=8.d0
  testmat(3,1)=3.d0
  testmat(3,2)=7.d0
  testmat(3,3)=11.d0
  testmat(3,4)=12.d0
  testmat(4,1)=4.d0
  testmat(4,2)=8.d0
  testmat(4,3)=12.d0
  testmat(4,4)=16.d0


  do j=1,testnum
     do i=1, testnum
        print*, testmat(j,i)
     end do
  end do

  StepSize=rhoMax/Num
  
  !Fill Matrix
!  do i =1,Num-1
!     do j=1,Num-1
!        if (i==j) then
!           a(j,i)=2.d0/(stepsize**2) + Potential(i*StepSize)
!        else if (j==i+1) then
!           a(j,i)=-1.d0/(stepsize**2)
!        else if (j==i-1) then
!           a(j,i)=-1.d0/(stepsize**2)
!        else
!           a(j,i)=0.d0
!        end if
!     end do      
!  end do

  do i=1,testnum
     do j=1,testnum
        if (i == j) then
           testr(j,i)=1.d0
        else
           testr(j,i)=0.d0
        end if
     end do
  end do
  
  
count = 0
max = 1.d0
tolerance=1.d-8

do while (abs(max) .gt. tolerance)
count=count+1
  !print*,count

  call findmax(testnum,testmat,maxj,maxi,max)
  !print*, max

  call rotate(testnum,testmat,testr,maxj,maxi)
end do

  do i=1,testnum
     do j=1,testnum
        print *, testmat(j,i)
     end do
  end do

  print*, 'vec'
  
  do j=1,testnum
     print*, testr(j,1)
  end do

!  open(12,file="Eigen.dat")
!  do i=1,Num-1
!     write(12,*), b(i)
! end do
!  close(12)
print*, count  
end program Jacobi

function Potential(x)
  implicit none
  real(8)::Potential,x,omega

  Potential = omega*omega*x*x + 1.d0/x
  return

end function Potential


Subroutine findmax(testnum,testmat,maxj,maxi,max)
  ! First find largest offdiagonal matrix element
  implicit none
  integer:: testnum, maxj,maxi,i,j
  real(8)::max
  real(8),dimension(testnum,testnum)::testmat
  
  max=0.d0

  do i=1,testnum
    do j=1,testnum
      if (j/=i) then
        if (abs(testmat(j,i)) .gt. abs(max)) then
          max=testmat(j,i)
          maxj=j
          maxi=i
        end if
      end if
    end do
 end do

 !print*, "max is inside ", max
 return
end subroutine findmax


subroutine rotate(testnum,testmat,testr,maxj,maxi)
  implicit none
  integer::i,j,maxi,maxj,testnum,q
  real(8), dimension(testnum,testnum)::testmat,testr
  real(8)::tau,t,s,c,tempii,tempjj,tempqj,tempqi,temprqi,temprqj
! calculate tau and tangent, sine, and cosine
  if (testmat(maxj,maxi) /=0) then
     tau = (testmat(maxi,maxi) - testmat(maxj,maxj))/(2 * testmat(maxj,maxi))
     if (tau .gt. 0.d0) then
        t = 1.d0/(tau + sqrt(1.d0+tau**2.d0))
     else
        t = -1.d0/(-tau + sqrt(1.d0+tau**2.d0))
     end if
     c = 1.d0/(sqrt(1.d0 + t**2.d0))
     s = t * c
  else
     c=1.d0
     s=0.d0
  end if 

  ! redefine matrix elements using temperary matrix elements
  tempii = testmat(maxi,maxi)
  tempjj = testmat(maxj,maxj)

  testmat(maxj,maxj) = tempjj*(c**2.d0) - 2.d0*testmat(maxj,maxi)*c*s + tempii*(s**2.d0)
  testmat(maxi,maxi) = tempjj*(s**2.d0) + 2.d0*testmat(maxj,maxi)*c*s + tempii*(c**2.d0)
  testmat(maxj,maxi) = 0.d0
  testmat(maxi,maxj) = 0.d0

  do q=1,testnum
     if (q /= maxj .and. q /=maxi) then
       tempqj=testmat(q,maxj)
       tempqi=testmat(q,maxi) 
       testmat(q,maxj)= tempqj*c - tempqi*s
       testmat(maxj,q)= testmat(q,maxj)
       testmat(q,maxi)= tempqi*c + tempqj*s
       testmat(maxi,q)= testmat(q,maxi)
    end if
    temprqj=testr(q,maxj)
    temprqi=testr(q,maxi)
    testr(q,maxj)=c*temprqj - s*temprqi
    testr(q,maxi)=c*temprqi + s*temprqj
  end do
  
  return
end subroutine rotate
