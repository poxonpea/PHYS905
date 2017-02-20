program Jacobi

  implicit none
  integer, parameter :: Num=300
  integer, parameter :: testnum=4
  real(8), parameter :: rhoMax= 10.d0
  real(8), dimension(Num-1,Num-1)::a
  real(8)::Stepsize,start,finish,potential,testmax,max,tau,c,s,t
  real(8)::tempii, tempjj, tempij, tempji, tolerance,test,tempqi,tempqj
  integer:: i,j,info,lwork,maxi,maxj,p,q,count
  real(8), dimension(Num-1)::b
  real(8), dimension((num-1)*(2+(num-1)/2)) :: work
  real(8), dimension(testnum,testnum)::testmat

!  do i=1,testnum
!    do j=1,testnum
!      testmat(j,i)=(1.d0+j)/i
!      print*,testmat(j,i)
!    end do
!  end do

  testmat(1,1)=1.d0
  testmat(2,1)=5.d0
  testmat(3,1)=9.d0
  testmat(4,1)=13.d0
  testmat(1,2)=2.d0
  testmat(2,2)=6.d0
  testmat(3,2)=10.d0
  testmat(4,2)=14.d0
  testmat(1,3)=3.d0
  testmat(2,3)=7.d0
  testmat(3,3)=11.d0
  testmat(4,3)=15.d0
  testmat(1,4)=4.d0
  testmat(2,4)=8.d0
  testmat(3,4)=12.d0
  testmat(4,4)=16.d0
  
  StepSize=rhoMax/Num
  lwork=(num-1)*(2+(num-1)/2)
  
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
!     b(i)=SourceValue(i*StepSize)*StepSize**2.d0
!  end do

count = 0
max = 1.d0
tolerance=1.d-10

do while (abs(max) .gt. tolerance)
count=count+1
  print*,count

  call findmax(testnum,testmat,maxj,maxi,max)
  print*, max


! calculate tau and tangent, sine, and cosine
  tau = (testmat(maxi,maxi) - testmat(maxj,maxj))/(2 * testmat(maxj,maxi))
!  print*,tau,maxi,maxj,testmat(maxi,maxi)

  if (testmat(maxj,maxi) /=0) then
     if (tau .gt. 0.d0) then
        t = 1.d0/(tau + sqrt(1.d0+tau**2.d0))
     else
        t = -1.d0/(-tau + sqrt(1.d0+tau**2.d0))
     end if
  else
     c=1.d0
     s=0.d0
  end if 
!  print *, 't is', t

  c = 1.d0/(sqrt(1.d0 + t**2.d0))
!  print*,"c", c
  s = t * c
!  print*,"s",s
  ! redefine matrix elements using temperary matrix elements
  tempii = testmat(maxi,maxi)
  tempjj = testmat(maxj,maxj)
  tempij = testmat(maxi,maxj)
  tempji = testmat(maxj,maxi)

  do q=1,testnum
     if (q /= maxj .and. q /=maxi) then
       tempqj=testmat(q,maxj)
       tempqi=testmat(q,maxi) 
       testmat(q,maxj)= tempqj*c - tempqi*s
       testmat(maxj,q)= testmat(q,maxj)
       testmat(q,maxi)= tempqi*c + tempqj*s
       testmat(maxi,q)= testmat(q,maxi)
    end if
  end do

  testmat(maxj,maxj) = tempjj*(c**2.d0) - 2.d0*tempji*c*s + tempii*(s**2.d0)
  testmat(maxi,maxi) = tempjj*(s**2.d0) + 2.d0*tempji*c*s + tempii*(c**2.d0)
  testmat(maxj,maxi) = 0.d0
  testmat(maxi,maxj) = 0.d0

  
  do i=1,testnum
     do j=1,testnum
        print *, testmat(j,i)
     end do
  end do

  maxj=0
  maxi=0
end do


!  open(12,file="Eigen.dat")
!  do i=1,Num-1
!     write(12,*), b(i)
! end do
!  close(12)
  
  end program Jacobi

function Potential(x)
  implicit none
  real(8)::Potential,x

  Potential = x**2
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

 print*, "max is inside ", max
 return
end subroutine findmax
!  print*, max
