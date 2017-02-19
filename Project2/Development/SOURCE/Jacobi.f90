program Jacobi

  implicit none
  integer, parameter :: Num=300
  integer, parameter :: testnum=3
  real(8), parameter :: rhoMax= 10.d0
  real(8), dimension(Num-1,Num-1)::a
  real(8)::Stepsize,start,finish,potential,testmax,max,tau,c,s,t
  real(8)::tempii, tempjj, tempij, tempji, tolerance,test
  integer:: i,j,info,lwork,maxi,maxj,p,q,count
  real(8), dimension(Num-1)::b
  real(8), dimension((num-1)*(2+(num-1)/2)) :: work
  real(8), dimension(3,3)::testmat

!  do i=1,testnum
!    do j=1,testnum
!      testmat(j,i)=(1.d0+j)/i
!      print*,testmat(j,i)
!    end do
!  end do

  testmat(1,1)=-2.d0
  testmat(2,1)=-2.d0
  testmat(3,1)=4.d0
  testmat(1,2)=-4.d0
  testmat(2,2)=1.d0
  testmat(3,2)=2.d0
  testmat(1,3)=2.d0
  testmat(2,3)=2.d0
  testmat(3,3)=5.d0
  
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

count = 0
test = 1
do while (abs(test) .gt. 1.d-10)
count=count+1
  print*,count

! First find largest offdiagonal matrix element
  max=0.d0

  do i=1,testnum
    do j=1,testnum
      if (j/=i) then
        testmax = testmat(j,i)
        if (abs(testmax) .ge. abs(max)) then
          max=testmax
          test = max
          !print *,"j",j,"i",i,"testmat(j,i)",testmat(j,i),"test",test,"max",max
          maxj=j
          maxi=i
!          print *, maxi, maxj
        end if
      end if
    end do
  end do
!  print*, max

! calculate tau and tangent, sine, and cosine
  tau = (testmat(maxi,maxi) - testmat(maxj,maxj))/(2 * testmat(maxj,maxi))
!  print*,tau,maxi,maxj,testmat(maxi,maxi)

  if (tau .gt. 0.d0) then
     t = 1.d0/(tau + sqrt(1.d0+tau**2))
  else
    t = 1.d0/(-tau + sqrt(1.d0+tau**2))
  end if

!  print *, 't is', t

  c = 1.d0/(sqrt(1.d0 + t**2))
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
      testmat(q,maxj)= testmat(q,maxj)*c - testmat(q,maxi)*s
      testmat(q,maxi)= testmat(q,maxi)*c + testmat(q,maxj)*s
    end if
  end do

  testmat(maxj,maxj) = tempjj*(c**2.d0) - 2.d0*tempji*c*s + tempii*(s**2.d0)
  testmat(maxi,maxi) = tempii*(s**2.d0) + 2.d0*tempji*c*s + tempjj*(c**2.d0)
  print*, "checking", maxj,maxi,tempjj, tempji, tempii,c,s, testmat(maxj,maxj)
  testmat(maxj,maxi) = 0.d0

  
  do i=1,testnum
     do j=1,testnum
        print *, testmat(j,i)
     end do
  end do

  maxj=0
  maxi=0
end do


  open(12,file="Eigen.dat")
  do i=1,Num-1
     write(12,*), b(i)
  end do
  close(12)
  
  end program Jacobi

function Potential(x)
  implicit none
  real(8)::Potential,x

  Potential = x**2
  return

end function Potential
