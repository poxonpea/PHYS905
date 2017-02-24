program Jacobi

  implicit none
  integer:: Num
  integer, parameter :: testnum=4
  real(8) ::rhoMax
  real(8), dimension(:,:),Allocatable::a,r
  real(8), dimension(:),Allocatable::eigenval
  real(8)::Stepsize,start,finish,potential,testmax,max,tolerance,test,omega,analyticsoln
  integer:: i,j,maxi,maxj,count,case,NumElectrons,verbose

  !Ask the user for information
  print *, "Do you want to run for a test matrix, or electrons in well? 1)test 2)electrons"
  read(*,*) case
  if (case == 2) then
     print*, "Are you interested in the case of one electron or two interacting electrons? 1) one 2) two"
     read(*,*) NumElectrons
     if (NumElectrons == 2) then
        print*,"what value do you want for omega?"
        read(*,*) omega
     else
        omega = 0.d0
     end if
  
     print*, "what is rho max?"
     read(*,*) rhoMax
     print*, "what is N?"
     read(*,*)Num
  end if
  print *, "Do you want me to print the orthogonality results? 1)yes 2)no"
  read(*,*) verbose 

  if (case==1) then
     Num = testnum + 1
  end if
  
  Allocate(a(Num-1,Num-1))
  Allocate(r(Num-1,Num-1))
  Allocate(eigenval(Num-1))
  
  ! Fill in test Matrix
  if (case==1) then
     a(1,1)=1.d0
     a(1,2)=2.d0
     a(1,3)=3.d0
     a(1,4)=4.d0
     a(2,1)=2.d0
     a(2,2)=6.d0
     a(2,3)=7.d0
     a(2,4)=8.d0
     a(3,1)=3.d0
     a(3,2)=7.d0
     a(3,3)=11.d0
     a(3,4)=12.d0
     a(4,1)=4.d0
     a(4,2)=8.d0
     a(4,3)=12.d0
     a(4,4)=16.d0
     print *, "The test matrix is "
     do j=1,testnum
        print*, a(j,1), a(j,2), a(j,3), a(j,4)
     end do
  end if
  
  if (case==2) then
     StepSize=rhoMax/Num
  
     !Fill in real Matrix
     do i =1,Num-1
        do j=1,Num-1
           if (i==j) then
              a(j,i)=2.d0/(stepsize**2) + Potential(i*StepSize,NumElectrons,omega)
           else if (j==i+1) then
              a(j,i)=-1.d0/(stepsize**2)
           else if (j==i-1) then
              a(j,i)=-1.d0/(stepsize**2)
           else
              a(j,i)=0.d0
           end if
        end do
     end do
  end if
  
  do i=1,Num-1
     do j=1,Num-1
        if (i == j) then
           r(j,i)=1.d0
        else
           r(j,i)=0.d0
        end if
     end do
  end do
  
  
count = 0
max = 1.d0
tolerance=1.d-8

call cpu_time(start)

do while (abs(max) .gt. tolerance)
  count=count+1

  !Find largest matrix element in array
  call findmax(Num-1,a,maxj,maxi,max)

  !Jacobi's algorithm to rotate largest element to 0
  call rotate(Num-1,a,r,maxj,maxi)

  !Test to check orthogonality of eigenvectors after every 10th iteration
  if (Mod(count,10) == 0) then
     call OrthoCheck(Num-1,r,maxi,maxj,count,verbose)
  end if
  
end do

call cpu_time(finish)
print*, 'cpu time',finish-start

max = 1000.d0

! Put eigenvalues in array and pick out the lowest value
  do i=1,Num-1
     eigenval(i)=a(i,i)
     if (eigenval(i) .lt. max)then
        max = eigenval(i)
        maxi = i
     end if   
  end do

  print*, "lowest eigenvalue",eigenval(maxi),maxi
  
  if (case == 1) then
     print*, "Eigenvectors from Matematica: 31.3397, 1.70691, 0.953375, 6.39959E-17"
     print*, "Eigenvectors from Jacobi: ", eigenval(1), eigenval(2), eigenval(3), eigenval(4)
  end if

  !write wavefunctions and eigenvalues to file

  open(12,file="wf.dat")
  do i=1,Num-1
     write(12,*), i*stepsize,r(i,maxi)
  end do
  close(12)

  open(13,file="Eigenvalue.dat/Subroutine")
  do i=1,Num-1
    write(13,*),eigenval(i)
  end do


  print*, "Number of Iterations", count
  
  Deallocate(a)
  Deallocate(r)
  Deallocate(eigenval)

end program Jacobi



function Potential(x,NumElectrons,omega)
  implicit none
  real(8)::Potential,x,omega
  integer:: NumElectrons

  if (NumElectrons == 1) then
     Potential = x*x
  else if (NumElectrons ==2) then
     Potential = omega*omega*x*x + 1.d0/x
  end if
  
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

 return
end subroutine findmax


subroutine rotate(testnum,testmat,testr,maxj,maxi)
  implicit none
  integer::i,j,maxi,maxj,testnum,q
  real(8),dimension(testnum,testnum)::testmat,testr
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

subroutine OrthoCheck(Num,r,maxi,maxj,count,verbose)
  implicit none
  integer:: Num,i,maxi,maxj,count,verbose
  real(8), dimension(Num,Num):: r
  real(8), dimension(Num)::col1,col2
  real(8) :: check

  do i=1,Num
     col1(i)=r(i,maxi)
     col2(i)=r(i,maxj)
  end do

! Calculate the dot product of two vectors
  check = dot_product(col1,col2)
  if (mod(count,10)==0) then
     if (verbose == 1) then
      print*,"Iteration",count,"Dot product of eigenvectors =", check
     end if
  end if
  if (abs(check) .gt. 1.d-8) then
     print*, "Error: Orthogonality violated"
     stop
  end if
end subroutine OrthoCheck
