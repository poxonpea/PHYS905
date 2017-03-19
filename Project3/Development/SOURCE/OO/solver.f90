module solver_class
  
  type solver
     integer::Bodies=2
     real(8),dimension(2) :: mass
     real(8),dimension(2,2) :: position
     real(8),dimension(2,2) :: velocity
  contains
    procedure ::relative_position
    procedure ::forces
    procedure ::calc_position

 end type solver

contains
  subroutine relative_position(system,Numbodies,relposition)
    implicit none

    class (solver),intent(in)::system
    integer :: i,j
    integer::NumBodies
    real(8),dimension(:,:,:),Allocatable :: relposition
    NumBodies = system%Bodies

    Allocate(relposition(NumBodies,Numbodies,3))
    
    do i=1,NumBodies
       do j =1,NumBodies
          if (i/=j) then
             relposition(j,i,1)=system%position(j,1)-system%position(i,1)
             relposition(j,i,2)=system%position(j,2)-system%position(i,2)
             relposition(j,i,3)=sqrt(relposition(j,i,1)**2.d0+relposition(j,i,2)**2.d0)
          end if
       end do
    end do
    
    return 
  end subroutine relative_position


  subroutine forces(system,Numbodies,relposition,relforce)
    implicit none

    class(solver),intent(in)::system
    integer :: j,i
    integer:: NumBodies
    real(8),dimension(Numbodies,2)::relforce
    real(8),dimension(Numbodies,Numbodies,3):: relposition
    real(8)::rrr,Fourpi2

    Fourpi2 = 4.d0*3.14159265*3.14159265

    do i=1,numbodies
       do j=1,numbodies
          if(j/=i) then
             rrr=(relposition(j,i,3)**3.d0)
             relforce(i,1) =relforce(i,1) - Fourpi2*system%mass(j)*relposition(i,j,1)/rrr
             relforce(i,2) =relforce(i,2) - Fourpi2*system%mass(j)*relposition(i,j,2)/rrr
          end if
       end do
    end do

  end subroutine forces


  subroutine calc_position(system,Numbodies,relforce,h)
    implicit none
    class(solver), intent(in)::system
    integer::Numbodies
    real(8),dimension(Numbodies,2)::relforce
    real(8)::h

    
end module solver_class
