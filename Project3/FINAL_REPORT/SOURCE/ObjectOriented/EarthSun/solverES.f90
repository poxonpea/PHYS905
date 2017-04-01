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
    procedure ::calc_velocities
    procedure ::kinetic_energy
    procedure ::potential_energy
    procedure ::angular_momentum

 end type solver

contains
  subroutine relative_position(system,numbodies,relposition)
    implicit none

    class (solver),intent(in)::system
    integer :: i,j
    integer::NumBodies
    real(8),dimension(Numbodies,Numbodies,3):: relposition
    
    do i=2,NumBodies
       do j =1,NumBodies
          if (i/=j) then
             relposition(j,i,1)=system%position(j,1)-system%position(i,1)
!             print*,'relx',i,j,relposition(j,i,1)
             relposition(j,i,2)=system%position(j,2)-system%position(i,2)
!             print*,'rely',i,j,system%position(j,2),system%position(i,2),relposition(j,i,2)
             relposition(j,i,3)=sqrt(relposition(j,i,1)*relposition(j,i,1)+relposition(j,i,2)*relposition(j,i,2))
!             print*,"inside",relposition(j,i,2),relposition(j,i,3)
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

    Fourpi2 = 4.d0*3.14*3.14

    do i=1,numbodies
       do j=1,2
          relforce(i,j)=0.d0
       end do
    end do
    
    
    do i=2,numbodies
       do j=1,numbodies
          if(j.ne.i) then
             rrr=(relposition(j,i,3)**3.d0)
             relforce(i,1) =relforce(i,1) - Fourpi2*system%mass(i)*relposition(j,i,1)/rrr
             relforce(i,2) =relforce(i,2) - Fourpi2*system%mass(i)*relposition(j,i,2)/rrr
!             print*,'position and force',i,j,relposition(j,i,2),relforce(i,2)
          end if
       end do
    end do

  end subroutine forces


  subroutine calc_position(system,Numbodies,relforce,h)
    implicit none
    class(solver), intent(inout)::system
    integer::Numbodies,i,j
    real(8),dimension(Numbodies,2)::relforce
    real(8)::h

    do i=2,numbodies
       system%position(i,1)=system%position(i,1)+h*system%velocity(i,1) - (h*h/2.d0)*relforce(i,1)
!       print*,system%position(i,1),system%velocity(i,1)
       system%position(i,2)=system%position(i,2)+h*system%velocity(i,2) - (h*h/2.d0)*relforce(i,2)
!       print*,"newy",relforce(2,2)
    end do

  end subroutine calc_position

  subroutine calc_velocities(system,numbodies,relforce,updatedforce,h)
    implicit none
    class(solver),intent(inout)::system
    integer::NumBodies,i,j
    real(8)::h
    real(8),dimension(Numbodies,2)::relforce
    real(8),dimension(Numbodies,2)::updatedforce

    do i=2, numbodies
       system%velocity(i,1)=system%velocity(i,1)-h*0.5*updatedforce(i,1)-h*0.5*relforce(i,1)
       system%velocity(i,2)=system%velocity(i,2)-h*0.5*updatedforce(i,2)-h*0.5*relforce(i,2)
!       print*,system%velocity(2,1)
    end do

  end subroutine calc_velocities
    
  subroutine kinetic_energy(system,numbodies,kinetic)
    implicit none
    class(solver), intent(in)::system
    integer::Numbodies,i
    real(8),dimension(numbodies+1)::kinetic
    real(8)::totalke

    totalke=0
    
    do i =1,numbodies
       kinetic(i)=(1.d0/2.d0)*system%mass(i)*(system%velocity(i,1)*system%velocity(i,1)+system%velocity(i,2)*system%velocity(i,2))
       totalke=totalke+kinetic(i)
    end do

    kinetic(numbodies+1)=totalke
  end subroutine kinetic_energy

  subroutine potential_energy(system, numbodies,relposition,potential)
    implicit none
    class(solver), intent(in)::system
    integer::Numbodies,i
    real(8),dimension(numbodies+1)::potential
    real(8),dimension(numbodies,numbodies,3)::relposition
    real(8)::totalpe

    totalpe=0

    do i=2,numbodies
       potential(i)=4*3.14*3.14*system%mass(i)/relposition(1,i,3)

       totalpe=totalpe+potential(i)
    end do

    potential(numbodies+1)=totalpe
  end subroutine potential_energy

  subroutine angular_momentum(system,numbodies,relposition,angular)
    implicit none
    class(solver),intent(in)::system
    integer::Numbodies,i
    real(8),dimension(numbodies+1)::angular
    real(8),dimension(numbodies,numbodies,3)::relposition
    real(8)::totalang

    totalang=0

    do i=2,numbodies
       angular(i)=system%mass(i)*relposition(1,i,3)*(sqrt(system%velocity(i,1)**2.d0 + system%velocity(i,2)**2.d0))

       totalang=totalang+angular(i)
    end do

    angular(numbodies+1)=totalang
  end subroutine angular_momentum
    
  
end module solver_class
