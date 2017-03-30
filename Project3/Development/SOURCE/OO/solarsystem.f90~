program threebody
!  use planet_class
  use solver_class

  implicit none

  integer,parameter::bodies=3
  real(8)::Tmax=12.0,h
  integer::Num_Steps=120,m
  real(8),dimension(bodies),parameter::masses(1)=1.d0
  real(8),dimension(3,2),parameter::position(1,1)=0.d0
  real(8),dimension(3,2),parameter::velocity(1,1)=1.d0
  real(8),dimension(bodies,bodies,3)::relposition,updaterel
  real(8),dimension(bodies,2)::relforce
  real(8),dimension(bodies,2)::updatedforce
  real(8),dimension(bodies+1)::kinetic,potential,angular

  type(solver)::jupiter=solver(3,masses,position,velocity)
  integer::numbodies
  
        
  jupiter%mass(1)=1.d0
  jupiter%mass(2)=3.d-6
  jupiter%mass(3)=9.5d-4
  jupiter%position(1,1)=0.0
  jupiter%position(1,2)=0.0
  jupiter%velocity(1,1)=0.0
  jupiter%velocity(1,2)=0.0
  jupiter%position(2,1)=-9.8825d-1
  jupiter%position(2,2)=8.49978d-2
  jupiter%velocity(2,1)=365.25*(-1.68024d-3)
  jupiter%velocity(2,2)=365.25*(-1.719988d-2)
  jupiter%position(3,1)=-5.23294d0
  jupiter%position(3,2)=-1.52515d0
  jupiter%velocity(3,1)=365.25*(2.022596d-3)
  jupiter%velocity(3,2)=365.25*(-6.88771645d-3)
  
  Numbodies=bodies

  open(3,file="Earth.dat")
  open(4,file="Jupiter.dat")
  open(5,file="energy.dat")
  open(7,file="momentum.dat")



  h=Tmax/Num_Steps
  do m=1,Num_Steps

     call relative_position(jupiter,numbodies,relposition)
!     print*,'exit rel_pos 1',relposition(2,3,1),relposition(1,2,1)
     call forces(jupiter,numbodies,relposition,relforce)
!     print*,'exit forces 1', relforce(2,1),relforce(3,1)
     call calc_position(jupiter,numbodies,relforce,h)
!     print*, 'exit calc_pos',jupiter%position(2,1),jupiter%position(2,2)
     call relative_position(jupiter,numbodies,updaterel)
!     print*, 'exit rel_pos 2',updaterel(2,3,1),updaterel(1,2,1)
     call forces(jupiter,numbodies,updaterel,updatedforce)
!     print*, 'exit forces 2'
     call calc_velocities(jupiter,numbodies,relforce,updatedforce,h)
!     print*, 'exit calc_vel'
     call kinetic_energy(jupiter,numbodies,kinetic)
     call potential_energy(jupiter,numbodies,updaterel,potential)
     call angular_momentum(jupiter,numbodies,updaterel,angular)
     
     write(5,*), kinetic(Numbodies+1), potential(numbodies+1), kinetic(Numbodies+1)+potential(numbodies+1)
     write(7,*), angular(Numbodies+1)
  
     write(3,*),jupiter%position(2,1),jupiter%position(2,2)
     write(4,*),jupiter%position(3,1),jupiter%position(3,2)

  end do

  close(3)
  close(4)
  close(5)
  close(7)
end program threebody

