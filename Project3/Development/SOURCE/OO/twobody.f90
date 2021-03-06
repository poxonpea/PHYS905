program twobody

  use solver_class

  implicit none

  integer,parameter::numbodies=2
  real(8)::Tmax=1.0,h
  integer::Num_Steps=10,m
  real(8),dimension(2)::masses
  real(8),dimension(2,2)::position
  real(8),dimension(2,2)::velocity
  real(8),dimension(numbodies,numbodies,2)::relposition,updaterel
  real(8),dimension(numbodies,2)::relforce
  real(8),dimension(numbodies,2)::updatedforce
  real(8),dimension(numbodies+1)::kinetic,potential,angular

  type(solver)::earthsun
  
  
  earthsun%mass(1)=1.d0
  earthsun%mass(2)=3.d-6
!  earthsun%position(2,1)=-0.988251
!  earthsun%position(2,2)=0.08499779
!  earthsun%velocity(2,1)=(365.25*(-0.00168024))
!  earthsun%velocity(2,2)=(365.25*(-0.01719988))
  earthsun%position(2,1)=1.0
  earthsun%position(2,2)=0.0
  earthsun%velocity(2,1)=0
  earthsun%velocity(2,2)=3.54
  earthsun%position(1,1)=0.d0
  earthsun%position(1,2)=0.d0
  earthsun%velocity(1,1)=0.d0
  earthsun%velocity(1,2)=0.d0

  open(3,file="OOVelVer.dat")
  open(4,file="energy.dat")
  open(5,file="momentum.dat")

  h=Tmax/Num_Steps
  do m=1,Num_Steps
     
     call relative_position(earthsun,numbodies,relposition)
     call forces(earthsun,numbodies,relposition,relforce)
     call calc_position(earthsun,numbodies,relforce,h)
     call relative_position(earthsun,numbodies,updaterel)
     call forces(earthsun,numbodies,updaterel,updatedforce)
     call calc_velocities(earthsun,numbodies,relforce,updatedforce,h)
     call kinetic_energy(earthsun,numbodies,kinetic)
     call potential_energy(earthsun,numbodies,updaterel,potential)
     call angular_momentum(earthsun,numbodies,updaterel,angular)
     
     write(4,*), kinetic(Numbodies+1), potential(numbodies+1)
     write(5,*), angular(Numbodies+1)
  
  write(3,*),earthsun%position(2,1), earthsun%position(2,2)

  end do

  close(3)
  close(4)
  close(5)
end program twobody
