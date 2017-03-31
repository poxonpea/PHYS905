program solarsystem
!  use planet_class
  use solver_class

  implicit none

  integer,parameter::Bodies=10
  real(8)::Tmax=250.0,h
  integer::Num_Steps=10000,m
  real(8),dimension(10)::masses
  real(8),dimension(10,2)::position
  real(8),dimension(10,2)::velocity
  real(8),dimension(bodies,bodies,3)::relposition,updaterel
  real(8),dimension(bodies,2)::relforce
  real(8),dimension(bodies,2)::updatedforce
  real(8),dimension(bodies+1)::kinetic,potential,angular

  type(solver)::solar
  integer,parameter::numbodies=10
  
   
  solar%mass(1)=1.d0
  solar%mass(2)=3.d-6
  solar%mass(3)=9.5d-4
  solar%mass(4)=3.227d-7 
  solar%mass(5)=2.4478d-6
  solar%mass(6)=2.85812d-4
  solar%mass(7)=1.66012d-7 
  solar%mass(8)=4.36576d-5
  solar%mass(9)=5.1503d-5
  solar%mass(10)=6.583d-9 
  solar%position(1,1)=3.1647d-3
  solar%position(1,2)=4.4307d-3
  solar%velocity(1,1)=365.25*(-3.37957d-6)
  solar%velocity(1,2)=365.25*(6.606862d-6)
  solar%position(2,1)=-9.8825d-1
  solar%position(2,2)=8.49978d-2
  solar%velocity(2,1)=365.25*(-1.68024d-3)
  solar%velocity(2,2)=365.25*(-1.719988d-2)
  solar%position(3,1)=-5.23294d0
  solar%position(3,2)=-1.52515d0
  solar%velocity(3,1)=365.25*(2.022596d-3)
  solar%velocity(3,2)=365.25*(-6.88771645d-3)
  solar%position(4,1)=7.78069d-1
  solar%position(4,2)=1.2797d0
  solar%velocity(4,1)=365.25*(-1.143145d-2)
  solar%velocity(4,2)=365.25*(8.4664712d-3)
  solar%position(5,1)=-7.02894d-1
  solar%position(5,2)=1.359581d-1
  solar%velocity(5,1)=365.25*(-3.8130624d-3)
  solar%velocity(5,2)=365.25*(-1.9968d-2)
  solar%position(6,1)=-1.48071d0
  solar%position(6,2)=-9.935855d0
  solar%velocity(6,1)=365.25*(5.212138d-3)
  solar%velocity(6,2)=365.25*(-8.3942195d-4)
  solar%position(7,1)=2.805339d-1
  solar%position(7,2)=1.7274317d-1
  solar%velocity(7,1)=365.25*(-2.01015d-2)
  solar%velocity(7,2)=365.25*(2.5290758d-2)
  solar%position(8,1)=1.822435d1
  solar%position(8,2)=8.083455d0
  solar%velocity(8,1)=365.25*(-1.623364d-3)
  solar%velocity(8,2)=365.25*(3.411947d-3)
  solar%position(9,1)=2.8412218d1
  solar%position(9,2)=-9.4680088d0
  solar%velocity(9,1)=365.25*(9.7114038d-4)
  solar%velocity(9,2)=365.25*(2.99682d-3)
  solar%position(10,1)=9.890335d0
  solar%position(10,2)=-3.177864195d1
  solar%velocity(10,1)=365.25*(3.06860367d-3)
  solar%velocity(10,2)=365.25*(2.905788d-4)
!  print*,numbodies

  open(3,file="Earth.dat")
  open(4,file="Jupiter.dat")
  open(5,file="energy.dat")
  open(7,file="momentum.dat")
  open(8,file="Mars.dat")
  open(9,file="Venus.dat")
  open(10,file="Saturn.dat")
  open(11,file="Mercury.dat")
  open(12,file="Uranus.dat")
  open(13,file="Neptune.dat")
  open(14,file="Pluto.dat")


  h=Tmax/Num_Steps
  do m=1,Num_Steps

     call relative_position(solar,numbodies,relposition)
!     print*,'numbodies out 1',numbodies
     call forces(solar,Numbodies,relposition,relforce)
!     print*,'numbodies out 2',numbodies
     call calc_position(solar,numbodies,relforce,h)
!     print*,'numbodies out 3',numbodies
     call relative_position(solar,numbodies,updaterel)
!     print*,'numbodies out 4',numbodies
     call forces(solar,numbodies,updaterel,updatedforce)
!     print*,'numbodies out 5',numbodies
     call calc_velocities(solar,numbodies,relforce,updatedforce,h)
!     print*,'numbodies out 6',numbodies
     call kinetic_energy(solar,numbodies,kinetic)
     call potential_energy(solar,numbodies,updaterel,potential)
     call angular_momentum(solar,numbodies,updaterel,angular)
     
     write(5,*), kinetic(Numbodies+1), potential(numbodies+1), kinetic(Numbodies+1)+potential(numbodies+1)
     write(7,*), angular(Numbodies+1)
  
     write(3,*),solar%position(2,1),solar%position(2,2)
     write(4,*),solar%position(3,1),solar%position(3,2)
     write(8,*),solar%position(4,1),solar%position(4,2)     
     write(9,*),solar%position(5,1),solar%position(5,2)
     write(10,*),solar%position(6,1),solar%position(6,2)
     write(11,*),solar%position(7,1),solar%position(7,2)
     write(12,*),solar%position(8,1),solar%position(8,2)
     write(13,*),solar%position(9,1),solar%position(9,2)
     write(14,*),solar%position(10,1),solar%position(10,2)

  end do

  
  close(3)
  close(4)
  close(5)
  close(7)
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
end program solarsystem
