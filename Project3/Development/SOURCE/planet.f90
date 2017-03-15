module planet_class
  implicit none

  type planet
     real(8) :: mass, potential, kinetic
     real(8), dimension(3)::position, velocity

     public::mass,potential,kinetic,position,velocity
     

end module planet_class

