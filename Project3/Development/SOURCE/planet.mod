module planet_class
  
  type planet
     real(8) :: pmass
     real(8) :: px
     real(8) :: py
     real(8) :: pvx
     real(8) :: pvy
     
  contains
    procedure, pass(this) ::setplanet
    procedure ::distance
    procedure ::GravitationalForce
    procedure ::Acceleration
    procedure ::KineticEnergy

 end type planet

contains
 include "planet.f90"
end module planet_class
