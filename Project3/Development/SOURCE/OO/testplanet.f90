
program test_ch2601
   use planet_class
   implicit none
   
   type (planet) :: mars = planet(100.d0,1,2,3,4)
   type (planet) :: earth = planet(100.d0,2,3,5,6)
   print *, ' get '
   print *, mars%pmass
   print*,mars%px
   print*,mars%py
   print*,mars%pvx
   print*,mars%pvy
   print*, distance(mars,earth)
   print*, gravitationalforce(mars,earth,1.d0)
   print*, Acceleration(mars,earth,1.d0)
   print*, KineticEnergy(earth)
 end program test_ch2601 
