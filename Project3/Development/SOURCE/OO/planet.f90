 module planet_class
  
  type, public :: planet
     real(8) :: pmass
     real(8) :: px
     real(8) :: py
     real(8) :: pvx
     real(8) :: pvy
     
  contains
    procedure ::distance
    procedure ::GravitationalForce
    procedure ::Acceleration
    procedure ::KineticEnergy
 end type planet

contains

!!$  subroutine setplanet(this,mass,x,y,vx,vy)
!!$    implicit none
!!$    class (planet), intent(inout) :: this
!!$    real(8),intent(in)::mass,x,y,vx,vy
!!$    this%pmass=mass
!!$    this%px=x
!!$    this%py=y
!!$    this%pvx=vx
!!$    this%pvy=vy
!!$ 
!!$  end subroutine setplanet

  function distance (this, other)
    implicit none
    class (planet),intent(in)::this, other
    real(8)::x1,x2,y1,y2,xx,yy,distance

    x1=this%px
    x2=other%px
    y1=this%py
    y2=other%py

    xx=x1-x2
    yy=y1-y2

    distance=sqrt(xx*xx+yy*yy)
    return
  end function distance

  function GravitationalForce(this, other, Gconst)
    implicit none
    class (planet),intent(in)::this,other
    real(8)::r,Gconst,GravitationalForce
    
    r=distance(this,other)
    if (r/=0.0) then
       GravitationalForce=( Gconst*this%pmass*other%pmass)/(r*r)
    end if
    return
  end function GravitationalForce
  
  function Acceleration(this,other,Gconst)
    implicit none
    class(planet),intent(in)::this,other
    real(8)::r,Gconst,Acceleration
    r=distance(this,other)
    if(r /= 0.0) then
      Acceleration=GravitationalForce(this,other,Gconst)/(this%pmass)
    end if
    return
  end function Acceleration

  function KineticEnergy(this)
    implicit none
    class(planet), intent(in)::this
    real(8)::KineticEnergy
    KineticEnergy=0.5*this%pmass*(this%pvx*this%pvx + this%pvy*this%pvy)
    return
  end function KineticEnergy

  end module planet_class

  
