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

  subroutine setplanet(this,mass,x,y,vx,vy)
    implicit none
    class (planet), intent(inout) :: this
    real(8),intent(in)::mass,x,y,vx,vy
    this%pmass=mass
    this%px=x
    this%py=y
    this%pvx=vx
    this%pvy=vy
 
  end subroutine setplanet

  function distance (p1,p2)
    implicit none
    class (planet),intent(in)::p1,p2
    real(8)::x1,x2,y1,y2,xx,yy,distance

    x1=p1%px
    x2=p2%px
    y1=p1%py
    y2=p2%py

    xx=x1-x2
    yy=y1-y2

    distance=sqrt(xx*xx+yy*yy)
    return
  end function distance

  function GravitationalForce(p1,p2,Gconst)
    implicit none
    class (planet),intent(in)::p1,p2
    real(8)::r,Gconst,GravitationalForce
    
    r=distance(p1,p2)
    if (r/=0.0) then
       GravitationalForce=( Gconst*p1%pmass*p2%pmass)/(r*r)
    end if
    return
  end function GravitationalForce
  
  function Acceleration(p1,p2,Gconst)
    implicit none
    class(planet),intent(in)::p1,p2
    real(8)::r,Gconst,Acceleration
    r=distance(p1,p2)
    if(r /= 0.0) then
      Acceleration=GravitationalForce(p1,p2,Gconst)/(p1%pmass)
    end if
    return
  end function Acceleration

  function KineticEnergy(p1)
    implicit none
    class(planet), intent(in)::p1
    real(8)::KineticEnergy
    KineticEnergy=0.5*p1%pmass*(p1%pvx*p1%pvx + p1%pvy*p1%pvy)
    return
  end function KineticEnergy

  end module planet_class

  