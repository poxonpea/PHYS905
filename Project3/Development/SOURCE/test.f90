module shape_module

  type shape_type
     integer :: x_ = 0
     integer :: y_ = 0

   contains
     procedure, pass(this) ::getx
     procedure, pass(this) ::gety
     procedure, pass(this) ::setx
     procedure, pass(this) ::sety
     procedure, pass(this) ::moveto
     procedure, pass(this) ::draw

  end type shape_type

contains
integer function getx(this)
   implicit none
   class (shape_type), intent (in) :: this
   getx = this%x_
 end function getx

 integer function gety(this)
   implicit none
   class (shape_type), intent (in) :: this
   gety = this%y_
 end function gety

 subroutine setx(this,x)
   implicit none 
   class (shape_type), intent (inout) :: this
   integer, intent (in) :: x
   this%x_ = x
 end subroutine setx
 
 subroutine sety(this,y)
   implicit none
   class (shape_type), intent (inout) :: this
   integer, intent (in) :: y
   this%y_ = y
 end subroutine sety

 subroutine moveto(this,newx, newy)
   implicit none
   class (shape_type), intent (inout) :: this
   integer, intent (in) :: newx
   integer, intent (in) :: newy
   this%x_ = newx
   this%y_ = newy
 end subroutine moveto

 subroutine draw(this)
   implicit none
   class (shape_type), intent (in) :: this
   print *, ' x = ', this%x_
   print *, ' y = ', this%y_
 end subroutine draw

end module shape_module

 program test_ch2601
   use shape_module
   implicit none
   type (shape_type) :: s1 = shape_type(10,20)
   integer :: x1 = 100
   integer :: y1 = 200
   print *, ' get '
   print *, s1%getx(), ' ', s1%gety()
   print *, ' draw '
   call s1%draw()
   print *, ' moveto '
   call s1%moveto (x1,y1)
   print *, ' draw '
   call s1%draw()
   print *, ' set '
   call s1%setx(99)
   call s1%sety(99)
   print *, ' draw'
   call s1%draw()
 end program test_ch2601 
