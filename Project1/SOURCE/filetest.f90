program filetest
	implicit none

	character (len=90) :: filename
	integer :: j

	do j = 31, 34
 		write (filename, '( "maltoLyo12per-reimage-set", I1, ".traj" )' )  j - 30
		OPEN(unit=j,status='old',file=filename)
		close (j)
	end do 
end program filetest