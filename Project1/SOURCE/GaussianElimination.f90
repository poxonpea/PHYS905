program GaussElimination
	implicit none
	real(8)::SourceValue
	real(8)::Fvalue,x,StepSize
	integer, parameter :: NumGridPoints = 10

	StepSize = 1.d0/NumGridPoints
	


!	x=2.d0
!	Fvalue = EquationValue(x)
!	print*, Fvalue


end program GaussElimination

function SourceValue(x)
	implicit none
	real(8)::SourceValue, x

	SourceValue = 100.d0 * Exp(10.d0 * x)
	return

end function EquationValue	