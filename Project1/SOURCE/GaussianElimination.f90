program GaussElimination
	implicit none
	! Declare Variables
	real(8)::SourceValue
	real(8)::x,StepSize,x0,e
	integer, parameter :: NumGridPoints = 10
	real(8), Dimension(NumGridPoints+1) :: u,d,f
	integer:: i

	!Initialize Variables
	StepSize = 1.d0/(NumGridPoints)
	u(1)=0
	u(NumGridPoints+1)=0
	d(2)=2
	d(NumGridPoints)=2
	e=-1.d0
	do i=1,NumGridPoints+1
		f(i)=SourceValue((i-1)*StepSize)*StepSize**2.d0
		!print *, (i-1)*StepSize, f(i)
	end do

	!Forward Substitution
	do i=3, NumGridPoints
		d(i)= 2 - (e * e)/ d(i-1)
		f(i)=f(i) - (e * f(i-1))/ d(i-1)
		print *, f(i)
	end do

	!Backward Substitution
	u(NumGridPoints) = f(NumGridPoints)/d(NumGridPoints)

	do i=NumGridPoints-1, 2, -1
		u(i)= (f(i) - e * u(i+1) )/ d(i)
	end do

	do i=1,NumGridPoints+1
		print *, (i-1)*StepSize, u(i)
	end do


end program GaussElimination

function SourceValue(x)
	implicit none
	real(8)::SourceValue, x

	SourceValue = 100.d0 * Exp(- 10.d0 * x)
	return

end function SourceValue	