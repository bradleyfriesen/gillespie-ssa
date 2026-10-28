module kinematics
use kind_parameters
implicit none


contains


pure function hill(x, K, n) result(h)
	real(dp), intent(in) :: x, K, n
	real(dp) :: h
	
	if (x == 0._dp) then
		if (n > 0._dp) then
			h = 0._dp
		else if (n < 0._dp) then
			h = 1._dp
		else
			h = 0.5_dp
		end if
	else
		h = 1._dp / (1._dp + (K/x)**n)
	end if
	
end function hill
end module
