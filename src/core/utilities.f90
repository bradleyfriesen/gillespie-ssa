module utilities
! Useful but miscellaneous
use kind_parameters
implicit none


contains


pure function relative_change(x, y) result(f)
	real(dp), intent(in) :: x, y
	real(dp) :: f
	f = abs((y - x) / (0.5_dp * (x + y)))
end function relative_change


end module utilities
