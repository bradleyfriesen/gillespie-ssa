module random_utils
use kind_parameters
implicit none


! Repeatable seed
integer, parameter :: seed_value = 12345

public :: initialize_random


contains


subroutine initialize_random()
	integer, allocatable :: seed(:)
	integer :: n

	call random_seed(size=n)
	allocate(seed(n))

	seed = seed_value
	call random_seed(put=seed)
	print *, seed

	deallocate(seed)
end subroutine initialize_random


subroutine random_exponential(rate, x)
	real(dp), intent(in) :: rate
	real(dp), intent(out) :: x
	real(dp) :: u

	call random_number(u)

	do while (u == 0._dp)
		call random_number(u)
	end do

	x = -log(u) / rate
	
end subroutine random_exponential

end module
