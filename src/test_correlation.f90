program cross_correlation
use kind_parameters
use correlation
implicit none

type(ring_buffer_t) :: buffer

integer :: i
integer :: state(2)
real(dp) :: time

call initialize_buffer(buffer, 5, 0.10_dp)

do i = 1, 100
	time = 0.01_dp * i
	state = [i, 100 - i]
	call add_observation(buffer, time, state)
end do

call print_buffer(buffer)

end program
