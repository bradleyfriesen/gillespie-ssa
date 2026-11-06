module correlation
! From the Gillespie algorithm, calculate correlations.
use kind_parameters
implicit none


type :: observation_t
	real(dp) :: time
	integer, allocatable :: state(:)
end type observation_t

type :: ring_buffer_t
	type(observation_t), allocatable :: observations(:)
	integer :: head = 0
	integer :: obs_count = 0
	real(dp) :: maximum_lag
end type ring_buffer_t

type :: correlation_t
	type(ring_buffer_t) :: buffer
	
	real(dp) :: maximum_lag
	real(dp) :: accumulation_time = 0.0_dp
	
	integer :: n
	real(dp), allocatable :: lag(:)
    
end type correlation_t


contains


subroutine initialize_buffer(buffer, initial_capacity, maximum_lag)
	type(ring_buffer_t), intent(out) :: buffer
	integer, intent(in) :: initial_capacity
	real(dp), intent(in) :: maximum_lag
	
	! Stop if we set up the program incorrectly.
	if (initial_capacity < 1) error stop "Invalid buffer capacity"
	if (maximum_lag < 0._dp) error stop "Invalid maximum lag"
	
	allocate(buffer%observations(initial_capacity))
	buffer%head = 0
	buffer%obs_count = 0
	buffer%maximum_lag = maximum_lag
	
end subroutine initialize_buffer


subroutine add_observation(buffer, time, state)
	type(ring_buffer_t), intent(inout) :: buffer
	real(dp), intent(in) :: time
	integer, intent(in) :: state(:)
	
	integer :: next
	
	if (buffer%obs_count > 0) then
		! Check time is increasing to ensure no head bugs
		if (time <= buffer%observations(buffer%head)%time) then
			error stop "Observation times must increase"
		end if
		
		! Retain only the history needed for maximum_lag,
		! plus the observation immediately preceding that history.
		call discard_old(buffer, time - buffer%maximum_lag)
	end if
	
	! Grow buffer if necessary.
	if (buffer%obs_count == size(buffer%observations)) then
		call grow_buffer(buffer)
	end if
	
	! Insert the new observation.
	next = modulo(buffer%head, size(buffer%observations)) + 1
	buffer%head = next
	buffer%obs_count = buffer%obs_count + 1

	buffer%observations(next)%time = time
	buffer%observations(next)%state = state
	
end subroutine add_observation


subroutine discard_old(buffer, cutoff)
! Discard old points by removing them from the useful obs_count
	type(ring_buffer_t), intent(inout) :: buffer
	real(dp), intent(in) :: cutoff

	integer :: oldest, second_oldest
	integer :: capacity

	capacity = size(buffer%observations)

	do while (buffer%obs_count > 1)
		oldest = oldest_index(buffer)
		second_oldest = modulo(oldest, capacity) + 1
		
		! "Forget" unneeded extra points by lowering obs_count
		if (buffer%observations(second_oldest)%time <= cutoff) then
			buffer%obs_count = buffer%obs_count - 1
		else
			exit
		end if
	end do
	
end subroutine discard_old


pure function oldest_index(buffer) result(idx)
! Find oldest index. One greater than head, wrapped with obs_count.
	type(ring_buffer_t), intent(in) :: buffer
	integer :: idx

	idx = modulo(buffer%head - buffer%obs_count, &
		size(buffer%observations)) + 1

end function oldest_index

	
subroutine grow_buffer(buffer)
! For continuous time, grow the buffer if necessary. In continuous time
! we record events, where the timing can vary. Thus we might need to
! add a few positiong in the buffer. Allocation is expensive, so try
! not to do it often by making a good size guess.
	type(ring_buffer_t), intent(inout) :: buffer
	
	type(observation_t), allocatable :: temp(:)
	integer :: old_capacity, new_capacity
	integer :: i, idx
	integer :: estimated_capacity
	real(dp) :: buffer_span, average_interval
	
	old_capacity = size(buffer%observations)
	
	buffer_span = buffer%observations(buffer%head)%time &
		- buffer%observations(oldest_index(buffer))%time
	
	! Estimate the number of observations needed to cover maximum_lag.
	if (buffer%obs_count >= 2 .AND. buffer_span > 0.0_dp) then
		average_interval = 1._dp * buffer_span / buffer%obs_count
	
		estimated_capacity = ceiling( &
			1.1_dp * buffer%maximum_lag / average_interval) + 2
	
		! If the estimate is insufficient, fall back to doubling.
		if (estimated_capacity > old_capacity) then
			new_capacity = estimated_capacity
		else
			new_capacity = 2 * old_capacity
		end if
	
	else
		new_capacity = 2 * old_capacity
	end if
	
	allocate(temp(new_capacity))
	
	! Copy observations into chronological order.
	do i = 1, buffer%obs_count
		idx = modulo(oldest_index(buffer) + i - 2, old_capacity) + 1
		temp(i) = buffer%observations(idx)
	end do
	
	call move_alloc(temp, buffer%observations)
	
	! The newest observation is now at the end of the copied data.
	buffer%head = buffer%obs_count
	
end subroutine grow_buffer


subroutine print_buffer(buffer)
	type(ring_buffer_t), intent(in) :: buffer

	integer :: i, idx, oldest

	! if (buffer%count == 0) return

	oldest = oldest_index(buffer)

	print *, "Capacity: ", size(buffer%observations)
	print *, "Observation Count:    ", buffer%obs_count
	print *, "------------------------------"

	do i = 1, buffer%obs_count
		idx = modulo(oldest + i - 2, size(buffer%observations)) + 1

		print *, idx, buffer%observations(idx)%time, &
				 buffer%observations(idx)%state
	end do
	
end subroutine print_buffer

end module correlation
