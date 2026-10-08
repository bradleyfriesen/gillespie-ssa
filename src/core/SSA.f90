module SSA

use kind_parameters
use model_interface
implicit none

type :: simulation_t
	real(dp) :: time
	real(dp) :: dt
	integer :: reaction
	integer, allocatable :: state(:)
	integer, allocatable :: reaction_count(:)
	real(dp), allocatable :: propensity(:)
end type


contains


subroutine initialize_SSA(sim, model, initial_state)
	type(simulation_t), intent(out) :: sim
	class(model_t), intent(in) :: model
	integer, intent(in) :: initial_state(:)
	
	
	allocate(sim%state(model%n_species))
	allocate(sim%reaction_count(model%n_reactions))
	allocate(sim%propensity(model%n_reactions))
	
	sim%time = 0._dp
	sim%dt = 0._dp
	sim%reaction = 0
	sim%state = initial_state
	sim%reaction_count = 0
	sim%propensity = 0
end subroutine


subroutine SSA_roll(sim, model)
	type(simulation_t), intent(inout) :: sim
	class(model_t), intent(in) :: model

	real(dp) :: r_total
	real(dp) :: r(model%n_reactions)
	real(dp) :: tau
	real(dp) :: u
	real(dp) :: threshold
	integer :: i
	
	call model%propensity(sim%state, r)
	sim%propensity = r
	
	r_total = sum(r)
	
	if (r_total <= 0.0_dp) then
		error stop "All reaction rates are zero"
	end if
	
	call random_number(u)
	u = max(u, tiny(u))
	
	sim%dt = -log(u) / r_total
	
	call random_number(u)
	threshold = u * r_total
	
	sim%reaction = model%n_reactions
	
	do i = 1, model%n_reactions
		threshold = threshold - r(i)
		
		if (threshold <= 0.0_dp) then
			sim%reaction = i
			exit
		end if
    end do
end subroutine SSA_roll


subroutine SSA_update(sim, model)
	type(simulation_t), intent(inout) :: sim
	class(model_t), intent(in) :: model
	
	sim%time = sim%time + sim%dt
	sim%state = sim%state + model%stoich(:, sim%reaction)
	sim%reaction_count(sim%reaction) = &
		sim%reaction_count(sim%reaction) + 1
end subroutine SSA_update


subroutine SSA_step(sim, model)
	type(simulation_t), intent(inout) :: sim
	class(model_t), intent(in) :: model
	
	real(dp) :: r_total
	real(dp) :: r(model%n_reactions)
	real(dp) :: tau
	real(dp) :: u
	real(dp) :: threshold
	integer :: reaction
	integer :: i
	
	call model%propensity(sim%state, r)
	
	r_total = sum(r)
	
	if (r_total <= 0.0_dp) then
		error stop "All reaction rates are zero"
	end if
	
	call random_number(u)
	u = max(u, tiny(u))
	
	tau = -log(u) / r_total
	
	call random_number(u)
	threshold = u * r_total
	
	reaction = model%n_reactions

	do i = 1, model%n_reactions
		threshold = threshold - r(i)    
		if (threshold <= 0.0_dp) then
			reaction = i
			exit
		end if
    end do
    
	sim%time = sim%time + tau
	sim%state = sim%state + model%stoich(:, reaction)
	sim%reaction_count(reaction) = sim%reaction_count(reaction) + 1
end subroutine SSA_step
	

end module
