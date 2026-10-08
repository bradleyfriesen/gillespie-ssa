program example
use kind_parameters
use command_line
use model_interface
use model_factory
use statistics
use SSA
implicit none

integer, parameter :: event_min = 10**6

class(model_t), allocatable :: model
type(simulation_t) :: sim
type(moments_t) :: moments

character(len=:), allocatable :: model_name
character(len=:), allocatable :: parameter_name

integer :: i
integer, allocatable :: initial_state(:)

! Command line arguments
call get_model_and_parameters(model_name, parameter_name)

! Make the model
call create_model(trim(model_name), trim(parameter_name), model)
! Output model parameters for human
call model%print_model()

! Set initial state and prepare the Gillespie algorithm
allocate(initial_state(model%n_species))
initial_state = 0
call initialize_SSA(sim, model, initial_state)

! Initialize statistics
call moments%initialize(model%n_species)

! Run the algorithm one step at a time
do while (minval(sim%reaction_count) < event_min)
	call SSA_roll(sim, model)
	call moments%update(sim%state, sim%dt)
	call SSA_update(sim, model)
end do

print *, "Simulation time: ", sim%time
print *, "State:           ", sim%state
print *, "Reaction counts: ", sim%reaction_count
print *, "Mean(s): ", moments%mean
print *, "Covariance matrix: ", moments%covariance()

end program
