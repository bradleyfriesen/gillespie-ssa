module oscillating_model
use kind_parameters
use model_interface
use kinematics, only:hill
implicit none

integer, parameter :: n_species = 3
integer, parameter :: n_reactions = 6

type :: params_t
	real(dp) :: lambda(3)
	real(dp) :: beta(3)
	real(dp) :: K(3)
	real(dp) :: n(3)
end type

type, extends(model_t) :: oscillating_t
	type(params_t) :: params
contains
	procedure :: propensity => oscillating_propensity
	procedure :: print_model => oscillating_print_model
end type


contains


subroutine initialize_oscillating(model, parameter_file)
	type(oscillating_t), intent(out) :: model
	character(len=*), intent(in) :: parameter_file
	real(dp), dimension(3) :: lambda, beta, K, n

	namelist /parameters/ lambda, beta, K, n

	model%n_species = n_species
	model%n_reactions = n_reactions

	allocate (model%stoich(model%n_species, model%n_reactions))

	model%stoich = reshape([ &
		1, 0, 0, &
		-1, 0, 0, &
		0, 1, 0, &
		0, -1, 0,&
		0, 0, 1, &
		0, 0, -1&
		], shape(model%stoich))

	open (unit=10, file=parameter_file, status="old", action="read")
	read (10, nml=parameters)
	close (10)
	
	model%params%lambda = lambda
	model%params%beta = beta
	model%params%K = K
	model%params%n = n

end subroutine initialize_oscillating


subroutine oscillating_propensity(self, state, r)
	class(oscillating_t), intent(in) :: self
	integer, intent(in) :: state(:)
	real(dp), intent(out) :: r(:)

	associate(x => real(state(:), dp), &
		lambda => self%params%lambda, &
		beta => self%params%beta, &
		K => self%params%K, &
		n => self%params%n)
		
		r(1) = lambda(1) * hill(x(3), K(1), n(1))
		r(2) = x(1) * beta(1)
		r(3) = lambda(2) * hill(x(1), K(2), n(2))
		r(4) = x(2) * beta(2)
		r(5) = lambda(3) * hill(x(2), K(3), n(3))
		r(6) = x(3) * beta(3)
	end associate

end subroutine oscillating_propensity


subroutine oscillating_print_model(self)
	class(oscillating_t), intent(in) :: self

	print *, "--------------------------------"
	print *, "Model: oscillating"
	print *, "Species:   ", self%n_species
	print *, "Reactions: ", self%n_reactions
	print *, "Stoichiometry:"
	print *, self%stoich
	print *, "Parameters:"
	print *, "  lambda = ", self%params%lambda
	print *, "  beta   = ", self%params%beta
	print *, "  K   = ", self%params%K
	print *, "  n   = ", self%params%n
	print *, "--------------------------------"

end subroutine oscillating_print_model

end module
