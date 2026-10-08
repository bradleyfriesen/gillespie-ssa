module homodimer_model
use kind_parameters
use model_interface
implicit none

integer, parameter :: n_species = 1
integer, parameter :: n_reactions = 2

type :: params_t
	real(dp) :: lambda
	real(dp) :: beta
end type

type, extends(model_t) :: homodimer_t
	type(params_t) :: params
contains
	procedure :: propensity => homodimer_propensity
	procedure :: print_model => homodimer_print_model
end type


contains


subroutine initialize_homodimer(model, parameter_file)
	type(homodimer_t), intent(out) :: model
	character(len=*), intent(in) :: parameter_file
	real(dp) :: lambda, beta

	namelist /parameters/ lambda, beta

	model%n_species = 1
	model%n_reactions = 2

	allocate (model%stoich(model%n_species, model%n_reactions))

	model%stoich(:, 1) = [1]
	model%stoich(:, 2) = [-2]

	open (unit=10, file=parameter_file, status="old", action="read")
	read (10, nml=parameters)
	close (10)
	
	model%params%lambda = lambda
	model%params%beta = beta

end subroutine initialize_homodimer


subroutine homodimer_propensity(self, state, r)
	class(homodimer_t), intent(in) :: self
	integer, intent(in) :: state(:)
	real(dp), intent(out) :: r(:)

	associate(x => real(state(1), dp), &
		lambda => self%params%lambda, &
		beta => self%params%beta)
		r(1) = lambda
		r(2) = x * (x-1._dp) * beta
	end associate

end subroutine homodimer_propensity


subroutine homodimer_print_model(self)
	class(homodimer_t), intent(in) :: self

	print *, "--------------------------------"
	print *, "Model: homodimer"
	print *, "Species:   ", self%n_species
	print *, "Reactions: ", self%n_reactions

	print *, "Stoichiometry:"
	print *, self%stoich

	print *, "Parameters:"
	print *, "  lambda = ", self%params%lambda
	print *, "  beta   = ", self%params%beta

	print *, "--------------------------------"

end subroutine homodimer_print_model


end module
