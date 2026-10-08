module birth_death_model
use kind_parameters
use model_interface
implicit none

integer, parameter :: n_species = 1
integer, parameter :: n_reactions = 2

type :: params_t
	real(dp) :: lambda
	real(dp) :: beta
end type

type, extends(model_t) :: birth_death_t
	type(params_t) :: params
contains
	procedure :: propensity => birth_death_propensity
	procedure :: print_model => birth_death_print_model
end type


contains


subroutine initialize_birth_death(model, parameter_file)
	type(birth_death_t), intent(out) :: model
	character(len=*), intent(in) :: parameter_file
	real(dp) :: lambda, beta

	namelist /parameters/ lambda, beta

	model%n_species = n_species
	model%n_reactions = n_reactions

	allocate (model%stoich(model%n_species, model%n_reactions))

	model%stoich(:, 1) = [1]
	model%stoich(:, 2) = [-1]

	open (unit=10, file=parameter_file, status="old", action="read")
	read (10, nml=parameters)
	close (10)
	
	model%params%lambda = lambda
	model%params%beta = beta
end subroutine initialize_birth_death


subroutine birth_death_propensity(self, state, r)
	class(birth_death_t), intent(in) :: self
	integer, intent(in) :: state(:)
	real(dp), intent(out) :: r(:)

	r(1) = self%params%lambda
	r(2) = real(state(1), dp)*self%params%beta

end subroutine birth_death_propensity


subroutine birth_death_print_model(self)
	class(birth_death_t), intent(in) :: self

	print *, "--------------------------------"
	print *, "Model: birth-death"
	print *, "Species:   ", self%n_species
	print *, "Reactions: ", self%n_reactions

	print *, "Stoichiometry:"
	print *, self%stoich

	print *, "Parameters:"
	print *, "  lambda = ", self%params%lambda
	print *, "  beta   = ", self%params%beta

	print *, "--------------------------------"
end subroutine birth_death_print_model


end module
