module brusselator_model
use kind_parameters
use model_interface
implicit none

integer, parameter :: n_species = 2
integer, parameter :: n_reactions = 4

type :: params_t
	real(dp) :: omega
	real(dp) :: b
	real(dp) :: a
end type

type, extends(model_t) :: brusselator_t
	type(params_t) :: params
contains
	procedure :: propensity => brusselator_propensity
	procedure :: print_model => brusselator_print_model
end type


contains


subroutine initialize_brusselator(model, parameter_file)
	type(brusselator_t), intent(out) :: model
	character(len=*), intent(in) :: parameter_file
	real(dp) :: omega, b, a

	namelist /parameters/ omega, b, a

	model%n_species = n_species
	model%n_reactions = n_reactions

	allocate (model%stoich(model%n_species, model%n_reactions))

	model%stoich = reshape([ &
		1, 0, &
		-1, 1, &
		1, -1, &
		-1, 0 & 
		], shape(model%stoich))

	open (unit=10, file=parameter_file, status="old", action="read")
	read (10, nml=parameters)
	close (10)
	
	model%params%omega = omega
	model%params%b = b
	model%params%a = a

end subroutine initialize_brusselator


subroutine brusselator_propensity(self, state, r)
	class(brusselator_t), intent(in) :: self
	integer, intent(in) :: state(:)
	real(dp), intent(out) :: r(:)

	associate(x => real(state(:), dp), &
		omega => self%params%omega, &
		b => self%params%b, &
		a => self%params%a)
		
		r(1) = omega
		r(2) = b * x(1)
		r(3) = a / omega**2 * x(1) * (x(1) - 1) * x(2)
		r(4) = x(1)
	end associate

end subroutine brusselator_propensity


subroutine brusselator_print_model(self)
	class(brusselator_t), intent(in) :: self

	print *, "--------------------------------"
	print *, "Model: brusselator"
	print *, "Species:   ", self%n_species
	print *, "Reactions: ", self%n_reactions
	print *, "Stoichiometry:"
	print *, self%stoich
	print *, "Parameters:"
	print *, "  omega = ", self%params%omega
	print *, "  b   = ", self%params%b
	print *, "  a   = ", self%params%a
	print *, "--------------------------------"

end subroutine brusselator_print_model

end module
