module gene_expression_model
use kind_parameters
use model_interface
implicit none

integer, parameter :: n_species = 2
integer, parameter :: n_reactions = 4

type :: params_t
	real(dp) :: lambda(2)
	real(dp) :: beta(2)
end type

type, extends(model_t) :: gene_expression_t
	type(params_t) :: params
contains
	procedure :: propensity => gene_expression_propensity
	procedure :: print_model => gene_expression_print_model
end type


contains


subroutine initialize_gene_expression(model, parameter_file)
	type(gene_expression_t), intent(out) :: model
	character(len=*), intent(in) :: parameter_file
	real(dp) :: lambda(2), beta(2)

	namelist /parameters/ lambda, beta

	model%n_species = n_species
	model%n_reactions = n_reactions

	allocate (model%stoich(model%n_species, model%n_reactions))

	model%stoich = reshape([ &
		1, 0, &
		-1, 0, &
		0, 1, &
		0, -1 &
		], shape(model%stoich))

	open (unit=10, file=parameter_file, status="old", action="read")
	read (10, nml=parameters)
	close (10)
	
	model%params%lambda = lambda
	model%params%beta = beta

end subroutine initialize_gene_expression


subroutine gene_expression_propensity(self, state, r)
	class(gene_expression_t), intent(in) :: self
	integer, intent(in) :: state(:)
	real(dp), intent(out) :: r(:)

	associate(x => real(state(:), dp), &
		lambda => self%params%lambda, &
		beta => self%params%beta)
		
		r(1) = lambda(1)
		r(2) = x(1) * beta(1)
		r(3) = lambda(2) * x(1)
		r(4) = x(2) * beta(2)
	end associate

end subroutine gene_expression_propensity


subroutine gene_expression_print_model(self)
	class(gene_expression_t), intent(in) :: self

	print *, "--------------------------------"
	print *, "Model: gene_expression"
	print *, "Species:   ", self%n_species
	print *, "Reactions: ", self%n_reactions

	print *, "Stoichiometry:"
	print *, self%stoich

	print *, "Parameters:"
	print *, "  lambda = ", self%params%lambda
	print *, "  beta   = ", self%params%beta

	print *, "--------------------------------"

end subroutine gene_expression_print_model

end module
