module model_interface

	use kind_parameters
	implicit none

	type, abstract :: model_t
		integer :: n_species
		integer :: n_reactions
		integer, allocatable :: stoich(:, :)
	contains
		procedure(propensity_interface), deferred :: propensity
		procedure(print_model_interface), deferred :: print_model
	end type model_t

	abstract interface
	
		subroutine propensity_interface(self, state, r)
			import :: model_t, dp
			
			class(model_t), intent(in) :: self
			integer, intent(in) :: state(:)
			real(dp), intent(out) :: r(:)
		end subroutine
		
		subroutine print_model_interface(self)
			import :: model_t
			class(model_t), intent(in) :: self
		end subroutine
		
	end interface

end module
