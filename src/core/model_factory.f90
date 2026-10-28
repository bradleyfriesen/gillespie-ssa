module model_factory
use model_interface
use birth_death_model
use brusselator_model
use gene_expression_model
use homodimer_model
use oscillating_model
implicit none


contains


subroutine list_models()
	print *, "Available models:"
	print *, "birth-death"
	print *, "brusselator"
	print *, "homodimer"
	print *, "gene-expression"
	print *, "oscillating"
end subroutine list_models


subroutine create_model(name, parameter_file, model)
	character(len=*), intent(in) :: name
	character(len=*), intent(in) :: parameter_file
	class(model_t), allocatable, intent(out) :: model

	select case (name)

	case("birth-death")
		allocate(birth_death_t :: model)
		select type(model)
		type is (birth_death_t)
			call initialize_birth_death(model, parameter_file)
		end select
		
	case("brusselator")
		allocate(brusselator_t :: model)
		select type(model)
		type is (brusselator_t)
			call initialize_brusselator(model, parameter_file)
		end select
		
	case("homodimer")
		allocate(homodimer_t :: model)
		select type(model)
		type is (homodimer_t)
			call initialize_homodimer(model, parameter_file)
		end select
		
	case("gene-expression")
		allocate(gene_expression_t :: model)
		select type(model)
		type is (gene_expression_t)
			call initialize_gene_expression(model, parameter_file)
		end select
		
	case("oscillating")
		allocate(oscillating_t :: model)
		select type(model)
		type is (oscillating_t)
			call initialize_oscillating(model, parameter_file)
		end select

	case default
		print *, "Unknown model: ", trim(name)
		call list_models()
		error stop

	end select

end subroutine create_model


end module
