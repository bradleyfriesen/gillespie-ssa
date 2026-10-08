module model_factory
use model_interface
use birth_death_model
use homodimer_model
implicit none


contains


subroutine list_models()
	print *, "Available models:"
	print *, "birth-death"
	print *, "homodimer"
end subroutine list_models


subroutine create_model(name, parameter_file, model)
	character(len=*), intent(in) :: name
	character(len=*), intent(in) :: parameter_file
	class(model_t), allocatable, intent(out) :: model

	select case (name)

	case ("birth-death")
		allocate (birth_death_t :: model)
		select type (model)
		type is (birth_death_t)
			call initialize_birth_death(model, parameter_file)
		end select
		
	case ("homodimer")
		allocate (homodimer_t :: model)
		select type (model)
		type is (homodimer_t)
			call initialize_homodimer(model, parameter_file)
		end select

	case default
		print *, "Unknown model: ", trim(name)
		call list_models()
		error stop

	end select

end subroutine create_model


end module
