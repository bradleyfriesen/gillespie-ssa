module command_line
use model_factory
implicit none


contains


subroutine get_model_and_parameters(model_name, parameter_file)
	character(len=:), allocatable, intent(out) :: model_name, parameter_file
	integer :: length
	
	if (command_argument_count() /= 2) then
		print *, "Usage:"
		print *, "./a.out MODEL PARAMSET"
		call list_models()
		error stop
	end if

	call get_command_argument(1, length=length)
	allocate(character(len=length) :: model_name)
	call get_command_argument(1, model_name)
	
	call get_command_argument(2, length=length)
	allocate(character(len=length) :: parameter_file)
	call get_command_argument(2, parameter_file)

end subroutine get_model_and_parameters

end module command_line
