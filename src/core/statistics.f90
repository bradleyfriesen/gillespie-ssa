module statistics
use kind_parameters
implicit none

type :: moments_t
	real(dp) :: total_time
	real(dp), allocatable :: mean(:)
	real(dp), allocatable :: m2(:,:)
contains
	procedure :: initialize => initialize_moments
	procedure :: update => update_moments
	procedure :: covariance => moments_covariance
end type moments_t


contains


subroutine initialize_moments(self, n_species)
	class(moments_t), intent(out) :: self
	integer, intent(in) :: n_species
	
	allocate(self%mean(n_species))
	allocate(self%m2(n_species, n_species))
	
	self%total_time = 0._dp
	self%mean = 0._dp
	self%m2 = 0._dp
	
end subroutine initialize_moments


subroutine update_moments(self, state, dt)
	class(moments_t), intent(inout) :: self
	integer, intent(in) :: state(:)
	real(dp), intent(in) :: dt

	real(dp) :: new_total_time
	real(dp), dimension(size(state)) :: delta, delta_new
	integer :: i, j

	new_total_time = self%total_time + dt
	
	delta = real(state, dp) - self%mean
	self%mean = self%mean + (dt / new_total_time) * delta
	delta_new = real(state, dp) - self%mean
	
	do i = 1, size(state)
		do j = 1, size(state)
			self%m2(i,j) = self%m2(i,j) + dt * delta(i) * delta_new(j)
		end do
	end do
	
	self%total_time = new_total_time
	
end subroutine update_moments


function moments_covariance(self) result(cov)
	class(moments_t), intent(in) :: self
	real(dp) :: cov(size(self%mean), size(self%mean))
	
	cov = self%m2 / self%total_time
end function moments_covariance


end module statistics
