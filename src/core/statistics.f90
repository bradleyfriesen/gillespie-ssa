module statistics
use kind_parameters
use utilities
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

type :: flux_statistics_t
	real(dp) :: total_time
	real(dp), allocatable :: mean_state(:)
	real(dp), allocatable :: mean_flux(:)
	real(dp), allocatable :: mean_creation(:)
	real(dp), allocatable :: mean_degradation(:)
	real(dp), allocatable :: m2_flux_state(:, :)
contains
	procedure :: initialize => initialize_flux_statistics
	procedure :: update => update_flux_statistics
	procedure :: covariance => flux_statistics_covariance
end type flux_statistics_t


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


subroutine initialize_flux_statistics(self, n_species)
	class(flux_statistics_t), intent(out) :: self
	integer, intent(in) :: n_species
	
	allocate(self%mean_state(n_species))
	allocate(self%mean_flux(n_species))
	allocate(self%mean_creation(n_species))
	allocate(self%mean_degradation(n_species))
	allocate(self%m2_flux_state(n_species, n_species))
	
	self%total_time = 0._dp
	self%mean_state = 0._dp
	self%mean_flux = 0._dp
	self%mean_creation = 0._dp
	self%mean_degradation = 0._dp
	self%m2_flux_state = 0._dp
	
end subroutine initialize_flux_statistics


subroutine update_flux_statistics(self, state, r, stoich, dt)
	class(flux_statistics_t), intent(inout) :: self
	integer, intent(in) :: state(:)
	real(dp), intent(in) :: r(:)
	integer, intent(in) :: stoich(:,:)
	real(dp), intent(in) :: dt
	
	real(dp) :: new_total_time
	real(dp), dimension(size(state)) :: delta_flux, delta_state, &
		delta_flux_new, delta_state_new, &
		creation, degradation, flux
	integer :: i, j
	
	new_total_time = self%total_time + dt
	
	! Calculate creation and degradation fluxes
	! Note the signs. Creation and Degradation are defined positive.
	creation = 0._dp
	degradation = 0._dp
	do i = 1, size(state)
		do j = 1, size(r)
			if (stoich(i,j) > 0) then
				creation(i) = creation(i) + stoich(i,j) * r(j)
			else if (stoich(i,j) < 0) then
				degradation(i) = degradation(i) - stoich(i,j) * r(j)
			end if
		end do
	end do
	flux = creation - degradation
	
	self%mean_creation = self%mean_creation + &
		(dt / new_total_time) * (creation - self%mean_creation)

	self%mean_degradation = self%mean_degradation + &
		(dt / new_total_time) * (degradation - self%mean_degradation)
	
	! Welford's online algorithm for covariances (update means first)
	delta_flux = flux - self%mean_flux
	delta_state = real(state, dp) - self%mean_state
	
	self%mean_state = self%mean_state + &
		(dt / new_total_time) * delta_state
	
	self%mean_flux = self%mean_flux + &
		(dt / new_total_time) * delta_flux
		
	delta_flux_new  = flux - self%mean_flux
	delta_state_new = real(state, dp) - self%mean_state
		
	self%total_time = new_total_time
	
	do i = 1, size(state)
		do j = 1, size(state)
			self%m2_flux_state(i,j) = self%m2_flux_state(i,j) + &
				dt * (-delta_flux(i)) * delta_state_new(j)
		end do
	end do
	
end subroutine update_flux_statistics


function flux_statistics_covariance(self) result(cov)
	class(flux_statistics_t), intent(in) :: self
	real(dp) :: cov(size(self%mean_flux), size(self%mean_flux))
	
	cov = self%m2_flux_state / self%total_time
end function flux_statistics_covariance


subroutine check_covariance_balance(flux, n_species, n_reactions, stoich, n_events)
	class(flux_statistics_t), intent(in) :: flux
	integer, intent(in) :: n_species, n_reactions, n_events(:), stoich(:,:)
	
	real(dp), parameter :: eps = tiny(eps)
	real(dp), dimension(n_species) :: flux_pm, tau
	real(dp), dimension(n_species, n_species) :: U, D, S, cov
	real(dp), dimension(n_reactions) :: r_avg
	real(dp) :: rel_change
	integer :: i, j, k
	
	flux_pm = (flux%mean_creation + flux%mean_degradation) / 2._dp
	cov = flux%covariance()
	r_avg = n_events / flux%total_time
	tau = 1._dp * flux%mean_state / flux%mean_degradation
	
	S = 0._dp
	! Calculate average step sizes
	do i = 1, n_species
	do j = 1, n_species
	do k = 1, n_reactions
		S(i, j) = S(i, j) &
			+ r_avg(k)*abs(stoich(i,k)) / sum(r_avg(:)*abs(stoich(i,:))) &
			* abs(stoich(j, k)) * sign(1._dp, real(stoich(i, k) * stoich(j, k), dp))
	end do
	end do
	end do
	
	! Calculate correlation/"drift" and diffusion matrix
	do i = 1, n_species
	do j = 1, n_species
		U(i, j) = cov(i, j) / flux_pm(i) / flux%mean_state(j) / tau(i)
		D(i, j) = S(i, j) / flux%mean_state(j) / tau(i) &
			+ S(j, i) / flux%mean_state(i) / tau(j)
	end do
	end do

	print *, "r_avg: ", r_avg
	print *, "tau: ", tau
	print *, "S: ", S
	print *, "U: ", U
	print *, "D: ", D
	
	! Output results
	write(*,*) "relative change between U(i,j)+U(j,i) and D(i,j), or U(i,j) and -U(i,j) if D(i,j)==0:"
	write(*,'(*(A20))') "D(ij)==0","i", "j", "Relative change", "U(ij)+U(ji)", "D(ij)"
	do i = 1, n_species
	do j = 1, n_species
		if (abs(D(i,j)) > eps) then
			rel_change = relative_change(U(i, j) + U(j, i), D(i, j))
			write(*,'(A20,2(I20),3(G20.12))') "FALSE", i, j, rel_change, U(i, j) + U(j, i), D(i, j)
		else
			rel_change = relative_change(U(i, j), -U(j, i))
			write(*,'(A20,2(I20),*(G20.12))') "TRUE", i, j, rel_change, U(i, j), -U(j, i)
		end if
	end do
	end do
	
end subroutine check_covariance_balance

end module statistics
