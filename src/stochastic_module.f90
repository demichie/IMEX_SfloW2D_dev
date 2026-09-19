!********************************************************************************
!> \brief Stochastic module
!
!> This moduel contains all the procedures for the stochastic variable
!
!> \date 11/06/2025
!> @author 
!> Zeno Geddo
!
!********************************************************************************
MODULE stochastic_module
  
  ! external variables
  USE parameters_2d, ONLY : wp, rheology_model, idx_stoch, idx_u, idx_v
  
  USE domain_2d, ONLY: domain_type
  USE state_2d, ONLY: state_type
  USE parameters_2d, ONLY : output_stoch_vars_flag, length_spatial_corr,        &
        stochastic_flag, stoch_transport_flag
  USE geometry_2d, ONLY : cell_size, comp_cells_x, comp_cells_y
  USE stochastic_random_2d, ONLY : initialize_stochastic_rng, gaussian_noise
     
  ! variables related to OU process
  
  !LOGICAL :: sym_noise_flag ! use symmetric noise or not

  !> Noise parameter:\n
  !> - 0    => the noise will be simmetric
  !> 1      => abs val is taken
  !> -1     => -abs val is taken
  !> .
  REAL(wp) :: sym_noise

  REAL(wp) :: tau_stochastic

  REAL(wp) :: std_max, std_min

  REAL(wp) :: std_slope_factor ! Fr_0_stochastic

  !> Velocity scale u_0 in the stochastic intensity law.
  REAL(wp) :: noise_activation_velocity
  
  REAL(wp) :: noise_pow_val !(|Z|^power)
  
  ! variables related to statistics
  REAL(wp) :: Z_min, Z_max, Z_mean, Z_std
  REAL(wp) :: percentiles(9) ! 5,10,20,30,50,70,80,90,95

  TYPE :: stochastic_workspace_type

     !> Stochastic field at cell centers
     REAL(wp), ALLOCATABLE :: Z(:,:)

     !> Possibly transformed fluctuation passed to the friction law. Z remains
     !> the Gaussian OU state transported through the conservative variable hZ.
     REAL(wp), ALLOCATABLE :: effective_Z(:,:)

     !> Spatial-correlation convolution kernel
     REAL(wp), ALLOCATABLE :: conv_kernel(:,:)

   CONTAINS

     PROCEDURE :: initialize => initialize_stochastic_workspace
     PROCEDURE :: finalize => finalize_stochastic_workspace
     PROCEDURE :: initialize_steady => getSteadyStateZ
     PROCEDURE :: generate_kernel => genConvolutionKernel
     PROCEDURE :: update => update_stochastic_variable
     PROCEDURE :: prepare_timestep => prepare_stochastic_timestep
     PROCEDURE :: refresh_effective => refresh_effective_stochastic_field

  END TYPE stochastic_workspace_type

CONTAINS

  SUBROUTINE initialize_stochastic_workspace(this)

    CLASS(stochastic_workspace_type), INTENT(INOUT) :: this

    IF ( ALLOCATED(this%Z) ) DEALLOCATE(this%Z)
    ALLOCATE(this%Z(comp_cells_x,comp_cells_y))
    this%Z = 0.0_wp

    IF ( ALLOCATED(this%effective_Z) ) DEALLOCATE(this%effective_Z)
    ALLOCATE(this%effective_Z(comp_cells_x,comp_cells_y))
    this%effective_Z = 0.0_wp

    IF (stochastic_flag) CALL initialize_stochastic_rng

  END SUBROUTINE initialize_stochastic_workspace

  SUBROUTINE finalize_stochastic_workspace(this)

    CLASS(stochastic_workspace_type), INTENT(INOUT) :: this

    IF ( ALLOCATED(this%Z) ) DEALLOCATE(this%Z)
    IF ( ALLOCATED(this%effective_Z) ) DEALLOCATE(this%effective_Z)
    IF ( ALLOCATED(this%conv_kernel) ) DEALLOCATE(this%conv_kernel)

  END SUBROUTINE finalize_stochastic_workspace
  
  REAL(wp) FUNCTION MeanFieldCorrection(g,h,Fr) ! working only in the case of mu(Fr)
    
    !> Compute the mean field correction 
    REAL(wp), INTENT(IN) :: g,h,Fr
    REAL(wp) :: factor, sFR, dsFR_dFr, spatial_corr
    
    ! Compute a normalization factor
    factor = (g**2._wp) / SQRT(g*h)
    
    ! Compute the intensity of the noise (std of the process)  
    sFR = std_max + ( std_min - std_max) * EXP(- Fr / std_slope_factor )
    
    ! Compute the derivative of the intensity of the noise
    dsFR_dFr = -(std_min - std_max) / std_slope_factor *                        &
         EXP(-Fr / std_slope_factor)
    
    ! Compute factor for spatial correlation
    spatial_corr = 2._wp * tau_stochastic ! ... will introduce correlation 
    
    ! Compute mean field correction
    MeanFieldCorrection = - factor * sFR * dsFR_dFr * spatial_corr

  END FUNCTION MeanFieldCorrection  


  SUBROUTINE getSteadyStateZ(this, state, domain)
  ! should find a better way to understand when the process is stable.
  ! Since sigma is varing in space and time, the process may never be stable.
    IMPLICIT none
    CLASS(stochastic_workspace_type), INTENT(INOUT) :: this
    CLASS(state_type), INTENT(INOUT) :: state
    CLASS(domain_type), INTENT(IN) :: domain
    INTEGER :: j, k , l 
    INTEGER :: i, n_iter
    REAL(wp) :: dt
    REAL(wp) :: t_steady
    
    ! Initialize all stochastic process to zero
    this%Z = 0.0_wp

    ! Define how many iteration to do for the burn in
    t_steady = 10_wp * tau_stochastic 
    dt = 0.1_wp * tau_stochastic
    n_iter = CEILING(t_steady/dt)
    ! WRITE(*,*) 'dt,n_inter',dt,n_iter

    ! Allocate and compute convolution kernel if needed (will be keept in memory)
    IF (length_spatial_corr > 0.0_wp) THEN
        CALL this%generate_kernel()
    END IF
    
    ! Compute n iterations to get to a steady state OU process
    DO i=1,n_iter
       CALL this%update(state, dt)
    END DO

    IF (stoch_transport_flag) THEN

       !$OMP PARALLEL DO private(j,k)

       DO l = 1,domain%solve_cells

          j = domain%j_cent(l)
          k = domain%k_cent(l)
          
          state%qp(idx_stoch,j,k) = this%Z(j,k)
          state%q(idx_stoch,j,k) = state%q(1,j,k) *                         &
               this%Z(j,k)
          
       END DO
       
       !$OMP END PARALLEL DO

       WRITE(*,*) 'end stochastic initialization'
       
    END IF

    RETURN
    
  END SUBROUTINE getSteadyStateZ 


  SUBROUTINE genConvolutionKernel(this)
      IMPLICIT none
      CLASS(stochastic_workspace_type), INTENT(INOUT) :: this
      INTEGER :: half_width, n_nodes_per_dim, x_index, y_index
      REAL(wp) :: bandwidth, x, y

      ! The thesis defines length_spatial_corr as l and the kernel bandwidth
      ! as l/2 (Eq. D.32).
      bandwidth = 0.5_wp * length_spatial_corr
      ! Four bandwidths on either side retain all but a negligible Gaussian tail.
      half_width = MAX(1, CEILING(4.0_wp * bandwidth / cell_size))
      n_nodes_per_dim = 2 * half_width + 1

      ! Allocate the output array
      IF (ALLOCATED(this%conv_kernel)) DEALLOCATE(this%conv_kernel)
      ALLOCATE(this%conv_kernel(n_nodes_per_dim,                             &
           n_nodes_per_dim))

      ! Compute 2D Gaussian values (at centers of cells)
      DO y_index = 1, n_nodes_per_dim 
        DO x_index = 1, n_nodes_per_dim
          x = cell_size * REAL(x_index - half_width - 1, wp)
          y = cell_size * REAL(y_index - half_width - 1, wp)
          this%conv_kernel(x_index, y_index) = evalGaussian2d(x, y)
        END DO
      END DO

      ! Renormalize values of the kernel
      this%conv_kernel = this%conv_kernel / SUM(this%conv_kernel)

  END SUBROUTINE genConvolutionKernel

  REAL(wp) FUNCTION evalGaussian2d(x, y)
      ! Sample the centered isotropic Gaussian kernel Q^(1/2).
      IMPLICIT none
      REAL(wp), INTENT(IN) :: x, y
      REAL(wp) :: s, PI

      s = 0.5_wp * length_spatial_corr
      PI = ACOS(-1.0_wp)
      evalGaussian2d = EXP(-0.5_wp * (((x/s)**2) + ((y/s)**2))) /             &
                        (2.0_wp * PI * s**2)

  END FUNCTION evalGaussian2d


  SUBROUTINE update_stochastic_variable(this, state, dt)
    ! UPDATE THE SOLUTION OF THE ORNSTEIN-UHLENBACK PROCESS USING EULER-MARUYAMA METHOD
    ! NOISE CAN BE TRANSPORTED, SOURCE TERM ADDED IN eval_mass_exchange_terms (IMPORTANT)
    IMPLICIT NONE
    CLASS(stochastic_workspace_type), INTENT(INOUT) :: this
    CLASS(state_type), INTENT(IN) :: state
    INTEGER :: j,k
    INTEGER :: noise_size
    REAL(wp), INTENT(IN) :: dt
    REAL(wp) :: sigma_noise
    REAL(wp):: noise(comp_cells_x, comp_cells_y)
    REAL(wp):: conv_result(comp_cells_x, comp_cells_y)

    ! Generate standard gaussian noise over entire domain-> N(0,1)
    noise_size = comp_cells_x*comp_cells_y
    noise = RESHAPE(gaussian_noise(noise_size),                              &
         [comp_cells_x, comp_cells_y])

    ! Convolve the gaussian noise to introduce spatial correlation if needed
    IF (length_spatial_corr > 0.0_wp) THEN
        ! (for convenience, the conv_kernel is generated only once before the burn in)
        CALL convolve_2d(noise, this%conv_kernel, conv_result)
        noise = conv_result
    END IF

    ! Loop over the entire grid to update stochastic process
    !$OMP PARALLEL
    !$OMP DO private(j,k,sigma_noise)
    DO k = 1,comp_cells_y
       DO j = 1,comp_cells_x  
          ! Update the Ornstein-Uhlenback process
          sigma_noise = getSigmaNoise(state,j,k)
          this%Z(j,k) = EulerMaruyamaScheme(                                &
               this%Z(j,k), dt, sigma_noise, noise(j,k))
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    ! Compute statistic in space at given time if needed as outputs
    IF (output_stoch_vars_flag) THEN
        CALL OUBasicStatsInSpaceAtGivenTime(this%Z, Z_min,                   &
             Z_mean, Z_max, Z_std)
        percentiles(:) = 0.0_wp ! Percentiles set to 0 to avoid the computations
        ! The percentiles shold not be computed at each iteration otherwise is too slow
        !CALL percentilesArrayAtGivenTime(Z, percentiles) ! it is very slow!!!!!!!!!!
    END IF
    
    CALL this%refresh_effective
   
  END SUBROUTINE update_stochastic_variable


  SUBROUTINE prepare_stochastic_timestep(this, state, dt)
    !> Apply the OU fractional step and initialize hZ for conservative transport.
    CLASS(stochastic_workspace_type), INTENT(INOUT) :: this
    CLASS(state_type), INTENT(INOUT) :: state
    REAL(wp), INTENT(IN) :: dt

    IF (stoch_transport_flag) THEN
       ! In wet cells, the OU process starts from the value transported during
       ! the previous timestep. Dry cells retain their independently evolving Z.
       WHERE (state%q(1,:,:) > EPSILON(1.0_wp))
          this%Z = state%q(idx_stoch,:,:) / state%q(1,:,:)
       END WHERE
    END IF

    CALL this%update(state, dt)

    IF (stoch_transport_flag) THEN
       WHERE (state%q(1,:,:) > EPSILON(1.0_wp))
          state%q(idx_stoch,:,:) = state%q(1,:,:) * this%Z
          state%qp(idx_stoch,:,:) = this%Z
       ELSEWHERE
          state%q(idx_stoch,:,:) = 0.0_wp
          state%qp(idx_stoch,:,:) = 0.0_wp
       END WHERE
    END IF

  END SUBROUTINE prepare_stochastic_timestep


  SUBROUTINE refresh_effective_stochastic_field(this)
    !> Derive the fluctuation used by friction without modifying the OU state.
    CLASS(stochastic_workspace_type), INTENT(INOUT) :: this

    IF (sym_noise > 0.0_wp) THEN
       this%effective_Z = ABS(this%Z)**noise_pow_val
    ELSEIF (sym_noise < 0.0_wp) THEN
       this%effective_Z = -(ABS(this%Z)**noise_pow_val)
    ELSE
       this%effective_Z = this%Z
    END IF

  END SUBROUTINE refresh_effective_stochastic_field


 REAL(wp) FUNCTION VelocitySquared(state,j,k)
  !> Compute |u|^2 at the cell center.
  IMPLICIT NONE
  CLASS(state_type), INTENT(IN) :: state
  INTEGER, INTENT(IN) :: j, k

  IF (state%q(1,j,k) > EPSILON(1.0_wp)) THEN
     ! qp already stores the physical velocities associated with q. Reusing
     ! them avoids repeating the full thermodynamic conversion in every cell.
     VelocitySquared = state%qp(idx_u,j,k)**2 + state%qp(idx_v,j,k)**2
  ELSE
     VelocitySquared = 0.0_wp
  END IF

END FUNCTION VelocitySquared


  REAL(wp) FUNCTION getSigmaNoise(state,j,k)
  ! Compute the intensity of the noise depending of the friction used
    IMPLICIT NONE
    CLASS(state_type), INTENT(IN) :: state
    INTEGER, INTENT(IN) :: j, k
    REAL(wp) :: velocity_squared

    IF (state%q(1,j,k) <= EPSILON(1.0_wp)) THEN
      ! Keep fluctuations available when flow enters a previously dry cell.
      getSigmaNoise = std_max
    ELSEIF ((rheology_model == 9) .OR. (rheology_model == 10)) THEN
      velocity_squared = VelocitySquared(state,j,k)
      getSigmaNoise = expFormNoise(velocity_squared)
    ELSE
      getSigmaNoise = std_max
    END IF  
    RETURN
  END FUNCTION getSigmaNoise

  REAL(wp) FUNCTION expFormNoise(val)
    ! Compute the bounded intensity of the noise of the stochastic process
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: val
    expFormNoise = std_max + ( std_min - std_max) * exp(- val / std_slope_factor )
  END FUNCTION expFormNoise

  REAL(wp) FUNCTION EulerMaruyamaScheme(Zij, dt, sigma, noise_ij)
    ! Update the Ornstein-Uhlenback process using EULER-MARUYAMA Method
    ! Z(j,k) = Z(j,k) - (dt / tau_stochastic) * Z(j,k) + &
    ! sigma_noise * SQRT(2.0_wp * (dt / tau_stochastic)) * noise(j,k)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: Zij, dt, sigma, noise_ij
    EulerMaruyamaScheme = Zij - (dt / tau_stochastic) * Zij +             &
            sigma * SQRT(2.0_wp * (dt / tau_stochastic)) * noise_ij
  END FUNCTION

subroutine convolve_2d(input_signal, kernel, result)
  !> Zero-padded two-dimensional convolution.
  !> All arrays follow the solver convention (x,y).
  implicit none
  real(wp), dimension(:,:), intent(in) :: input_signal, kernel
  real(wp), dimension(:,:), intent(out) :: result
  integer :: len_kernel_x, len_kernel_y, len_inp_sig_x,                 &
  len_inp_sig_y
  integer :: half_len_k_x, half_len_k_y
  integer :: left_padding, right_padding, south_padding,                &
  north_padding
  integer :: inx_x, inx_y, idx_centre_x, idx_centre_y,                  &
  idx_k_x, idx_k_y, shift_x, shift_y
  real(wp), dimension(:,:), allocatable :: padded_signal
  real(wp) :: sum_at_position
  
  ! Get shape inp signal (should be moved outside)
  len_inp_sig_x = size(input_signal, 1)
  len_inp_sig_y = size(input_signal, 2)
  ! Get shape kernel (should be moved outside)
  len_kernel_x = size(kernel, 1)
  len_kernel_y = size(kernel, 2)
  half_len_k_x = (len_kernel_x - 1) / 2
  half_len_k_y = (len_kernel_y - 1) / 2

  ! Calculate the required padding for each dimension
  ! If padding is odd, one more zero is added to the right side.
  left_padding = half_len_k_x
  right_padding = (len_kernel_x - 1) - left_padding
  south_padding = half_len_k_y
  north_padding = (len_kernel_y - 1) - south_padding

  ! Allocate array for padded inp signal
  allocate(padded_signal(len_inp_sig_x + left_padding + right_padding, &
   len_inp_sig_y + south_padding + north_padding))
  
  ! Pad the input signal with zeros
  padded_signal = 0._wp
  padded_signal(left_padding + 1:left_padding + len_inp_sig_x,          &
  south_padding + 1:south_padding + len_inp_sig_y) = input_signal
  
  ! Initialize to zero the result
  result = 0._wp

  ! Loop over all inp signal cells
  do inx_x = 1, len_inp_sig_x
    do inx_y = 1, len_inp_sig_y
      idx_centre_x = left_padding + inx_x
      idx_centre_y = south_padding + inx_y
      sum_at_position = 0._wp

      ! Loop to convolve values around the current cell
      do idx_k_x = 1, len_kernel_x
        do idx_k_y = 1, len_kernel_y
          ! Convolution offset: for kernel index k the flipped shift is
          ! -(k - 1 - half_len_k) = -k + half_len_k + 1. The previous +2
          ! overran the padded array by one, which is what tripped the
          ! bound checks below on an identity kernel.
          shift_x = -idx_k_x + half_len_k_x + 1
          shift_y = -idx_k_y + half_len_k_y + 1
          ! sum_at_position = 0.1_wp

          ! here there may be a bug in the indexing (segmentation fault-invalid memory reference) 
          if ((idx_centre_y + shift_y .LT. 1) .or.  (idx_centre_y + shift_y .GT. size(padded_signal,2))) THEN
            
            print*,"problem y idx"
            print*,shape(kernel)
            print*,shape(input_signal)
            print*,south_padding, north_padding
            print*,left_padding, right_padding
            print*,shape(padded_signal)
            print*, inx_x, inx_y, idx_k_x, idx_k_y
            print*,idx_centre_y + shift_y
            STOP
          END if
          if ((idx_centre_x + shift_x .LT. 1) .or.  (idx_centre_x + shift_x .GT. size(padded_signal,1))) THEN
            print*,"problem x idx"
            print*,shape(kernel)
            print*,shape(input_signal)
            print*,south_padding, north_padding
            print*,left_padding, right_padding
            print*,shape(padded_signal)
            print*, inx_x, inx_y, idx_k_x, idx_k_y
            print*,idx_centre_x + shift_x
            STOP  
          END if  

          sum_at_position = sum_at_position +                          &
          padded_signal(idx_centre_x + shift_x, idx_centre_y + shift_y)&
          * kernel(idx_k_x, idx_k_y)
        end do
      end do

      result(inx_x, inx_y) = sum_at_position
    end do
  end do

  deallocate(padded_signal)
end subroutine convolve_2d


subroutine OUBasicStatsInSpaceAtGivenTime(OUSolution, min_val, mean_val, max_val, std_dev)
  !> Computes statistics (in space) of a given array representing the solution (at a given time) of an Ornstein-Uhlenbeck process.
  !> Contains the entire domain, where many values may be zero !  
  !> This subroutine calculates the minimum, maximum, mean, standard deviation for the provided array.
  !>
  !> INPUT:
  !> - OUSolution: Array representing the time series of an Ornstein-Uhlenbeck process.
  !>
  !> OUTPUT:
  !> - min_val: Minimum value of the array.
  !> - max_val: Maximum value of the array.
  !> - mean_val: Mean value of the array.
  !> - std_dev: Standard deviation of the array.

  implicit none
  real(wp), dimension(:,:), intent(in) :: OUSolution
  real(wp), intent(out) :: min_val, mean_val, max_val, std_dev 
  real(wp) :: sumSol, sum_sq
  integer :: n

  ! Get array size
  n = size(OUSolution)

  ! Check for empty array
  if (n == 0) then
    write(*, *) "Error: Empty array passed to statsOUAtGivenTime"
    stop
  end if

  ! Calculate min and max
  min_val = minval(OUSolution)
  max_val = maxval(OUSolution)

  ! Calculate mean
  sumSol = sum(OUSolution)
  mean_val = sumSol / real(n)

  ! Calculate standard deviation
  sum_sq = sum((OUSolution - mean_val)**2)
  std_dev = sqrt(sum_sq / real(n))

  !WRITE(*,*) 'min_val',min_val
  !WRITE(*,*) 'max_val',max_val
  !WRITE(*,*) 'mean_val',mean_val
end subroutine OUBasicStatsInSpaceAtGivenTime


subroutine percentilesArrayAtGivenTime(Array2d, percentiles)
!> Compute the percentile of the OU solutions at actual time
!> Percentiles found using the quicksort algoritm
!> Should add a flag to conside only the cells where h>0
  implicit none
  real(wp), dimension(:, :), intent(in) :: Array2d
  real(wp), dimension(9), intent(out) :: percentiles
  real(wp), dimension(:), allocatable :: flatArray
  integer :: n

  ! Get total number of elements
  n = size(Array2d) 

  ! Check for empty array
  if (n == 0) then
    write(*, *) "Error: Empty array passed to statsOUAtGivenTime"
    stop
  end if

  ! Flatten the 2D array
  allocate(flatArray(n))
  flatArray = reshape(Array2d, shape(flatArray))

  ! Sort Flattened array
  call quickSort(flatArray, 1, n)

  ! Calculate percentiles
  percentiles(1) = flatArray(percentileIndex(n, 5._wp))
  percentiles(2) = flatArray(percentileIndex(n, 10._wp))
  percentiles(3) = flatArray(percentileIndex(n, 20._wp))
  percentiles(4) = flatArray(percentileIndex(n, 30._wp))
  percentiles(5) = flatArray(percentileIndex(n, 50._wp))
  percentiles(6) = flatArray(percentileIndex(n, 70._wp))
  percentiles(7) = flatArray(percentileIndex(n, 80._wp))
  percentiles(8) = flatArray(percentileIndex(n, 90._wp))
  percentiles(9) = flatArray(percentileIndex(n, 95._wp))

  ! Deallocate temporary array
  deallocate(flatArray)

end subroutine percentilesArrayAtGivenTime


subroutine quickSort(Array1d, start_idx, end_idx)
!> Sort the given array using the QuickSort algorithm
!> Run-time complexity : 
!> 1) Best case O(n log(n))
!> 2) Average case O(n log(n))
!> 3) Worst case O(n**2)
    implicit none
    real(wp), intent(in out) :: Array1d(:)
    integer, intent(in) :: start_idx, end_idx ! idx array
    integer :: pivot ! partition index

    ! Check if there are more than one element in the array
    if (start_idx < end_idx) then 
        pivot = partition(Array1d, start_idx, end_idx)
        ! Recursively apply quickSort to the left and right subarrays
        call quickSort(Array1d, start_idx, pivot - 1)
        call quickSort(Array1d, pivot + 1, end_idx)
    end if
end subroutine quickSort


integer function partition(Array1d, start_idx, end_idx)
!> Helper function of the quickSort algorithm
    implicit none
    real(wp), intent(in out) :: Array1d(:)
    integer, intent(in) :: start_idx, end_idx
    integer :: i, j
    ! pivot and temp hold elements of the REAL array Array1d; declaring
    ! them INTEGER truncated every value that passed through them
    real(wp) :: pivot, temp

    ! Get the pivot element 
    pivot = Array1d(end_idx)
    ! Initialize the index of the smaller element
    i = start_idx - 1
    ! Loop through the array and rearrange elements based on the pivot value
    do j = start_idx, end_idx - 1
        if (Array1d(j) <= pivot) then
            ! Swap arr(i+1) and arr(j)
            i = i + 1
            temp = Array1d(i)
            Array1d(i) = Array1d(j)
            Array1d(j) = temp
        end if
    end do

    ! Swap arr(i+1) and arr(high) to place the pivot in its correct position
    i = i + 1
    temp = Array1d(i)
    Array1d(i) = Array1d(end_idx)
    Array1d(end_idx) = temp
    partition = i
end function partition


integer function percentileIndex(size_arr, p)
!> Find the index of the percentile p in array of size n
  integer, intent(in) :: size_arr
  real(wp), intent(in) :: p
  percentileIndex = max(1, min(size_arr, ceiling(real(size_arr) * p / 100._wp)))
end function percentileIndex
  

END MODULE stochastic_module
