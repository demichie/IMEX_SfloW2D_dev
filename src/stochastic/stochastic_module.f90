!********************************************************************************
!> \brief Stochastic module
!
!> This module contains all the procedures for the stochastic variable
!
!> \date 11/06/2025
!> @author 
!> Zeno Geddo
!
!>
!> Separates the Gaussian Ornstein-Uhlenbeck state Z, its conservative transport, and the
!> transformed effective_Z used by friction. Random samples are generated outside parallel cell
!> updates.
!********************************************************************************

MODULE stochastic_module
  
  ! external variables
  USE parameters_2d, ONLY : wp, rheology_model, idx_stoch, idx_u, idx_v
  
  USE domain_2d, ONLY: domain_type
  USE state_2d, ONLY: state_type
  USE parameters_2d, ONLY : length_spatial_corr,        &
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
  

  !> \brief Gaussian OU state, friction-only transformation and correlation kernel.
  !> \details Z and effective_Z are cell maps; only Z is conservatively transported
  !>          and checkpointed. Friction receives the separately derived effective_Z.
  TYPE :: stochastic_workspace_type

     !> Stochastic field at cell centers
     REAL(wp), ALLOCATABLE :: Z(:,:)

     !> Possibly transformed fluctuation passed to the friction law. Z remains
     !> the Gaussian OU state transported through the conservative variable M*Z.
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

  !> \brief Allocate OU/friction fields and seed the generator when stochastic physics is enabled.
  !>
  !> \param[in,out] this Gaussian OU field, effective-friction field and correlation kernel.

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

  !> \brief Release OU, effective-friction and convolution-kernel fields.
  !>
  !> \param[in,out] this Gaussian OU field, effective-friction field and correlation kernel.

  SUBROUTINE finalize_stochastic_workspace(this)

    CLASS(stochastic_workspace_type), INTENT(INOUT) :: this

    IF ( ALLOCATED(this%Z) ) DEALLOCATE(this%Z)
    IF ( ALLOCATED(this%effective_Z) ) DEALLOCATE(this%effective_Z)
    IF ( ALLOCATED(this%conv_kernel) ) DEALLOCATE(this%conv_kernel)

  END SUBROUTINE finalize_stochastic_workspace

  !> \brief Burn in the OU field and initialize the transported stochastic state when configured.
  !>
  !> \param[in,out] this Gaussian OU field, effective-friction field and correlation kernel.
  !> \param[in,out] state Current flow state; stochastic q/qp entries are initialized after burn-in
  !>                      when transport is enabled.
  !> \param[in] domain Active-cell/face lists and reconstruction halo for this simulation.
  !>
  !> \note This is a finite burn-in, not a proof of statistical stationarity when the
  !>       velocity-dependent noise intensity changes.

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

  !> \brief Build the normalized spatial-correlation convolution kernel.
  !>
  !> \param[in,out] this Gaussian OU field, effective-friction field and correlation kernel.
  !>
  !> \note Reads length_spatial_corr and cell_size; writes this%conv_kernel for subsequent spatial
  !>       filtering of Gaussian samples.

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

  !> \brief Sample the centered Gaussian used to construct the spatial-correlation kernel.
  !>
  !> \param[in] x X displacement from the correlation-kernel centre [m].
  !> \param[in] y Y displacement from the correlation-kernel centre [m].
  !> \return Gaussian kernel value at the requested spatial offset.
  !>
  !> \note Uses standard deviation length_spatial_corr/2; called for positive correlation length.

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

  !> \brief Advance the OU field by Euler-Maruyama using velocity-dependent noise intensity.
  !>
  !> \param[in,out] this Gaussian OU field, effective-friction field and correlation kernel.
  !> \param[in] state Simulation cell states and accumulated diagnostics.
  !> \param[in] dt Time increment [s].

  SUBROUTINE update_stochastic_variable(this, state, dt)
    ! OU evolution is a fractional step distinct from conservative transport.
    IMPLICIT NONE
    CLASS(stochastic_workspace_type), INTENT(INOUT) :: this
    CLASS(state_type), INTENT(IN) :: state
    INTEGER :: j,k
    INTEGER :: noise_size
    REAL(wp), INTENT(IN) :: dt
    REAL(wp) :: sigma_noise
    REAL(wp):: noise(comp_cells_x, comp_cells_y)
    REAL(wp):: conv_result(comp_cells_x, comp_cells_y)

    ! Draw random numbers serially before the OpenMP update: thread scheduling
    ! must not change the random stream associated with a fixed seed.
    ! Generate standard Gaussian noise over the entire domain: N(0,1).
    noise_size = comp_cells_x*comp_cells_y
    noise = RESHAPE(gaussian_noise(noise_size),                              &
         [comp_cells_x, comp_cells_y])

    ! Convolve the gaussian noise to introduce spatial correlation if needed
    IF (length_spatial_corr > 0.0_wp) THEN
        ! (for convenience, the conv_kernel is generated only once before the burn in)
        CALL convolve_2d(noise, this%conv_kernel, conv_result)
        noise = conv_result
    END IF

    ! Evolve dry cells too, so an OU state exists when they become wet later.
    ! Parallelism here changes neither random-number generation nor ownership.
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

    
    CALL this%refresh_effective
   
  END SUBROUTINE update_stochastic_variable

  !> \brief Apply the OU fractional step and initialize conservative stochastic transport.
  !>
  !> \param[in,out] this Gaussian OU field, effective-friction field and correlation kernel.
  !> \param[in,out] state Current flow state; transported stochastic q/qp entries are updated after
  !>                      the OU step.
  !> \param[in] dt Time increment [s].

  SUBROUTINE prepare_stochastic_timestep(this, state, dt)
    ! Apply the OU fractional step and initialize M*Z for conservative transport.
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

  !> \brief Derive the friction fluctuation from Z without changing the Gaussian OU state.
  !>
  !> \param[in,out] this Gaussian OU field, effective-friction field and correlation kernel.
  !>
  !> \note The transformed fluctuation is for friction only; conservative transport and restart
  !>       retain the Gaussian OU state.

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

 !> \brief Return squared cell speed from the current physical velocities.
 !>
 !> \param[in] state Simulation cell states and accumulated diagnostics.
 !> \param[in] j X index of the current computational cell.
 !> \param[in] k Y index of the current computational cell.
 !> \return Squared speed [m^2 s^-2], or zero when the cell mass is negligible.

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

  !> \brief Evaluate the configured noise intensity for a wet or dry cell.
  !>
  !> \param[in] state Simulation cell states and accumulated diagnostics.
  !> \param[in] j X index of the current computational cell.
  !> \param[in] k Y index of the current computational cell.
  !> \return Local OU noise intensity; dry cells use std_max.

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

  !> \brief Evaluate the exponential transition between minimum and maximum noise intensities.
  !>
  !> \param[in] val Squared velocity entering the noise-intensity transition [m^2 s^-2].
  !> \return Noise intensity std_max+(std_min-std_max)*exp(-val/std_slope_factor).

  REAL(wp) FUNCTION expFormNoise(val)
    ! Compute the bounded intensity of the noise of the stochastic process
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: val
    expFormNoise = std_max + ( std_min - std_max) * exp(- val / std_slope_factor )
  END FUNCTION expFormNoise

  !> \brief Apply one Euler-Maruyama step to a scalar Ornstein-Uhlenbeck state.
  !>
  !> \param[in] Zij Current Gaussian OU state before this stochastic update.
  !> \param[in] dt Time increment [s].
  !> \param[in] sigma Local OU noise intensity.
  !> \param[in] noise_ij Standard-normal forcing sample for the current cell.
  !> \return Updated Gaussian OU state after damping and stochastic forcing.

  REAL(wp) FUNCTION EulerMaruyamaScheme(Zij, dt, sigma, noise_ij)
    ! Update the Ornstein-Uhlenback process using EULER-MARUYAMA Method
    ! Z(j,k) = Z(j,k) - (dt / tau_stochastic) * Z(j,k) + &
    ! sigma_noise * SQRT(2.0_wp * (dt / tau_stochastic)) * noise(j,k)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: Zij, dt, sigma, noise_ij
    EulerMaruyamaScheme = Zij - (dt / tau_stochastic) * Zij +             &
            sigma * SQRT(2.0_wp * (dt / tau_stochastic)) * noise_ij
  END FUNCTION

!> \brief Convolve a cell-centered signal using a centred kernel and zero padding.
!>
!> \param[in] input_signal Cell-centered input field indexed as (x,y).
!> \param[in] kernel Centred convolution weights indexed as (x-offset,y-offset).
!> \param[out] result Convolved field with the same shape as input_signal.

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
  

END MODULE stochastic_module
