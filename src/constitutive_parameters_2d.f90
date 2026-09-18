!********************************************************************************
!> \brief Constitutive model parameters
!********************************************************************************
MODULE constitutive_parameters_2d

  USE parameters_2d, ONLY : wp

  IMPLICIT NONE

  !> flag to activate air entrainment
  LOGICAL :: entrainment_flag

  !> gravitational acceleration
  REAL(wp) :: grav
  REAL(wp) :: inv_grav


  !> Optional collective-settling drag-law coefficients
  REAL(wp) :: A_drag
  REAL(wp) :: B_drag
  LOGICAL :: collective_settling_flag

  !> drag coefficients (Voellmy-Salm model)
  REAL(wp) :: mu
  REAL(wp) :: xi
  REAL(wp) :: xi_temp

  !> friction coefficients (function Coulomb: mu(Fr) and mu:(U))
  REAL(wp) :: mu_0
  REAL(wp) :: mu_inf
  REAL(wp) :: Fr_0
  REAL(wp) :: U_w

  !> mu(I) rheology parameters
  REAL(wp) :: I_0 !< reference inertial number
  REAL(wp) :: mu_2 !< friction at high inertial number above which flow accelerates
  REAL(wp) :: mu_s !< static friction coefficient
  REAL(wp) :: muI_inf !< friction at high I to avoid plateau (Barker et al. 2017)
  REAL(wp) :: I_transition = 4.0d-3 !< transition inertial number

   !> f_inhibit
   REAL(wp) :: alpha_trans !< solid volume fraction at which f_inhibit starts decreasing
   INTEGER :: N_inh = 2 !< 2N+1 determines order of the smoothstep function for f_inhibit
   CHARACTER(LEN=8) :: f_inhibit_mode = 'OFF'  ! options: 'OFF','STATIC','DYNAMIC'

  !> drag coefficients (B&W model)
  REAL(wp) :: friction_factor

  !> drag coefficients (plastic model)
  REAL(wp) :: tau

  !> evironment temperature [K]
  REAL(wp) :: T_env

  !> reference temperature [K]
  REAL(wp) :: T_ref

  !> reference kinematic viscosity [m2/s]
  REAL(wp) :: nu_ref

  !> viscosity parameter [K-1] (b in Table 1 Costa & Macedonio, 2005)
  REAL(wp) :: visc_par

  !> yield strength for lava rheology [kg m-1 s-2] (Eq.4 Kelfoun & Varga, 2015)
  REAL(wp) :: tau0

  !> velocity boundary layer fraction of total thickness
  REAL(wp) :: emme

  !> specific heat [J kg-1 K-1]
  REAL(wp) :: c_p

  !> coefficient for the convective term in the temperature equation for lava
  REAL(wp) :: convective_term_coeff

  !> atmospheric heat trasnfer coefficient [W m-2 K-1] (lambda in C&M, 2005)
  REAL(wp) :: atm_heat_transf_coeff

  !> fractional area of the exposed inner core (f in C&M, 2005)
  REAL(wp) :: exp_area_fract

  !> coefficient for the radiative term in the temperature equation for lava
  REAL(wp) :: radiative_term_coeff

  !> Stephan-Boltzmann constant [W m-2 K-4]
  REAL(wp), PARAMETER :: SBconst = 5.67E-8_wp

  !> emissivity (eps in Costa & Macedonio, 2005)
  REAL(wp) :: emissivity

  !> thermal boundary layer fraction of total thickness
  REAL(wp) :: enne

  !> temperature of lava-ground interface
  REAL(wp) :: T_ground

  !> thermal conductivity [W m-1 K-1] (k in Costa & Macedonio, 2005)
  REAL(wp) :: thermal_conductivity

  !--- START Lahars rheology model parameters

  !> 1st param for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL(wp) :: alpha2    ! (units: kg m-1 s-2)

  !> 2nd param for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL(wp) :: beta2     ! (units: nondimensional)

  !> ratio between reference value from input and computed values from eq.
  REAL(wp) :: alpha1_coeff ! (units: nondimensional )

  !> 2nd param for fluid viscosity empirical relationship (O'Brian et al, 1993)
  REAL(wp) :: beta1     ! (units: nondimensional, input parameter)

  !> Empirical resistance parameter (dimensionless, input parameter)
  REAL(wp) :: Kappa

  !> Mannings roughness coefficient ( units: T L^(-1/3) )
  REAL(wp) :: n_td
  REAL(wp) :: n_td2

  !--- END Lahars rheology model parameters

  !> Specific heat of carrier phase (gas or liquid)
  REAL(wp) :: sp_heat_c  ! ( initialized from input)

  !> Specific gas constant of gas mixture (units: J kg-1 K-1)

  !> Density of carrier phase in substrate ( units: kg m-3 )
  REAL(wp) :: rho_c_sub

  !> Ambient density of air ( units: kg m-3 )
  REAL(wp) :: rho_a_amb

  !> Specific heat of air (units: J K-1 kg-1)
  REAL(wp) :: sp_heat_a

  !> Specific gas constant of air (units: J kg-1 K-1)
  REAL(wp) :: sp_gas_const_a

  !> Kinematic viscosity of air (units: m2 s-1)
  REAL(wp) :: kin_visc_a

  !> Specific heat of additional gas (units: J K-1 kg-1)
  REAL(wp), ALLOCATABLE :: sp_heat_g(:)

  !> Specific gas constant of additional gas (units: J kg-1 K-1)
  REAL(wp), ALLOCATABLE :: sp_gas_const_g(:)

  !> Kinematic viscosity of liquid (units: m2 s-1)
  REAL(wp) :: kin_visc_l

  !> Kinematic viscosity of carrier phase (units: m2 s-1)
  REAL(wp) :: kin_visc_c

  !> Reference temperature for Sutherland’s law (units: K)
  REAL(wp) :: Tref_Suth

  !> Reference viscosity for Sutherland's law (units: kg m-1 s-1)
  REAL(wp) :: muRef_Suth

  !> Sutherland constant for Sutherland’s law (units: K)
  REAL(wp) :: S_mu

  !> Temperature of ambient air (units: K)
  REAL(wp) :: T_ambient

  !> Density of sediments ( units: kg m-3 )
  REAL(wp), ALLOCATABLE :: rho_s(:)

  !> Reciprocal of density of sediments ( units: kg m-3 )
  REAL(wp), ALLOCATABLE :: inv_rho_s(:)

  !> Diameter of sediments ( units: m )
  REAL(wp), ALLOCATABLE :: diam_s(:)

  !> Sphericity of sediments (dimensionless, default 1.0)
  REAL(wp), ALLOCATABLE :: sphericity_s(:)

  !> Specific heat of solids (units: J K-1 kg-1)
  REAL(wp), ALLOCATABLE :: sp_heat_s(:)

  !> Flag to determine if sedimentation is active
  LOGICAL :: settling_flag

  !> Minimum volume fraction of solids in the flow
  REAL(wp) :: alphastot_min

  !> erosion model coefficient  (units: m-1 )
  REAL(wp), ALLOCATABLE :: erosion_coeff

  !> water loss_rate (unit: m s-1)
  REAL(wp), ALLOCATABLE :: loss_rate

  !> erodible substrate solid relative volume fractions
  REAL(wp), ALLOCATABLE :: erodible_fract(:)

  !> erodible substrate porosity (we assume filled by continous phase)
  REAL(wp) :: erodible_porosity

  !> coefficient to compute (eroded/deposited) volume of continuous phase
  !> from volume of solid
  REAL(wp) :: coeff_porosity

  !> temperature of erodible substrate (units: K)
  REAL(wp) :: T_erodible

  !> ambient pressure (units: Pa)
  REAL(wp) :: pres

  !> reciprocal of ambient pressure (units: Pa)
  REAL(wp) :: inv_pres

  !> liquid density (units: kg m-3)
  REAL(wp) :: rho_l

  !> reciprocal of liquid density (units: kg m-3)
  REAL(wp) :: inv_rho_l

  !> Sepcific heat of liquid (units: J K-1 kg-1)
  REAL(wp) :: sp_heat_l

  !> Fraction of heat lost by particles producing steam
  REAL(wp) :: gamma_steam

  !> Von Karman constant used for the diagnostic Rouse number
  REAL(wp) :: vonK

  !> Hydraulic permeability (units: m2)
  REAL(wp) :: hydraulic_permeability
  !> Flag to activate dynamic permeability
  LOGICAL :: dynamic_permeability_flag

  !> Maximum solid packing fraction
  REAL(wp) :: maximum_solid_packing

  ! store combined pascal coefficients: coeff(n) = pascal1(n) * pascal2(n)
  REAL(wp), ALLOCATABLE :: pascal_coeff(:)   ! index with 0..N_inh
  LOGICAL :: pascal_coeff_precomputed

CONTAINS

  !> Function that calculates the pascal coefficients for the smooth transition of f_inhibit
  SUBROUTINE precompute_pascal_coefficient(N_order)
    !! Computes the two Pascal-like coefficients for each n = 0..N
    !! Returns coeff(2, N+1)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N_order
    INTEGER :: n, k
    REAL(wp) :: pascal1, pascal2

    ! allocate or reallocate pascal_coeff with 0-based indices
   IF ( ALLOCATED(pascal_coeff) ) THEN
      IF ( LBOUND(pascal_coeff,1) /= 0 .OR. UBOUND(pascal_coeff,1) /= N_order ) THEN
         DEALLOCATE(pascal_coeff)
         ALLOCATE(pascal_coeff(0:N_order))
      END IF
   ELSE
      ALLOCATE(pascal_coeff(0:N_order))
   END IF

    ! Validate input
    IF (N_order < 0) THEN
        WRITE(*,*) 'ERROR: N must be >= 0. N =', N_order
        STOP
    END IF

    ! Loop over n = 0..N and compute the two coefficients
    DO n = 0, N_order
        ! ----- First coefficient: pascal_triangle(-N-1, n) -----
        pascal1 = 1.0_wp
        IF (n > 0) THEN
            DO k = 0, n-1
               pascal1 = pascal1 * ( -REAL(N_order,kind=wp) - 1.0_wp - &
                              REAL(k, kind=wp) ) / REAL(k+1,kind=wp)
            END DO
        END IF
        ! ----- Second coefficient: pascal_triangle(2N+1, N-n) -----
        pascal2 = 1.0_wp
        IF (N_order - n > 0) THEN
            DO k = 0, N_order - n - 1
               pascal2 = pascal2 * (2.0_wp*REAL(N_order,kind=wp)+1.0_wp- &
                              REAL(k,kind=wp) ) / REAL(k+1,kind=wp)
            END DO
        END IF

        pascal_coeff(n) = pascal1 * pascal2
    END DO

      pascal_coeff_precomputed = .TRUE.

 END SUBROUTINE precompute_pascal_coefficient

END MODULE constitutive_parameters_2d
