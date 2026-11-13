!********************************************************************************
!> \brief Shared material and physical parameters.
!>
!> Stores configuration read during setup that is shared across constitutive
!> models. Rheology-specific coefficients belong to their model modules.
!********************************************************************************

MODULE constitutive_parameters_2d

    USE parameters_2d, ONLY : wp

    IMPLICIT NONE

    PRIVATE

    ! Public constitutive-model configuration. Keep this list explicit so that
    ! implementation details and names imported from other modules are not
    ! re-exported accidentally.
    PUBLIC :: entrainment_flag
    PUBLIC :: grav, inv_grav
    PUBLIC :: A_drag, B_drag, collective_settling_flag
    PUBLIC :: SBconst
    PUBLIC :: rho_c_sub, rho_a_amb
    PUBLIC :: sp_heat_a, sp_gas_const_a, kin_visc_a
    PUBLIC :: sp_heat_g, sp_gas_const_g
    PUBLIC :: kin_visc_l, kin_visc_c, Tref_Suth, muRef_Suth, S_mu, T_ambient
    PUBLIC :: rho_s, inv_rho_s, diam_s, sphericity_s, sp_heat_s
    PUBLIC :: settling_flag, alphastot_min
    PUBLIC :: erosion_coeff, loss_rate, erodible_fract, erodible_porosity
    PUBLIC :: coeff_porosity, T_erodible
    PUBLIC :: pres, inv_pres, rho_l, inv_rho_l, sp_heat_l, gamma_steam, vonK
    PUBLIC :: maximum_solid_packing

    !> flag to activate air entrainment
    LOGICAL :: entrainment_flag

    !> gravitational acceleration
    REAL(wp) :: grav
    REAL(wp) :: inv_grav


    !> Optional collective-settling drag-law coefficients
    REAL(wp) :: A_drag
    REAL(wp) :: B_drag
    LOGICAL :: collective_settling_flag
    !> Stefan-Boltzmann constant [W m-2 K-4].
    REAL(wp), PARAMETER :: SBconst = 5.67E-8_wp


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

    !> Maximum solid packing fraction
    REAL(wp) :: maximum_solid_packing

END MODULE constitutive_parameters_2d
