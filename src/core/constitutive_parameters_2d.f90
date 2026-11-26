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
    PUBLIC :: grav, inv_grav
    PUBLIC :: SBconst
    PUBLIC :: rho_a_amb
    PUBLIC :: sp_heat_a, sp_gas_const_a, kin_visc_a
    PUBLIC :: sp_heat_g, sp_gas_const_g
    PUBLIC :: kin_visc_l, kin_visc_c, Tref_Suth, muRef_Suth, S_mu, T_ambient
    PUBLIC :: alphastot_min
    PUBLIC :: pres, inv_pres, rho_l, inv_rho_l, sp_heat_l, vonK

    !> gravitational acceleration
    REAL(wp) :: grav
    REAL(wp) :: inv_grav


    !> Stefan-Boltzmann constant [W m-2 K-4].
    REAL(wp), PARAMETER :: SBconst = 5.67E-8_wp


    !> Specific gas constant of gas mixture (units: J kg-1 K-1)

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

    !> Minimum volume fraction of solids in the flow
    REAL(wp) :: alphastot_min

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

    !> Von Karman constant used for the diagnostic Rouse number
    REAL(wp) :: vonK


END MODULE constitutive_parameters_2d
