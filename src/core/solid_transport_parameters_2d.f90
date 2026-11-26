!********************************************************************************
!> @brief Canonical SHARED properties of solid species and packing.
!> @details All layers use the same solid-species catalog. Derived density
!>          reciprocals are computed from rho_s during input initialization.
!>          Process-specific drag and deposition laws are configured separately.
!********************************************************************************
MODULE solid_transport_parameters_2d
    USE parameters_2d, ONLY : wp
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: rho_s, inv_rho_s, diam_s, sphericity_s, sp_heat_s
    PUBLIC :: maximum_solid_packing

    !> Density of each solid species [kg m-3].
    REAL(wp), ALLOCATABLE :: rho_s(:)

    !> Derived reciprocal of the species density [m3 kg-1]; never read from input.
    REAL(wp), ALLOCATABLE :: inv_rho_s(:)

    !> Diameter of each solid species [m].
    REAL(wp), ALLOCATABLE :: diam_s(:)

    !> Species sphericity [-]; defaults are set by the input initialization.
    REAL(wp), ALLOCATABLE :: sphericity_s(:)

    !> Specific heat of each solid species [J kg-1 K-1].
    REAL(wp), ALLOCATABLE :: sp_heat_s(:)

    !> Shared maximum physically admissible total solid volume fraction [-].
    REAL(wp) :: maximum_solid_packing
END MODULE solid_transport_parameters_2d
