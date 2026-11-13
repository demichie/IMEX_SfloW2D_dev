!> @brief Canonical SHARED terrain configuration and derived substrate properties.
!> @details One substrate is used by all material layers. Dynamic bed inventories
!>          and topographic geometry remain owned by the simulation state/geometry.
MODULE substrate_parameters_2d
    USE parameters_2d, ONLY : wp
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: substrate_config_type, substrate, erodible_fract, coeff_porosity

    !> Shared terrain and erodible-bed configuration (no layer-specific copies).
    TYPE :: substrate_config_type
        REAL(wp) :: T_ground=-1.0_wp              !< Ground temperature [K].
        REAL(wp) :: erodible_porosity=-1.0_wp     !< Pore volume fraction of the bed.
        CHARACTER(LEN=40) :: topography_file='topography_dem.asc'
        CHARACTER(LEN=40) :: erodible_file=''
        REAL(wp) :: initial_erodible_thickness=-1.0_wp
        LOGICAL :: topo_change_flag=.FALSE.
        LOGICAL :: erodible_deposit_flag=.FALSE.
    END TYPE substrate_config_type

    TYPE(substrate_config_type), TARGET :: substrate

    !> Canonical class fractions are kept outside the derived type because gfortran
    !> NAMELIST input cannot read a derived object with ALLOCATABLE components.
    !> The array is allocated to n_solid before SUBSTRATE_PARAMETERS is read.
    REAL(wp), ALLOCATABLE :: erodible_fract(:)

    !> Derived values, not independent input parameters.
    REAL(wp) :: coeff_porosity=0.0_wp
END MODULE substrate_parameters_2d
