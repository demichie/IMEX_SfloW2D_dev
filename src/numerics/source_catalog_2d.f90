MODULE source_catalog_2d
  USE parameters_2d, ONLY : rheology_id_voellmy, rheology_id_plastic, &
       rheology_id_lava, rheology_id_lahar, rheology_id_empirical, &
       rheology_id_quadratic, rheology_id_froude, rheology_id_lucas, &
       rheology_id_mui, rheology_id_mui_voellmy
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: source_term_type, SOURCE_CATALOG, SOURCE_TERM_COUNT
  PUBLIC :: source_term_descriptor_by_name, source_catalog_names_unique
  PUBLIC :: rheology_source_name, rheology_has_role, canonical_rheology_catalog_unique
  INTEGER, PARAMETER, PUBLIC :: SOURCE_METHOD_EXPLICIT=1, SOURCE_METHOD_NEWTON=2, SOURCE_METHOD_ANALYTIC=3
  INTEGER, PARAMETER, PUBLIC :: SOURCE_CADENCE_STAGE=1, SOURCE_CADENCE_STEP=2
  INTEGER, PARAMETER, PUBLIC :: SOURCE_PHASE_RK_IMPLICIT=1, SOURCE_PHASE_RK_EXPLICIT=2, SOURCE_PHASE_POST_RK=3
  INTEGER, PARAMETER, PUBLIC :: SOURCE_SOLVE_NEWTON_BLOCK=1, SOURCE_SOLVE_ANALYTIC_VELOCITY=2
  INTEGER, PARAMETER, PUBLIC :: SOURCE_SOLVE_ANALYTIC_TEMPERATURE=3, SOURCE_SOLVE_ANALYTIC_PORE_PRESSURE=4
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_IMPLICIT=1, SOURCE_ROLE_SEMI_IMPLICIT=2
  TYPE :: source_term_type
    CHARACTER(LEN=64) :: name = ''
    CHARACTER(LEN=32) :: target = ''
    LOGICAL :: has_explicit = .FALSE.
    LOGICAL :: has_newton = .FALSE.
    LOGICAL :: has_analytic = .FALSE.
  END TYPE source_term_type
  INTEGER, PARAMETER :: SOURCE_TERM_COUNT=24
  TYPE(source_term_type), PARAMETER :: SOURCE_CATALOG(SOURCE_TERM_COUNT)=[ &
       source_term_type('rheology_voellmy','velocity',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_plastic','velocity',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_lava','velocity',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_lahar','velocity',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_empirical','velocity',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_quadratic','velocity',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_froude','velocity',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_lucas','velocity',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_mui','velocity',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_mui_voellmy','velocity',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('thermal_cooling','temperature',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('pore_pressure_dissipation','pore_pressure',.FALSE.,.TRUE.,.TRUE.), &
       source_term_type('bottom_circular_injection','explicit',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('bottom_fissural_injection','explicit',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('radial_injection','explicit',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('geometry_curvature','explicit',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('boundary_injection','explicit',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('mass_exchange_erosion','mass_exchange',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('mass_exchange_deposition','mass_exchange',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('mass_exchange_carrier_liquid_loss','mass_exchange',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('mass_exchange_carrier_deposition_loss','mass_exchange',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('mass_exchange_carrier_darcy_loss','mass_exchange',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('mass_exchange_entrainment','mass_exchange',.TRUE.,.FALSE.,.FALSE.), &
       source_term_type('mass_exchange_vaporization','mass_exchange',.TRUE.,.FALSE.,.FALSE.) &
       ]
CONTAINS
  PURE FUNCTION source_term_descriptor_by_name(name) RESULT(term)
    CHARACTER(*), INTENT(IN) :: name
    TYPE(source_term_type) :: term
    INTEGER :: i
    term=source_term_type()
    DO i=1,SOURCE_TERM_COUNT
      IF (TRIM(name)==TRIM(SOURCE_CATALOG(i)%name)) THEN
        term=SOURCE_CATALOG(i)
        RETURN
      END IF
    END DO
  END FUNCTION source_term_descriptor_by_name
  PURE LOGICAL FUNCTION source_catalog_names_unique() RESULT(ok)
    INTEGER :: i,j
    ok=.FALSE.
    DO i=1,SOURCE_TERM_COUNT
      IF (LEN_TRIM(SOURCE_CATALOG(i)%name)==0) RETURN
      DO j=i+1,SOURCE_TERM_COUNT
        IF (SOURCE_CATALOG(i)%name==SOURCE_CATALOG(j)%name) RETURN
      END DO
    END DO
    ok=.TRUE.
  END FUNCTION source_catalog_names_unique
  PURE FUNCTION rheology_source_name(model_id) RESULT(name)
    INTEGER, INTENT(IN) :: model_id
    CHARACTER(LEN=64) :: name
    name=''
    SELECT CASE(model_id)
    CASE(rheology_id_voellmy); name='rheology_voellmy'
    CASE(rheology_id_plastic); name='rheology_plastic'
    CASE(rheology_id_lava); name='rheology_lava'
    CASE(rheology_id_lahar); name='rheology_lahar'
    CASE(rheology_id_empirical); name='rheology_empirical'
    CASE(rheology_id_quadratic); name='rheology_quadratic'
    CASE(rheology_id_froude); name='rheology_froude'
    CASE(rheology_id_lucas); name='rheology_lucas'
    CASE(rheology_id_mui); name='rheology_mui'
    CASE(rheology_id_mui_voellmy); name='rheology_mui_voellmy'
    END SELECT
  END FUNCTION rheology_source_name
  PURE LOGICAL FUNCTION canonical_rheology_catalog_unique() RESULT(ok)
    INTEGER :: i,j
    ok=.FALSE.
    DO i=1,10
      DO j=i+1,10
        IF (SOURCE_CATALOG(i)%name==SOURCE_CATALOG(j)%name) RETURN
      END DO
    END DO
    ok=.TRUE.
  END FUNCTION canonical_rheology_catalog_unique
  PURE LOGICAL FUNCTION rheology_has_role(model_id,role) RESULT(found)
    INTEGER, INTENT(IN) :: model_id,role
    found=.FALSE.
    SELECT CASE(role)
    CASE(SOURCE_ROLE_IMPLICIT)
      SELECT CASE(model_id)
      CASE(rheology_id_voellmy,rheology_id_plastic,rheology_id_lava, &
           rheology_id_lahar,rheology_id_empirical,rheology_id_quadratic,rheology_id_mui_voellmy)
        found=.TRUE.
      END SELECT
    CASE(SOURCE_ROLE_SEMI_IMPLICIT)
      SELECT CASE(model_id)
      CASE(rheology_id_voellmy,rheology_id_lava,rheology_id_lahar, &
           rheology_id_froude,rheology_id_lucas,rheology_id_mui,rheology_id_mui_voellmy)
        found=.TRUE.
      END SELECT
    END SELECT
  END FUNCTION rheology_has_role
END MODULE source_catalog_2d
