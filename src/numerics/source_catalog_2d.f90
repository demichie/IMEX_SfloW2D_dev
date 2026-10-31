MODULE source_catalog_2d
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: source_term_type, SOURCE_CATALOG, SOURCE_TERM_COUNT
  PUBLIC :: source_term_descriptor_by_name, source_catalog_names_unique
  PUBLIC :: canonical_rheology_source_name, rheology_name_has_role, canonical_rheology_catalog_unique
  INTEGER, PARAMETER, PUBLIC :: SOURCE_METHOD_EXPLICIT=1, SOURCE_METHOD_NEWTON=2, SOURCE_METHOD_ANALYTIC=3
  INTEGER, PARAMETER, PUBLIC :: SOURCE_CADENCE_STAGE=1, SOURCE_CADENCE_STEP=2
  INTEGER, PARAMETER, PUBLIC :: SOURCE_SOLVE_NEWTON_BLOCK=1, SOURCE_SOLVE_ANALYTIC_VELOCITY=2
  INTEGER, PARAMETER, PUBLIC :: SOURCE_SOLVE_ANALYTIC_TEMPERATURE=3, SOURCE_SOLVE_ANALYTIC_PORE_PRESSURE=4
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_IMPLICIT=1, SOURCE_ROLE_SEMI_IMPLICIT=2
  ! Equation-group membership is configured independently of source identity.
  TYPE :: source_term_type
    CHARACTER(LEN=64) :: name = ''
    LOGICAL :: has_explicit = .FALSE.
    LOGICAL :: has_newton = .FALSE.
    LOGICAL :: has_analytic = .FALSE.
    INTEGER :: default_method = SOURCE_METHOD_NEWTON
    INTEGER :: default_cadence = SOURCE_CADENCE_STAGE
    ! Explicit per-entry permissions; absent permissions never enable a new source.
    LOGICAL :: allows_gas = .FALSE.
    LOGICAL :: allows_liquid = .FALSE.
    LOGICAL :: allows_lava = .FALSE.
  END TYPE source_term_type
  ! Default method/cadence are explicit in each catalog entry; constants
  ! SOURCE_METHOD_* remain generic method IDs also used by the solvers.
  INTEGER, PARAMETER :: SOURCE_TERM_COUNT=24
  TYPE(source_term_type), PARAMETER :: SOURCE_CATALOG(SOURCE_TERM_COUNT)=[ &
       source_term_type('rheology_voellmy',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('rheology_plastic',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('rheology_lava',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.FALSE., allows_liquid=.FALSE., allows_lava=.TRUE.), &
       source_term_type('rheology_lahar',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.FALSE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('rheology_empirical',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('rheology_quadratic',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('rheology_froude',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('rheology_lucas',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('rheology_mui',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('rheology_mui_voellmy',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('rheology_newmodel',.FALSE.,.TRUE.,.FALSE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('thermal_cooling',.TRUE.,.TRUE.,.TRUE., &
             default_method=SOURCE_METHOD_NEWTON, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.TRUE.), &
       source_term_type('pore_pressure_dissipation',.FALSE.,.FALSE.,.TRUE., &
             default_method=SOURCE_METHOD_ANALYTIC, &
             default_cadence=SOURCE_CADENCE_STEP, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('bottom_circular_injection',.TRUE.,.FALSE.,.FALSE., &
            default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('bottom_fissural_injection',.TRUE.,.FALSE.,.FALSE., &
            default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('radial_injection',.TRUE.,.FALSE.,.FALSE., &
             default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('geometry_curvature',.TRUE.,.FALSE.,.FALSE., &
             default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('boundary_injection',.TRUE.,.FALSE.,.FALSE., &
             default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STAGE, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('mass_exchange_lava_deposition',.TRUE.,.FALSE.,.FALSE., &
            default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STEP, &
             allows_gas=.FALSE., allows_liquid=.FALSE., allows_lava=.TRUE.), &
       source_term_type('mass_exchange_erosion',.TRUE.,.FALSE.,.FALSE., &
            default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STEP, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('mass_exchange_deposition',.TRUE.,.FALSE.,.FALSE., &
            default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STEP, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('mass_exchange_carrier_liquid_loss',.TRUE.,.FALSE.,.FALSE., &
            default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STEP, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.), &
       source_term_type('mass_exchange_entrainment',.TRUE.,.FALSE.,.FALSE., &
            default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STEP, &
             allows_gas=.TRUE., allows_liquid=.FALSE., allows_lava=.FALSE.), &
       source_term_type('mass_exchange_vaporization',.TRUE.,.FALSE.,.FALSE., &
            default_method=SOURCE_METHOD_EXPLICIT, &
             default_cadence=SOURCE_CADENCE_STEP, &
             allows_gas=.TRUE., allows_liquid=.TRUE., allows_lava=.FALSE.) &
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
  PURE FUNCTION canonical_rheology_source_name(model_name) RESULT(name)
    CHARACTER(*), INTENT(IN) :: model_name
    CHARACTER(LEN=64) :: name
    INTEGER :: i
    name=''
    IF (TRIM(model_name)=='none') RETURN
    DO i=1,SOURCE_TERM_COUNT
      IF (INDEX(TRIM(SOURCE_CATALOG(i)%name),'rheology_')/=1) CYCLE
      IF (TRIM(model_name)==TRIM(SOURCE_CATALOG(i)%name(10:))) THEN
        name=SOURCE_CATALOG(i)%name
        RETURN
      END IF
    END DO
  END FUNCTION canonical_rheology_source_name
  PURE LOGICAL FUNCTION canonical_rheology_catalog_unique() RESULT(ok)
    INTEGER :: i,j
    ok=.FALSE.
    DO i=1,SOURCE_TERM_COUNT
      IF (INDEX(TRIM(SOURCE_CATALOG(i)%name),'rheology_')/=1) CYCLE
      DO j=i+1,SOURCE_TERM_COUNT
        IF (INDEX(TRIM(SOURCE_CATALOG(j)%name),'rheology_')/=1) CYCLE
        IF (SOURCE_CATALOG(i)%name==SOURCE_CATALOG(j)%name) RETURN
      END DO
    END DO
    ok=.TRUE.
  END FUNCTION canonical_rheology_catalog_unique
  !> Check a rheology's actual SI/IM contribution by canonical source name.
  !> This is independent of the legacy numeric model ID.
  PURE LOGICAL FUNCTION rheology_name_has_role(source_name,role) RESULT(found)
    CHARACTER(*), INTENT(IN) :: source_name
    INTEGER, INTENT(IN) :: role
    found=.FALSE.
    SELECT CASE(role)
    CASE(SOURCE_ROLE_IMPLICIT)
      SELECT CASE(TRIM(source_name))
      CASE('rheology_voellmy','rheology_plastic','rheology_lava', &
           'rheology_lahar','rheology_empirical','rheology_quadratic', &
           'rheology_mui_voellmy','rheology_newmodel')
        found=.TRUE.
      END SELECT
    CASE(SOURCE_ROLE_SEMI_IMPLICIT)
      SELECT CASE(TRIM(source_name))
      CASE('rheology_voellmy','rheology_lava','rheology_lahar', &
           'rheology_froude','rheology_lucas','rheology_mui', &
           'rheology_mui_voellmy')
        found=.TRUE.
      END SELECT
    END SELECT
  END FUNCTION rheology_name_has_role
END MODULE source_catalog_2d
