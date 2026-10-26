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
  ! Integer IDs remain only for equation-block targets, not source catalog identities.
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TARGET_VELOCITY=1, SOURCE_TARGET_TEMPERATURE=2
  INTEGER, PARAMETER, PUBLIC :: SOURCE_TARGET_PORE_PRESSURE=3
  INTEGER, PARAMETER, PUBLIC :: SOURCE_ROLE_IMPLICIT=1, SOURCE_ROLE_SEMI_IMPLICIT=2
  TYPE :: source_term_type
    CHARACTER(LEN=64) :: name = ''
    CHARACTER(LEN=32) :: target = ''
    LOGICAL :: has_explicit = .FALSE.
    LOGICAL :: has_newton = .FALSE.
    LOGICAL :: has_analytic = .FALSE.
    INTEGER :: default_method = SOURCE_METHOD_NEWTON
    INTEGER :: default_cadence = SOURCE_CADENCE_STAGE
  END TYPE source_term_type
  INTEGER, PARAMETER :: SOURCE_TERM_COUNT=24
  TYPE(source_term_type), PARAMETER :: SOURCE_CATALOG(SOURCE_TERM_COUNT)=[ &
       source_term_type('rheology_voellmy','velocity',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_plastic','velocity',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_lava','velocity',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_lahar','velocity',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_empirical','velocity',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_quadratic','velocity',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_froude','velocity',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_lucas','velocity',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_mui','velocity',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_mui_voellmy','velocity',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('rheology_newmodel','velocity',.TRUE.,.TRUE.,.FALSE.), &
       source_term_type('thermal_cooling','temperature',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('pore_pressure_dissipation','pore_pressure',.TRUE.,.TRUE.,.TRUE.), &
       source_term_type('bottom_circular_injection','explicit',.TRUE.,.FALSE.,.FALSE., &
            SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STAGE), &
       source_term_type('bottom_fissural_injection','explicit',.TRUE.,.FALSE.,.FALSE., &
            SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STAGE), &
       source_term_type('radial_injection','explicit',.TRUE.,.FALSE.,.FALSE.,SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STAGE), &
       source_term_type('geometry_curvature','explicit',.TRUE.,.FALSE.,.FALSE.,SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STAGE), &
       source_term_type('boundary_injection','explicit',.TRUE.,.FALSE.,.FALSE.,SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STAGE), &
       source_term_type('mass_exchange_erosion','mass_exchange',.TRUE.,.FALSE.,.FALSE., &
            SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STEP), &
       source_term_type('mass_exchange_deposition','mass_exchange',.TRUE.,.FALSE.,.FALSE., &
            SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STEP), &
       source_term_type('mass_exchange_carrier_liquid_loss','mass_exchange',.TRUE.,.FALSE.,.FALSE., &
            SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STEP), &
       source_term_type('mass_exchange_carrier_darcy_loss','mass_exchange',.TRUE.,.FALSE.,.FALSE., &
            SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STEP), &
       source_term_type('mass_exchange_entrainment','mass_exchange',.TRUE.,.FALSE.,.FALSE., &
            SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STEP), &
       source_term_type('mass_exchange_vaporization','mass_exchange',.TRUE.,.FALSE.,.FALSE., &
            SOURCE_METHOD_EXPLICIT,SOURCE_CADENCE_STEP) &
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
