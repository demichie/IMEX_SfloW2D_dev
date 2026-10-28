!********************************************************************************
!> \brief Translate configured source activation into validated D17 metadata.
!> \details Runs once from physical source settings. Owns no
!>          physical state and performs no rate evaluation or temporal update.
!>          Per-cell wetness, packing and window gates remain in existing kernels.
!********************************************************************************
MODULE source_configuration_2d
  USE source_catalog_2d
  USE source_registry_2d, ONLY : source_registry_type, source_group_type
  USE equation_metadata_2d, ONLY : equation_partition_type, &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON
  USE source_input_2d
  USE mass_exchange_lava_deposition_2d, ONLY : lava_deposition_flag
  USE parameters_2d, ONLY : wp, n_eqns, rheology_flag, rheology_model, &
        thermal_cooling_flag, &
       pore_pressure_flag, curvature_term_flag, bottom_radial_source_flag, &
       bottom_fissural_source_flag, radial_source_flag, lateral_source_flag, &
       liquid_flag, gas_flag, gas_loss_flag, liquid_vaporization_flag, &
       idx_pore
  USE constitutive_parameters_2d, ONLY : erosion_coeff, settling_flag, &
       coeff_porosity, entrainment_flag
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_configured_sources
CONTAINS
  !> Identify active physical sources exactly once for both input syntaxes.
  !> Thermal activation is determined by the physical cooling flag.
  SUBROUTINE collect_active_source_names(include_thermal,active_names,error)
    LOGICAL, INTENT(IN) :: include_thermal
    CHARACTER(LEN=64), ALLOCATABLE, INTENT(OUT) :: active_names(:)
    CHARACTER(*), INTENT(OUT) :: error
    ALLOCATE(active_names(0))
    IF (rheology_flag) THEN
      IF (LEN_TRIM(canonical_rheology_source_name(rheology_model))>0) &
           active_names=[CHARACTER(LEN=64) :: active_names,canonical_rheology_source_name(rheology_model)]
      IF (LEN_TRIM(canonical_rheology_source_name(rheology_model))==0 .AND. &
           TRIM(rheology_model)/='none') THEN
        error='D17: unsupported configured rheology'; RETURN
      END IF
    END IF
    IF (include_thermal) active_names=[CHARACTER(LEN=64) :: active_names, &
         'thermal_cooling']
    IF (pore_pressure_flag) active_names=[CHARACTER(LEN=64) :: active_names,'pore_pressure_dissipation']
    IF (curvature_term_flag) active_names=[CHARACTER(LEN=64) :: active_names,'geometry_curvature']
    IF (bottom_radial_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'bottom_circular_injection']
    IF (bottom_fissural_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'bottom_fissural_injection']
    IF (radial_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'radial_injection']
    IF (lateral_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'boundary_injection']
    IF (ALLOCATED(erosion_coeff)) THEN
      IF (erosion_coeff>0.0_wp) &
           active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_erosion']
    END IF
    IF (lava_deposition_flag) THEN
      IF (.NOT.liquid_flag .OR. TRIM(rheology_model)/='lava') THEN
        error='Lava deposition requires liquid lava rheology'; RETURN
      END IF
      active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_lava_deposition']
    END IF
    IF (settling_flag) THEN
      active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_deposition']
    END IF
    IF (liquid_flag) active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_carrier_liquid_loss']
    IF (pore_pressure_flag .AND. gas_loss_flag) &
         active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_carrier_darcy_loss']
    IF (entrainment_flag) active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_entrainment']
    IF (gas_flag .AND. liquid_vaporization_flag) &
         active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_vaporization']
    error=''
  END SUBROUTINE collect_active_source_names

  PURE INTEGER FUNCTION group_method_rank(group_method) RESULT(rank)
    INTEGER, INTENT(IN) :: group_method
    rank=1
    IF (group_method==SOURCE_METHOD_EXPLICIT) rank=2
  END FUNCTION group_method_rank

  !> \brief Build simulation-owned source metadata from current validated input.
  !> \param[in,out] registry Destination, cleared if translation fails.
  !> \param[in] partition Existing compact equation partition, never modified.
  !> \param[out] error Empty on success, otherwise a pre-simulation diagnostic.
  !> \note This phase records possible activation; zero per-cell rates stay valid.
  !> Construct final physical groups once, deriving the equation partition.
  !> No temporary registry or subsequent source reassignment is constructed.
  SUBROUTINE initialize_configured_sources(registry,partition,error)
    TYPE(source_registry_type), INTENT(INOUT) :: registry
    TYPE(equation_partition_type), INTENT(INOUT) :: partition
    CHARACTER(*), INTENT(OUT) :: error
    CHARACTER(LEN=64), ALLOCATABLE :: active_names(:)
    TYPE(source_group_type), ALLOCATABLE :: stages(:),steps(:)
    TYPE(source_group_type) :: group,tmp
    TYPE(source_term_type) :: term
    CHARACTER(LEN=16) :: requested_method,requested_cadence
    INTEGER :: i,j,k,id,method,cadence,pore,eq,compact,n
    INTEGER :: target
    INTEGER :: membership(n_eqns)
    INTEGER, ALLOCATABLE :: partition_methods(:)
    INTEGER, ALLOCATABLE :: stage_order(:)
    INTEGER :: swap_order
    LOGICAL :: mask(n_eqns)
    CALL registry%finalize()
    error=''
    IF (n_eqns<4) THEN
      error='D18: invalid equation layout'; RETURN
    END IF
    CALL collect_active_source_names(thermal_cooling_flag,active_names,error)
    IF (LEN_TRIM(error)>0) RETURN
    pore=0
    IF (pore_pressure_flag) pore=idx_pore
    ALLOCATE(stages(0),steps(0),stage_order(0))
    DO i=1,SIZE(active_names)
      term=source_term_descriptor_by_name(active_names(i))
      IF (LEN_TRIM(term%name)==0) THEN
        error='D18: unknown active source'; RETURN
      END IF
      method=term%default_method
      cadence=term%default_cadence
      id=0
      requested_method=''; requested_cadence=''
      SELECT CASE(TRIM(term%target))
      CASE('velocity')
        target=1
        requested_method=rheology_method
        requested_cadence=rheology_cadence
        id=rheology_group
      CASE('temperature')
        target=2
        requested_method=thermal_method
        requested_cadence=thermal_cadence
        id=thermal_group
      CASE('pore_pressure')
        target=3
        requested_method=pore_pressure_method
        requested_cadence=pore_pressure_cadence
        id=pore_pressure_group
      CASE('explicit')
        target=0
        SELECT CASE(TRIM(active_names(i)))
        CASE('bottom_circular_injection')
          requested_method=bottom_circular_method
          requested_cadence=bottom_circular_cadence
          id=bottom_circular_group
        CASE('bottom_fissural_injection')
          requested_method=bottom_fissural_method
          requested_cadence=bottom_fissural_cadence
          id=bottom_fissural_group
        CASE('radial_injection')
          requested_method=radial_injection_method
          requested_cadence=radial_injection_cadence
          id=radial_injection_group
        CASE('geometry_curvature')
          requested_method=geometry_curvature_method
          requested_cadence=geometry_curvature_cadence
          id=geometry_curvature_group
        CASE('boundary_injection')
          requested_method=boundary_injection_method
          requested_cadence=boundary_injection_cadence
          id=boundary_injection_group
        END SELECT
      CASE('mass_exchange')
        target=-1
        requested_method=mass_exchange_method
        requested_cadence=mass_exchange_cadence
        id=mass_exchange_group
      CASE DEFAULT
        target=0
      END SELECT
      IF (id<0) THEN
        error='D18: GROUP must be nonnegative'; RETURN
      END IF
      IF (LEN_TRIM(requested_method)>0) THEN
        SELECT CASE(TRIM(requested_method))
        CASE('EXPLICIT'); method=SOURCE_METHOD_EXPLICIT
        CASE('NEWTON'); method=SOURCE_METHOD_NEWTON
        CASE('ANALYTIC'); method=SOURCE_METHOD_ANALYTIC
        CASE DEFAULT
          error='D18: invalid source METHOD'; RETURN
        END SELECT
      END IF
      IF (LEN_TRIM(requested_cadence)>0) THEN
        SELECT CASE(TRIM(requested_cadence))
        CASE('STAGE'); cadence=SOURCE_CADENCE_STAGE
        CASE('STEP'); cadence=SOURCE_CADENCE_STEP
        CASE DEFAULT
          error='D18: invalid source CADENCE'; RETURN
        END SELECT
      END IF
      IF (method==SOURCE_METHOD_EXPLICIT .AND. .NOT.term%has_explicit) THEN
        error='D18: explicit interface unavailable for source'; RETURN
      END IF
      IF (method==SOURCE_METHOD_ANALYTIC .AND. .NOT.term%has_analytic) THEN
        error='D18: analytic interface unavailable for source'; RETURN
      END IF
      IF (method==SOURCE_METHOD_NEWTON .AND. .NOT.term%has_newton) THEN
        error='D18: Newton interface unavailable for source'; RETURN
      END IF
      IF (target>0 .AND. method==SOURCE_METHOD_ANALYTIC .AND. &
           .NOT.term%has_analytic) THEN
        error='D18: analytic interface unavailable for source'; RETURN
      END IF
      IF (target>0 .AND. method==SOURCE_METHOD_NEWTON .AND. &
           .NOT.term%has_newton) THEN
        error='D18: Newton interface unavailable for source'; RETURN
      END IF
      k=0
      IF (cadence==SOURCE_CADENCE_STAGE) THEN
        DO j=1,SIZE(stages)
          IF (stages(j)%method==method .AND. stages(j)%id==id) k=j
        END DO
      ELSE
        DO j=1,SIZE(steps)
          IF (steps(j)%method==method .AND. steps(j)%id==id) k=j
        END DO
      END IF
      ! A numeric ID identifies one STEP operator, independent of method.
      ! At STAGE, only EXPLICIT and NEWTON may share an ID: the IMEX
      ! tableau determines their distinct numerical contributions.
      IF (cadence==SOURCE_CADENCE_STEP) THEN
        DO j=1,SIZE(steps)
          IF (steps(j)%id==id .AND. steps(j)%method/=method) THEN
            error='D18: STEP GROUP cannot mix METHOD values'
            RETURN
          END IF
        END DO
      ELSE
        DO j=1,SIZE(stages)
          IF (stages(j)%id/=id .OR. stages(j)%method==method) CYCLE
          IF (method==SOURCE_METHOD_ANALYTIC .OR. &
               stages(j)%method==SOURCE_METHOD_ANALYTIC) THEN
            error='D18: ANALYTIC STAGE GROUP must be exclusive'
            RETURN
          END IF
        END DO
      END IF
      IF (k/=0 .AND. method==SOURCE_METHOD_ANALYTIC) THEN
        error='D18: ANALYTIC GROUP must contain a single source'
        RETURN
      END IF
      IF (k==0) THEN
        group%id=id
        group%method=method; group%cadence=cadence
        group%source_names=[CHARACTER(LEN=64)::active_names(i)]
        group%equation_ids=[INTEGER::]
        IF (method==SOURCE_METHOD_EXPLICIT) THEN
          group%equation_ids=[(n,n=1,n_eqns)]
        ELSE
          SELECT CASE(target)
          CASE(1); group%equation_ids=[2,3]
          CASE(2); group%equation_ids=[4]
          CASE(3); group%equation_ids=[pore]
          END SELECT
        END IF
        IF (cadence==SOURCE_CADENCE_STAGE) THEN
          stages=[stages,group]
          stage_order=[stage_order,id]
        ELSE
          steps=[steps,group]
        END IF
      ELSE IF (cadence==SOURCE_CADENCE_STAGE) THEN
        stages(k)%source_names=[CHARACTER(LEN=64)::stages(k)%source_names,active_names(i)]
        IF (method/=SOURCE_METHOD_EXPLICIT) THEN
          SELECT CASE(target)
          CASE(1); stages(k)%equation_ids=[stages(k)%equation_ids,2,3]
          CASE(2); stages(k)%equation_ids=[stages(k)%equation_ids,4]
          CASE(3); stages(k)%equation_ids=[stages(k)%equation_ids,pore]
          END SELECT
        END IF
      ELSE
        steps(k)%source_names=[CHARACTER(LEN=64)::steps(k)%source_names,active_names(i)]
        IF (method/=SOURCE_METHOD_EXPLICIT) THEN
          SELECT CASE(target)
          CASE(1); steps(k)%equation_ids=[steps(k)%equation_ids,2,3]
          CASE(2); steps(k)%equation_ids=[steps(k)%equation_ids,4]
          CASE(3); steps(k)%equation_ids=[steps(k)%equation_ids,pore]
          END SELECT
        END IF
      END IF
    END DO
    ! STEP is a true post-RK operator sequence ordered by GROUP.
    ! Sources inside a single group use the same input state.
    DO i=1,SIZE(steps)
      DO j=i+1,SIZE(steps)
        IF (steps(j)%id<steps(i)%id) THEN
          tmp=steps(i); steps(i)=steps(j); steps(j)=tmp
        END IF
      END DO
    END DO
    ! Keep an empty explicit stage group required by the RK stage interface.
    IF (.NOT.ANY([(stages(j)%method==SOURCE_METHOD_EXPLICIT,j=1,SIZE(stages))])) THEN
      group%id=0
      group%method=SOURCE_METHOD_EXPLICIT
      group%cadence=SOURCE_CADENCE_STAGE
      group%source_names=[CHARACTER(LEN=64)::]
      group%equation_ids=[(n,n=1,n_eqns)]
      stages=[stages,group]
      stage_order=[stage_order,0]
    END IF
    ! Keep only the empty stage group required by the RK interface.
    ! Runtime array order derives from method class and GROUP; no extra state.
    DO i=1,SIZE(stages)
      DO j=i+1,SIZE(stages)
        ! Ordering priority is fixed during group construction.
        ! No query of the pre-existing equation partition in this sort.
        IF (stage_order(j)<stage_order(i) .OR. &
             (stage_order(j)==stage_order(i) .AND. &
              group_method_rank(stages(j)%method)<group_method_rank(stages(i)%method))) THEN
          tmp=stages(i); stages(i)=stages(j); stages(j)=tmp
          swap_order=stage_order(i)
          stage_order(i)=stage_order(j)
          stage_order(j)=swap_order
        END IF
      END DO
    END DO
    ALLOCATE(partition_methods(SIZE(stages)))
    mask=.FALSE.; membership=0; partition_methods=0; compact=0
    DO i=1,SIZE(stages)
      IF (stages(i)%method==SOURCE_METHOD_EXPLICIT) CYCLE
      compact=compact+1
      partition_methods(compact)=IM_GROUP_METHOD_NEWTON
      IF (stages(i)%method==SOURCE_METHOD_ANALYTIC) &
           partition_methods(compact)=IM_GROUP_METHOD_ANALYTIC
      DO j=1,SIZE(stages(i)%equation_ids)
        eq=stages(i)%equation_ids(j)
        IF (membership(eq)>0) THEN
          error='D18: equation belongs to multiple implicit groups'; RETURN
        END IF
        membership(eq)=compact; mask(eq)=.TRUE.
      END DO
    END DO
    CALL partition%initialize(mask,membership,partition_methods(1:compact))
    CALL registry%initialize(n_eqns,pore,stages,steps,error)
  END SUBROUTINE initialize_configured_sources

END MODULE source_configuration_2d
