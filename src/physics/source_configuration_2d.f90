!********************************************************************************
!> \brief Translate configured source activation into validated D17 metadata.
!> \details Runs once after legacy equation partition initialization. Owns no
!>          physical state and performs no rate evaluation or temporal update.
!>          Per-cell wetness, packing and window gates remain in existing kernels.
!********************************************************************************
MODULE source_configuration_2d
  USE source_catalog_2d
  USE source_registry_2d, ONLY : source_registry_type, source_group_type
  USE equation_metadata_2d, ONLY : equation_partition_type, resolve_imex_group_method, &
       IM_GROUP_METHOD_ANALYTIC, IM_GROUP_METHOD_NEWTON
  USE source_input_2d
  USE parameters_2d, ONLY : wp, n_eqns, rheology_flag, rheology_model_id, &
       rheology_id_none, rheology_id_lava, &
       pore_pressure_flag, curvature_term_flag, bottom_radial_source_flag, &
       bottom_fissural_source_flag, radial_source_flag, lateral_source_flag, &
       liquid_flag, gas_flag, gas_loss_flag, liquid_vaporization_flag, &
       imex_velocity_group, imex_energy_group, imex_pore_pressure_group, imex_group_method, idx_pore
  USE constitutive_parameters_2d, ONLY : erosion_coeff, settling_flag, &
       coeff_porosity, entrainment_flag
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_configured_sources
  PUBLIC :: initialize_native_source_configuration
CONTAINS
  !> \brief Build simulation-owned source metadata from current validated input.
  !> \param[in,out] registry Destination, cleared if translation fails.
  !> \param[in] partition Existing compact equation partition, never modified.
  !> \param[out] error Empty on success, otherwise a pre-simulation diagnostic.
  !> \note This phase records possible activation; zero per-cell rates stay valid.
  SUBROUTINE initialize_configured_sources(registry,partition,error)
    TYPE(source_registry_type), INTENT(INOUT) :: registry
    TYPE(equation_partition_type), INTENT(IN) :: partition
    CHARACTER(*), INTENT(OUT) :: error
    INTEGER, ALLOCATABLE :: terms(:)
    INTEGER :: methods(3),i,pore
    CHARACTER(LEN=16) :: method
    CALL registry%finalize()
    error=''
    IF (.NOT.partition%is_initialized()) THEN
      error='D17: source configuration requires initialized equation metadata'; RETURN
    END IF
    IF (partition%n_equations<4) THEN
      error='D17: source layout requires at least four equations'; RETURN
    END IF
    terms=[INTEGER::]
    IF (rheology_flag) THEN
      terms=rheology_source_terms(rheology_model_id)
      IF (SIZE(terms)==0 .AND. rheology_model_id/=rheology_id_none) THEN
        error='D17: unsupported configured rheology'; RETURN
      END IF
    END IF
    IF (partition%implicit(4)) terms=[terms,SOURCE_TERM_THERMAL_RADIATIVE_COOLING, &
         SOURCE_TERM_THERMAL_CONVECTIVE_COOLING]
    pore=0
    IF (pore_pressure_flag) THEN
      pore=idx_pore
      terms=[terms,SOURCE_TERM_PORE_PRESSURE_DISSIPATION]
    END IF
    IF (curvature_term_flag) terms=[terms,SOURCE_TERM_GEOMETRY_CURVATURE]
    IF (bottom_radial_source_flag) terms=[terms,SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION]
    IF (bottom_fissural_source_flag) terms=[terms,SOURCE_TERM_BOTTOM_FISSURAL_INJECTION]
    IF (radial_source_flag) terms=[terms,SOURCE_TERM_RADIAL_INJECTION]
    IF (lateral_source_flag) terms=[terms,SOURCE_TERM_BOUNDARY_INJECTION]
    IF (erosion_coeff>0.0_wp) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_EROSION]
    IF (settling_flag) THEN
      terms=[terms,SOURCE_TERM_MASS_EXCHANGE_DEPOSITION]
      IF (ABS(coeff_porosity)>0.0_wp) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS]
    END IF
    IF (liquid_flag) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS]
    IF (pore_pressure_flag .AND. gas_loss_flag) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS]
    IF (entrainment_flag) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT]
    IF (gas_flag .AND. liquid_vaporization_flag) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION]
    DO i=1,3
      method=imex_group_method(i)
      CALL resolve_imex_group_method(method,methods(i))
    END DO
    CALL registry%initialize_legacy(partition, &
         [imex_velocity_group,imex_energy_group,imex_pore_pressure_group], &
         methods,pore,terms,error)
  END SUBROUTINE initialize_configured_sources

  !> \brief Translate the native namelist into an equation partition and registry.
  !> \param[in,out] registry Registry destination, cleared on error.
  !> \param[in,out] partition Equation partition destination, finalized on error.
  !> \param[out] error Empty on success, otherwise a startup diagnostic.
  SUBROUTINE initialize_native_source_configuration(registry,partition,error)
    TYPE(source_registry_type), INTENT(INOUT) :: registry
    TYPE(equation_partition_type), INTENT(INOUT) :: partition
    CHARACTER(*), INTENT(OUT) :: error
    TYPE(source_group_type), ALLOCATABLE :: stages(:),steps(:)
    TYPE(source_group_type) :: group
    INTEGER, ALLOCATABLE :: parsed_terms(:,:),ordered(:),targets(:),equations(:)
    INTEGER, ALLOCATABLE :: group_terms(:),active_terms(:),law_terms(:)
    INTEGER :: cadence(SOURCE_INPUT_MAX_GROUPS),method(SOURCE_INPUT_MAX_GROUPS)
    INTEGER :: phase(SOURCE_INPUT_MAX_GROUPS),solver_count(SOURCE_INPUT_MAX_GROUPS)
    INTEGER :: parsed_count(SOURCE_INPUT_MAX_GROUPS)
    INTEGER :: target_count(1,SOURCE_INPUT_MAX_GROUPS)
    INTEGER :: target_ids(SOURCE_INPUT_MAX_TARGETS,1, &
         SOURCE_INPUT_MAX_GROUPS)
    INTEGER :: target_group(3),methods(SOURCE_INPUT_MAX_GROUPS),membership(n_eqns)
    LOGICAL :: implicit_mask(n_eqns),seen_terms(SOURCE_TERM_COUNT),expected_targets(3)
    INTEGER :: i,j,k,g,s,t,slot,compact_group,stage_count,step_count,term_id,target
    CHARACTER(LEN=256) :: validation_error

    CALL registry%finalize()
    CALL partition%finalize()
    error=''
    IF (.NOT.native_source_configuration_present) THEN
      error='D17: native source configuration was not supplied'; RETURN
    END IF
    IF (source_group_count<1 .OR. source_group_count>SOURCE_INPUT_MAX_GROUPS) THEN
      error='D17: SOURCE_GROUP_COUNT is outside the supported range'; RETURN
    END IF
    IF (n_eqns<3 .OR. (rheology_model_id==rheology_id_lava .AND. n_eqns<4)) THEN
      error='D17: native source configuration requires a valid equation layout'; RETURN
    END IF
    IF (pore_pressure_flag .AND. (idx_pore<1 .OR. idx_pore>n_eqns .OR. &
         idx_pore==2 .OR. idx_pore==3 .OR. &
         (rheology_model_id==rheology_id_lava .AND. idx_pore==4))) THEN
      error='D17: pore-pressure equation index is outside the configured layout'; RETURN
    END IF

    ALLOCATE(parsed_terms(SOURCE_INPUT_MAX_TERMS,source_group_count))
    parsed_terms=0; parsed_count=0; cadence=0; method=0; phase=0; solver_count=0
    target_count=0; target_ids=0
    seen_terms=.FALSE.
    DO g=1,source_group_count
      error='D17: source group IDs and execution orders must be positive'
      IF (source_group_id(g)<1 .OR. source_group_order(g)<1) RETURN
      DO i=1,g-1
        IF (source_group_cadence(g)==source_group_cadence(i) .AND. &
             source_group_id(g)==source_group_id(i)) THEN
          error='D17: duplicate source group ID within a cadence'; RETURN
        END IF
      END DO
      cadence(g)=native_cadence(source_group_cadence(g))
      method(g)=native_method(source_group_method(g))
      error='D17: unsupported source group cadence or method'
      IF (cadence(g)==0 .OR. method(g)==0) RETURN
      SELECT CASE(cadence(g))
      CASE(SOURCE_CADENCE_STAGE)
        IF (method(g)==SOURCE_METHOD_EXPLICIT) THEN
          phase(g)=SOURCE_PHASE_RK_EXPLICIT
        ELSE
          phase(g)=SOURCE_PHASE_RK_IMPLICIT
        END IF
      CASE(SOURCE_CADENCE_STEP)
        phase(g)=SOURCE_PHASE_POST_RK
        IF (method(g)/=SOURCE_METHOD_EXPLICIT) THEN
          error='D17: post-RK source groups must be explicit'; RETURN
        END IF
      END SELECT
      DO i=1,g-1
        IF (phase(g)==phase(i) .AND. cadence(g)==cadence(i) .AND. &
             source_group_order(g)==source_group_order(i)) THEN
          error='D17: duplicate execution order within a phase'; RETURN
        END IF
      END DO
      IF (source_group_term_count(g)<0 .OR. &
           source_group_term_count(g)>SOURCE_INPUT_MAX_TERMS) THEN
        error='D17: invalid source term count'; RETURN
      END IF
      IF (source_group_term_count(g)==0 .AND. method(g)==SOURCE_METHOD_EXPLICIT) THEN
        error='D17: explicit source groups must declare at least one term'; RETURN
      END IF
      DO i=1,source_group_term_count(g)
        ! One canonical rheology name expands to the complete legacy law.
        IF (rheology_flag .AND. &
             LEN_TRIM(rheology_source_name(rheology_model_id))>0 .AND. &
             TRIM(uppercase_source_text(source_group_terms(i,g)))== &
             TRIM(uppercase_source_text(rheology_source_name(rheology_model_id)))) THEN
          law_terms=rheology_source_terms(rheology_model_id)
          error='D17: unsupported complete rheology source'
          IF (SIZE(law_terms)==0) RETURN
        ELSE
          term_id=native_term_id(source_group_terms(i,g))
          error='D17: unknown, dependent or unschedulable source term name'
          IF (term_id==0) RETURN
          law_terms=[term_id]
        END IF
        DO t=1,SIZE(law_terms)
          term_id=law_terms(t)
          error='D17: source group expansion exceeds supported catalog capacity'
          IF (parsed_count(g)>=SOURCE_INPUT_MAX_TERMS) RETURN
          error='D17: a source term may be registered only once globally'
          IF (seen_terms(term_id)) RETURN
          parsed_count(g)=parsed_count(g)+1
          seen_terms(term_id)=.TRUE.
          parsed_terms(parsed_count(g),g)=term_id
        END DO
      END DO
      ! An implicit group has a single target list, not nested solvers.
      error='D17: invalid group target count'
      IF (source_group_target_count(g)<0 .OR. &
           source_group_target_count(g)>SOURCE_INPUT_MAX_TARGETS) RETURN
      IF (method(g)==SOURCE_METHOD_EXPLICIT) THEN
        error='D17: explicit source groups cannot declare group targets'
        IF (source_group_target_count(g)/=0) RETURN
        IF (ANY(source_group_targets(:,g)/='')) RETURN
        solver_count(g)=0
      ELSE
        error='D17: Newton and analytic groups require group targets'
        IF (source_group_target_count(g)<1) RETURN
        solver_count(g)=1
        target_count(1,g)=source_group_target_count(g)
        DO t=1,target_count(1,g)
          target=native_target(source_group_targets(t,g))
          error='D17: unknown or duplicate group target'
          IF (target==0) RETURN
          IF (ANY(target_ids(:,1,g)==target)) RETURN
          target_ids(t,1,g)=target
        END DO
        IF (method(g)==SOURCE_METHOD_ANALYTIC) THEN
          error='D17: analytic groups must have exactly one target'
          IF (target_count(1,g)/=1) RETURN
        END IF
      END IF
    END DO

    active_terms=[INTEGER::]
    IF (rheology_flag) THEN
      active_terms=rheology_source_terms(rheology_model_id)
      IF (SIZE(active_terms)==0 .AND. rheology_model_id/=rheology_id_none) THEN
        error='D17: unsupported configured rheology'; RETURN
      END IF
    END IF
    IF (rheology_model_id==rheology_id_lava) active_terms=[active_terms, &
         SOURCE_TERM_THERMAL_RADIATIVE_COOLING,SOURCE_TERM_THERMAL_CONVECTIVE_COOLING]
    IF (pore_pressure_flag) active_terms=[active_terms,SOURCE_TERM_PORE_PRESSURE_DISSIPATION]
    IF (curvature_term_flag) active_terms=[active_terms,SOURCE_TERM_GEOMETRY_CURVATURE]
    IF (bottom_radial_source_flag) active_terms=[active_terms,SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION]
    IF (bottom_fissural_source_flag) active_terms=[active_terms,SOURCE_TERM_BOTTOM_FISSURAL_INJECTION]
    IF (radial_source_flag) active_terms=[active_terms,SOURCE_TERM_RADIAL_INJECTION]
    IF (lateral_source_flag) active_terms=[active_terms,SOURCE_TERM_BOUNDARY_INJECTION]
    IF (erosion_coeff>0.0_wp) active_terms=[active_terms,SOURCE_TERM_MASS_EXCHANGE_EROSION]
    IF (settling_flag) THEN
      active_terms=[active_terms,SOURCE_TERM_MASS_EXCHANGE_DEPOSITION]
      IF (ABS(coeff_porosity)>0.0_wp) &
           active_terms=[active_terms,SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS]
    END IF
    IF (liquid_flag) active_terms=[active_terms,SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS]
    IF (pore_pressure_flag .AND. gas_loss_flag) &
         active_terms=[active_terms,SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS]
    IF (entrainment_flag) active_terms=[active_terms,SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT]
    IF (gas_flag .AND. liquid_vaporization_flag) &
         active_terms=[active_terms,SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION]
    error='D17: native schedule must register every configured active term exactly once'
    IF (COUNT(seen_terms)/=SIZE(active_terms)) RETURN
    DO i=1,SIZE(active_terms)
      IF (.NOT.seen_terms(active_terms(i))) RETURN
    END DO

    target_group=0
    expected_targets=[.TRUE.,rheology_model_id==rheology_id_lava,pore_pressure_flag]
    DO g=1,source_group_count
      IF (cadence(g)/=SOURCE_CADENCE_STAGE .OR. method(g)==SOURCE_METHOD_EXPLICIT) CYCLE
      DO s=1,solver_count(g)
        DO t=1,target_count(s,g)
          target=target_ids(t,s,g)
          IF (.NOT.expected_targets(target)) THEN
            error='D17: source solver target is not implicit in the configured model'; RETURN
          END IF
          IF (target_group(target)/=0) THEN
            error='D17: implicit target is assigned to more than one source group'; RETURN
          END IF
          target_group(target)=g
        END DO
      END DO
    END DO
    DO target=1,3
      IF (expected_targets(target) .NEQV. (target_group(target)/=0)) THEN
        error='D17: native schedule is missing a required implicit target'; RETURN
      END IF
    END DO

    implicit_mask=.FALSE.; membership=0
    implicit_mask(2:3)=.TRUE.
    IF (expected_targets(2)) implicit_mask(4)=.TRUE.
    IF (expected_targets(3)) implicit_mask(idx_pore)=.TRUE.
    ordered=[(g,g=1,source_group_count)]
    DO i=2,source_group_count
      slot=ordered(i); j=i-1
      DO WHILE (j>=1)
        g=ordered(j)
        IF (phase(g)<phase(slot)) EXIT
        IF (phase(g)==phase(slot) .AND. source_group_order(g)<=source_group_order(slot)) EXIT
        ordered(j+1)=ordered(j); j=j-1
      END DO
      ordered(j+1)=slot
    END DO

    compact_group=0
    DO i=1,source_group_count
      g=ordered(i)
      IF (phase(g)/=SOURCE_PHASE_RK_IMPLICIT) CYCLE
      compact_group=compact_group+1
      methods(compact_group)=IM_GROUP_METHOD_NEWTON
      IF (method(g)==SOURCE_METHOD_ANALYTIC) methods(compact_group)=IM_GROUP_METHOD_ANALYTIC
      ! Exactly one solver is admitted for each non-explicit group.
      DO t=1,target_count(1,g)
        SELECT CASE(target_ids(t,1,g))
        CASE(SOURCE_TARGET_VELOCITY); membership(2:3)=compact_group
        CASE(SOURCE_TARGET_TEMPERATURE); membership(4)=compact_group
        CASE(SOURCE_TARGET_PORE_PRESSURE); membership(idx_pore)=compact_group
        END SELECT
      END DO
    END DO
    CALL partition%initialize(implicit_mask,membership,methods(1:compact_group))

    stage_count=COUNT(cadence(1:source_group_count)==SOURCE_CADENCE_STAGE)
    step_count=source_group_count-stage_count
    ALLOCATE(stages(stage_count),steps(step_count))
    j=0; k=0
    DO i=1,source_group_count
      g=ordered(i)
      group%id=source_group_id(g); group%legacy_id=0
      group%cadence=cadence(g); group%method=method(g)
      group%phase=phase(g); group%execution_order=source_group_order(g)
      group_terms=parsed_terms(1:parsed_count(g),g)
      group%term_ids=group_terms
      group%equation_ids=[INTEGER::]
      IF (phase(g)==SOURCE_PHASE_RK_IMPLICIT) THEN
        equations=[INTEGER::]
        targets=target_ids(1:target_count(1,g),1,g)
        DO t=1,SIZE(targets)
          SELECT CASE(targets(t))
          CASE(SOURCE_TARGET_VELOCITY); equations=[equations,2,3]
          CASE(SOURCE_TARGET_TEMPERATURE); equations=[equations,4]
          CASE(SOURCE_TARGET_PORE_PRESSURE); equations=[equations,idx_pore]
          END SELECT
        END DO
        group%equation_ids=equations
      ELSE
        group%equation_ids=[(t,t=1,n_eqns)]
      END IF
      IF (cadence(g)==SOURCE_CADENCE_STAGE) THEN
        j=j+1; stages(j)=group
      ELSE
        k=k+1; steps(k)=group
      END IF
    END DO
    CALL registry%initialize(n_eqns,MERGE(idx_pore,0,pore_pressure_flag), &
         stages,steps,validation_error)
    IF (LEN_TRIM(validation_error)>0) THEN
      CALL partition%finalize()
      error=TRIM(validation_error)
      RETURN
    END IF
    error=''
  END SUBROUTINE initialize_native_source_configuration

  !> \brief Map a native cadence name to its catalog identifier.
  !> \param[in] name Input cadence name.
  !> \return Cadence identifier, or zero when unsupported.
  PURE INTEGER FUNCTION native_cadence(name)
    CHARACTER(*), INTENT(IN) :: name
    SELECT CASE(TRIM(uppercase_source_text(name)))
    CASE('STAGE'); native_cadence=SOURCE_CADENCE_STAGE
    CASE('STEP'); native_cadence=SOURCE_CADENCE_STEP
    CASE DEFAULT; native_cadence=0
    END SELECT
  END FUNCTION native_cadence

  !> \brief Map a native solver-method name to its catalog identifier.
  !> \param[in] name Input method name.
  !> \return Method identifier, or zero when unsupported.
  PURE INTEGER FUNCTION native_method(name)
    CHARACTER(*), INTENT(IN) :: name
    SELECT CASE(TRIM(uppercase_source_text(name)))
    CASE('EXPLICIT'); native_method=SOURCE_METHOD_EXPLICIT
    CASE('NEWTON'); native_method=SOURCE_METHOD_NEWTON
    CASE('ANALYTIC'); native_method=SOURCE_METHOD_ANALYTIC
    CASE DEFAULT; native_method=0
    END SELECT
  END FUNCTION native_method

  !> \brief Map a native solver-target name to its catalog identifier.
  !> \param[in] name Input target name.
  !> \return Target identifier, or zero when unsupported.
  PURE INTEGER FUNCTION native_target(name)
    CHARACTER(*), INTENT(IN) :: name
    SELECT CASE(TRIM(uppercase_source_text(name)))
    CASE('VELOCITY'); native_target=SOURCE_TARGET_VELOCITY
    CASE('TEMPERATURE'); native_target=SOURCE_TARGET_TEMPERATURE
    CASE('PORE_PRESSURE'); native_target=SOURCE_TARGET_PORE_PRESSURE
    CASE DEFAULT; native_target=0
    END SELECT
  END FUNCTION native_target

  !> \brief Resolve a native source-term name to its schedulable catalog identifier.
  !> \param[in] name Input physical source-term name.
  !> \return Schedulable term identifier, or zero when unsupported.
  PURE INTEGER FUNCTION native_term_id(name)
    CHARACTER(*), INTENT(IN) :: name
    TYPE(source_term_type) :: term
    INTEGER :: id
    native_term_id=0
    DO id=1,SOURCE_TERM_COUNT
      term=source_term_descriptor(id)
      IF (TRIM(uppercase_source_text(name))==TRIM(term%name)) THEN
        IF (term%schedulable) native_term_id=id
        RETURN
      END IF
    END DO
  END FUNCTION native_term_id
END MODULE source_configuration_2d
