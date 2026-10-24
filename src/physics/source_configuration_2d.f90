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
  USE parameters_2d, ONLY : wp, n_eqns, rheology_flag, rheology_model, &
        thermal_cooling_flag, &
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
    CHARACTER(LEN=64), ALLOCATABLE :: active_names(:)
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
    ALLOCATE(active_names(0))
    IF (rheology_flag) THEN
      IF (LEN_TRIM(canonical_rheology_source_name(rheology_model))>0) &
           active_names=[CHARACTER(LEN=64) :: active_names,canonical_rheology_source_name(rheology_model)]
      IF (LEN_TRIM(canonical_rheology_source_name(rheology_model))==0 .AND. &
           TRIM(rheology_model)/='none') THEN
        error='D17: unsupported configured rheology'; RETURN
      END IF
    END IF
    IF (partition%implicit(4)) active_names=[CHARACTER(LEN=64) :: active_names,'thermal_cooling']
    pore=0
    IF (pore_pressure_flag) THEN
      pore=idx_pore
      active_names=[CHARACTER(LEN=64) :: active_names,'pore_pressure_dissipation']
    END IF
    IF (curvature_term_flag) active_names=[CHARACTER(LEN=64) :: active_names,'geometry_curvature']
    IF (bottom_radial_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'bottom_circular_injection']
    IF (bottom_fissural_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'bottom_fissural_injection']
    IF (radial_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'radial_injection']
    IF (lateral_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'boundary_injection']
    IF (erosion_coeff>0.0_wp) active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_erosion']
    IF (settling_flag) THEN
      active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_deposition']
    END IF
    IF (liquid_flag) active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_carrier_liquid_loss']
    IF (pore_pressure_flag .AND. gas_loss_flag) &
         active_names=[CHARACTER(LEN=64) :: active_names, &
         'mass_exchange_carrier_darcy_loss']
    IF (entrainment_flag) active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_entrainment']
    IF (gas_flag .AND. liquid_vaporization_flag) &
         active_names=[CHARACTER(LEN=64) :: active_names, &
         'mass_exchange_vaporization']
    DO i=1,3
      method=imex_group_method(i)
      CALL resolve_imex_group_method(method,methods(i))
    END DO
    CALL registry%initialize_legacy(partition, &
         [imex_velocity_group,imex_energy_group,imex_pore_pressure_group], &
         methods,pore,active_names,error)
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
    TYPE(source_term_type) :: component
    INTEGER, ALLOCATABLE :: ordered(:),targets(:),equations(:)
    CHARACTER(LEN=64), ALLOCATABLE :: parsed_names(:,:)
    CHARACTER(LEN=64), ALLOCATABLE :: active_names(:)
    INTEGER :: cadence(SOURCE_INPUT_MAX_GROUPS),method(SOURCE_INPUT_MAX_GROUPS)
    INTEGER :: phase(SOURCE_INPUT_MAX_GROUPS),solver_count(SOURCE_INPUT_MAX_GROUPS)
    INTEGER :: parsed_count(SOURCE_INPUT_MAX_GROUPS)
    INTEGER :: target_count(1,SOURCE_INPUT_MAX_GROUPS)
    INTEGER :: target_ids(SOURCE_INPUT_MAX_TARGETS,1, &
         SOURCE_INPUT_MAX_GROUPS)
    INTEGER :: target_group(3),methods(SOURCE_INPUT_MAX_GROUPS),membership(n_eqns)
    LOGICAL :: implicit_mask(n_eqns),expected_targets(3)
    CHARACTER(LEN=64) :: seen_names(SOURCE_INPUT_MAX_GROUPS*SOURCE_INPUT_MAX_TERMS)
    INTEGER :: n_seen
    INTEGER :: i,j,k,g,s,t,slot,compact_group,stage_count,step_count,target
    CHARACTER(LEN=64) :: source_name
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
    IF (n_eqns<3 .OR. (thermal_cooling_flag .AND. n_eqns<4)) THEN
      error='D17: native source configuration requires a valid equation layout'; RETURN
    END IF
    IF (pore_pressure_flag .AND. (idx_pore<1 .OR. idx_pore>n_eqns .OR. &
         idx_pore==2 .OR. idx_pore==3 .OR. &
         (thermal_cooling_flag .AND. idx_pore==4))) THEN
      error='D17: pore-pressure equation index is outside the configured layout'; RETURN
    END IF

    ALLOCATE(parsed_names(SOURCE_INPUT_MAX_TERMS,source_group_count))
    parsed_names=''; parsed_count=0; cadence=0; method=0; phase=0; solver_count=0
    target_count=0; target_ids=0
    seen_names=''; n_seen=0
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
        ! A named rheological law is one schedulable physical source.
        source_name=TRIM(lowercase_source_text(source_group_terms(i,g)))
        IF (rheology_flag .AND. &
             LEN_TRIM(canonical_rheology_source_name(rheology_model))>0 .AND. &
             source_name==TRIM(canonical_rheology_source_name(rheology_model))) THEN
          source_name=canonical_rheology_source_name(rheology_model)
        ELSE
          component=source_term_descriptor_by_name(source_name)
          error='D17: unknown, dependent or unschedulable source term name'
          IF (LEN_TRIM(component%name)==0) RETURN
          source_name=component%name
        END IF
        error='D17: source group exceeds supported catalog capacity'
        IF (parsed_count(g)>=SOURCE_INPUT_MAX_TERMS) RETURN
        error='D17: a source term may be registered only once globally'
        IF (ANY(seen_names(1:n_seen)==source_name)) RETURN
        parsed_count(g)=parsed_count(g)+1
        n_seen=n_seen+1
        seen_names(n_seen)=source_name
        parsed_names(parsed_count(g),g)=source_name
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

    ALLOCATE(active_names(0))
    IF (rheology_flag) THEN
      IF (LEN_TRIM(canonical_rheology_source_name(rheology_model))>0) &
           active_names=[CHARACTER(LEN=64) :: active_names,canonical_rheology_source_name(rheology_model)]
      IF (LEN_TRIM(canonical_rheology_source_name(rheology_model))==0 .AND. &
           TRIM(rheology_model)/='none') THEN
        error='D17: unsupported configured rheology'; RETURN
      END IF
    END IF
    IF (thermal_cooling_flag) active_names=[CHARACTER(LEN=64) :: active_names, &
         'thermal_cooling']
    IF (pore_pressure_flag) active_names=[CHARACTER(LEN=64) :: active_names,'pore_pressure_dissipation']
    IF (curvature_term_flag) active_names=[CHARACTER(LEN=64) :: active_names,'geometry_curvature']
    IF (bottom_radial_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'bottom_circular_injection']
    IF (bottom_fissural_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'bottom_fissural_injection']
    IF (radial_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'radial_injection']
    IF (lateral_source_flag) active_names=[CHARACTER(LEN=64) :: active_names,'boundary_injection']
    IF (erosion_coeff>0.0_wp) active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_erosion']
    IF (settling_flag) THEN
      active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_deposition']
    END IF
    IF (liquid_flag) active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_carrier_liquid_loss']
    IF (pore_pressure_flag .AND. gas_loss_flag) &
         active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_carrier_darcy_loss']
    IF (entrainment_flag) active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_entrainment']
    IF (gas_flag .AND. liquid_vaporization_flag) &
         active_names=[CHARACTER(LEN=64) :: active_names,'mass_exchange_vaporization']
    error='D17: native schedule must register every configured active term exactly once'
    IF (n_seen/=SIZE(active_names)) RETURN
    DO i=1,SIZE(active_names)
      IF (.NOT.ANY(seen_names(1:n_seen)==active_names(i))) RETURN
    END DO

    target_group=0
    expected_targets=[.TRUE.,thermal_cooling_flag,pore_pressure_flag]
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
      IF (ALLOCATED(group%source_names)) DEALLOCATE(group%source_names)
      ALLOCATE(group%source_names(parsed_count(g)))
      group%source_names=parsed_names(1:parsed_count(g),g)
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

  !> Normalize user-supplied physical-source names without accepting legacy aliases.
  PURE FUNCTION lowercase_source_text(value) RESULT(lower)
    CHARACTER(*), INTENT(IN) :: value
    CHARACTER(LEN(value)) :: lower
    INTEGER :: i,code
    lower=value
    DO i=1,LEN(value)
      code=IACHAR(value(i:i))
      IF (code>=IACHAR('A') .AND. code<=IACHAR('Z')) &
           lower(i:i)=ACHAR(code+IACHAR('a')-IACHAR('A'))
    END DO
  END FUNCTION lowercase_source_text

END MODULE source_configuration_2d
