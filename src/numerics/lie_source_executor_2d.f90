!********************************************************************************
!> \brief Execute ordered first-order Lie physical source groups.
!> \details Ordered groups advance each cell independently while geometry
!>          remains frozen. Physical rates are assembled by source_terms_2d;
!>          accepted bed changes are projected through geometry_2d.
!********************************************************************************

MODULE lie_source_executor_2d
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp, n_eqns, n_vars, n_solid, mass_exchange_rates_type, idx_totMassEqn,                   &
        idx_solidEqn_first, idx_solidEqn_last
    USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
    USE domain_2d, ONLY : domain_type
    USE source_registry_2d, ONLY : source_registry_type, source_group_type
    USE source_catalog_2d, ONLY : SOURCE_ID_GEOMETRY_CURVATURE, SOURCE_ID_EXCHANGE_DEPOSITION,                         &
        SOURCE_ID_EXCHANGE_LAVA_DEPOSITION, SOURCE_ID_EXCHANGE_VAPORIZATION
    USE state_fields, ONLY : state_type
    USE state_types, ONLY : physical_state_real, physical_state_complex, &
        initialize_physical_state_real, initialize_physical_state_complex
    USE mass_exchange_2d, ONLY : mass_exchange_scratch_type

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: execute_lie_source_groups,initialize_lie_source_executor,finalize_lie_source_executor,                   &
        report_lie_source_profile,lie_source_executor_type
    !> @brief Optional cumulative counters for Lie group profiling.
    TYPE :: group_profile_type
        INTEGER :: id=0,method=0,calls=0
        REAL(wp) :: setup=0.0_wp,evaluate=0.0_wp,commit=0.0_wp
        REAL(wp) :: elapsed=0.0_wp,physical=0.0_wp
        INTEGER :: evaluated_cells=0
    END TYPE group_profile_type
    !> @brief Stable Lie configuration and reusable per-simulation work arrays.
    !> @details Initialized once; cell-indexed buffers are sized for the full
    !>          grid. Detailed failure diagnostics are disabled by default.
    TYPE :: lie_source_executor_type
        TYPE(source_group_type), ALLOCATABLE :: groups(:)
        LOGICAL, ALLOCATABLE :: exchange_group(:),dry_noop_group(:)
        TYPE(source_group_type), ALLOCATABLE :: non_exchange_groups(:)
        TYPE(physical_state_real), ALLOCATABLE :: source_physical_scratch(:)
        ! Per-worker complex-step workspace persists across cells and timesteps.
        TYPE(physical_state_complex), ALLOCATABLE :: source_complex_scratch(:)
        TYPE(mass_exchange_scratch_type), ALLOCATABLE :: source_exchange_scratch(:)
        TYPE(mass_exchange_rates_type), ALLOCATABLE :: source_transfers(:)
        LOGICAL :: has_exchange_group=.FALSE.,has_bed_group=.FALSE.,has_newton_group=.FALSE.
        LOGICAL :: all_groups_dry_noop=.FALSE.
        LOGICAL :: diagnostics_enabled=.FALSE.,newton_witness=.FALSE.,initialized=.FALSE.
        LOGICAL :: reject_active=.FALSE.
        INTEGER :: bed_group_id=-1,reject_group=-1
        INTEGER, ALLOCATABLE :: failure_code(:),failure_group(:),failure_detail(:),failure_iterations(:)
        REAL(wp), ALLOCATABLE :: failure_mass(:),failure_components(:)
        REAL(wp), ALLOCATABLE :: failure_reference(:,:),failure_rate(:,:),failure_candidate(:,:)
        REAL(wp), ALLOCATABLE :: injection_mass(:),injection_solids(:),pre_buoyancy_mass(:)
        REAL(wp), ALLOCATABLE :: pre_buoyancy_solids(:),buoyancy_reduced_gravity(:)
        INTEGER, ALLOCATABLE :: newton_iterations(:,:)
        REAL(wp), ALLOCATABLE :: bed_rate_cell(:,:),bed_rate_vertex(:,:),bed_candidate(:,:)
    END TYPE lie_source_executor_type
    TYPE(group_profile_type), ALLOCATABLE, SAVE :: group_profile(:)
    LOGICAL, SAVE :: profiling_initialized=.FALSE.,profiling_enabled=.FALSE.
    REAL(wp), SAVE :: profile_cell_sweep=0.0_wp,profile_final_commit=0.0_wp
    REAL(wp), SAVE :: profile_thread_cell=0.0_wp

CONTAINS

    !> @brief Prepare group metadata and allocate all executor work arrays.
    !> @details Called once after source registration and before the first step.
    !> @param[in,out] executor Persistent Lie configuration and work arrays.
    !> @param[in] source_registry Validated ordered physical-source groups.
    SUBROUTINE initialize_lie_source_executor(executor,source_registry)
        USE source_catalog_2d, ONLY : SOURCE_METHOD_ANALYTIC,SOURCE_METHOD_EXPLICIT,SOURCE_METHOD_NEWTON
        USE parameters_2d, ONLY : topo_change_flag
        USE constitutive_parameters_2d, ONLY : N_inh,f_inhibit_mode
        USE geometry_2d, ONLY : comp_interfaces_x,comp_interfaces_y
        TYPE(lie_source_executor_type), INTENT(INOUT) :: executor
        TYPE(source_registry_type), INTENT(IN) :: source_registry
        INTEGER :: group_index,env_status,parse_status,ncells,source_id
        CHARACTER(LEN=24) :: env_value

        IF (.NOT.source_registry%is_initialized()) CALL fatal_error('D19: uninitialized Lie source registry')
        IF (executor%initialized) CALL fatal_error('Lie source executor already initialized')
        ! Cache the groups and classify their methods once per simulation.
        executor%groups=source_registry%groups()
        ALLOCATE(executor%exchange_group(SIZE(executor%groups)),                                                       &
            executor%dry_noop_group(SIZE(executor%groups)),                                                           &
            executor%non_exchange_groups(SIZE(executor%groups)))
        executor%exchange_group=.FALSE.
        executor%dry_noop_group=.FALSE.
        DO group_index=1,SIZE(executor%groups)
            IF (.NOT.ALLOCATED(executor%groups(group_index)%source_names))                                        &
                CALL fatal_error('D20: Lie source group has no allocated source list')
            IF (SIZE(executor%groups(group_index)%source_names)==0) CYCLE
            SELECT CASE(executor%groups(group_index)%method)
            CASE(SOURCE_METHOD_EXPLICIT,SOURCE_METHOD_ANALYTIC,SOURCE_METHOD_NEWTON)
                CONTINUE
            CASE DEFAULT
                CALL fatal_error('D20: invalid Lie source-group method before OpenMP evaluation')
            END SELECT
            IF (ANY(executor%groups(group_index)%source_names=='pore_pressure_dissipation')                            &
                .OR. ANY(executor%groups(group_index)%source_names=='pore_pressure_darcy')) THEN
                IF (TRIM(ADJUSTL(f_inhibit_mode))=='STATIC' .OR. TRIM(ADJUSTL(f_inhibit_mode))=='DYNAMIC') THEN
                    IF (N_inh<0) CALL fatal_error('D20: negative Darcy inhibition order')
                END IF
            END IF
            ! Resolved catalogue IDs 18:23 are all mass-exchange processes.
            executor%exchange_group(group_index)=                                                                      &
                IAND(executor%groups(group_index)%source_mask,                                                        &
                ISHFT(IBITS(NOT(0),0,SOURCE_ID_EXCHANGE_VAPORIZATION-SOURCE_ID_EXCHANGE_LAVA_DEPOSITION+1),           &
                SOURCE_ID_EXCHANGE_LAVA_DEPOSITION-1))/=0
            executor%non_exchange_groups(group_index)=executor%groups(group_index)
            IF (executor%exchange_group(group_index)) THEN
                DO source_id=SOURCE_ID_EXCHANGE_LAVA_DEPOSITION,SOURCE_ID_EXCHANGE_VAPORIZATION
                    executor%non_exchange_groups(group_index)%source_mask=                                            &
                        IBCLR(executor%non_exchange_groups(group_index)%source_mask,source_id-1)
                END DO
            END IF
            executor%dry_noop_group(group_index)=                                                                     &
                IAND(executor%groups(group_index)%source_mask,                                                        &
                NOT(IOR(IBSET(0,SOURCE_ID_GEOMETRY_CURVATURE-1),                                                       &
                IBSET(0,SOURCE_ID_EXCHANGE_DEPOSITION-1))))==0
            IF (executor%exchange_group(group_index)) THEN
                IF (executor%groups(group_index)%method/=SOURCE_METHOD_EXPLICIT)                                   &
                    CALL fatal_error('D19: mass exchange requires explicit group method')
                executor%has_exchange_group=.TRUE.
                IF (topo_change_flag) executor%bed_group_id=executor%groups(group_index)%id
            END IF
            IF (executor%groups(group_index)%method==SOURCE_METHOD_NEWTON) executor%has_newton_group=.TRUE.
        END DO
        executor%all_groups_dry_noop=ALL(executor%dry_noop_group .OR.                                   &
            [(SIZE(executor%groups(group_index)%source_names)==0,group_index=1,SIZE(executor%groups))])
        executor%has_bed_group=executor%has_exchange_group .AND. topo_change_flag
        ! Read diagnostic and D20 test settings only at startup.
        env_value=''
        CALL GET_ENVIRONMENT_VARIABLE('IMEX_LIE_DIAGNOSTICS',env_value)
        executor%diagnostics_enabled=TRIM(env_value)=='1'
        env_value=''
        CALL GET_ENVIRONMENT_VARIABLE('IMEX_NEWTON_WITNESS',env_value)
        executor%newton_witness=TRIM(env_value)=='1'
        env_value=''
        CALL GET_ENVIRONMENT_VARIABLE('IMEX_LIE_TEST_REJECT_AFTER_GROUP',env_value,STATUS=env_status)
        IF (env_status==0 .AND. LEN_TRIM(env_value)>0) THEN
            READ(env_value,*,IOSTAT=parse_status) executor%reject_group
            IF (parse_status/=0 .OR. executor%reject_group<1)                                                       &
                CALL fatal_error('Invalid IMEX_LIE_TEST_REJECT_AFTER_GROUP')
        END IF
        executor%reject_active=ANY(executor%groups%id==executor%reject_group)
        ! Allocate for the maximum worklist size; only enabled diagnostics
        ! receive their detailed per-cell arrays.
        ncells=comp_cells_x*comp_cells_y
        ALLOCATE(executor%failure_code(ncells))
        IF (executor%diagnostics_enabled) THEN
            ALLOCATE(executor%failure_group(ncells),executor%failure_detail(ncells))
            ALLOCATE(executor%failure_iterations(ncells))
            ALLOCATE(executor%failure_mass(ncells),executor%failure_components(ncells))
            ALLOCATE(executor%failure_reference(n_eqns,ncells),executor%failure_rate(n_eqns,ncells))
            ALLOCATE(executor%failure_candidate(n_eqns,ncells))
            ALLOCATE(executor%injection_mass(ncells),executor%injection_solids(ncells))
            ALLOCATE(executor%pre_buoyancy_mass(ncells),executor%pre_buoyancy_solids(ncells))
            ALLOCATE(executor%buoyancy_reduced_gravity(ncells))
        END IF
        IF (executor%newton_witness                                                                                    &
            .AND. executor%has_newton_group) ALLOCATE(executor%newton_iterations(ncells,SIZE(executor%groups)))
        IF (executor%has_bed_group) THEN
            ALLOCATE(executor%bed_rate_cell(comp_cells_x,comp_cells_y))
            ALLOCATE(executor%bed_rate_vertex(comp_interfaces_x,comp_interfaces_y))
            ALLOCATE(executor%bed_candidate(comp_interfaces_x,comp_interfaces_y))
        END IF
        IF (.NOT.profiling_initialized) THEN
            env_value=''
            CALL GET_ENVIRONMENT_VARIABLE('IMEX_PROFILE_LIE',env_value,STATUS=env_status)
            profiling_enabled=env_status==0 .AND. TRIM(env_value)=='1'
            profiling_initialized=.TRUE.
            IF (profiling_enabled) THEN
                ALLOCATE(group_profile(SIZE(executor%groups)))
                DO group_index=1,SIZE(executor%groups)
                    group_profile(group_index)%id=executor%groups(group_index)%id
                    group_profile(group_index)%method=executor%groups(group_index)%method
                END DO
            END IF
        END IF
        executor%initialized=.TRUE.
    END SUBROUTINE initialize_lie_source_executor

    !> @brief Release cached groups, work arrays and diagnostic storage.
    !> @param[in,out] executor Lie executor to reset.
    SUBROUTINE finalize_lie_source_executor(executor)
        TYPE(lie_source_executor_type), INTENT(INOUT) :: executor
        IF (ALLOCATED(executor%groups)) DEALLOCATE(executor%groups)
        IF (ALLOCATED(executor%exchange_group)) DEALLOCATE(executor%exchange_group)
        IF (ALLOCATED(executor%dry_noop_group)) DEALLOCATE(executor%dry_noop_group)
        IF (ALLOCATED(executor%non_exchange_groups)) DEALLOCATE(executor%non_exchange_groups)
        IF (ALLOCATED(executor%source_physical_scratch)) DEALLOCATE(executor%source_physical_scratch)
        IF (ALLOCATED(executor%source_complex_scratch)) DEALLOCATE(executor%source_complex_scratch)
        IF (ALLOCATED(executor%source_exchange_scratch)) DEALLOCATE(executor%source_exchange_scratch)
        IF (ALLOCATED(executor%source_transfers)) DEALLOCATE(executor%source_transfers)
        IF (ALLOCATED(executor%failure_code)) DEALLOCATE(executor%failure_code)
        IF (ALLOCATED(executor%failure_group)) DEALLOCATE(executor%failure_group)
        IF (ALLOCATED(executor%failure_detail)) DEALLOCATE(executor%failure_detail)
        IF (ALLOCATED(executor%failure_iterations)) DEALLOCATE(executor%failure_iterations)
        IF (ALLOCATED(executor%failure_mass)) DEALLOCATE(executor%failure_mass)
        IF (ALLOCATED(executor%failure_components)) DEALLOCATE(executor%failure_components)
        IF (ALLOCATED(executor%failure_reference)) DEALLOCATE(executor%failure_reference)
        IF (ALLOCATED(executor%failure_rate)) DEALLOCATE(executor%failure_rate)
        IF (ALLOCATED(executor%failure_candidate)) DEALLOCATE(executor%failure_candidate)
        IF (ALLOCATED(executor%injection_mass)) DEALLOCATE(executor%injection_mass)
        IF (ALLOCATED(executor%injection_solids)) DEALLOCATE(executor%injection_solids)
        IF (ALLOCATED(executor%pre_buoyancy_mass)) DEALLOCATE(executor%pre_buoyancy_mass)
        IF (ALLOCATED(executor%pre_buoyancy_solids)) DEALLOCATE(executor%pre_buoyancy_solids)
        IF (ALLOCATED(executor%buoyancy_reduced_gravity)) DEALLOCATE(executor%buoyancy_reduced_gravity)
        IF (ALLOCATED(executor%newton_iterations)) DEALLOCATE(executor%newton_iterations)
        IF (ALLOCATED(executor%bed_rate_cell)) DEALLOCATE(executor%bed_rate_cell)
        IF (ALLOCATED(executor%bed_rate_vertex)) DEALLOCATE(executor%bed_rate_vertex)
        IF (ALLOCATED(executor%bed_candidate)) DEALLOCATE(executor%bed_candidate)
        executor%has_exchange_group=.FALSE.
        executor%has_bed_group=.FALSE.
        executor%has_newton_group=.FALSE.
        executor%all_groups_dry_noop=.FALSE.
        executor%diagnostics_enabled=.FALSE.
        executor%newton_witness=.FALSE.
        executor%reject_active=.FALSE.
        executor%reject_group=-1
        executor%bed_group_id=-1
        executor%initialized=.FALSE.
    END SUBROUTINE finalize_lie_source_executor

    !> @brief Apply ordered physical-source groups to the trial state.
    !> @details Each OpenMP iteration processes one cell and visits the groups
    !>          sequentially. Group zero reads the pre-spatial snapshot while
    !>          incrementing the post-spatial trial state. Geometry remains
    !>          frozen until all cells succeed; failure triggers caller rollback.
    !> @param[in,out] executor Configured Lie groups and reusable work arrays.
    !> @param[in,out] state Owns substrate inventories and aggregate fields.
    !> @param[in,out] q Conservative cell states (component,x-cell,y-cell).
    !> @param[in,out] qp Primitive cell states (component,x-cell,y-cell).
    !> @param[in] time_step Positive Lie substep duration [s].
    !> @param[in] domain Active cell list and per-cell coordinate lookup.
    !> @param[in] stochastic_values Optional stochastic cell samples.
    !> @param[in] step_start_time Optional beginning time for time-dependent sources [s].
    !> @param[in] q_original Pre-spatial conservative snapshot (read-only).
    !> @param[in] qp_original Pre-spatial primitive snapshot (read-only).
    !> @param[out] step_failed True when candidate group state is rejected.
    !> @param[out] failure_message Optional diagnostic describing rejected group.
    SUBROUTINE execute_lie_source_groups(executor,state,q,qp,time_step,domain,stochastic_values,                       &
        step_start_time,q_original,qp_original,step_failed,failure_message)
        USE geometry_2d, ONLY : B_prime_x_geom,B_prime_y_geom,                                                         &
            B_second_xx_geom,B_second_xy_geom,B_second_yy_geom,grav_coeff,                                             &
            project_cell_field_to_vertices,assign_topography_vertices,B_vertex, comp_interfaces_x,comp_interfaces_y,   &
            B_zone,cell_source_fractions,cell_fissure_fractions,                                              &
            cell_arc_perim,cell_arc_n_x,cell_arc_n_y,dx,dy
        USE source_terms_2d, ONLY : assemble_lie_explicit_transaction, eval_source_analytic,source_cell_context_type
        USE source_catalog_2d, ONLY : SOURCE_METHOD_ANALYTIC,SOURCE_METHOD_EXPLICIT, SOURCE_METHOD_NEWTON
        USE nonlinear_solver_2d, ONLY : solve_backward_euler_group
        USE state_conversions, ONLY : conservative_to_primitive
        USE state_physical, ONLY : build_physical_state_from_qp
        USE state_fields, ONLY : initialize_dry_physical_cell,conservative_mass_admissible,                               &
            physical_cell_state_finite,finite_real_value
        USE parameters_2d, ONLY : topo_change_flag,erodible_deposit_flag,sea_mass_removal_flag,idx_h
        !$ USE omp_lib, ONLY : omp_get_max_threads,omp_get_thread_num
        TYPE(lie_source_executor_type), INTENT(INOUT) :: executor
        CLASS(state_type), INTENT(INOUT) :: state
        TYPE(source_cell_context_type) :: cell_context
        REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(INOUT) :: qp(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: time_step
        CLASS(domain_type), INTENT(IN) :: domain
        REAL(wp), INTENT(IN), OPTIONAL :: stochastic_values(:,:)
        REAL(wp), INTENT(IN), OPTIONAL :: step_start_time
        REAL(wp), INTENT(IN) :: q_original(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: qp_original(n_vars,comp_cells_x,comp_cells_y)
        LOGICAL, INTENT(OUT) :: step_failed
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: failure_message
        REAL(wp) :: qpentry(n_vars),q_candidate(n_vars),qp_candidate(n_vars)
        REAL(wp) :: q_local(n_vars),qp_local(n_vars)
        REAL(wp) :: deposit_local(n_solid),erosion_local(n_solid),erodible_local(n_solid),lava_local
        REAL(wp) :: deposit_candidate(n_solid),erosion_candidate(n_solid)
        REAL(wp) :: erodible_candidate(n_solid),lava_candidate
        REAL(wp) :: total_rate(n_eqns)
        REAL(wp) :: analytic_rate(n_eqns),analytic_input(n_vars),source_time,noise
        LOGICAL :: cell_has_exchange,ok,failed,bed_accepted,conversion_ok,cell_updated
        INTEGER :: group_index,cell,j,k,status,iterations,linear_info,first_failure
        INTEGER :: cell_error,cell_detail,cell_iterations,cell_group
        REAL(wp) :: cell_mass,cell_components
        REAL(wp) :: red_grav
        ! Profiling counters are thread-owned: no atomics inside the cell loop.
        REAL(wp), ALLOCATABLE :: thread_group(:,:),thread_physics(:,:),thread_cell(:)
        INTEGER, ALLOCATABLE :: thread_evaluations(:,:)
        INTEGER :: nthreads,tid
        INTEGER(KIND=8) :: cell_tick,group_tick,physical_tick,profile_tock,physical_tick_rate

        CHARACTER(LEN=96) :: context
        CHARACTER(LEN=256) :: failure_reason
        INTEGER :: tick,tock,clock_rate

        IF (.NOT.executor%initialized) CALL fatal_error('Lie source executor not initialized')
        ! Reset only per-step results and accumulators; never rebuild metadata.
        step_failed=.FALSE.
        failure_reason=''
        IF (PRESENT(failure_message)) failure_message=''
        source_time=0.0_wp
        IF (PRESENT(step_start_time)) source_time=step_start_time
        IF (time_step<=0.0_wp) CALL fatal_error('D19: nonpositive time interval')
        IF (profiling_enabled) THEN
            CALL SYSTEM_CLOCK(tick,clock_rate)
            DO group_index=1,SIZE(executor%groups)
                IF (SIZE(executor%groups(group_index)%source_names)==0) CYCLE
                group_profile(group_index)%calls=group_profile(group_index)%calls+1
            END DO
        END IF
        IF (executor%has_bed_group) executor%bed_rate_cell=0.0_wp
        IF (executor%diagnostics_enabled) THEN
            executor%injection_mass(1:domain%solve_cells)=0.0_wp
            executor%injection_solids(1:domain%solve_cells)=0.0_wp
            executor%pre_buoyancy_mass(1:domain%solve_cells)=0.0_wp
            executor%pre_buoyancy_solids(1:domain%solve_cells)=0.0_wp
            executor%buoyancy_reduced_gravity(1:domain%solve_cells)=0.0_wp
        END IF
        IF (ALLOCATED(executor%newton_iterations)) executor%newton_iterations(1:domain%solve_cells,:)=0
        IF (profiling_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            IF (SIZE(group_profile)>0)                                                                                 &
                group_profile(1)%setup=group_profile(1)%setup+REAL(tock-tick,wp)/REAL(clock_rate,wp)
            tick=tock
        END IF
        ! Allocate one reusable physical-state record per available worker.
        nthreads=1
        !$ nthreads=omp_get_max_threads()
        IF (ALLOCATED(executor%source_physical_scratch)) THEN
            IF (SIZE(executor%source_physical_scratch)<nthreads)                                                       &
                DEALLOCATE(executor%source_physical_scratch)
        END IF
        IF (.NOT.ALLOCATED(executor%source_physical_scratch))                                                         &
            ALLOCATE(executor%source_physical_scratch(nthreads))
        IF (executor%has_newton_group) THEN
            IF (ALLOCATED(executor%source_complex_scratch)) THEN
                IF (SIZE(executor%source_complex_scratch)<nthreads)                                                   &
                    DEALLOCATE(executor%source_complex_scratch)
            END IF
            IF (.NOT.ALLOCATED(executor%source_complex_scratch))                                                       &
                ALLOCATE(executor%source_complex_scratch(nthreads))
        END IF
        ! Set up composition storage once per OpenMP worker, outside the hot loops.
        DO tid=1,nthreads
            IF (.NOT.executor%source_physical_scratch(tid)%layout_valid) &
                CALL initialize_physical_state_real(executor%source_physical_scratch(tid))
            IF (executor%has_newton_group) THEN
                IF (.NOT.executor%source_complex_scratch(tid)%layout_valid) &
                    CALL initialize_physical_state_complex(executor%source_complex_scratch(tid))
            END IF
        END DO
        IF (ALLOCATED(executor%source_exchange_scratch)) THEN
            IF (SIZE(executor%source_exchange_scratch)<nthreads) DEALLOCATE(executor%source_exchange_scratch)
        END IF
        IF (.NOT.ALLOCATED(executor%source_exchange_scratch))                                                         &
            ALLOCATE(executor%source_exchange_scratch(nthreads))
        IF (ALLOCATED(executor%source_transfers)) THEN
            IF (SIZE(executor%source_transfers)<nthreads) DEALLOCATE(executor%source_transfers)
        END IF
        IF (.NOT.ALLOCATED(executor%source_transfers)) ALLOCATE(executor%source_transfers(nthreads))
        IF (profiling_enabled) THEN
            ALLOCATE(thread_group(SIZE(executor%groups),nthreads))
            ALLOCATE(thread_physics(SIZE(executor%groups),nthreads))
            ALLOCATE(thread_evaluations(SIZE(executor%groups),nthreads))
            ALLOCATE(thread_cell(nthreads))
            thread_group=0.0_wp
            thread_physics=0.0_wp
            thread_evaluations=0
            thread_cell=0.0_wp
        END IF
        ! Independent cells share a fixed geometry. Reduce the first failed
        ! worklist index without a separate serial scan.
        first_failure=HUGE(0)
        !$OMP PARALLEL DEFAULT(SHARED)                                                                              &
        !$OMP& PRIVATE(cell,j,k,group_index,qpentry,q_candidate,qp_candidate,                                          &
        !$OMP& q_local,qp_local,deposit_local,erosion_local,erodible_local,lava_local,                                 &
        !$OMP& deposit_candidate,erosion_candidate,erodible_candidate,lava_candidate,total_rate,                       &
        !$OMP& analytic_rate,analytic_input,noise,cell_context,                                           &
        !$OMP& cell_has_exchange,ok,failed,conversion_ok,status,iterations,linear_info,context,                       &
        !$OMP& cell_error,cell_updated,cell_detail,cell_iterations,cell_group,cell_mass,cell_components,              &
        !$OMP& red_grav,tid,                                                                                         &
        !$OMP& cell_tick,group_tick,                                                                                 &
        !$OMP& physical_tick,profile_tock,physical_tick_rate)                                   &
        !$OMP& REDUCTION(MIN:first_failure)
        tid=1
        !$ tid=omp_get_thread_num()+1
        !$OMP DO
        cell_loop: DO cell=1,domain%solve_cells
            IF (profiling_enabled) CALL SYSTEM_CLOCK(cell_tick,physical_tick_rate)
            cell_error=0
            IF (executor%diagnostics_enabled) THEN
                cell_detail=0
                cell_iterations=0
                cell_group=0
                cell_mass=0.0_wp
                cell_components=0.0_wp
            END IF
            ! Load the candidate state and optional substrate inventories.
            j=domain%j_cent(cell)
            k=domain%k_cent(cell)
            ! Group zero uses the pre-spatial snapshot: never evaluate it over absorbing cells.
            IF (sea_mass_removal_flag .AND. B_zone(j,k)==1) CYCLE cell_loop
            ! No configured group can create mass in a dry cell in this case.
            IF (executor%all_groups_dry_noop) THEN
                IF (q(idx_totMassEqn,j,k)<=0.0_wp .AND. q_original(idx_totMassEqn,j,k)<=0.0_wp) CYCLE cell_loop
            END IF
            q_local=q(:,j,k)
            qp_local=qp(:,j,k)
            IF (executor%has_exchange_group) THEN
                deposit_local=state%deposit(j,k,1:n_solid)
                erosion_local=state%erosion(j,k,1:n_solid)
                erodible_local=state%erodible(1:n_solid,j,k)
                lava_local=state%lava_deposit(j,k)
            END IF
            ! Capture cell-local geometry and stochastic forcing for all groups.
            cell_context%Bprime_x=B_prime_x_geom(j,k)
            cell_context%Bprime_y=B_prime_y_geom(j,k)
            cell_context%Bsecond_xx=B_second_xx_geom(j,k)
            cell_context%Bsecond_xy=B_second_xy_geom(j,k)
            cell_context%Bsecond_yy=B_second_yy_geom(j,k)
            cell_context%grav_coeff=grav_coeff(j,k)
            noise=0.0_wp
            IF (PRESENT(stochastic_values)) noise=stochastic_values(j,k)
            cell_context%Zij=noise

            ! The Lie groups run in order within each OpenMP cell iteration.
            group_loop: DO group_index=1,SIZE(executor%groups)
                IF (SIZE(executor%groups(group_index)%source_names)==0) CYCLE group_loop
                ! Group zero reads the pre-spatial primitive state; subsequent
                ! groups read the previous group's accepted local state.
                IF (executor%groups(group_index)%id==0) THEN
                    qpentry=qp_original(:,j,k)
                ELSE
                    qpentry=qp_local
                END IF
                ! Curvature and deposition cannot create flow in a dry input cell.
                ! Do not skip groups containing injection or other dry-cell sources.
                IF (executor%dry_noop_group(group_index)) THEN
                    IF (qpentry(idx_h)<=0.0_wp) CYCLE group_loop
                END IF
                IF (profiling_enabled) THEN
                    CALL SYSTEM_CLOCK(group_tick)
                    thread_evaluations(group_index,tid)=thread_evaluations(group_index,tid)+1
                END IF

                context='Lie source group'
                IF (executor%exchange_group(group_index)) THEN
                    IF (ANY(erodible_local<0.0_wp) .OR. ANY(deposit_local<0.0_wp) .OR.                             &
                        ANY(erosion_local<0.0_wp) .OR. lava_local<0.0_wp) THEN
                        cell_error=3
                        EXIT group_loop
                    END IF
                END IF

                SELECT CASE(executor%groups(group_index)%method)
                CASE(SOURCE_METHOD_EXPLICIT)
                    ! The current recipient provides strict conservation bounds.
                    ! The pre-spatial inventory separately caps withdrawals by
                    ! mass present before any source group (legacy Lie semantics).
                    ! Newly injected solids cannot settle within the same step.
                    IF (profiling_enabled) CALL SYSTEM_CLOCK(physical_tick)
                    CALL assemble_lie_explicit_transaction(qpentry, q_local(1:n_eqns), time_step, source_time,           &
                        executor%groups(group_index), erodible_local, total_rate,                                      &
                        executor%source_transfers(tid),cell_has_exchange,ok,                                          &
                        noise, B_zone(j,k), cell_context, cell_source_fractions(j,k),                              &
                        cell_fissure_fractions(j,k,:),cell_arc_perim(j,k),cell_arc_n_x(j,k),                         &
                        cell_arc_n_y(j,k),dx*dy,                                                                     &
                        executor%non_exchange_groups(group_index),executor%exchange_group(group_index),                &
                        executor%source_physical_scratch(tid),executor%source_exchange_scratch(tid),                &
                        q_original(1:n_eqns,j,k))
                    IF (profiling_enabled) THEN
                        CALL SYSTEM_CLOCK(profile_tock)
                        thread_physics(group_index,tid)=thread_physics(group_index,tid)+                            &
                            REAL(profile_tock-physical_tick,wp)/REAL(physical_tick_rate,wp)
                    END IF

                    IF (cell_has_exchange .NEQV. executor%exchange_group(group_index)) THEN
                        cell_error=4
                        EXIT group_loop
                    END IF

                    IF (.NOT.ok) THEN
                        cell_error=5
                        EXIT group_loop
                    END IF
                    ! Stage the entire local fluid/substrate transfer before committing it.
                    q_candidate=q_local
                    q_candidate(1:n_eqns)=q_candidate(1:n_eqns)+time_step*total_rate

                    IF (.NOT.conservative_mass_admissible(q_candidate(1:n_eqns))) THEN
                        cell_error=6
                        IF (executor%diagnostics_enabled) THEN
                            cell_group=executor%groups(group_index)%id
                            cell_mass=q_candidate(idx_totMassEqn)
                            cell_components=SUM(q_candidate(idx_solidEqn_first:idx_solidEqn_last))
                            executor%failure_reference(:,cell)=q_local(1:n_eqns)
                            executor%failure_rate(:,cell)=total_rate
                            executor%failure_candidate(:,cell)=q_candidate(1:n_eqns)
                        END IF
                        EXIT group_loop
                    END IF

                    IF (executor%exchange_group(group_index)) THEN
                        deposit_candidate=deposit_local+time_step*executor%source_transfers(tid)%solid_deposition_rate
                        erosion_candidate=erosion_local+time_step*executor%source_transfers(tid)%solid_erosion_rate
                        erodible_candidate=erodible_local-time_step*executor%source_transfers(tid)%solid_erosion_rate
                        IF (erodible_deposit_flag) erodible_candidate=                                                 &
                            erodible_candidate+time_step*executor%source_transfers(tid)%solid_deposition_rate
                        lava_candidate=lava_local+time_step*executor%source_transfers(tid)%lava_deposit_rate
                        IF (ANY(erodible_candidate<0.0_wp) .OR. ANY(deposit_candidate<0.0_wp) .OR.                 &
                            ANY(erosion_candidate<0.0_wp) .OR. lava_candidate<0.0_wp) THEN
                            cell_error=7
                            EXIT group_loop
                        END IF
                    END IF

                    ! No-op conservative groups leave the synchronized primitive
                    ! representation valid. Exchange inventories may still change.
                    IF (ANY(q_candidate/=q_local)) THEN
                        IF (q_candidate(idx_totMassEqn)>0.0_wp) THEN
                            CALL conservative_to_primitive(q_candidate,qp_candidate,accepted=conversion_ok,             &
                                numerical_reconstruction=.TRUE.)
                            IF (.NOT.conversion_ok) THEN
                                cell_error=8
                                EXIT group_loop
                            END IF
                        ELSE
                            CALL initialize_dry_physical_cell(qp_candidate)
                        END IF
                        qp_local=qp_candidate
                    END IF
                    ! Publish to cell-local state; the trial array is written once after all groups.
                    q_local=q_candidate
                    IF (executor%diagnostics_enabled .AND. executor%groups(group_index)%id==0) THEN
                        IF (cell_arc_perim(j,k)>0.0_wp) THEN
                            executor%injection_mass(cell)=q_local(idx_totMassEqn)
                            executor%injection_solids(cell)=SUM(q_local(idx_solidEqn_first:idx_solidEqn_last))
                        END IF
                    END IF

                    IF (executor%exchange_group(group_index)) THEN
                        deposit_local=deposit_candidate
                        erosion_local=erosion_candidate
                        erodible_local=erodible_candidate
                        lava_local=lava_candidate
                        IF (topo_change_flag) executor%bed_rate_cell(j,k)=executor%bed_rate_cell(j,k)+                &
                            executor%source_transfers(tid)%dB_dt
                    END IF
                CASE(SOURCE_METHOD_ANALYTIC)
                    ! Group 0 evaluates on the immutable pre-spatial state;
                    ! its increment is applied to the post-spatial group candidate.
                    analytic_input=q_local
                    IF (executor%groups(group_index)%id==0) analytic_input=q_original(:,j,k)
                    IF (analytic_input(idx_totMassEqn)<=0.0_wp) CYCLE
                    IF (profiling_enabled) CALL SYSTEM_CLOCK(physical_tick)
                    CALL eval_source_analytic(executor%groups(group_index), analytic_input, cell_context,          &
                        time_step, analytic_rate, status, executor%source_physical_scratch(tid))
                    IF (profiling_enabled) THEN
                        CALL SYSTEM_CLOCK(profile_tock)
                        thread_physics(group_index,tid)=thread_physics(group_index,tid)+                            &
                            REAL(profile_tock-physical_tick,wp)/REAL(physical_tick_rate,wp)
                    END IF

                    IF (status/=0) THEN
                        cell_error=9
                        IF (executor%diagnostics_enabled) cell_detail=status
                        EXIT group_loop
                    END IF

                    q_candidate=q_local
                    q_candidate(1:n_eqns)=q_candidate(1:n_eqns)+time_step*analytic_rate
                CASE(SOURCE_METHOD_NEWTON)
                    ! Solve the coupled implicit update on a cell-local candidate.
                    IF (q(idx_totMassEqn,j,k)<=0.0_wp) CYCLE
                    q_candidate=q_local
                    IF (profiling_enabled) CALL SYSTEM_CLOCK(physical_tick)
                    CALL solve_backward_euler_group(executor%groups(group_index), q_candidate, time_step,              &
                        cell_context%Bprime_x, cell_context%Bprime_y, cell_context%grav_coeff, noise,      &
                        iterations, ok, linear_info, failed, executor%source_complex_scratch(tid),                    &
                        executor%source_physical_scratch(tid))
                    IF (profiling_enabled) THEN
                        CALL SYSTEM_CLOCK(profile_tock)
                        thread_physics(group_index,tid)=thread_physics(group_index,tid)+                            &
                            REAL(profile_tock-physical_tick,wp)/REAL(physical_tick_rate,wp)
                    END IF
                    IF (executor%newton_witness                                                                        &
                        .AND. ALLOCATED(executor%newton_iterations))                                                   &
                        executor%newton_iterations(cell,group_index)=iterations

                    IF (.NOT.ok .OR. failed .OR. linear_info/=0) THEN
                        cell_error=11
                        IF (executor%diagnostics_enabled) THEN
                            cell_detail=linear_info
                            cell_iterations=iterations
                        END IF
                        EXIT group_loop
                    END IF
                CASE DEFAULT
                    cell_error=12
                    EXIT group_loop

                END SELECT
                ! Explicit candidates were already validated and converted above.
                IF (executor%groups(group_index)%method==SOURCE_METHOD_EXPLICIT) THEN
                    IF (profiling_enabled) THEN
                        CALL SYSTEM_CLOCK(profile_tock)
                        thread_group(group_index,tid)=thread_group(group_index,tid)+                              &
                            REAL(profile_tock-group_tick,wp)/REAL(physical_tick_rate,wp)
                    END IF
                    CYCLE group_loop
                END IF
                ! Only check admissibility: never clip the accepted source update.
                IF (.NOT.conservative_mass_admissible(q_candidate(1:n_eqns))) THEN
                    cell_error=13
                    EXIT group_loop
                END IF

                ! Analytic and Newton candidates only invalidate primitive values
                ! when the conservative vector actually changes.
                IF (ANY(q_candidate/=q_local)) THEN
                    IF (q_candidate(idx_totMassEqn)>0.0_wp) THEN
                        CALL conservative_to_primitive(q_candidate,qp_candidate,accepted=conversion_ok,                &
                            numerical_reconstruction=.TRUE.)
                        IF (.NOT.conversion_ok) THEN
                            cell_error=14
                            EXIT group_loop
                        END IF
                    ELSE
                        CALL initialize_dry_physical_cell(qp_candidate)
                    END IF
                    qp_local=qp_candidate
                END IF

                q_local=q_candidate
                IF (profiling_enabled) THEN
                    CALL SYSTEM_CLOCK(profile_tock)
                    thread_group(group_index,tid)=thread_group(group_index,tid)+                                  &
                        REAL(profile_tock-group_tick,wp)/REAL(physical_tick_rate,wp)
                END IF

            END DO group_loop

            ! Restore the v2 buoyancy removal after all mass-exchange groups.
            ! This is a physical sink, not a numerical limiter on the source rates.
            IF (cell_error==0 .AND. q_local(idx_totMassEqn)>0.0_wp) THEN
                CALL build_physical_state_from_qp(qp_local,executor%source_physical_scratch(tid),                 &
                    accepted=conversion_ok)
                IF (.NOT.conversion_ok) THEN
                    cell_error=15
                ELSE
                    red_grav=executor%source_physical_scratch(tid)%red_grav
                END IF
                IF (executor%diagnostics_enabled .AND. cell_arc_perim(j,k)>0.0_wp) THEN
                    executor%pre_buoyancy_mass(cell)=q_local(idx_totMassEqn)
                    executor%pre_buoyancy_solids(cell)=SUM(q_local(idx_solidEqn_first:idx_solidEqn_last))
                    executor%buoyancy_reduced_gravity(cell)=red_grav
                END IF
                IF (cell_error==0 .AND. red_grav<=0.0_wp) THEN
                    q_local=0.0_wp
                    CALL initialize_dry_physical_cell(qp_local)
                END IF
            END IF
            ! Check finiteness only for cells actually changed by the complete
            ! spatial + Lie update. Test once, after all groups and buoyancy,
            ! before publishing the local fluid/substrate candidate.
            IF (cell_error==0) THEN
                cell_updated=ANY(q_local/=q_original(:,j,k))
                IF (.NOT.cell_updated) cell_updated=ANY(qp_local/=qp_original(:,j,k))
                IF (executor%has_exchange_group .AND. .NOT.cell_updated) THEN
                    cell_updated=ANY(deposit_local/=state%deposit(j,k,1:n_solid)) .OR.                              &
                        ANY(erosion_local/=state%erosion(j,k,1:n_solid)) .OR.                                        &
                        ANY(erodible_local/=state%erodible(1:n_solid,j,k)) .OR.                                      &
                        lava_local/=state%lava_deposit(j,k)
                END IF
                IF (executor%has_bed_group .AND. .NOT.cell_updated)                                                  &
                    cell_updated=executor%bed_rate_cell(j,k)/=0.0_wp
                IF (cell_updated) THEN
                    IF (.NOT.ALL(finite_real_value(q_local)) .OR. .NOT.physical_cell_state_finite(qp_local))         &
                        cell_error=2
                    IF (executor%has_exchange_group .AND. cell_error==0) THEN
                        IF (.NOT.ALL(finite_real_value(deposit_local)) .OR.                                           &
                            .NOT.ALL(finite_real_value(erosion_local)) .OR.                                           &
                            .NOT.ALL(finite_real_value(erodible_local)) .OR.                                          &
                            .NOT.finite_real_value(lava_local)) cell_error=3
                    END IF
                END IF
            END IF
            ! Publish this cell only after the final numerical check.
            IF (cell_error==0) THEN
                q(:,j,k)=q_local
                qp(:,j,k)=qp_local
                IF (executor%has_exchange_group) THEN
                    state%deposit(j,k,1:n_solid)=deposit_local
                    state%erosion(j,k,1:n_solid)=erosion_local
                    state%erodible(1:n_solid,j,k)=erodible_local
                    state%lava_deposit(j,k)=lava_local
                END IF
            END IF
            IF (cell_error/=0) THEN
                executor%failure_code(cell)=cell_error
                IF (executor%diagnostics_enabled) THEN
                    executor%failure_group(cell)=cell_group
                    executor%failure_detail(cell)=cell_detail
                    executor%failure_iterations(cell)=cell_iterations
                    executor%failure_mass(cell)=cell_mass
                    executor%failure_components(cell)=cell_components
                END IF
                first_failure=MIN(first_failure,cell)
            END IF
            IF (profiling_enabled) THEN
                CALL SYSTEM_CLOCK(profile_tock)
                thread_cell(tid)=thread_cell(tid)+REAL(profile_tock-cell_tick,wp)/REAL(physical_tick_rate,wp)
            END IF
        END DO cell_loop
        !$OMP END DO
        !$OMP END PARALLEL
        IF (profiling_enabled) THEN
            DO group_index=1,SIZE(executor%groups)
                group_profile(group_index)%elapsed=group_profile(group_index)%elapsed+SUM(thread_group(group_index,:))
                group_profile(group_index)%physical=group_profile(group_index)%physical+                             &
                    SUM(thread_physics(group_index,:))
                group_profile(group_index)%evaluated_cells=group_profile(group_index)%evaluated_cells+                &
                    SUM(thread_evaluations(group_index,:))
            END DO
            profile_thread_cell=profile_thread_cell+SUM(thread_cell)
            DEALLOCATE(thread_group,thread_physics,thread_evaluations,thread_cell)
        END IF
        IF (profiling_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            ! With cell-major traversal this is a *whole-sweep* wall-clock time,
            ! not an individual group evaluation. Report it separately.
            profile_cell_sweep=profile_cell_sweep+REAL(tock-tick,wp)/REAL(clock_rate,wp)
            tick=tock
        END IF
        ! One witness per timestep, gathered after the OpenMP cell sweep.
        ! Compare injection, deposition and final buoyancy deletion.
        IF (executor%diagnostics_enabled) THEN
            WRITE(*,'(A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
                'D19_SOURCE_AUDIT t=',source_time, &
                ' injection_mass=',SUM(executor%injection_mass(1:domain%solve_cells)), &
                ' injection_solids=',SUM(executor%injection_solids(1:domain%solve_cells)), &
                ' after_exchange_mass=',SUM(executor%pre_buoyancy_mass(1:domain%solve_cells))
            WRITE(*,'(A,ES14.6,A,I0,A,I0)') &
                'D19_SOURCE_AUDIT after_exchange_solids=', &
                SUM(executor%pre_buoyancy_solids(1:domain%solve_cells)), &
                ' injection_cells=',COUNT(executor%injection_mass(1:domain%solve_cells)>0.0_wp), &
                ' buoyancy_removed_cells=',COUNT(executor%pre_buoyancy_mass(1:domain%solve_cells)>0.0_wp .AND. &
                executor%buoyancy_reduced_gravity(1:domain%solve_cells)<=0.0_wp)
        END IF
        ! Build a human-readable failure message only when a cell was rejected.
        IF (first_failure/=HUGE(0)) THEN
            step_failed=.TRUE.
            SELECT CASE(executor%failure_code(first_failure))
            CASE(1)
                failure_reason='D19: non-finite source-group geometry'
            CASE(2)
                failure_reason='D19: non-finite source-group physical input'
            CASE(3)
                failure_reason='D19: invalid substrate inventory before mass exchange'
            CASE(4)
                failure_reason='D20: inconsistent cell exchange group metadata'
            CASE(5)
                failure_reason='D19: physical source transaction rejected'
            CASE(7)
                failure_reason='D19: invalid candidate substrate inventory'
            CASE(8)
                failure_reason='D19: explicit candidate conversion failed'
            CASE(10)
                failure_reason='D19: analytic source returned non-finite rate'
            CASE(12)
                failure_reason='D19: invalid group method'
            CASE(13)
                failure_reason='D19: inadmissible conservative source group state'
            CASE(14)
                failure_reason='D19: source group conversion failed'
            CASE(6)
                IF (executor%diagnostics_enabled) THEN
                    j=domain%j_cent(first_failure)
                    k=domain%k_cent(first_failure)
                    WRITE(failure_reason,'(A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4)')                               &
                        'D19: inadmissible source candidate group=',executor%failure_group(first_failure),           &
                        ' cell=',j,',',k,' mass=',executor%failure_mass(first_failure),                                &
                        ' components=',executor%failure_components(first_failure)
                    WRITE(*,'(A,2(I0,1X),A,I0)') 'D19_BUDGET cell=',j,k, &
                        'group=',executor%failure_group(first_failure)
                    WRITE(*,'(A,ES24.16)') 'D19_BUDGET dt=',time_step
                    WRITE(*,'(A,*(ES24.16,1X))') 'D19_BUDGET recipient=', &
                        executor%failure_reference(:,first_failure)
                    WRITE(*,'(A,*(ES24.16,1X))') 'D19_BUDGET scaled_rate=', &
                        executor%failure_rate(:,first_failure)
                    WRITE(*,'(A,*(ES24.16,1X))') 'D19_BUDGET candidate=', &
                        executor%failure_candidate(:,first_failure)
                    WRITE(*,'(A,L1)') 'D19_BUDGET reconstruction_match=', &
                        ALL(executor%failure_candidate(:,first_failure)== &
                        executor%failure_reference(:,first_failure)+ &
                        time_step*executor%failure_rate(:,first_failure))
                ELSE
                    failure_reason='D19: inadmissible source candidate'
                END IF
            CASE(9)
                failure_reason='D19: analytic source group failed'
                IF (executor%diagnostics_enabled) WRITE(failure_reason,'(A,I0)')                            &
                    'D19: analytic group status ',executor%failure_detail(first_failure)
            CASE(11)
                failure_reason='D19: Newton source group failed'
                IF (executor%diagnostics_enabled) WRITE(failure_reason,'(A,I0,A,I0)')                       &
                    'D19: Newton failed, linear_info=',executor%failure_detail(first_failure),                      &
                    ', iterations=',executor%failure_iterations(first_failure)
            CASE DEFAULT
                failure_reason='D19: unknown cell failure'
            END SELECT
            IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
            RETURN
        END IF
        ! The optional Newton witness is reported after the parallel sweep.
        IF (ALLOCATED(executor%newton_iterations)) THEN
            DO group_index=1,SIZE(executor%groups)
                IF (executor%groups(group_index)%method/=SOURCE_METHOD_NEWTON) CYCLE
                WRITE(*,'(A,I0,A,I0,A,I0)') 'D20_NEWTON_WITNESS group=',                                               &
                    executor%groups(group_index)%id,' solved_cells=',                                                  &
                    COUNT(executor%newton_iterations(1:domain%solve_cells,group_index)>0),                             &
                    ' total_iterations=',SUM(executor%newton_iterations(1:domain%solve_cells,group_index))
            END DO
        END IF
        ! Geometry is the same for every source group in the step. Project
        ! all mass-exchange bed rates only once, after the cell sweep.
        IF (executor%has_bed_group) THEN
            CALL project_cell_field_to_vertices(executor%bed_rate_cell,executor%bed_rate_vertex)
            executor%bed_candidate=B_vertex+time_step*executor%bed_rate_vertex
            IF (executor%diagnostics_enabled .OR. executor%reject_active) THEN
                WRITE(*,'(A,I0,A,ES16.8)') 'D20_BED_WITNESS group=',executor%bed_group_id,                             &
                    ' max_abs_delta=',MAXVAL(ABS(executor%bed_candidate-B_vertex))
            END IF
            CALL assign_topography_vertices(executor%bed_candidate,bed_accepted)
            IF (.NOT.bed_accepted) THEN
                CALL reject_cell('D19: invalid refreshed geometry')
                IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
                RETURN
            END IF
        END IF
        ! Refresh derived inventories only after all cells have been accepted.
        IF (executor%has_exchange_group) CALL state%refresh_substrate_totals()
        IF (profiling_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            profile_final_commit=profile_final_commit+REAL(tock-tick,wp)/REAL(clock_rate,wp)
        END IF
        ! Inject after the deferred geometry refresh, so the D20 rollback
        ! witness still verifies restoration of genuinely changed geometry.
        IF (executor%reject_active) THEN
            step_failed=.TRUE.
            failure_reason='D20_INJECTED_LATE_LIE_FAILURE'
            IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
            RETURN
        END IF

        CONTAINS

        ! Return a recoverable cell/group rejection without leaving the executor
        ! before returning the source-group status.
        SUBROUTINE reject_cell(message)
            CHARACTER(*), INTENT(IN) :: message
            failure_reason=message
            step_failed=.TRUE.
        END SUBROUTINE reject_cell
    END SUBROUTINE execute_lie_source_groups

    !> @brief Report cell-major sweep wall time and per-group summed thread times.
    !> @details Group times are summed thread CPU-like times, not wall-clock durations.
    SUBROUTINE report_lie_source_profile()
        INTEGER :: n
        IF (.NOT.profiling_initialized) RETURN
        IF (.NOT.profiling_enabled) RETURN
        WRITE(*,'(A)') 'LIE_GROUP_PROFILE_BEGIN'
        DO n=1,SIZE(group_profile)
            WRITE(*,'(A,I0,A,I0,A,I0)') 'LIE_GROUP_PROFILE id=',group_profile(n)%id,                                   &
                ' method=',group_profile(n)%method, ' calls=',group_profile(n)%calls
        END DO
        WRITE(*,'(A)') 'LIE_GROUP_PROFILE_END'
        WRITE(*,'(A)') 'LIE_SOURCE_DETAIL_BEGIN'
        DO n=1,SIZE(group_profile)
            WRITE(*,'(A,I0,A,I0,A,ES15.7,A,ES15.7,A,ES15.7)')                                     &
                'LIE_SOURCE_DETAIL id=',group_profile(n)%id,                                                     &
                ' evaluated_cells=',group_profile(n)%evaluated_cells,                                             &
                ' thread_seconds=',group_profile(n)%elapsed,                                                      &
                ' physics_thread_seconds=',group_profile(n)%physical,                                              &
                ' overhead_thread_seconds=',MAX(0.0_wp,group_profile(n)%elapsed-group_profile(n)%physical)
        END DO
        WRITE(*,'(A,ES15.7)') 'LIE_SOURCE_DETAIL cell_thread_seconds=',profile_thread_cell
        WRITE(*,'(A,ES15.7)') 'LIE_SOURCE_DETAIL outside_groups_thread_seconds=',                                  &
            MAX(0.0_wp,profile_thread_cell-SUM(group_profile%elapsed))
        WRITE(*,'(A)') 'LIE_SOURCE_DETAIL_END'
        WRITE(*,'(A,ES15.7)') 'LIE_CELL_SWEEP_PROFILE seconds=',profile_cell_sweep
        WRITE(*,'(A,ES15.7)') 'LIE_GEOMETRY_COMMIT_PROFILE seconds=',profile_final_commit
    END SUBROUTINE report_lie_source_profile

END MODULE lie_source_executor_2d
