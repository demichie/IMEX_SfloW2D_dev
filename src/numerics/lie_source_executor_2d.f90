!********************************************************************************
!> \brief Execute ordered first-order Lie physical source groups.
!> \details Coordinates simultaneous cell-local rate evaluation, physical
!>          limitation and atomic fluid/substrate transfer. Rates come from source_terms_2d,
!>          while bed projection and geometry refresh remain in geometry_2d.
!********************************************************************************

MODULE lie_source_executor_2d
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp, n_eqns, n_vars, n_solid, mass_exchange_rates_type,                                   &
        idx_totMassEqn, idx_solidEqn_first, idx_solidEqn_last
    USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
    USE domain_2d, ONLY : domain_type
    USE source_registry_2d, ONLY : source_registry_type, source_group_type
    USE state_2d, ONLY : state_type

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: execute_lie_source_groups,report_lie_source_profile
    TYPE :: group_profile_type
        INTEGER :: id=0,method=0,calls=0
        REAL(wp) :: setup=0.0_wp,evaluate=0.0_wp,commit=0.0_wp
    END TYPE group_profile_type
    TYPE(group_profile_type), ALLOCATABLE, SAVE :: group_profile(:)
    LOGICAL, SAVE :: profiling_initialized=.FALSE.,profiling_enabled=.FALSE.
    REAL(wp), SAVE :: profile_cell_sweep=0.0_wp,profile_final_commit=0.0_wp

CONTAINS

    !> @brief Evaluate and commit each ordered Lie physical-source group.
    !> @details Cell-major parallel sweep with frozen geometry. In group zero
    !>          evaluate source rates using the pre-spatial
    !>          snapshots but apply their limited increments to the post-FV
    !>          candidate. Later groups read the previously committed group
    !>          state; an inadmissible candidate rejects the group transaction.
    !> @param[in,out] state Owns cell-independent substrate inventories and aggregates.
    !> @param[in,out] q Conservative cell states (component,x-cell,y-cell).
    !> @param[in,out] qp Primitive cell states (component,x-cell,y-cell).
    !> @param[in] time_step Positive Lie substep duration [s].
    !> @param[in] domain Active cell list and per-cell coordinate lookup.
    !> @param[in] source_registry Initialized ordered Lie group definitions.
    !> @param[in] stochastic_values Optional stochastic cell samples.
    !> @param[in] step_start_time Optional beginning time for time-dependent sources [s].
    !> @param[in] q_original Pre-spatial conservative snapshot (read-only).
    !> @param[in] qp_original Pre-spatial primitive snapshot (read-only).
    !> @param[out] step_failed True when candidate group state is rejected.
    !> @param[out] failure_message Optional diagnostic describing rejected group.
    SUBROUTINE execute_lie_source_groups(state, q, qp, time_step, domain, source_registry, stochastic_values,          &
        step_start_time, q_original, qp_original, step_failed, failure_message)
    USE geometry_2d, ONLY : B_prime_x_geom,B_prime_y_geom,                                                             &
            B_second_xx_geom,B_second_xy_geom,B_second_yy_geom,grav_coeff,                                             &
            project_cell_field_to_vertices,assign_topography_vertices,B_vertex,                                        &
            comp_interfaces_x,comp_interfaces_y, B_zone, cell_source_fractions,                                        &
            cell_fissure_fractions,cell_arc_perim,cell_arc_n_x,cell_arc_n_y,dx,dy
    USE source_terms_2d, ONLY : assemble_lie_explicit_transaction,                                                     &
            eval_source_analytic,analytic_term_context_type
    USE source_catalog_2d, ONLY : SOURCE_METHOD_ANALYTIC,SOURCE_METHOD_EXPLICIT,                                       &
            SOURCE_METHOD_NEWTON
    USE nonlinear_solver_2d, ONLY : solve_backward_euler_group
    USE state_conversion_2d, ONLY : qc_to_qp
    USE state_2d, ONLY : initialize_dry_physical_cell,conservative_mass_admissible,                                    &
            physical_cell_state_finite
    USE parameters_2d, ONLY : topo_change_flag,erodible_deposit_flag
    USE constitutive_parameters_2d, ONLY : N_inh,f_inhibit_mode
        CLASS(state_type), INTENT(INOUT) :: state
        TYPE(source_group_type), ALLOCATABLE :: groups(:)
        TYPE(mass_exchange_rates_type) :: transfers
        TYPE(analytic_term_context_type) :: analytic_context
        REAL(wp), INTENT(INOUT) :: q(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(INOUT) :: qp(n_vars+2,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: time_step
        CLASS(domain_type), INTENT(IN) :: domain
        CLASS(source_registry_type), INTENT(IN) :: source_registry
        REAL(wp), INTENT(IN), OPTIONAL :: stochastic_values(:,:)
        REAL(wp), INTENT(IN), OPTIONAL :: step_start_time
        REAL(wp), INTENT(IN) :: q_original(n_vars,comp_cells_x,comp_cells_y)
        REAL(wp), INTENT(IN) :: qp_original(n_vars+2,comp_cells_x,comp_cells_y)
        LOGICAL, INTENT(OUT) :: step_failed
        CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: failure_message
        REAL(wp) :: qpentry(n_vars+2),q_candidate(n_vars),qp_candidate(n_vars+2)
        REAL(wp) :: q_local(n_vars),qp_local(n_vars+2)
        REAL(wp) :: deposit_local(n_solid),erosion_local(n_solid),erodible_local(n_solid),lava_local
        REAL(wp) :: deposit_candidate(n_solid),erosion_candidate(n_solid)
        REAL(wp) :: erodible_candidate(n_solid),lava_candidate
        INTEGER, ALLOCATABLE :: newton_iterations(:,:)
        REAL(wp), ALLOCATABLE :: bed_rate_cell(:,:),bed_rate_vertex(:,:),bed_candidate(:,:)
        REAL(wp), ALLOCATABLE :: staged_bed_rate(:)
        REAL(wp) :: total_rate(n_eqns)
        REAL(wp) :: analytic_rate(n_eqns),analytic_input(n_vars),pressure,source_time,noise
        LOGICAL :: with_exchange,cell_has_exchange,ok,failed,bed_accepted,conversion_ok,has_bed_group
        LOGICAL, ALLOCATABLE :: exchange_group(:)
        INTEGER :: group_index,cell,j,k,status,iterations,linear_info,first_failure,bed_group_id
        CHARACTER(LEN=256), ALLOCATABLE :: cell_failure(:)

        CHARACTER(LEN=96) :: context
        CHARACTER(LEN=256) :: failure_reason
        CHARACTER(LEN=24) :: newton_witness_env
        CHARACTER(LEN=24) :: reject_env
        INTEGER :: reject_group,env_status,parse_status
        INTEGER :: tick,tock,clock_rate
        CHARACTER(LEN=8) :: profile_env

        step_failed=.FALSE.
        failure_reason=''
        reject_group=-1
        reject_env=''
        CALL GET_ENVIRONMENT_VARIABLE('IMEX_LIE_TEST_REJECT_AFTER_GROUP', reject_env,STATUS=env_status)

        IF (env_status==0 .AND. LEN_TRIM(reject_env)>0) THEN
            READ(reject_env,*,IOSTAT=parse_status) reject_group
            IF (parse_status/=0 .OR. reject_group<1)                                                                   &
                CALL fatal_error('Invalid IMEX_LIE_TEST_REJECT_AFTER_GROUP')

        END IF

        IF (PRESENT(failure_message)) failure_message=''
        IF (.NOT.source_registry%is_initialized())                                                                     &
            CALL fatal_error('D19: uninitialized Lie source registry')
        groups=source_registry%groups()
        IF (.NOT.profiling_initialized) THEN
            profile_env=''
            CALL GET_ENVIRONMENT_VARIABLE('IMEX_PROFILE_LIE',profile_env,STATUS=env_status)
            profiling_enabled=env_status==0 .AND. TRIM(profile_env)=='1'
            profiling_initialized=.TRUE.
            IF (profiling_enabled) THEN
                ALLOCATE(group_profile(SIZE(groups)))
                DO group_index=1,SIZE(groups)
                    group_profile(group_index)%id=groups(group_index)%id
                    group_profile(group_index)%method=groups(group_index)%method
                END DO
            END IF
        END IF
        IF (profiling_enabled) THEN
            IF (SIZE(group_profile)/=SIZE(groups)) CALL fatal_error('Lie profiling group count changed')
        END IF
        source_time=0.0_wp
        IF (PRESENT(step_start_time)) source_time=step_start_time
        IF (.NOT.ieee_is_finite(time_step)) CALL fatal_error('D19: non-finite Lie timestep')
        IF (time_step<=0.0_wp) CALL fatal_error('D19: nonpositive time interval')

        ! The ordered Lie composition is cell-local. Shared geometry remains
        ! frozen until the entire trial source sweep has succeeded.
        ALLOCATE(exchange_group(SIZE(groups)))
        exchange_group=.FALSE.
        bed_group_id=-1
        has_bed_group=.FALSE.
        IF (profiling_enabled) CALL SYSTEM_CLOCK(tick,clock_rate)
        DO group_index=1,SIZE(groups)
            IF (.NOT.ALLOCATED(groups(group_index)%source_names))                                                      &
                CALL fatal_error('D20: Lie source group has no allocated source list')
            IF (SIZE(groups(group_index)%source_names)==0) CYCLE
            SELECT CASE(groups(group_index)%method)
            CASE(SOURCE_METHOD_EXPLICIT,SOURCE_METHOD_ANALYTIC,SOURCE_METHOD_NEWTON)
                CONTINUE
            CASE DEFAULT
                CALL fatal_error('D20: invalid Lie source-group method before OpenMP evaluation')
            END SELECT
            IF (ANY(groups(group_index)%source_names=='pore_pressure_dissipation') .OR.                                &
                ANY(groups(group_index)%source_names=='pore_pressure_darcy')) THEN
                IF (TRIM(ADJUSTL(f_inhibit_mode))=='STATIC' .OR.                                                       &
                    TRIM(ADJUSTL(f_inhibit_mode))=='DYNAMIC') THEN
                    IF (N_inh<0) CALL fatal_error('D20: negative Darcy inhibition order')
                END IF
            END IF
            exchange_group(group_index)=ANY(INDEX(groups(group_index)%source_names,'mass_exchange_')==1)
            IF (exchange_group(group_index)) THEN
                IF (groups(group_index)%method/=SOURCE_METHOD_EXPLICIT)                                                &
                    CALL fatal_error('D19: mass exchange requires explicit group method')
                IF (topo_change_flag) THEN
                    has_bed_group=.TRUE.
                    bed_group_id=groups(group_index)%id
                END IF
            END IF
            IF (profiling_enabled) group_profile(group_index)%calls=group_profile(group_index)%calls+1
        END DO
        ALLOCATE(cell_failure(domain%solve_cells))
        cell_failure=''
        IF (has_bed_group) THEN
            ALLOCATE(staged_bed_rate(domain%solve_cells))
            staged_bed_rate=0.0_wp
            ALLOCATE(bed_rate_cell(comp_cells_x,comp_cells_y))
            bed_rate_cell=0.0_wp
        END IF
        newton_witness_env=''
        CALL GET_ENVIRONMENT_VARIABLE('IMEX_NEWTON_WITNESS',newton_witness_env)
        IF (TRIM(newton_witness_env)=='1') THEN
            IF (ANY(groups%method==SOURCE_METHOD_NEWTON)) THEN
                ALLOCATE(newton_iterations(domain%solve_cells,SIZE(groups)))
                newton_iterations=0
            END IF
        END IF
        IF (profiling_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            group_profile(1)%setup=group_profile(1)%setup+REAL(tock-tick,wp)/REAL(clock_rate,wp)
            tick=tock
        END IF
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP& PRIVATE(cell,j,k,group_index,with_exchange,qpentry,q_candidate,qp_candidate, &
        !$OMP& q_local,qp_local,deposit_local,erosion_local,erodible_local,lava_local, &
        !$OMP& deposit_candidate,erosion_candidate,erodible_candidate,lava_candidate,total_rate, &
        !$OMP& analytic_rate,analytic_input,pressure,noise,analytic_context, &
        !$OMP& transfers,cell_has_exchange,ok,failed,conversion_ok,status,iterations,linear_info,context)
        cell_loop: DO cell=1,domain%solve_cells
            j=domain%j_cent(cell)
            k=domain%k_cent(cell)
            q_local=q(:,j,k)
            qp_local=qp(:,j,k)
            deposit_local=state%deposit(j,k,1:n_solid)
            erosion_local=state%erosion(j,k,1:n_solid)
            erodible_local=state%erodible(1:n_solid,j,k)
            lava_local=state%lava_deposit(j,k)
            group_loop: DO group_index=1,SIZE(groups)
                IF (SIZE(groups(group_index)%source_names)==0) CYCLE group_loop
                with_exchange=exchange_group(group_index)
                IF (groups(group_index)%id==0) THEN
                    qpentry=qp_original(:,j,k)
                ELSE
                    qpentry=qp_local

                END IF

                context='Lie source group'
                analytic_context%Bprime_x=B_prime_x_geom(j,k)
                analytic_context%Bprime_y=B_prime_y_geom(j,k)
                analytic_context%Bsecond_xx=B_second_xx_geom(j,k)
                analytic_context%Bsecond_xy=B_second_xy_geom(j,k)
                analytic_context%Bsecond_yy=B_second_yy_geom(j,k)
                analytic_context%grav_coeff=grav_coeff(j,k)
                noise=0.0_wp
                IF (PRESENT(stochastic_values)) noise=stochastic_values(j,k)
                analytic_context%Zij=noise

                ! Validate constitutive inputs before dispatch. Model/method errors
                ! remain configuration errors, not recoverable timestep rejections.
                IF (.NOT.physical_cell_state_finite(qpentry) .OR.                                                      &
                    .NOT.ieee_is_finite(noise)) THEN
                    cell_failure(cell)='D19: non-finite source-group physical input'
                    EXIT group_loop

                END IF

                IF (.NOT.ieee_is_finite(B_prime_x_geom(j,k)) .OR.                                                      &
                    .NOT.ieee_is_finite(B_prime_y_geom(j,k)) .OR.                                                      &
                    .NOT.ieee_is_finite(B_second_xx_geom(j,k)) .OR.                                                    &
                    .NOT.ieee_is_finite(B_second_xy_geom(j,k)) .OR.                                                    &
                    .NOT.ieee_is_finite(B_second_yy_geom(j,k)) .OR.                                                    &
                    .NOT.ieee_is_finite(grav_coeff(j,k))) THEN
                    cell_failure(cell)='D19: non-finite source-group geometry'
                    EXIT group_loop

                END IF

                IF (with_exchange) THEN
                    ! All immobilized and erodible inventories are nonnegative.
                    ! Validate the group input before computing any exchange rate.
                    IF (.NOT.ALL(ieee_is_finite(erodible_local)) .OR.                                   &
                        .NOT.ALL(ieee_is_finite(deposit_local)) .OR.                                    &
                        .NOT.ALL(ieee_is_finite(erosion_local)) .OR.                                    &
                        .NOT.ieee_is_finite(lava_local) .OR.                                              &
                        ANY(erodible_local<0.0_wp) .OR.                                                 &
                        ANY(deposit_local<0.0_wp) .OR.                                                  &
                        ANY(erosion_local<0.0_wp) .OR.                                                  &
                        lava_local<0.0_wp) THEN
                        cell_failure(cell)='D19: invalid substrate inventory before mass exchange'
                        EXIT group_loop

                    END IF

                END IF

                SELECT CASE(groups(group_index)%method)
                CASE(SOURCE_METHOD_EXPLICIT)
                    CALL assemble_lie_explicit_transaction(qpentry, q(1:n_eqns,j,k), time_step, source_time,           &
                        groups(group_index), erodible_local, total_rate, transfers, cell_has_exchange,  &
                        ok, noise, B_zone(j,k), analytic_context, cell_source_fractions(j,k),                          &
                        cell_fissure_fractions(j,k,:), cell_arc_perim(j,k), cell_arc_n_x(j,k), cell_arc_n_y(j,k), dx*dy)

                    IF (cell_has_exchange .NEQV. with_exchange) THEN
                        cell_failure(cell)='D20: inconsistent cell exchange group metadata'
                        EXIT group_loop

                    END IF

                    IF (.NOT.ok) THEN
                        cell_failure(cell)='D19: physical source transaction rejected'
                        EXIT group_loop

                    END IF
                    ! Stage the entire local fluid/substrate transfer before committing it.
                    q_candidate=q_local
                    q_candidate(1:n_eqns)=q_candidate(1:n_eqns)+time_step*total_rate

                    IF (.NOT.conservative_mass_admissible(q_candidate(1:n_eqns))) THEN
                        WRITE(cell_failure(cell),'(A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4)')                                 &
                            'D19: inadmissible source candidate group=',groups(group_index)%id,                        &
                            ' cell=',j,',',k,' mass=',q_candidate(idx_totMassEqn),                                     &
                            ' components=',SUM(q_candidate(idx_solidEqn_first:idx_solidEqn_last))
                        EXIT group_loop

                    END IF

                    IF (with_exchange) THEN
                        deposit_candidate=deposit_local+                                                &
                            time_step*transfers%solid_deposition_rate
                        erosion_candidate=erosion_local+                                                &
                            time_step*transfers%solid_erosion_rate
                        erodible_candidate=erodible_local-                                              &
                            time_step*transfers%solid_erosion_rate
                        IF (erodible_deposit_flag) erodible_candidate=                                                 &
                            erodible_candidate+time_step*transfers%solid_deposition_rate
                        lava_candidate=lava_local+time_step*transfers%lava_deposit_rate
                        IF (.NOT.ALL(ieee_is_finite(deposit_candidate)) .OR.                                           &
                            .NOT.ALL(ieee_is_finite(erosion_candidate)) .OR.                                           &
                            .NOT.ALL(ieee_is_finite(erodible_candidate)) .OR.                                          &
                            .NOT.ieee_is_finite(lava_candidate) .OR.                                                   &
                            ANY(erodible_candidate<0.0_wp) .OR.                                                        &
                            ANY(deposit_candidate<0.0_wp) .OR.                                                         &
                            ANY(erosion_candidate<0.0_wp) .OR.                                                         &
                            lava_candidate<0.0_wp) THEN
                            cell_failure(cell)='D19: invalid candidate substrate inventory'
                            EXIT group_loop

                        END IF

                    END IF

                    IF (q_candidate(idx_totMassEqn)>0.0_wp) THEN
                        CALL qc_to_qp(q_candidate,qp_candidate,pressure,'Lie explicit group candidate',conversion_ok)

                        IF (.NOT.conversion_ok) THEN
                            cell_failure(cell)='D19: explicit candidate conversion failed'
                            EXIT group_loop

                        END IF
                    ELSE
                        CALL initialize_dry_physical_cell(qp_candidate)

                    END IF
                    ! Publish to cell-local state; the trial array is written once after all groups.
                    q_local=q_candidate
                    qp_local=qp_candidate

                    IF (with_exchange) THEN
                        deposit_local=deposit_candidate
                        erosion_local=erosion_candidate
                        erodible_local=erodible_candidate
                        lava_local=lava_candidate
                        IF (topo_change_flag) staged_bed_rate(cell)=staged_bed_rate(cell)+transfers%dB_dt

                    END IF
                CASE(SOURCE_METHOD_ANALYTIC)
                    ! Group 0 evaluates on the immutable pre-spatial state;
                    ! its increment is applied to the post-spatial group candidate.
                    analytic_input=q_local
                    IF (groups(group_index)%id==0) analytic_input=q_original(:,j,k)
                    IF (analytic_input(idx_totMassEqn)<=0.0_wp) CYCLE
                    CALL eval_source_analytic(groups(group_index), analytic_input, analytic_context, time_step,        &
                        analytic_rate, status)

                    IF (status/=0) THEN
                        WRITE(cell_failure(cell),'(A,I0)') 'D19: analytic group status ',status
                        EXIT group_loop

                    END IF

                    IF (.NOT.ALL(ieee_is_finite(analytic_rate))) THEN
                        cell_failure(cell)='D19: analytic source returned non-finite rate'
                        EXIT group_loop

                    END IF
                    q_candidate=q_local
                    q_candidate(1:n_eqns)=q_candidate(1:n_eqns)+time_step*analytic_rate
                CASE(SOURCE_METHOD_NEWTON)
                    IF (q(idx_totMassEqn,j,k)<=0.0_wp) CYCLE
                    q_candidate=q_local
                    CALL solve_backward_euler_group(groups(group_index), q_candidate, time_step,                      &
                        B_prime_x_geom(j,k), B_prime_y_geom(j,k), grav_coeff(j,k), noise, iterations, ok, linear_info, &
                        failed)
                    IF (ALLOCATED(newton_iterations)) newton_iterations(cell,group_index)=iterations

                    IF (.NOT.ok .OR. failed .OR. linear_info/=0) THEN
                        WRITE(cell_failure(cell),'(A,I0,A,I0)')                                                        &
                            'D19: Newton failed, linear_info=',linear_info,', iterations=',iterations
                        EXIT group_loop

                    END IF
                CASE DEFAULT
                    cell_failure(cell)='D19: invalid group method'
                    EXIT group_loop

                END SELECT
                ! Explicit candidates were already validated and converted above.
                IF (groups(group_index)%method==SOURCE_METHOD_EXPLICIT) CYCLE group_loop
                ! Only check admissibility: never clip the accepted source update.
                IF (.NOT.conservative_mass_admissible(q_candidate(1:n_eqns))) THEN
                    cell_failure(cell)='D19: inadmissible conservative source group state'
                    EXIT group_loop

                END IF

                IF (q_candidate(idx_totMassEqn)>0.0_wp) THEN
                    CALL qc_to_qp(q_candidate,qp_candidate,pressure,TRIM(context),conversion_ok)

                    IF (.NOT.conversion_ok) THEN
                        cell_failure(cell)='D19: source group conversion failed'
                        EXIT group_loop

                    END IF
                ELSE
                    CALL initialize_dry_physical_cell(qp_candidate)

                END IF

                q_local=q_candidate
                qp_local=qp_candidate


            END DO group_loop
            IF (LEN_TRIM(cell_failure(cell))==0) THEN
                q(:,j,k)=q_local
                qp(:,j,k)=qp_local
                state%deposit(j,k,1:n_solid)=deposit_local
                state%erosion(j,k,1:n_solid)=erosion_local
                state%erodible(1:n_solid,j,k)=erodible_local
                state%lava_deposit(j,k)=lava_local
            END IF
        END DO cell_loop
        !$OMP END PARALLEL DO
        IF (profiling_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            ! With cell-major traversal this is a *whole-sweep* wall-clock time,
            ! not an individual group evaluation. Report it separately.
            profile_cell_sweep=profile_cell_sweep+REAL(tock-tick,wp)/REAL(clock_rate,wp)
            tick=tock
        END IF
        first_failure=0
        DO cell=1,domain%solve_cells
            IF (LEN_TRIM(cell_failure(cell))>0) THEN
                first_failure=cell
                EXIT
            END IF
        END DO
        IF (first_failure/=0) THEN
            step_failed=.TRUE.
            failure_reason=cell_failure(first_failure)
            IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
            RETURN
        END IF
        IF (ALLOCATED(newton_iterations)) THEN
            DO group_index=1,SIZE(groups)
                IF (groups(group_index)%method/=SOURCE_METHOD_NEWTON) CYCLE
                WRITE(*,'(A,I0,A,I0,A,I0)') 'D20_NEWTON_WITNESS group=',                                               &
                    groups(group_index)%id,' solved_cells=',COUNT(newton_iterations(:,group_index)>0),                 &
                    ' total_iterations=',SUM(newton_iterations(:,group_index))
            END DO
            DEALLOCATE(newton_iterations)
        END IF
        ! Geometry is the same for every source group in the step. Project
        ! all mass-exchange bed rates only once, after the cell sweep.
        IF (has_bed_group) THEN
            ALLOCATE(bed_rate_vertex(comp_interfaces_x,comp_interfaces_y))
            CALL record_group_bed_rates
            bed_accepted=ALL(ieee_is_finite(bed_rate_cell))
            IF (bed_accepted) THEN
                CALL project_cell_field_to_vertices(bed_rate_cell,bed_rate_vertex)
                bed_accepted=ALL(ieee_is_finite(bed_rate_vertex))
            END IF
            IF (bed_accepted) THEN
                bed_candidate=B_vertex+time_step*bed_rate_vertex
                IF (.NOT.ALL(ieee_is_finite(bed_candidate))) bed_accepted=.FALSE.
            END IF
            IF (.NOT.bed_accepted) THEN
                CALL reject_cell('D19: invalid candidate projected bed')
                IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
                RETURN
            END IF
            WRITE(*,'(A,I0,A,ES16.8)') 'D20_BED_WITNESS group=',bed_group_id,                                           &
                ' max_abs_delta=',MAXVAL(ABS(bed_candidate-B_vertex))
            CALL assign_topography_vertices(bed_candidate,bed_accepted)
            IF (.NOT.bed_accepted) THEN
                CALL reject_cell('D19: invalid refreshed geometry')
                IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
                RETURN
            END IF
        END IF
        IF (ANY(exchange_group)) CALL state%refresh_substrate_totals()
        IF (profiling_enabled) THEN
            CALL SYSTEM_CLOCK(tock)
            profile_final_commit=profile_final_commit+REAL(tock-tick,wp)/REAL(clock_rate,wp)
        END IF
        ! Inject after the deferred geometry refresh, so the D20 rollback
        ! witness still verifies restoration of genuinely changed geometry.
        IF (ANY(groups%id==reject_group)) THEN
            step_failed=.TRUE.
            failure_reason='D20_INJECTED_LATE_LIE_FAILURE'
            IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
            RETURN
        END IF

    CONTAINS
    SUBROUTINE record_group_bed_rates
        INTEGER :: n,j_cell,k_cell

        DO n=1,domain%solve_cells
            j_cell=domain%j_cent(n)
            k_cell=domain%k_cent(n)
            bed_rate_cell(j_cell,k_cell)=staged_bed_rate(n)

        END DO
    END SUBROUTINE record_group_bed_rates

        ! Return a recoverable cell/group rejection without leaving the executor
        ! before returning the source-group status.
    SUBROUTINE reject_cell(message)
        CHARACTER(*), INTENT(IN) :: message
        failure_reason=message
        step_failed=.TRUE.
    END SUBROUTINE reject_cell
    END SUBROUTINE execute_lie_source_groups

    !> @brief Emit group call counts and cell-major sweep/commit wall times.
    !> @details Separate group wall times are not defined for the cell-major sweep.
    SUBROUTINE report_lie_source_profile()
        INTEGER :: n
        IF (.NOT.profiling_initialized) RETURN
        IF (.NOT.profiling_enabled) RETURN
        WRITE(*,'(A)') 'LIE_GROUP_PROFILE_BEGIN'
        DO n=1,SIZE(group_profile)
            WRITE(*,'(A,I0,A,I0,A,I0)') &
                'LIE_GROUP_PROFILE id=',group_profile(n)%id, &
                ' method=',group_profile(n)%method, &
                ' calls=',group_profile(n)%calls
        END DO
        WRITE(*,'(A)') 'LIE_GROUP_PROFILE_END'
        WRITE(*,'(A,ES15.7)') 'LIE_CELL_SWEEP_PROFILE seconds=',profile_cell_sweep
        WRITE(*,'(A,ES15.7)') 'LIE_GEOMETRY_COMMIT_PROFILE seconds=',profile_final_commit
    END SUBROUTINE report_lie_source_profile

END MODULE lie_source_executor_2d
