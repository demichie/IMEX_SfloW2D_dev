!********************************************************************************
!> \brief Execute ordered first-order Lie physical source groups.
!> \details Coordinates simultaneous cell-local rate evaluation, physical
!>          limitation and atomic fluid/substrate transfer. Rates come from source_terms_2d,
!>          while bed projection and geometry refresh remain in geometry_2d.
!********************************************************************************

MODULE lie_source_executor_2d
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    USE diagnostics_2d, ONLY : fatal_error
    USE parameters_2d, ONLY : wp, n_eqns, n_vars, n_solid, mass_exchange_rates_type, &
        idx_totMassEqn, idx_solidEqn_first, idx_solidEqn_last
    USE geometry_2d, ONLY : comp_cells_x, comp_cells_y
    USE domain_2d, ONLY : domain_type
    USE source_registry_2d, ONLY : source_registry_type, source_group_type
    USE state_2d, ONLY : state_type

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: execute_lie_source_groups

CONTAINS

    !> @brief Evaluate and commit each ordered Lie physical-source group.
    !> @details In group zero evaluate source rates using the pre-spatial
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
    SUBROUTINE execute_lie_source_groups(state,q,qp,time_step,domain,source_registry, &
            stochastic_values,step_start_time,q_original,qp_original,step_failed,failure_message)
        USE geometry_2d, ONLY : B_prime_x_geom,B_prime_y_geom, &
            B_second_xx_geom,B_second_xy_geom,B_second_yy_geom,grav_coeff, &
            project_cell_field_to_vertices,assign_topography_vertices,B_vertex, &
            comp_interfaces_x,comp_interfaces_y, B_zone, cell_source_fractions, &
            cell_fissure_fractions,cell_arc_perim,cell_arc_n_x,cell_arc_n_y,dx,dy
        USE source_terms_2d, ONLY : assemble_lie_explicit_transaction, &
            eval_source_analytic,analytic_term_context_type
        USE source_catalog_2d, ONLY : SOURCE_METHOD_ANALYTIC,SOURCE_METHOD_EXPLICIT, &
            SOURCE_METHOD_NEWTON
        USE nonlinear_solver_2d, ONLY : solve_backward_euler_group
        USE state_conversion_2d, ONLY : qc_to_qp
        USE state_2d, ONLY : initialize_dry_physical_cell,conservative_mass_admissible, &
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
        REAL(wp), ALLOCATABLE :: q_group(:,:,:),qp_group(:,:,:)
        REAL(wp) :: deposit_candidate(n_solid),erosion_candidate(n_solid)
        REAL(wp) :: erodible_candidate(n_solid),lava_candidate
        INTEGER, ALLOCATABLE :: group_j(:),group_k(:),newton_iterations(:)
        REAL(wp), ALLOCATABLE :: bed_rate_cell(:,:),bed_rate_vertex(:,:),bed_candidate(:,:)
        REAL(wp), ALLOCATABLE :: staged_deposit(:,:),staged_erosion(:,:), &
            staged_erodible(:,:),staged_lava(:),staged_bed_rate(:)
        REAL(wp) :: total_rate(n_eqns)
        REAL(wp) :: analytic_rate(n_eqns),analytic_input(n_vars),pressure,source_time,noise
        LOGICAL :: with_exchange,cell_has_exchange,ok,failed,bed_accepted,conversion_ok
        INTEGER :: group_index,cell,j,k,status,iterations,linear_info,first_failure
        CHARACTER(LEN=256), ALLOCATABLE :: cell_failure(:)

        CHARACTER(LEN=96) :: context
        CHARACTER(LEN=256) :: failure_reason

        step_failed=.FALSE.
        failure_reason=''
        IF (PRESENT(failure_message)) failure_message=''
        IF (.NOT.source_registry%is_initialized()) &
            CALL fatal_error('D19: uninitialized Lie source registry')
        groups=source_registry%groups()
        source_time=0.0_wp
        IF (PRESENT(step_start_time)) source_time=step_start_time
        IF (.NOT.ieee_is_finite(time_step)) CALL fatal_error('D19: non-finite Lie timestep')
        IF (time_step<=0.0_wp) CALL fatal_error('D19: nonpositive time interval')

        DO group_index=1,SIZE(groups)
            q_group=q
            qp_group=qp
            ! Validate configuration before spawning workers: fatal dispatch errors
            ! must not originate from a concurrent cell evaluation.
            IF (.NOT.ALLOCATED(groups(group_index)%source_names)) &
                CALL fatal_error('D20: Lie source group has no allocated source list')
            IF (SIZE(groups(group_index)%source_names)==0) CYCLE
            SELECT CASE(groups(group_index)%method)
            CASE(SOURCE_METHOD_EXPLICIT,SOURCE_METHOD_ANALYTIC,SOURCE_METHOD_NEWTON)
                CONTINUE
            CASE DEFAULT
                CALL fatal_error('D20: invalid Lie source-group method before OpenMP evaluation')
            END SELECT
            ! Darcy's inhibition polynomial must have a valid configured order
            ! before any worker evaluates its cell-local pressure source.
            IF (ANY(groups(group_index)%source_names=='pore_pressure_dissipation') .OR. &
                ANY(groups(group_index)%source_names=='pore_pressure_darcy')) THEN
                IF (TRIM(ADJUSTL(f_inhibit_mode))=='STATIC' .OR. &
                    TRIM(ADJUSTL(f_inhibit_mode))=='DYNAMIC') THEN
                    IF (N_inh<0) CALL fatal_error('D20: negative Darcy inhibition order')
                END IF
            END IF
            with_exchange=ANY(INDEX(groups(group_index)%source_names,'mass_exchange_')==1)
            IF (with_exchange) THEN
                IF (groups(group_index)%method/=SOURCE_METHOD_EXPLICIT) &
                    CALL fatal_error('D19: mass exchange requires explicit group method')
                ALLOCATE(group_j(domain%solve_cells),group_k(domain%solve_cells))
                ALLOCATE(staged_deposit(n_solid,domain%solve_cells), &
                    staged_erosion(n_solid,domain%solve_cells), &
                    staged_erodible(n_solid,domain%solve_cells), &
                    staged_lava(domain%solve_cells),staged_bed_rate(domain%solve_cells))
                IF (topo_change_flag) THEN
                    ALLOCATE(bed_rate_cell(comp_cells_x,comp_cells_y))
                    bed_rate_cell=0.0_wp
                END IF
            END IF

            ALLOCATE(cell_failure(domain%solve_cells))
            cell_failure=''
            IF (groups(group_index)%method==SOURCE_METHOD_NEWTON) THEN
                ALLOCATE(newton_iterations(domain%solve_cells))
                newton_iterations=0
            END IF
            ! q and qp remain immutable group-entry states until the serial commit.
            ! Only q_group/qp_group candidates and per-cell staging are modified.
            ! Each active cell owns distinct candidate slices and staging slots.
            !$OMP PARALLEL DO DEFAULT(SHARED) &
            !$OMP& PRIVATE(cell,j,k,qpentry,q_candidate,qp_candidate,deposit_candidate, &
            !$OMP& erosion_candidate,erodible_candidate,lava_candidate,total_rate, &
            !$OMP& analytic_rate,analytic_input,pressure,noise,analytic_context, &
            !$OMP& transfers,cell_has_exchange,ok,failed,conversion_ok,status,iterations,linear_info,context)
            cell_loop: DO cell=1,domain%solve_cells
                j=domain%j_cent(cell)
                k=domain%k_cent(cell)
                IF (groups(group_index)%id==0) THEN
                    qpentry=qp_original(:,j,k)
                ELSE
                    qpentry=qp(:,j,k)
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
                IF (.NOT.physical_cell_state_finite(qpentry) .OR. &
                    .NOT.ieee_is_finite(noise)) THEN
                    cell_failure(cell)='D19: non-finite source-group physical input'
                    CYCLE cell_loop
                END IF
                IF (.NOT.ieee_is_finite(B_prime_x_geom(j,k)) .OR. &
                    .NOT.ieee_is_finite(B_prime_y_geom(j,k)) .OR. &
                    .NOT.ieee_is_finite(B_second_xx_geom(j,k)) .OR. &
                    .NOT.ieee_is_finite(B_second_xy_geom(j,k)) .OR. &
                    .NOT.ieee_is_finite(B_second_yy_geom(j,k)) .OR. &
                    .NOT.ieee_is_finite(grav_coeff(j,k))) THEN
                    cell_failure(cell)='D19: non-finite source-group geometry'
                    CYCLE cell_loop
                END IF
                IF (with_exchange) THEN
                    ! All immobilized and erodible inventories are nonnegative.
                    ! Validate the group input before computing any exchange rate.
                    IF (.NOT.ALL(ieee_is_finite(state%erodible(1:n_solid,j,k))) .OR. &
                        .NOT.ALL(ieee_is_finite(state%deposit(j,k,1:n_solid))) .OR. &
                        .NOT.ALL(ieee_is_finite(state%erosion(j,k,1:n_solid))) .OR. &
                        .NOT.ieee_is_finite(state%lava_deposit(j,k)) .OR. &
                        ANY(state%erodible(1:n_solid,j,k)<0.0_wp) .OR. &
                        ANY(state%deposit(j,k,1:n_solid)<0.0_wp) .OR. &
                        ANY(state%erosion(j,k,1:n_solid)<0.0_wp) .OR. &
                        state%lava_deposit(j,k)<0.0_wp) THEN
                        cell_failure(cell)='D19: invalid substrate inventory before mass exchange'
                        CYCLE cell_loop
                    END IF
                END IF

                SELECT CASE(groups(group_index)%method)
                CASE(SOURCE_METHOD_EXPLICIT)
                    CALL assemble_lie_explicit_transaction(qpentry,q(1:n_eqns,j,k), &
                        time_step,source_time,groups(group_index),state%erodible(1:n_solid,j,k), &
                        total_rate,transfers,cell_has_exchange,ok,noise,B_zone(j,k),analytic_context, &
                        cell_source_fractions(j,k),cell_fissure_fractions(j,k,:), &
                        cell_arc_perim(j,k),cell_arc_n_x(j,k),cell_arc_n_y(j,k),dx*dy)
                    IF (cell_has_exchange .NEQV. with_exchange) THEN
                        cell_failure(cell)='D20: inconsistent cell exchange group metadata'
                        CYCLE cell_loop
                    END IF
                    IF (.NOT.ok) THEN
                        cell_failure(cell)='D19: physical source transaction rejected'
                        CYCLE cell_loop
                    END IF
                    ! Stage the entire local fluid/substrate transfer before committing it.
                    q_candidate=q(:,j,k)
                    q_candidate(1:n_eqns)=q_candidate(1:n_eqns)+time_step*total_rate
                    IF (.NOT.conservative_mass_admissible(q_candidate(1:n_eqns))) THEN
                        WRITE(cell_failure(cell),'(A,I0,A,I0,A,I0,A,ES12.4,A,ES12.4)') &
                            'D19: inadmissible source candidate group=',groups(group_index)%id, &
                            ' cell=',j,',',k,' mass=',q_candidate(idx_totMassEqn), &
                            ' components=',SUM(q_candidate(idx_solidEqn_first:idx_solidEqn_last))
                        CYCLE cell_loop
                    END IF
                    IF (with_exchange) THEN
                        deposit_candidate=state%deposit(j,k,1:n_solid)+ &
                            time_step*transfers%solid_deposition_rate
                        erosion_candidate=state%erosion(j,k,1:n_solid)+ &
                            time_step*transfers%solid_erosion_rate
                        erodible_candidate=state%erodible(1:n_solid,j,k)- &
                            time_step*transfers%solid_erosion_rate
                        IF (erodible_deposit_flag) erodible_candidate= &
                            erodible_candidate+time_step*transfers%solid_deposition_rate
                        lava_candidate=state%lava_deposit(j,k)+time_step*transfers%lava_deposit_rate
                        IF (.NOT.ALL(ieee_is_finite(deposit_candidate)) .OR. &
                            .NOT.ALL(ieee_is_finite(erosion_candidate)) .OR. &
                            .NOT.ALL(ieee_is_finite(erodible_candidate)) .OR. &
                            .NOT.ieee_is_finite(lava_candidate) .OR. &
                            ANY(erodible_candidate<0.0_wp) .OR. &
                            ANY(deposit_candidate<0.0_wp) .OR. &
                            ANY(erosion_candidate<0.0_wp) .OR. &
                            lava_candidate<0.0_wp) THEN
                            cell_failure(cell)='D19: invalid candidate substrate inventory'
                            CYCLE cell_loop
                        END IF
                    END IF
                    IF (q_candidate(idx_totMassEqn)>0.0_wp) THEN
                        CALL qc_to_qp(q_candidate,qp_candidate,pressure,'Lie explicit group candidate',conversion_ok)
                        IF (.NOT.conversion_ok) THEN
                            cell_failure(cell)='D19: explicit candidate conversion failed'
                            CYCLE cell_loop
                        END IF
                    ELSE
                        CALL initialize_dry_physical_cell(qp_candidate)
                    END IF
                    ! Commit only after the complete cell-local candidate is admissible.
                    q_group(:,j,k)=q_candidate
                    qp_group(:,j,k)=qp_candidate
                    IF (with_exchange) THEN
                        group_j(cell)=j
                        group_k(cell)=k
                        staged_deposit(:,cell)=deposit_candidate
                        staged_erosion(:,cell)=erosion_candidate
                        staged_erodible(:,cell)=erodible_candidate
                        staged_lava(cell)=lava_candidate
                        staged_bed_rate(cell)=transfers%dB_dt
                    END IF
                CASE(SOURCE_METHOD_ANALYTIC)
                    ! Group 0 evaluates on the immutable pre-spatial state;
                    ! its increment is applied to the post-spatial group candidate.
                    analytic_input=q(:,j,k)
                    IF (groups(group_index)%id==0) analytic_input=q_original(:,j,k)
                    IF (analytic_input(idx_totMassEqn)<=0.0_wp) CYCLE
                    CALL eval_source_analytic(groups(group_index),analytic_input, &
                        analytic_context,time_step,analytic_rate,status)
                    IF (status/=0) THEN
                        WRITE(cell_failure(cell),'(A,I0)') 'D19: analytic group status ',status
                        CYCLE cell_loop
                    END IF
                    IF (.NOT.ALL(ieee_is_finite(analytic_rate))) THEN
                        cell_failure(cell)='D19: analytic source returned non-finite rate'
                        CYCLE cell_loop
                    END IF
                    ! The analytic method must update the group candidate only.
                    ! Publishing q here would bypass the atomic group commit.
                    q_group(1:n_eqns,j,k)=q_group(1:n_eqns,j,k)+time_step*analytic_rate
                CASE(SOURCE_METHOD_NEWTON)
                    IF (q_group(idx_totMassEqn,j,k)<=0.0_wp) CYCLE
                    CALL solve_backward_euler_group(groups(group_index), &
                        q_group(:,j,k),time_step,B_prime_x_geom(j,k),B_prime_y_geom(j,k), &
                        grav_coeff(j,k),noise,iterations,ok,linear_info,failed)
                    newton_iterations(cell)=iterations
                    IF (.NOT.ok .OR. failed .OR. linear_info/=0) THEN
                        WRITE(cell_failure(cell),'(A,I0,A,I0)') &
                            'D19: Newton failed, linear_info=',linear_info,', iterations=',iterations
                        CYCLE cell_loop
                    END IF
                CASE DEFAULT
                    cell_failure(cell)='D19: invalid group method'
                    CYCLE cell_loop
                END SELECT
                ! Explicit candidates were already validated and converted above.
                IF (groups(group_index)%method==SOURCE_METHOD_EXPLICIT) CYCLE cell_loop
                ! Only check admissibility: never clip the accepted source update.
                IF (.NOT.conservative_mass_admissible(q_group(1:n_eqns,j,k))) THEN
                    cell_failure(cell)='D19: inadmissible conservative source group state'
                    CYCLE cell_loop
                END IF
                IF (q_group(idx_totMassEqn,j,k)>0.0_wp) THEN
                    CALL qc_to_qp(q_group(:,j,k),qp_group(:,j,k),pressure,TRIM(context),conversion_ok)
                    IF (.NOT.conversion_ok) THEN
                        cell_failure(cell)='D19: source group conversion failed'
                        CYCLE cell_loop
                    END IF
                ELSE
                    CALL initialize_dry_physical_cell(qp_group(:,j,k))
                END IF
            END DO cell_loop
            !$OMP END PARALLEL DO
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
            END IF
            DEALLOCATE(cell_failure)
            IF (ALLOCATED(newton_iterations)) THEN
                WRITE(*,'(A,I0,A,I0,A,I0)') 'D20_NEWTON_WITNESS group=', &
                    groups(group_index)%id,' solved_cells=',COUNT(newton_iterations>0), &
                    ' total_iterations=',SUM(newton_iterations)
                DEALLOCATE(newton_iterations)
            END IF
            IF (step_failed) THEN
                IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
                RETURN
            END IF
            ! Publish inventory changes only after all cells in the group succeed.
            IF (with_exchange .AND. topo_change_flag) THEN
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
            END IF
            IF (with_exchange) THEN
                ! Accept the geometric projection before publishing inventory changes.
                ! A rejected geometric refresh therefore leaves both inventories
                ! and the conserved/primitive group state unmodified.
                IF (topo_change_flag) THEN
                    WRITE(*,'(A,I0,A,ES16.8)') 'D20_BED_WITNESS group=', &
                        groups(group_index)%id,' max_abs_delta=', &
                        MAXVAL(ABS(bed_candidate-B_vertex))
                    CALL assign_topography_vertices(bed_candidate,bed_accepted)
                    IF (.NOT.bed_accepted) THEN
                        CALL reject_cell('D19: invalid refreshed geometry')
                        IF (PRESENT(failure_message)) failure_message=TRIM(failure_reason)
                        RETURN
                    END IF
                    DEALLOCATE(bed_rate_cell,bed_rate_vertex,bed_candidate)
                END IF
                DO cell=1,domain%solve_cells
                    j=group_j(cell)
                    k=group_k(cell)
                    state%deposit(j,k,1:n_solid)=staged_deposit(:,cell)
                    state%erosion(j,k,1:n_solid)=staged_erosion(:,cell)
                    state%erodible(1:n_solid,j,k)=staged_erodible(:,cell)
                    state%lava_deposit(j,k)=staged_lava(cell)
                END DO
                CALL state%refresh_substrate_totals()
                DEALLOCATE(group_j,group_k,staged_deposit,staged_erosion, &
                    staged_erodible,staged_lava,staged_bed_rate)
            END IF
            q=q_group
            qp=qp_group
        END DO
    CONTAINS
        SUBROUTINE record_group_bed_rates
        INTEGER :: n,j_cell,k_cell
        DO n=1,domain%solve_cells
            j_cell=group_j(n)
            k_cell=group_k(n)
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

END MODULE lie_source_executor_2d
