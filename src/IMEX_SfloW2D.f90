!> @mainpage IMEX-SfloW2D - Finite-volume solver for geophysical flows
!> @brief Open-source Fortran solver for one- and two-layer depth-averaged flows.
!>
!> IMEX-SfloW2D solves two-dimensional depth-averaged mixture balance equations
!> for gravity-driven flows over complex terrain, including granular, lahar,
!> pyroclastic-current and lava configurations. Material composition, sensible
!> thermal content, pore-pressure variables and optional stochastic transport
!> are included according to the selected physical configuration.
!>
!> @section model Model and numerical formulation
!> The spatial discretization is a hydrostatic-path, path-conservative
!> central-upwind finite-volume method (HP-PCCU). The implemented timestep
!> advances the spatial residual first and then applies ordered Lie physical
!> source groups. Explicit, semi-implicit, Newton/complex-step implicit and
!> analytic treatments are selected according to the supported source method.
!> The solver uses an active-cell workset with OpenMP parallelism.
!>
!> @section layers Supported material layers
!> N_LAYERS=1 and N_LAYERS=2 use the same execution path. In two-layer runs,
!> the layers have independently configured physical laws and coupled spatial
!> pressure and interlayer transfer/drag sources. Terrain and substrate
!> inventories are shared; the conservative and primitive state arrays use
!> q/qp(component,layer,x,y). No environment activation switch is needed.
!>
!> @section limitations Scientific applicability
!> The coupled two-layer characteristic spectrum is not unconditionally real.
!> Density inversion and sufficiently strong interlayer shear can produce
!> complex wave speeds. The total-depth CFL bound does not establish
!> hyperbolicity or quantitative validity of arbitrary physical closures.
!> See the two-layer hyperbolicity and physical-closure documents for limits.
!>
!> @section output Input, output and continuation
!> The runtime reads IMEX_SfloW2D.inp and its configured initial/terrain data.
!> It writes configured physical outputs (including NetCDF), diagnostics and
!> a canonical binary restart checkpoint for the selected layer count.
!>
!> @section resources Documentation and project resources
!> - Project: https://github.com/demichie/IMEX_SfloW2D_dev
!> - User guide: https://demichie.github.io/IMEX_SfloW2D_dev/
!> - Fortran API: https://demichie.github.io/IMEX_SfloW2D_dev/doxygen/
!> - Architecture: docs/TWO_LAYER_ARCHITECTURE_AND_ROADMAP.md
!> - Hyperbolicity: docs/TWO_LAYER_HYPERBOLICITY.md
!> - Physical closures: docs/TL13B_C_PHYSICAL_CLOSURES.md
!>
!> @author Mattia de' Michieli Vitturi (original solver author).

!> @brief Initialize, advance and write one depth-averaged simulation.

PROGRAM IMEX_SfloW2D

    USE, intrinsic :: iso_fortran_env
    USE, intrinsic :: ieee_arithmetic

    USE constitutive_parameters_2d, ONLY : T_ambient


    USE geometry_2d, ONLY : init_grid
    USE geometry_2d, ONLY : init_source
    USE geometry_2d, ONLY : refresh_topography_geometry

    USE geometry_2d, ONLY : dx,dy
    ! USE geometry_2d, ONLY : comp_cells_x,comp_cells_y

    USE init_2d, ONLY : collapsing_volume
    USE init_2d, ONLY : init_empty, release_initialization_fields

    USE inpout_2d, ONLY : init_param
    USE inpout_2d, ONLY : read_param
    USE inpout_2d, ONLY : update_param
    USE inpout_2d, ONLY : output_solution
    USE inpout_2d, ONLY : output_runout
    USE inpout_2d, ONLY : open_runout_outputs
    USE inpout_2d, ONLY : validate_runout_outputs
    USE inpout_2d, ONLY : read_solution
    USE inpout_2d, ONLY : close_units
    USE inpout_2d, ONLY : n_probes
    USE inpout_2d, ONLY : output_probes

    USE inpout_2d, ONLY : init_netcdf_output
    USE inpout_2d, ONLY : close_netcdf
    USE inpout_2d, ONLY : restart, restart_file
    USE inpout_2d, ONLY : read_restart_file, write_restart_file
    USE inpout_2d, ONLY : validate_restart_output_manifest
    USE inpout_2d, ONLY : write_param_backup
    USE output_paths, ONLY : output_file_path
    USE inpout_2d, ONLY : output_runout_flag
    USE inpout_2d, ONLY : output_cons_flag
    USE inpout_2d, ONLY : output_esri_flag
    USE inpout_2d, ONLY : output_netcdf_flag

    USE simulation_2d, ONLY : simulation

    USE parameters_2d, ONLY : wp, idx_h, idx_T, idx_u, idx_v, idx_totMassEqn

    USE parameters_2d, ONLY : t_start
    USE parameters_2d, ONLY : t_end
    USE parameters_2d, ONLY : dt_output
    USE parameters_2d, ONLY : dt0
    USE substrate_parameters_2d, ONLY : substrate
    USE parameters_2d, ONLY : verbose_level
    USE parameters_2d, ONLY : n_solid
    USE parameters_2d, ONLY : n_vars,n_layers
    USE rheology_parameters_2d, ONLY : rheology_l1,rheology_l2
    USE parameters_2d, ONLY : radial_source_flag
    USE parameters_2d, ONLY : lateral_source_flag, bottom_radial_source_flag
    USE parameters_2d, ONLY : collapsing_volume_flag
    USE parameters_2d, ONLY : serial_flag


    USE parameters_2d, ONLY : n_thickness_levels , n_dyn_pres_levels ,                                                 &
        thickness_levels , dyn_pres_levels

    USE state_conversions, ONLY : conservative_to_primitive, dynamic_pressure_from_synchronized_state
    USE diagnostics_2d, ONLY : fatal_error


    USE OMP_LIB

    IMPLICIT NONE

    REAL(wp) :: t1 , t2 , t3

    REAL(wp) :: rate
    INTEGER :: st1 , st2 , st3 , cr , cm

    LOGICAL :: stop_flag
    LOGICAL :: thickness_exceeds, pressure_exceeds
    LOGICAL :: stop_flag_old

    INTEGER :: j,k,l,layer

    INTEGER :: i_pdyn_lev , i_thk_lev , i_table

    INTEGER n_threads

    LOGICAL :: use_openmp = .false.

    !> Dynamic pressure
    REAL(wp) :: p_dyn
    LOGICAL :: conversion_ok,conversion_failed

    REAL(wp) :: mod_vel , mod_vel2, r_u, r_v

    INTEGER :: len_fname
    CHARACTER(LEN=4) :: ext
    LOGICAL :: is_binary_restart
    LOGICAL :: runout_exists, mass_exists
    LOGICAL :: timestep_diagnostics_enabled
    INTEGER :: timestep_count
    REAL(wp) :: cfl_dt,dt_pos,cfl_max_a_x,cfl_max_a_y
    CHARACTER(LEN=8) :: timestep_diagnostics_setting

    WRITE(*,*) '---------------------'
    WRITE(*,*) 'IMEX_SfloW2D 2.0'
    WRITE(*,*) '---------------------'

    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    rate = DBLE(cr)

    CALL cpu_time(t1)
    CALL system_clock (st1)

    CALL init_param

    CALL read_param
    timestep_diagnostics_setting=''
    CALL GET_ENVIRONMENT_VARIABLE('IMEX_TIMESTEP_DIAGNOSTICS',timestep_diagnostics_setting)
    timestep_diagnostics_enabled=TRIM(timestep_diagnostics_setting)=='1'
    timestep_count=0

    runout_exists = .FALSE.
    mass_exists = .FALSE.

!$ use_openmp = .true.
!$ print *, "OpenMP program"

    IF ( .NOT. use_openmp) THEN

        PRINT *, "Non-OpenMP simulation"

    ELSE

        ! Check the serial flag
        IF (serial_flag) THEN
            CALL omp_set_num_threads(1)  ! Use only 1 thread
            WRITE(*,*) 'Serial run'

        END IF

!$    n_threads = omp_get_max_threads()
!$    CALL OMP_SET_NUM_THREADS(n_threads)
        WRITE(*,*) 'Parallel run: number of threads used',n_threads

    END IF

    CALL init_grid

    CALL simulation%initialize

    is_binary_restart = .FALSE.
    simulation%runtime%t = t_start
    simulation%runtime%t_output = t_start + dt_output
    simulation%runtime%t_runout = 0.0_wp
    simulation%runtime%t_probes = 0.0_wp

    IF ( restart ) THEN

        ! Calculate string length and extract the last 4 characters
        len_fname = LEN_TRIM(restart_file)

        IF (len_fname >= 4) THEN
            ext = restart_file(len_fname-3:len_fname)
        ELSE
            ext = '    '

        END IF

        IF ( ext .EQ. '.bin' ) THEN
            ! ---------------------------------------------------------
            ! CASE 1: BINARY RESTART (Exact)
            ! ---------------------------------------------------------
            WRITE(*,*) 'Binary file (.bin) detected: Performing full restart.'
            CALL read_restart_file(TRIM(restart_file), simulation%runtime, simulation%stochastic, simulation%state, &
                simulation%domain)
            CALL validate_restart_output_manifest(simulation%runtime)

            is_binary_restart = .TRUE.

            ! A checkpoint restores each OU state. Regenerate only derived
            ! correlation kernels and friction maps for active layer slots.
            DO layer=1,n_layers
                IF (layer==1) THEN
                    IF (.NOT.rheology_l1%stochastic_rheology_flag .AND. &
                        .NOT.simulation%stochastic(layer)%config%stochastic_transport_flag) CYCLE
                ELSE
                    IF (.NOT.rheology_l2%stochastic_rheology_flag .AND. &
                        .NOT.simulation%stochastic(layer)%config%stochastic_transport_flag) CYCLE
                END IF
                IF (simulation%stochastic(layer)%config%length_spatial_corr>0.0_wp) &
                    CALL simulation%stochastic(layer)%generate_kernel()
                CALL simulation%stochastic(layer)%refresh_effective()
            END DO

        ELSE
            ! ---------------------------------------------------------
            ! CASE 2: CLASSIC INITIALIZATION (Python / ASCII)
            ! ---------------------------------------------------------
            WRITE(*,*) 'Standard format detected: Executing read_solution.'
            CALL read_solution(simulation%state, simulation%domain)

        END IF

    ELSE

        IF ( collapsing_volume_flag ) THEN

            CALL collapsing_volume(simulation%state)

        ELSE

            CALL init_empty(simulation%state)

        END IF

    END IF

    CALL release_initialization_fields

    IF ( output_netcdf_flag ) CALL init_netcdf_output

    IF ( radial_source_flag .OR. lateral_source_flag ) CALL init_source

    IF (n_layers==2) THEN
        CALL simulation%domain%check_solve(simulation%state%q(:,1,:,:), &
            simulation%runtime%t,.TRUE.,q_other=simulation%state%q(:,2,:,:))
    ELSE
        CALL simulation%domain%check_solve(simulation%state%q(:,1,:,:), &
            simulation%runtime%t,.TRUE.)
    END IF
    CALL simulation%spatial_operators(1)%prepare_hp_workset(simulation%domain)

    ! For a new run or a classic restart, initialize the stochastic field only
    ! after check_solve has populated solve_cells, j_cent and k_cent. A binary
    ! restart already contains Z and only needs its convolution kernel above.
    IF (.NOT.is_binary_restart) THEN
        DO layer=1,n_layers
            IF (layer==1) THEN
                IF (.NOT.rheology_l1%stochastic_rheology_flag .AND. &
                    .NOT.simulation%stochastic(layer)%config%stochastic_transport_flag) CYCLE
            ELSE
                IF (.NOT.rheology_l2%stochastic_rheology_flag .AND. &
                    .NOT.simulation%stochastic(layer)%config%stochastic_transport_flag) CYCLE
            END IF
            CALL simulation%stochastic(layer)%initialize_steady( &
                simulation%state,simulation%domain, &
                simulation%stochastic(layer)%config%stochastic_transport_flag)
        END DO
    END IF

    IF ( substrate%topo_change_flag ) CALL refresh_topography_geometry

    IF ( verbose_level .GE. 0 ) THEN

        WRITE(*,*)
        WRITE(*,*) '******** START COMPUTATION *********'
        WRITE(*,*)

    END IF

    IF ( is_binary_restart ) THEN
        ! The current step and its history were restored by read_restart_file.
        WRITE(*,*) 'Restarting with saved timestep dt =', simulation%runtime%dt
    ELSE
        ! If it is a new simulation or a legacy restart, we use dt0 from the input
        simulation%runtime%dt = dt0
        simulation%runtime%dt_old = dt0
        simulation%runtime%dt_old_old = dt0

    END IF

    stop_flag = .FALSE.

    IF (is_binary_restart) THEN
        stop_flag = simulation%runtime%stop_requested
        IF (.NOT. stop_flag) simulation%runtime%steady_deadline = t_end
    ELSE
        simulation%runtime%steady_deadline = t_end

    END IF

    IF ( .NOT. is_binary_restart ) simulation%state%vuln_table = .FALSE.
    conversion_failed=.FALSE.

    !$OMP PARALLEL DO private(j,k,p_dyn,i_table,i_thk_lev,i_pdyn_lev,mod_vel2, &
    !$OMP & thickness_exceeds,pressure_exceeds, &
    !$OMP & mod_vel,r_u,r_v,conversion_ok) REDUCTION(.OR.:conversion_failed)

    DO l = 1,simulation%domain%solve_cells

        j = simulation%domain%j_cent(l)
        k = simulation%domain%k_cent(l)

        IF ( simulation%state%q(idx_totMassEqn,1,j,k) .GT. 0.0_wp ) THEN

            CALL conservative_to_primitive(simulation%state%q(1:n_vars,1,j,k), &
                simulation%state%qp(1:n_vars,1,j,k),accepted=conversion_ok,p_dyn=p_dyn)
            IF (.NOT.conversion_ok) conversion_failed=.TRUE.

            IF ( .NOT. is_binary_restart )                                                                             &
                    simulation%state%hmax(j,k) = simulation%state%qp(idx_h,1,j,k)

            r_u = simulation%state%qp(idx_u,1,j,k)
            r_v = simulation%state%qp(idx_v,1,j,k)

            mod_vel2 = r_u**2 + r_v**2
            mod_vel = SQRT( mod_vel2 )

            IF ( ( .NOT. is_binary_restart ) .AND.                                                                     &
                    ( simulation%state%qp(idx_h,1,j,k) .GT. 0.001_wp ) ) THEN

                simulation%state%pdynmax(j,k) = p_dyn
                simulation%state%mod_vel_max(j,k) = mod_vel

            END IF

            IF ( .NOT. is_binary_restart ) THEN

                i_table = 0

                DO i_thk_lev=1,n_thickness_levels

                    thickness_exceeds =                                                                                &
                            ( simulation%state%qp(idx_h,1,j,k) .GE.                                                      &
                            thickness_levels(i_thk_lev) )

                    DO i_pdyn_lev=1,n_dyn_pres_levels

                        pressure_exceeds =                                                                             &
                                ( p_dyn .GE. dyn_pres_levels(i_pdyn_lev) )

                        i_table = i_table + 1

                        simulation%state%vuln_table(i_table,j,k) =                                                     &
                                ( thickness_exceeds .AND.                                                              &
                                pressure_exceeds )

                    END DO

                END DO

            END IF

        ELSE

            simulation%state%qp(1:n_vars,1,j,k) = 0.0_wp
            simulation%state%qp(idx_T,1,j,k) = T_ambient

            IF ( .NOT. is_binary_restart ) THEN
                simulation%state%hmax(j,k) = 0.0_wp
                simulation%state%pdynmax(j,k) = 0.0_wp
                simulation%state%mod_vel_max(j,k) = 0.0_wp

            END IF

        END IF

    END DO

    !$OMP END PARALLEL DO
    IF (conversion_failed) CALL fatal_error('Initial diagnostic state conversion rejected')

    IF (is_binary_restart) THEN
        CALL validate_runout_outputs(simulation%runtime, runout_exists, mass_exists)
        CALL write_param_backup()

    END IF

    CALL open_runout_outputs(simulation%runtime, is_binary_restart, runout_exists, mass_exists)

    IF ( output_runout_flag .AND. (.NOT. is_binary_restart) )                                                          &
            CALL output_runout(simulation%runtime%t, stop_flag, simulation%state, simulation%runtime)

    IF ( output_cons_flag .OR. output_esri_flag .OR.                                                                   &
        output_netcdf_flag ) THEN
        IF ( .NOT. is_binary_restart ) CALL output_solution(                                                           &
        simulation%runtime%t, simulation%state, simulation%runtime)

    END IF

    IF ( ( n_probes .GT. 0 ) .AND. (.NOT. is_binary_restart) )                                                         &
            CALL output_probes(simulation%runtime%t, simulation%state, simulation%runtime)

    IF ( SUM(simulation%state%q(idx_totMassEqn,1,:,:)) .EQ. 0.0_wp )                                                     &
        simulation%runtime%steady_deadline = t_end

    IF ( verbose_level .GE. 0 ) THEN

        WRITE(*,FMT="(A3,F11.4,A5,F9.5,A9,ES11.3E3,A11,ES11.3E3,A9,ES11.3E3,A15,                                       &
        &ES11.3E3,A11,ES11.3E3)")                                                                                      &
            't =',simulation%runtime%t,'dt =',simulation%runtime%dt,                                                   &
            ' mass = ',dx*dy*SUM(simulation%state%q(idx_totMassEqn,1,:,:)) ,                                             &
            ' volume = ',dx*dy*SUM(simulation%state%qp(idx_h,1,:,:)) ,                                                   &
            ' area = ',dx*dy*COUNT(simulation%state%q(idx_totMassEqn,1,:,:).GT.1.D-5) ,                                  &
            ' solid mass = ',dx*dy*SUM(simulation%state%q(5:4+n_solid,1,:,:)),                                           &
            ' runout = ',simulation%runtime%runout%last_distance

    END IF

    CALL cpu_time(t2)
    CALL system_clock (st2)

    DO WHILE ( ( simulation%runtime%t .LT. t_end ) .AND.                                                               &
            ( simulation%runtime%t .LT. simulation%runtime%steady_deadline ) )

        CALL update_param

        IF (n_layers==2) THEN
            CALL simulation%domain%check_solve(simulation%state%q(:,1,:,:), &
                simulation%runtime%t,.FALSE.,q_other=simulation%state%q(:,2,:,:))
        ELSE
            CALL simulation%domain%check_solve(simulation%state%q(:,1,:,:), &
                simulation%runtime%t,.FALSE.)
        END IF

        CALL simulation%spatial_operators(1)%prepare_hp_workset(simulation%domain)
        IF (n_layers==2) CALL simulation%spatial_operators(2)%prepare_hp_workset(simulation%domain)

        ! Assemble the spatial residual once, then reduce CFL and mixture-mass
        ! positivity on the active cells before applying any timestep sources.
        ! Keep the current CFL boundary-trace override unchanged. Lateral Dirichlet
        ! boundary injection is spatial and is not registered as a Lie source.
        CALL simulation%time_integration%prepare_spatial(simulation%state,simulation%runtime%t, &
            simulation%domain,simulation%spatial_operators,cfl_dt,dt_pos,cfl_max_a_x,cfl_max_a_y, &
            boundary_source_enabled=.FALSE., &
            diagnose_max_x_face=timestep_diagnostics_enabled .AND. &
                (timestep_count<20 .OR. MOD(timestep_count+1,100)==0))
        simulation%runtime%dt=MIN(cfl_dt,dt_pos)

        IF ( t_end - simulation%runtime%t_output < 1.0E-7_WP )                                                         &
            simulation%runtime%t_output = t_end
        IF ( t_end - simulation%runtime%t_runout < 1.0E-7_WP )                                                         &
            simulation%runtime%t_runout = t_end
        IF ( t_end - simulation%runtime%t_probes < 1.0E-7_WP )                                                         &
            simulation%runtime%t_probes = t_end

        IF ( simulation%runtime%t+simulation%runtime%dt .GT. t_end )                                                   &
                simulation%runtime%dt = t_end - simulation%runtime%t
        IF ( simulation%runtime%t+simulation%runtime%dt .GT.                                                           &
                simulation%runtime%t_output )                                                                          &
                simulation%runtime%dt = simulation%runtime%t_output - simulation%runtime%t

        IF ( output_runout_flag ) THEN

            IF ( simulation%runtime%t+simulation%runtime%dt .GT.                                                       &
                    simulation%runtime%t_runout )                                                                      &
                    simulation%runtime%dt = simulation%runtime%t_runout - simulation%runtime%t

        END IF

        IF ( n_probes .GT. 0 ) THEN

            IF ( simulation%runtime%t+simulation%runtime%dt .GT.                                                       &
                    simulation%runtime%t_probes )                                                                      &
                    simulation%runtime%dt = simulation%runtime%t_probes - simulation%runtime%t

        END IF

        simulation%runtime%dt = MIN(simulation%runtime%dt,1.1_wp * 0.5_wp *                                            &
                ( simulation%runtime%dt_old + simulation%runtime%dt_old_old ) )

        timestep_count=timestep_count+1
        IF (timestep_diagnostics_enabled .AND. &
            (timestep_count<=20 .OR. MOD(timestep_count,100)==0)) THEN
            WRITE(*,'(A,I0,5(A,ES24.16))') 'TIMESTEP_DIAG step=',timestep_count, &
                ' t=',simulation%runtime%t,' cfl_dt=',cfl_dt, &
                ' accepted_dt=',simulation%runtime%dt, &
                ' max_a_x=',cfl_max_a_x,' max_a_y=',cfl_max_a_y
            WRITE(*,'(A,I0,A,ES24.16)') 'DT_POS_DIAG step=',timestep_count,' dt_pos=',dt_pos
        END IF
        IF (.NOT.ieee_is_finite(simulation%runtime%dt) .OR. &
            simulation%runtime%dt<=0.0_wp .OR. &
            simulation%runtime%t+simulation%runtime%dt<=simulation%runtime%t) THEN
            WRITE(*,'(A,3(A,ES24.16))') 'INVALID_TIMESTEP', &
                ' t=',simulation%runtime%t,' cfl_dt=',cfl_dt, &
                ' accepted_dt=',simulation%runtime%dt
            WRITE(*,'(A,ES24.16)') 'INVALID_TIMESTEP dt_pos=',dt_pos
            CALL fatal_error('Selected timestep does not advance simulation time')
        END IF
        simulation%runtime%dt_old_old = simulation%runtime%dt_old
        simulation%runtime%dt_old = simulation%runtime%dt

        ! The unique pre-step spatial residual remains unchanged until Lie
        ! consumes it. All time-dependent source averages use the final dt.
        CALL simulation%advance_timestep(pre_spatial_primitives_ready=.TRUE., &
            spatial_residual_ready=.TRUE.)

        ! advance_lie has already published a synchronized q/qp pair.
        ! Only the diagnostic dynamic pressure is needed; never reconvert qp.
        !$OMP PARALLEL DO private(j,k,p_dyn,i_table,i_thk_lev,i_pdyn_lev,mod_vel2, &
        !$OMP & thickness_exceeds,pressure_exceeds,mod_vel,r_u,r_v)


        DO l = 1,simulation%domain%solve_cells

            j = simulation%domain%j_cent(l)
            k = simulation%domain%k_cent(l)

            IF ( simulation%state%q(idx_totMassEqn,1,j,k) .GT. 0.0_wp ) THEN

                p_dyn=dynamic_pressure_from_synchronized_state(                                                     &
                    simulation%state%q(idx_totMassEqn,1,j,k),simulation%state%qp(:,1,j,k))

                simulation%state%hmax(j,k) = MAX(simulation%state%hmax(j,k),                                           &
                        simulation%state%qp(idx_h,1,j,k))

                r_u = simulation%state%qp(idx_u,1,j,k)
                r_v = simulation%state%qp(idx_v,1,j,k)

                mod_vel2 = r_u**2 + r_v**2
                mod_vel = SQRT( mod_vel2 )

                IF ( simulation%state%qp(idx_h,1,j,k) .GT. 0.001_wp ) THEN

                    simulation%state%pdynmax(j,k) =                                                                    &
                            MAX(simulation%state%pdynmax(j,k), p_dyn)
                    simulation%state%mod_vel_max(j,k) =                                                                &
                            MAX(simulation%state%mod_vel_max(j,k), mod_vel)

                END IF

                i_table = 0

                DO i_thk_lev=1,n_thickness_levels

                    thickness_exceeds =                                                                                &
                            ( simulation%state%qp(idx_h,1,j,k) .GE.                                                      &
                            thickness_levels(i_thk_lev) )

                    DO i_pdyn_lev=1,n_dyn_pres_levels

                        pressure_exceeds =                                                                             &
                                ( p_dyn .GE. dyn_pres_levels(i_pdyn_lev) )

                        i_table = i_table + 1

                        simulation%state%vuln_table(i_table,j,k) =                                                     &
                                simulation%state%vuln_table(i_table,j,k) .OR.                                          &
                                ( thickness_exceeds .AND.                                                              &
                                pressure_exceeds )

                    END DO

                END DO

            ELSE

                simulation%state%qp(1:n_vars,1,j,k) = 0.0_wp
                simulation%state%qp(idx_T,1,j,k) = T_ambient

            END IF

        END DO

        !$OMP END PARALLEL DO
        ! Conversion acceptance was already checked during the Lie transaction.

        IF ( verbose_level .GE. 0 ) THEN

            WRITE(*,FMT="(A3,F11.4,A5,F9.5,A9,ES11.3E3,A11,ES11.3E3,A9,ES11.3E3,A15,                                   &
            &ES11.3E3,A11,ES11.3E3)")                                                                                  &
                't =',simulation%runtime%t,'dt =',simulation%runtime%dt,                                               &
                ' mass = ',dx*dy*SUM(simulation%state%q(idx_totMassEqn,1,:,:)) ,                                         &
                ' volume = ',dx*dy*SUM(simulation%state%qp(idx_h,1,:,:)) ,                                               &
                ' area = ',dx*dy*COUNT(simulation%state%q(idx_totMassEqn,1,:,:).GT.1.D-7) ,                              &
                ' solid mass = ',dx*dy*SUM(                                                                            &
                simulation%state%q(5:4+n_solid,1,:,:)),                                                                  &
                ' runout = ',simulation%runtime%runout%last_distance

        END IF

        IF ( output_runout_flag ) THEN

            IF ( ( simulation%runtime%t .GE. simulation%runtime%t_runout ) .OR.                                        &
                    ( simulation%runtime%t .GE.                                                                        &
                        simulation%runtime%steady_deadline ) ) THEN

                stop_flag_old = stop_flag

                IF ( output_runout_flag ) CALL output_runout(                                                          &
                        simulation%runtime%t, stop_flag, simulation%state, simulation%runtime)

                IF ( ( stop_flag ) .AND. (.NOT.stop_flag_old) ) THEN

                    simulation%runtime%steady_deadline =                                                               &
                            MIN(t_end, simulation%runtime%t_output)

                END IF

            END IF

        END IF

        IF ( n_probes .GT. 0 ) THEN

            IF ( simulation%runtime%t .GE. simulation%runtime%t_probes )                                               &
                    CALL output_probes(simulation%runtime%t, simulation%state, simulation%runtime)

        END IF

        IF ( ( simulation%runtime%t .GE. simulation%runtime%t_output ) .OR.                                            &
                ( simulation%runtime%t .GE. t_end ) ) THEN

            CALL cpu_time(t3)
            CALL system_clock(st3)

            IF ( verbose_level .GE. 0 ) THEN

                WRITE(*,*) 'Time taken by iterations is',t3-t2,'seconds'
                WRITE(*,*) 'Elapsed real time = ', DBLE( st3 - st2 ) / rate,'seconds'

            END IF

            IF ( output_cons_flag .OR. output_esri_flag .OR. output_netcdf_flag ) THEN

                CALL output_solution(simulation%runtime%t, simulation%state, simulation%runtime)
                CALL write_restart_file(output_file_path('restart', 'restart.bin'), simulation%runtime,                &
                    simulation%stochastic, simulation%state, simulation%domain)

            END IF

        END IF

    END DO

    CALL simulation%finalize

    IF ( output_netcdf_flag ) CALL close_netcdf

    CALL close_units

    CALL cpu_time(t3)
    CALL system_clock(st3)

    WRITE(*,*) 'Total time taken by the code is',t3-t1,'seconds'
    WRITE(*,*) 'Total elapsed real time is', DBLE( st3 - st1 ) / rate,'seconds'

END PROGRAM IMEX_SfloW2D
