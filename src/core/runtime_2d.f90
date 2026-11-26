!********************************************************************************
!> \brief Runtime state ownership
!
!> This module owns the mutable time state of a simulation.
!>
!> Groups the current simulation time and timestep history independently of physical state arrays.
!> Binary restarts preserve this history together with the output schedule.
!********************************************************************************

MODULE runtime_2d

    USE parameters_2d, ONLY : wp

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: runtime_state_type
    PUBLIC :: runout_state_type

    !> \brief Restartable runout history, independent of formatted output files.
    !> Origins and initial extents refer to the original initial condition, not
    !> the wet footprint at the restart time. CSV tails guard append operations.
    TYPE :: runout_state_type
        !> True after runout reference geometry and history are initialized.
        LOGICAL :: initialized = .FALSE.
        !> Original runout-reference center coordinates [m] in x and y.
        REAL(wp) :: origin_x = 0.0_wp, origin_y = 0.0_wp
        !> Initial radial extent relative to the runout reference [m].
        REAL(wp) :: initial_distance = 0.0_wp
        !> Initial extents in the x and y directions [m].
        REAL(wp) :: initial_distance_x = 0.0_wp, initial_distance_y = 0.0_wp
        !> Previous center-of-mass position in x and y [m].
        REAL(wp) :: previous_center_x = 0.0_wp, previous_center_y = 0.0_wp
        !> Time of the last written runout record [s].
        REAL(wp) :: last_time = 0.0_wp
        !> Last and maximum recorded runout distances [m].
        REAL(wp) :: last_distance = -1.0_wp, maximum_distance = 0.0_wp
        !> Previous formatted runout and mass-center rows for restart-safe appending.
        CHARACTER(LEN=128) :: runout_row = '', mass_center_row = ''
    END TYPE runout_state_type

    !> \brief Mutable clock, output schedule, timestep history and runout state for one run.
    TYPE :: runtime_state_type

        !> Current simulation time
        REAL(wp) :: t = 0.0_wp

        !> Current time step
        REAL(wp) :: dt = 0.0_wp

        !> Time steps used by the temporal smoothing limiter
        REAL(wp) :: dt_old = 0.0_wp
        !> Second-previous accepted timestep for temporal smoothing [s].
        REAL(wp) :: dt_old_old = 0.0_wp

        !> Next output times, restored with the binary restart state.
        REAL(wp) :: t_output = 0.0_wp
        !> Next scheduled runout-record time [s].
        REAL(wp) :: t_runout = 0.0_wp
        !> Next scheduled probe-output time [s].
        REAL(wp) :: t_probes = 0.0_wp

        !> Runout/steady-stop state saved even when no CSV files are available.
        TYPE(runout_state_type) :: runout
        !> Whether runout-based stopping has been requested.
        LOGICAL :: stop_requested = .FALSE.
        !> Latest simulation time permitted by the steady-stop criterion [s].
        REAL(wp) :: steady_deadline = 0.0_wp

    END TYPE runtime_state_type

END MODULE runtime_2d
