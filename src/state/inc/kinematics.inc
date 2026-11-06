!> @file kinematics.inc
!> @brief Compute local velocity magnitudes and bed-curvature acceleration.
    state%w=0.0_wp
    IF (slope_correction_flag) state%w=state%u*bx+state%v*by
    state%mod_vel_hor=SQRT(state%u**2+state%v**2)
    state%mod_vel2=state%u**2+state%v**2+state%w**2
    state%mod_vel=SQRT(state%mod_vel2)
    state%curvature_acceleration=0.0_wp
    IF (curvature_term_flag) state%curvature_acceleration=state%u**2*bxx+                                              &
        2.0_wp*state%u*state%v*bxy+state%v**2*byy
