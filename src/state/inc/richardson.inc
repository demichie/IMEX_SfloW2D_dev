!> @file richardson.inc
!> @brief Compute the depth-averaged Richardson number with a large-value cap.
    speed2=u**2+v**2
    ri=0.0_wp
    IF (REAL(speed2,wp)>0.0_wp) THEN
        IF (REAL(red_grav*h/1.0E15_wp,wp)>REAL(speed2,wp)) THEN
            ri=1.0E15_wp
        ELSE
            ri=red_grav*h/speed2
        END IF
    END IF
