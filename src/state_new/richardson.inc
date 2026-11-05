    speed2=u**2+v**2
    ri=0.0_wp
    IF (REAL(speed2,wp)>0.0_wp) THEN
        ri=red_grav*h/speed2
        IF (REAL(ri,wp)>1.0E15_wp) ri=CMPLX(1.0E15_wp,0.0_wp,wp)
    END IF
