! Quantitative D17.20 test through the production rheology dispatcher.
PROGRAM test_d17_production_rheology
  USE parameters_2d, ONLY : wp
  USE rheology_2d, ONLY : eval_rheology_implicit,eval_rheology_semi_implicit
  IMPLICIT NONE
  REAL(wp), PARAMETER :: drag=0.1_wp, epsilon_cs=1.0E-24_wp
  REAL(wp) :: density_real,velocity(2),reference(2),jacobian(2,2)
  REAL(wp) :: si(2),solid(1),step
  COMPLEX(wp) :: density,vel(2),rate(2),rate_reference(2)
  INTEGER :: trial,j
  solid=0.0_wp
  DO trial=1,4
    SELECT CASE(trial)
    CASE(1)
      density_real=1000.0_wp;velocity=[2.0_wp,-3.0_wp]
    CASE(2)
      density_real=2500.0_wp;velocity=[-0.001_wp,0.004_wp]
    CASE(3)
      density_real=850.0_wp;velocity=[0.0_wp,0.0_wp]
    CASE DEFAULT
      density_real=1450.0_wp;velocity=[-8.0_wp,-6.0_wp]
    END SELECT
    density=CMPLX(density_real,0.0_wp,wp)
    vel=CMPLX(velocity,0.0_wp,wp)
    CALL implicit(vel,rate)
    reference=-drag*density_real*velocity
    CALL compare(REAL(rate,wp),reference,'production dispatch rate')
    CALL eval_rheology_semi_implicit('rheology_newmodel',2.0_wp, &
         velocity(1),velocity(2),SQRT(SUM(velocity**2)),SQRT(SUM(velocity**2)), &
         density_real,8.0_wp,9.81_wp,solid,.FALSE., &
         0.0_wp,0.0_wp,0.0_wp,.FALSE.,0.0_wp,.FALSE.,0.0_wp,si)
    CALL compare(si,[0.0_wp,0.0_wp],'production SI zero')
    DO j=1,2
      vel=CMPLX(velocity,0.0_wp,wp)
      vel(j)=vel(j)+CMPLX(0.0_wp,epsilon_cs,wp)
      CALL implicit(vel,rate)
      jacobian(:,j)=AIMAG(rate)/epsilon_cs
      reference=0.0_wp
      reference(j)=-drag*density_real
      CALL compare(jacobian(:,j),reference,'production complex-step')
    END DO
    ! Compare the production momentum rate with the independent implicit
    ! momentum update q_new=q_old/(1+drag*dt), for several stage intervals.
    step=0.2_wp*REAL(trial,wp)
    vel=CMPLX(velocity/(1.0_wp+drag*step),0.0_wp,wp)
    CALL implicit(vel,rate)
    rate_reference=CMPLX(density_real*(REAL(vel,wp)-velocity)/step,0.0_wp,wp)
    CALL compare(REAL(rate,wp),REAL(rate_reference,wp),'implicit balance')
  END DO
  PRINT *, 'PASS D17 production newmodel: SI, rate, complex-step and implicit balance'
CONTAINS
  SUBROUTINE implicit(v,r)
    COMPLEX(wp), INTENT(IN) :: v(2)
    COMPLEX(wp), INTENT(OUT) :: r(2)
    COMPLEX(wp) :: speed,horizontal,squared
    horizontal=SQRT(v(1)**2+v(2)**2)
    speed=horizontal
    squared=speed**2
    CALL eval_rheology_implicit('rheology_newmodel',CMPLX(2.0_wp,0.0_wp,wp), &
         CMPLX(300.0_wp,0.0_wp,wp),density,1.0_wp/density, &
         v(1),v(2),speed,horizontal,squared, &
         CMPLX(8.0_wp,0.0_wp,wp),CMPLX(solid,0.0_wp,wp), &
         .TRUE.,.FALSE.,0.0_wp,r)
  END SUBROUTINE implicit
  SUBROUTINE compare(actual,expected,label)
    REAL(wp), INTENT(IN) :: actual(2),expected(2)
    CHARACTER(*), INTENT(IN) :: label
    IF (ANY(ABS(actual-expected)>3.0E-11_wp*MAX(1.0_wp,MAXVAL(ABS(expected))))) THEN
      PRINT *, 'FAIL: ',label,actual,expected
      ERROR STOP 1
    END IF
  END SUBROUTINE compare
END PROGRAM test_d17_production_rheology
