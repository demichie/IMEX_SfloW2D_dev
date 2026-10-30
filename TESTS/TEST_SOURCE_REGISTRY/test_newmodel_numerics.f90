! D17.20 independent quantitative fixture for the Newton-only extension law.
! Compare against the exact linear drag rate, Jacobian and backward-Euler map.
PROGRAM test_newmodel_numerics
  USE parameters_2d, ONLY : wp
  USE source_catalog_2d, ONLY : source_term_descriptor_by_name, &
       canonical_rheology_source_name, source_term_type, rheology_name_has_role, &
       SOURCE_ROLE_IMPLICIT, SOURCE_ROLE_SEMI_IMPLICIT
  USE state_conversion_2d, ONLY : physical_state_type, complex_physical_state_type
  USE rheology_newmodel_2d, ONLY : eval_rheology_newmodel_implicit, &
       eval_rheology_newmodel_semi_implicit
  IMPLICIT NONE
  TYPE(source_term_type) :: term
  TYPE(physical_state_type) :: rs
  TYPE(complex_physical_state_type) :: cs
  REAL(wp), PARAMETER :: k=0.1_wp, eps=1.0E-24_wp
  REAL(wp) :: rho, u(2), q_before(2), q_after(2), expected(2)
  REAL(wp) :: dt, analytic_jac(2,2), cs_jac(2,2), residual(2), &
       si_x,si_y, det, inverse_jac(2,2)
  COMPLEX(wp) :: uc(2), rate(2), density
  INTEGER :: i,j,case_number
  term=source_term_descriptor_by_name('rheology_newmodel')
  IF (TRIM(term%name)/='rheology_newmodel') ERROR STOP 'newmodel missing'
  IF (TRIM(canonical_rheology_source_name('newmodel'))/=TRIM(term%name)) &
       ERROR STOP 'noncanonical newmodel lookup'
  IF (.NOT.term%has_newton .OR. term%has_explicit .OR. term%has_analytic) &
       ERROR STOP 'newmodel capabilities incorrect'
  IF (.NOT.rheology_name_has_role(term%name,SOURCE_ROLE_IMPLICIT)) &
       ERROR STOP 'newmodel implicit component not selected'
  IF (rheology_name_has_role(term%name,SOURCE_ROLE_SEMI_IMPLICIT)) &
       ERROR STOP 'newmodel SI is intentionally zero'
  CALL eval_rheology_newmodel_semi_implicit(rs,si_x,si_y)
  IF (ABS(si_x)+ABS(si_y)>TINY(1.0_wp)) ERROR STOP 'newmodel SI nonzero'
  DO case_number=1,4
    SELECT CASE(case_number)
    CASE(1)
      rho=1000.0_wp;u=[3.0_wp,-4.0_wp];dt=0.75_wp
    CASE(2)
      rho=1800.0_wp;u=[0.0_wp,0.0_wp];dt=0.25_wp
    CASE(3)
      rho=850.0_wp;u=[-1.0E-5_wp,2.0E-5_wp];dt=10.0_wp
    CASE(4)
      rho=2750.0_wp;u=[-8.0_wp,-6.0_wp];dt=0.001_wp
    END SELECT
    density=CMPLX(rho,0.0_wp,wp)
    uc=CMPLX(u,0.0_wp,wp)
    cs%u=uc(1)
    cs%v=uc(2)
    cs%rho_m=density
    CALL eval_rheology_newmodel_implicit(cs,rate(1),rate(2))
    expected=-k*rho*u
    CALL close_vector(REAL(rate,wp),expected,'implicit rate')
    analytic_jac=0.0_wp
    analytic_jac(1,1)=-k*rho
    analytic_jac(2,2)=-k*rho
    DO j=1,2
      uc=CMPLX(u,0.0_wp,wp)
      uc(j)=uc(j)+CMPLX(0.0_wp,eps,wp)
      cs%u=uc(1)
    cs%v=uc(2)
    cs%rho_m=density
    CALL eval_rheology_newmodel_implicit(cs,rate(1),rate(2))
      cs_jac(:,j)=AIMAG(rate)/eps
    END DO
    DO j=1,2
      CALL close_vector(cs_jac(:,j),analytic_jac(:,j),'complex-step column')
    END DO
    ! Independent solution of the conservative 2x2 implicit momentum system:
    ! (I - dt * dS/dq) q_after = q_before, with q = rho * u.
    q_before=rho*u
    q_after=q_before
    DO i=1,5
      uc=CMPLX(q_after/rho,0.0_wp,wp)
      cs%u=uc(1)
    cs%v=uc(2)
    cs%rho_m=density
    CALL eval_rheology_newmodel_implicit(cs,rate(1),rate(2))
      residual=q_after-q_before-dt*REAL(rate,wp)
      ! Jacobian of the momentum residual (not of the velocity rate).
      inverse_jac=0.0_wp
      det=(1.0_wp+k*dt)**2
      inverse_jac(1,1)=(1.0_wp+k*dt)/det
      inverse_jac(2,2)=(1.0_wp+k*dt)/det
      q_after=q_after-MATMUL(inverse_jac,residual)
    END DO
    expected=q_before/(1.0_wp+k*dt)
    CALL close_vector(q_after,expected,'Newton 2x2 versus exact step')
    uc=CMPLX(q_after/rho,0.0_wp,wp)
    cs%u=uc(1)
    cs%v=uc(2)
    cs%rho_m=density
    CALL eval_rheology_newmodel_implicit(cs,rate(1),rate(2))
    residual=q_after-q_before-dt*REAL(rate,wp)
    CALL close_vector(residual,[0.0_wp,0.0_wp],'Newton final residual')
  END DO
  WRITE(*,'(A)') 'PASS D17.20 newmodel rate, SI zero, complex-step 2x2, Newton exact map (4 cases)'
CONTAINS
  SUBROUTINE close_vector(actual,reference,label)
    REAL(wp), INTENT(IN) :: actual(2),reference(2)
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp) :: tolerance
    tolerance=3.0E-11_wp*MAX(1.0_wp,MAXVAL(ABS(reference)))
    IF (ANY(ABS(actual-reference)>tolerance)) THEN
      WRITE(*,*) 'FAIL ',label,actual,reference
      ERROR STOP 1
    END IF
  END SUBROUTINE close_vector
END PROGRAM test_newmodel_numerics
