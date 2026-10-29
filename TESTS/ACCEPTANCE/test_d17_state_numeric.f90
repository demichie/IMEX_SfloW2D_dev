PROGRAM test_d17_state_numeric
  USE parameters_2d, ONLY : wp, n_vars, n_solid, n_add_gas, idx_h, idx_hu, idx_hv, idx_T, idx_u, idx_v, &
       idx_solid_first, idx_solid_last, idx_add_gas_first, idx_add_gas_last, &
       gas_flag, liquid_flag
  USE constitutive_parameters_2d, ONLY : T_ambient
  USE state_2d, ONLY : finalize_imex_cell_state, initialize_dry_physical_cell
  IMPLICIT NONE
  REAL(wp) :: q(4), qp(6), previous(4)
  n_vars=4
  n_solid=0
  n_add_gas=0
  idx_solid_first=5
  idx_solid_last=4
  idx_add_gas_first=5
  idx_add_gas_last=4
  idx_h=1
  idx_hu=2
  idx_hv=3
  idx_T=4
  idx_u=5
  idx_v=6
  gas_flag=.FALSE.
  liquid_flag=.FALSE.
  T_ambient=293.15_wp
  qp=-77.0_wp
  CALL initialize_dry_physical_cell(qp)
  CALL expect('dry initializer', ALL(qp(1:3)==0.0_wp) .AND. &
       qp(4)==T_ambient .AND. ALL(qp(5:6)==0.0_wp))
  previous=0.0_wp
  q=0.0_wp
  qp=-77.0_wp
  CALL finalize_imex_cell_state(q,qp,previous,.FALSE.)
  CALL expect('zero mass remains zero',ALL(q==0.0_wp))
  CALL expect('zero mass has zero velocities',ALL(qp(2:3)==0.0_wp) &
       .AND. ALL(qp(5:6)==0.0_wp))
  q=[-1.0e-8_wp,3.0_wp,4.0_wp,100.0_wp]
  qp=-77.0_wp
  CALL finalize_imex_cell_state(q,qp,previous,.FALSE.)
  CALL expect('small negative mass is zeroed',ALL(q==0.0_wp))
  q=0.0_wp
  qp=-77.0_wp
  CALL finalize_imex_cell_state(q,qp,previous,.TRUE.)
  CALL expect('nodata conservative mask',ALL(q==0.0_wp))
  PRINT *, 'PASS: D17 numeric dry-cell finalization and initialization'
CONTAINS
  SUBROUTINE expect(label,ok)
    CHARACTER(*),INTENT(IN)::label
    LOGICAL,INTENT(IN)::ok
    IF (.NOT.ok) THEN
       PRINT *, 'FAIL: ',label
       ERROR STOP 1
    END IF
  END SUBROUTINE expect
END PROGRAM test_d17_state_numeric
