! D17.20: quantitative production mass-exchange balance regression.
PROGRAM test_d17_mass_balance
  USE parameters_2d
  USE constitutive_parameters_2d
  USE mass_exchange_balance_2d, ONLY : assemble_mass_exchange_balance
  IMPLICIT NONE
  REAL(wp) :: er(2),dep(2),rates(6),topo
  REAL(wp) :: solid_erosion,solid_deposition,carrier_erosion,carrier_loss
  REAL(wp) :: u,v,temp,carrier_rho,carrier_heat,air,bed_temp
  REAL(wp) :: expected(6),expected_topo,dt
  INTEGER :: scenario
  n_solid=2;n_eqns=6;n_vars=6;n_add_gas=0
  idx_solidEqn_first=5;idx_solidEqn_last=6
  idx_stochEqn=0;idx_poreEqn=0
  pore_pressure_flag=.FALSE.;stoch_transport_flag=.FALSE.
  liquid_vaporization_flag=.FALSE.;gas_flag=.FALSE.
  rho_a_amb=1.2_wp;sp_heat_a=1005.0_wp;T_ambient=295.0_wp
  rho_c_sub=1100.0_wp;erodible_porosity=0.35_wp
  ALLOCATE(rho_s(2),sp_heat_s(2))
  rho_s=[2100.0_wp,2700.0_wp]
  sp_heat_s=[850.0_wp,1050.0_wp]
  u=3.0_wp;v=-4.0_wp;temp=600.0_wp;bed_temp=310.0_wp
  carrier_rho=1020.0_wp;carrier_heat=3900.0_wp
  dt=0.4_wp
  DO scenario=1,4
    SELECT CASE(scenario)
    CASE(1)
      er=[0.0_wp,0.0_wp];dep=[0.01_wp,0.02_wp]
      carrier_erosion=0.0_wp;carrier_loss=0.002_wp;air=0.0_wp
    CASE(2)
      er=[0.03_wp,0.01_wp];dep=0.0_wp
      carrier_erosion=0.015_wp;carrier_loss=0.0_wp;air=0.002_wp
    CASE(3)
      er=[0.008_wp,0.004_wp];dep=[0.006_wp,0.012_wp]
      carrier_erosion=0.004_wp;carrier_loss=0.003_wp;air=0.001_wp
    CASE DEFAULT
      er=0.0_wp;dep=0.0_wp
      carrier_erosion=0.0_wp;carrier_loss=0.0_wp;air=0.0_wp
    END SELECT
    CALL assemble_mass_exchange_balance(0,u,v,temp,carrier_rho,carrier_heat, &
         bed_temp,air,0.0_wp,0.0_wp,er,dep,carrier_erosion,carrier_loss, &
         rates,topo)
    solid_erosion=SUM(rho_s*er)
    solid_deposition=SUM(rho_s*dep)
    expected=0.0_wp
    expected(1)=solid_erosion-solid_deposition+rho_c_sub*carrier_erosion &
         -carrier_rho*carrier_loss+rho_a_amb*air
    expected(2)=-u*(solid_deposition+carrier_rho*carrier_loss)
    expected(3)=-v*(solid_deposition+carrier_rho*carrier_loss)
    expected(4)=bed_temp*(SUM(rho_s*sp_heat_s*er)+ &
         rho_c_sub*carrier_heat*carrier_erosion) &
         -temp*(SUM(rho_s*sp_heat_s*dep)+carrier_rho*carrier_heat*carrier_loss) &
         +T_ambient*sp_heat_a*rho_a_amb*air
    expected(5:6)=rho_s*(er-dep)
    expected_topo=(SUM(dep)-SUM(er))/(1.0_wp-erodible_porosity)
    CALL compare(rates,expected,'mass/momentum/energy/solids')
    IF (ABS(topo-expected_topo)>1.0E-12_wp*MAX(1.0_wp,ABS(expected_topo))) &
         ERROR STOP 'bed volume balance'
    ! The bath and flow inventory exchange must close for solids.
    IF (ANY(ABS(dt*rates(5:6)+dt*rho_s*(dep-er))>1.0E-12_wp)) &
         ERROR STOP 'solid flow-bed exchange'
    IF (scenario==4 .AND. ANY(ABS(rates)>0.0_wp)) ERROR STOP 'zero exchange'
  END DO
  PRINT *, 'PASS D17 production mass-exchange conservation (4 scenarios)'
CONTAINS
  SUBROUTINE compare(actual,reference,label)
    REAL(wp),INTENT(IN)::actual(6),reference(6)
    CHARACTER(*),INTENT(IN)::label
    IF (ANY(ABS(actual-reference)>2.0E-12_wp*MAX(1.0_wp,MAXVAL(ABS(reference))))) THEN
      PRINT *, 'FAIL ',label,actual,reference
      ERROR STOP 1
    END IF
  END SUBROUTINE compare
END PROGRAM test_d17_mass_balance
