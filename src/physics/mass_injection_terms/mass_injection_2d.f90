!> Assemble complete conservative rates from volumetric injection processes.
MODULE mass_injection_2d
  USE parameters_2d, ONLY : wp, n_eqns
  USE state_conversion_2d, ONLY : physical_state_type
  USE mass_bottom_circular_injection_2d, ONLY : eval_bottom_circular_injection_explicit
  USE mass_bottom_fissural_injection_2d, ONLY : eval_bottom_fissural_injection_explicit
  USE mass_radial_injection_2d, ONLY : eval_radial_injection_explicit
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: eval_mass_injection_explicit
CONTAINS
  SUBROUTINE eval_mass_injection_explicit(state,time,cell_fract,cell_fissure_fract, &
       lat_arc_perim,lat_n_x,lat_n_y,cell_area,use_circular,use_fissural,use_radial,rate)
    TYPE(physical_state_type), INTENT(IN) :: state
    REAL(wp), INTENT(IN) :: time,cell_fract,cell_fissure_fract(:)
    REAL(wp), INTENT(IN) :: lat_arc_perim,lat_n_x,lat_n_y,cell_area
    LOGICAL, INTENT(IN) :: use_circular,use_fissural,use_radial
    REAL(wp), INTENT(OUT) :: rate(n_eqns)
    REAL(wp) :: process_rate(n_eqns)
    rate=0.0_wp
    IF (use_circular) THEN
      CALL eval_bottom_circular_injection_explicit(state,time,cell_fract,process_rate)
      rate=rate+process_rate
    END IF
    IF (use_fissural) THEN
      CALL eval_bottom_fissural_injection_explicit(state,time,cell_fissure_fract,process_rate)
      rate=rate+process_rate
    END IF
    IF (use_radial) THEN
      CALL eval_radial_injection_explicit(state,time,lat_arc_perim,lat_n_x,lat_n_y, &
           cell_area,process_rate)
      rate=rate+process_rate
    END IF
  END SUBROUTINE eval_mass_injection_explicit
END MODULE mass_injection_2d
