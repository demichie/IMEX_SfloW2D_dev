PROGRAM test_model_layout_reject

  USE model_layout_2d, ONLY : model_layout_type

  IMPLICIT NONE

  TYPE(model_layout_type) :: layout

  ! This call must terminate with the explicit unsupported-layer diagnostic for three layers.
  CALL layout%initialize(3, 4, 4, 0, 0, 0, 0)

  ERROR STOP 'N_LAYERS=3 was accepted unexpectedly'

END PROGRAM test_model_layout_reject
