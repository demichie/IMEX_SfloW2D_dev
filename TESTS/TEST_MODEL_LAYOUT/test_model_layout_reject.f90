PROGRAM test_model_layout_reject

  USE model_layout_2d, ONLY : model_layout_type

  IMPLICIT NONE

  TYPE(model_layout_type) :: layout

  ! This call must terminate with the explicit unsupported-layer diagnostic.
  CALL layout%initialize(2, 4, 4, 0, 0, 0, 0, .FALSE.)

  ERROR STOP 'N_LAYERS=2 was accepted unexpectedly'

END PROGRAM test_model_layout_reject
