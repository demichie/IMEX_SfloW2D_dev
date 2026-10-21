PROGRAM test_rheology_model_names

  USE parameters_2d, ONLY : rheology_id_none, rheology_id_voellmy,             &
       rheology_id_plastic, rheology_id_lava, rheology_id_lahar,               &
       rheology_id_empirical, rheology_id_quadratic, rheology_id_froude,       &
       rheology_id_lucas, rheology_id_mui, rheology_id_mui_voellmy
  USE rheology_2d, ONLY : resolve_rheology_model, rheology_name_from_id

  IMPLICIT NONE

  INTEGER, PARAMETER :: model_count = 11
  CHARACTER(LEN=32), PARAMETER :: names(model_count) = [CHARACTER(LEN=32) ::  &
       'none', 'voellmy', 'plastic', 'lava', 'lahar', 'empirical',            &
       'quadratic', 'froude', 'lucas', 'mui', 'mui_voellmy']
  INTEGER, PARAMETER :: model_ids(model_count) = [                            &
       rheology_id_none, rheology_id_voellmy, rheology_id_plastic,            &
       rheology_id_lava, rheology_id_lahar, rheology_id_empirical,            &
       rheology_id_quadratic, rheology_id_froude, rheology_id_lucas,          &
       rheology_id_mui, rheology_id_mui_voellmy]

  CHARACTER(LEN=32) :: name
  INTEGER :: model, model_id

  DO model = 1, model_count
     name = names(model)
     IF (model .EQ. 2) name = 'VOELLMY'
     CALL resolve_rheology_model(name, model_id)
     IF (model_id .NE. model_ids(model)) ERROR STOP 'incorrect rheology model ID'
     IF (TRIM(name) .NE. TRIM(names(model))) ERROR STOP 'name not canonicalized'
     IF (TRIM(rheology_name_from_id(model_id)) .NE. TRIM(names(model)))       &
          ERROR STOP 'incorrect reverse rheology name'
  END DO

  name = 'unknown-law'
  CALL resolve_rheology_model(name, model_id)
  IF (model_id .NE. -1) ERROR STOP 'unknown rheology name was accepted'

END PROGRAM test_rheology_model_names
