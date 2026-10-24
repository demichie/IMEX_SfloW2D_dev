PROGRAM test_rheology_model_names
  USE source_catalog_2d, ONLY : canonical_rheology_source_name, &
       source_term_descriptor_by_name, source_term_type
  IMPLICIT NONE
  CHARACTER(LEN=32), PARAMETER :: names(12)=[CHARACTER(LEN=32) :: &
       'none','voellmy','plastic','lava','lahar','empirical', &
       'quadratic','froude','lucas','mui','mui_voellmy','newmodel']
  TYPE(source_term_type) :: descriptor
  CHARACTER(LEN=64) :: source_name
  INTEGER :: model
  DO model=1,SIZE(names)
    source_name=canonical_rheology_source_name(names(model))
    IF (TRIM(names(model))=='none') THEN
      IF (LEN_TRIM(source_name)/=0) ERROR STOP 'none is not disabled'
    ELSE
      IF (TRIM(source_name)/='rheology_'//TRIM(names(model))) &
           ERROR STOP 'incorrect canonical rheology name'
      descriptor=source_term_descriptor_by_name(source_name)
      IF (descriptor%name/=source_name) ERROR STOP 'missing canonical source'
      IF (TRIM(names(model))=='newmodel') THEN
        IF (.NOT.descriptor%has_newton .OR. descriptor%has_explicit .OR. &
             descriptor%has_analytic) ERROR STOP 'incorrect newmodel capabilities'
      END IF
    END IF
  END DO
  IF (LEN_TRIM(canonical_rheology_source_name('unknown-law'))/=0) &
       ERROR STOP 'unknown rheology accepted'
END PROGRAM test_rheology_model_names
