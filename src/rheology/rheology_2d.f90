!********************************************************************************
!> \brief Rheology model selection and naming.
!********************************************************************************

MODULE rheology_2d

  USE parameters_2d, ONLY : rheology_id_none, rheology_id_voellmy, rheology_id_plastic, &
       rheology_id_lava, rheology_id_lahar, rheology_id_empirical, rheology_id_quadratic,  &
       rheology_id_froude, rheology_id_lucas, rheology_id_mui, rheology_id_mui_voellmy

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: resolve_rheology_model, rheology_name_from_id

CONTAINS

  !******************************************************************************
  !> \brief Resolve a user-supplied rheology name to its internal model ID.
  !>
  !> \param[in] name Input model name.
  !> \param[out] model_id Internal model ID, or -1 for an unsupported name.
  !******************************************************************************

  SUBROUTINE resolve_rheology_model( name, model_id )

    CHARACTER(LEN=*), INTENT(INOUT) :: name
    INTEGER, INTENT(OUT) :: model_id

    CHARACTER(LEN=LEN(name)) :: canonical_name
    INTEGER :: character_index, character_code

    canonical_name = ADJUSTL(name)
    DO character_index = 1, LEN_TRIM(canonical_name)
       character_code = IACHAR(canonical_name(character_index:character_index))
       IF ( character_code .GE. IACHAR('A') .AND. &
            character_code .LE. IACHAR('Z') ) THEN
          canonical_name(character_index:character_index) =                  &
               ACHAR(character_code + IACHAR('a') - IACHAR('A'))
       END IF
    END DO
    name = TRIM(canonical_name)

    SELECT CASE (TRIM(name))
    CASE ('none')
       model_id = rheology_id_none
    CASE ('voellmy')
       model_id = rheology_id_voellmy
    CASE ('plastic')
       model_id = rheology_id_plastic
    CASE ('lava')
       model_id = rheology_id_lava
    CASE ('lahar')
       model_id = rheology_id_lahar
    CASE ('empirical')
       model_id = rheology_id_empirical
    CASE ('quadratic')
       model_id = rheology_id_quadratic
    CASE ('froude')
       model_id = rheology_id_froude
    CASE ('lucas')
       model_id = rheology_id_lucas
    CASE ('mui')
       model_id = rheology_id_mui
    CASE ('mui_voellmy')
       model_id = rheology_id_mui_voellmy
    CASE DEFAULT
       model_id = -1
    END SELECT

  END SUBROUTINE resolve_rheology_model


  !******************************************************************************
  !> \brief Return the canonical input name corresponding to an internal ID.
  !>
  !> \param[in] model_id Internal rheology model ID.
  !> \return Canonical name, or 'unknown' for an unsupported ID.
  !******************************************************************************

  PURE FUNCTION rheology_name_from_id( model_id ) RESULT(name)

    INTEGER, INTENT(IN) :: model_id
    CHARACTER(LEN=32) :: name

    SELECT CASE (model_id)
    CASE (rheology_id_none)
       name = 'none'
    CASE (rheology_id_voellmy)
       name = 'voellmy'
    CASE (rheology_id_plastic)
       name = 'plastic'
    CASE (rheology_id_lava)
       name = 'lava'
    CASE (rheology_id_lahar)
       name = 'lahar'
    CASE (rheology_id_empirical)
       name = 'empirical'
    CASE (rheology_id_quadratic)
       name = 'quadratic'
    CASE (rheology_id_froude)
       name = 'froude'
    CASE (rheology_id_lucas)
       name = 'lucas'
    CASE (rheology_id_mui)
       name = 'mui'
    CASE (rheology_id_mui_voellmy)
       name = 'mui_voellmy'
    CASE DEFAULT
       name = 'unknown'
    END SELECT

  END FUNCTION rheology_name_from_id

END MODULE rheology_2d
