!> @brief Cell-local geometrical contacts for one or two depth-averaged material layers.
!> @details Contact classification uses the existing canonical wet/dry thickness tolerance.
!>          This module adds no exchange law and does not activate two-layer execution.
MODULE layer_contacts_2d
    USE parameters_2d, ONLY : wp, dry_thickness_tolerance
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: layer_contacts_type, evaluate_layer_contacts

    !> @brief Exposure of each layer to atmosphere and terrain, and L1/L2 interface contact.
    !> @details Only active column layers can be exposed; the two slots are L1 and L2.
    TYPE :: layer_contacts_type
        !> True for the uppermost wet layer of a material column.
        LOGICAL :: atmosphere_exposed(2)=.FALSE.
        !> True for the lowermost wet layer of a material column.
        LOGICAL :: terrain_exposed(2)=.FALSE.
        !> True only when both L1 and L2 are wet; cold-start transfers are not prohibited.
        LOGICAL :: interlayer_contact=.FALSE.
    END TYPE layer_contacts_type

CONTAINS

    !> @brief Determine instantaneous column contacts from candidate layer thicknesses.
    !> @details Accept one or two nonnegative finite thicknesses (m), ordered bottom-up.
    !>          A layer is wet iff h > dry_thickness_tolerance. Invalid inputs return
    !>          valid=.FALSE. and all contacts false, without modifying any physical state.
    !>          The caller must use current candidate thicknesses, not published old state.
    !> @param[in] thickness Candidate physical layer thicknesses (m), L1 then optional L2.
    !> @param[out] contacts Exposure/contact classification; both slots are always defined.
    !> @param[out] valid True only for 1-2 finite nonnegative candidate thicknesses.
    PURE SUBROUTINE evaluate_layer_contacts(thickness,contacts,valid)
        REAL(wp), INTENT(IN) :: thickness(:)
        TYPE(layer_contacts_type), INTENT(OUT) :: contacts
        LOGICAL, INTENT(OUT) :: valid
        LOGICAL :: wet(2)
        INTEGER :: n_layers

        contacts=layer_contacts_type()
        valid=.FALSE.
        n_layers=SIZE(thickness)

        IF (n_layers/=1 .AND. n_layers/=2) RETURN
        IF (.NOT.ALL(ieee_is_finite(thickness))) RETURN
        IF (ANY(thickness<0.0_wp)) RETURN

        valid=.TRUE.
        wet=.FALSE.
        wet(1:n_layers)=thickness>dry_thickness_tolerance

        IF (wet(1)) THEN
            contacts%terrain_exposed(1)=.TRUE.
            IF (n_layers==1) contacts%atmosphere_exposed(1)=.TRUE.
        END IF

        IF (n_layers==2) THEN
            IF (wet(2)) THEN
                contacts%atmosphere_exposed(2)=.TRUE.
                IF (.NOT.wet(1)) contacts%terrain_exposed(2)=.TRUE.
            ELSE
                IF (wet(1)) contacts%atmosphere_exposed(1)=.TRUE.
            END IF

            contacts%interlayer_contact=wet(1) .AND. wet(2)
        END IF
    END SUBROUTINE evaluate_layer_contacts

END MODULE layer_contacts_2d
