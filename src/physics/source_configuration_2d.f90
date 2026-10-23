!********************************************************************************
!> \brief Translate configured source activation into validated D17 metadata.
!> \details Runs once after legacy equation partition initialization. Owns no
!>          physical state and performs no rate evaluation or temporal update.
!>          Per-cell wetness, packing and window gates remain in existing kernels.
!********************************************************************************
MODULE source_configuration_2d
  USE source_catalog_2d
  USE source_registry_2d, ONLY : source_registry_type
  USE equation_metadata_2d, ONLY : equation_partition_type, resolve_imex_group_method
  USE parameters_2d, ONLY : wp, rheology_flag, rheology_model_id, rheology_id_none, &
       pore_pressure_flag, curvature_term_flag, bottom_radial_source_flag, &
       bottom_fissural_source_flag, radial_source_flag, lateral_source_flag, &
       liquid_flag, gas_flag, gas_loss_flag, liquid_vaporization_flag, &
       imex_velocity_group, imex_energy_group, imex_pore_pressure_group, imex_group_method, idx_pore
  USE constitutive_parameters_2d, ONLY : erosion_coeff, settling_flag, &
       coeff_porosity, entrainment_flag
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_configured_sources
CONTAINS
  !> \brief Build simulation-owned source metadata from current validated input.
  !> \param[in,out] registry Destination, cleared if translation fails.
  !> \param[in] partition Existing compact equation partition, never modified.
  !> \param[out] error Empty on success, otherwise a pre-simulation diagnostic.
  !> \note This phase records possible activation; zero per-cell rates stay valid.
  SUBROUTINE initialize_configured_sources(registry,partition,error)
    TYPE(source_registry_type), INTENT(INOUT) :: registry
    TYPE(equation_partition_type), INTENT(IN) :: partition
    CHARACTER(*), INTENT(OUT) :: error
    INTEGER, ALLOCATABLE :: terms(:)
    INTEGER :: methods(3),i,pore
    CHARACTER(LEN=16) :: method
    CALL registry%finalize()
    error=''
    IF (.NOT.partition%is_initialized()) THEN
      error='D17: source configuration requires initialized equation metadata'; RETURN
    END IF
    IF (partition%n_equations<4) THEN
      error='D17: source layout requires at least four equations'; RETURN
    END IF
    terms=[INTEGER::]
    IF (rheology_flag) THEN
      terms=rheology_source_terms(rheology_model_id)
      IF (SIZE(terms)==0 .AND. rheology_model_id/=rheology_id_none) THEN
        error='D17: unsupported configured rheology'; RETURN
      END IF
    END IF
    IF (partition%implicit(4)) terms=[terms,SOURCE_TERM_THERMAL_RADIATIVE_COOLING, &
         SOURCE_TERM_THERMAL_CONVECTIVE_COOLING]
    pore=0
    IF (pore_pressure_flag) THEN
      pore=idx_pore
      terms=[terms,SOURCE_TERM_PORE_PRESSURE_DISSIPATION]
    END IF
    IF (curvature_term_flag) terms=[terms,SOURCE_TERM_GEOMETRY_CURVATURE]
    IF (bottom_radial_source_flag) terms=[terms,SOURCE_TERM_BOTTOM_CIRCULAR_INJECTION]
    IF (bottom_fissural_source_flag) terms=[terms,SOURCE_TERM_BOTTOM_FISSURAL_INJECTION]
    IF (radial_source_flag) terms=[terms,SOURCE_TERM_RADIAL_INJECTION]
    IF (lateral_source_flag) terms=[terms,SOURCE_TERM_BOUNDARY_INJECTION]
    IF (erosion_coeff>0.0_wp) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_EROSION]
    IF (settling_flag) THEN
      terms=[terms,SOURCE_TERM_MASS_EXCHANGE_DEPOSITION]
      IF (ABS(coeff_porosity)>0.0_wp) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_CARRIER_DEPOSITION_LOSS]
    END IF
    IF (liquid_flag) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_CARRIER_LIQUID_LOSS]
    IF (pore_pressure_flag .AND. gas_loss_flag) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_CARRIER_DARCY_LOSS]
    IF (entrainment_flag) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_ENTRAINMENT]
    IF (gas_flag .AND. liquid_vaporization_flag) terms=[terms,SOURCE_TERM_MASS_EXCHANGE_VAPORIZATION]
    DO i=1,3
      method=imex_group_method(i)
      CALL resolve_imex_group_method(method,methods(i))
    END DO
    CALL registry%initialize_legacy(partition, &
         [imex_velocity_group,imex_energy_group,imex_pore_pressure_group], &
         methods,pore,terms,error)
  END SUBROUTINE initialize_configured_sources
END MODULE source_configuration_2d
