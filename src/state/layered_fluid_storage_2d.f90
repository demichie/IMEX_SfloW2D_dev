!> @brief Rank-four fluid storage for the staged two-layer state migration.
!> @details Layout is (component, layer, x-cell, y-cell). Cell-local kernel
!>          arguments use the contiguous section (:,layer,j,k).
!>          This type is not yet owned by state_type; integration into
!>          production state and all consumers belongs to TL01/TL02.
MODULE layered_fluid_storage_2d
    USE, INTRINSIC :: iso_fortran_env, ONLY : real64
    USE model_layout_2d, ONLY : model_layout_type
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: layered_fluid_storage_type

    TYPE :: layered_fluid_storage_type
        REAL(real64), ALLOCATABLE :: q(:,:,:,:)
        REAL(real64), ALLOCATABLE :: qp(:,:,:,:)
    CONTAINS
        PROCEDURE :: initialize => initialize_layered_fluid_storage
        PROCEDURE :: finalize => finalize_layered_fluid_storage
    END TYPE layered_fluid_storage_type

CONTAINS

    SUBROUTINE initialize_layered_fluid_storage(this,layout,nx,ny)
        CLASS(layered_fluid_storage_type), INTENT(INOUT) :: this
        TYPE(model_layout_type), INTENT(IN) :: layout
        INTEGER, INTENT(IN) :: nx,ny
        IF (.NOT.layout%is_initialized()) &
            ERROR STOP 'Layered fluid storage requires an initialized layout'
        IF (nx<1 .OR. ny<1) &
            ERROR STOP 'Layered fluid storage requires positive grid extents'
        CALL this%finalize()
        ALLOCATE(this%q(layout%n_variables,layout%n_layers,nx,ny))
        ALLOCATE(this%qp(layout%n_physical_variables,layout%n_layers,nx,ny))
        this%q=0.0_real64
        this%qp=0.0_real64
    END SUBROUTINE initialize_layered_fluid_storage

    SUBROUTINE finalize_layered_fluid_storage(this)
        CLASS(layered_fluid_storage_type), INTENT(INOUT) :: this
        IF (ALLOCATED(this%q)) DEALLOCATE(this%q)
        IF (ALLOCATED(this%qp)) DEALLOCATE(this%qp)
    END SUBROUTINE finalize_layered_fluid_storage
END MODULE layered_fluid_storage_2d
