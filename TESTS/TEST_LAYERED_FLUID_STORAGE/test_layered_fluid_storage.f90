PROGRAM test_layered_fluid_storage
    USE, INTRINSIC :: iso_fortran_env, ONLY : real64
    USE, INTRINSIC :: iso_c_binding, ONLY : c_loc, c_associated
    USE model_layout_2d, ONLY : model_layout_type
    USE layered_fluid_storage_2d, ONLY : layered_fluid_storage_type
    IMPLICIT NONE
    TYPE(model_layout_type) :: layout
    TYPE(layered_fluid_storage_type) :: storage
    INTEGER :: i,j,k
    REAL(real64) :: cell(6)
    CALL layout%initialize(1,6,6,1,1,0,0)
    CALL storage%initialize(layout,4,3)
    IF (ANY(SHAPE(storage%q)/=[6,1,4,3])) ERROR STOP 'wrong q shape'
    IF (ANY(SHAPE(storage%qp)/=[6,1,4,3])) ERROR STOP 'wrong qp shape'
    DO k=1,3
        DO j=1,4
            storage%q(:,1,j,k)=[(REAL(i+10*j+100*k,real64),i=1,6)]
            storage%qp(:,1,j,k)=storage%q(:,1,j,k)*2.0_real64
            IF (.NOT.IS_CONTIGUOUS(storage%q(:,1,j,k))) &
                ERROR STOP 'non-contiguous q cell'
            IF (.NOT.IS_CONTIGUOUS(storage%qp(:,1,j,k))) &
                ERROR STOP 'non-contiguous qp cell'
            cell=storage%q(:,1,j,k)
            IF (ANY(cell/=[(REAL(i+10*j+100*k,real64),i=1,6)])) &
                ERROR STOP 'cell slice changed values'
        END DO
    END DO
    ! A layer-local kernel must be able to update a cell through
    ! the section without touching adjacent cells or layers.
    CALL update_cell(storage%q(:,1,2,2))
    IF (ANY(storage%q(:,1,2,2)/=-1.0_real64)) ERROR STOP 'cell update failed'
    IF (ANY(storage%q(:,1,3,2)/=[(REAL(i+30+200,real64),i=1,6)])) &
        ERROR STOP 'neighbor was modified'
    CALL storage%finalize()
    IF (ALLOCATED(storage%q) .OR. ALLOCATED(storage%qp)) &
        ERROR STOP 'storage finalize failed'
    CALL storage%initialize(layout,2,2)
    IF (ANY(SHAPE(storage%q)/=[6,1,2,2])) ERROR STOP 'storage reuse failed'
    CALL storage%finalize()
    PRINT *, 'PASS: rank-four single-layer fluid storage and cell slices'
CONTAINS
    SUBROUTINE update_cell(q_cell)
        REAL(real64), INTENT(INOUT) :: q_cell(:)
        q_cell=-1.0_real64
    END SUBROUTINE update_cell
END PROGRAM test_layered_fluid_storage
