!********************************************************************************
!> \brief Initial solution
!
!> This module contains the variables and the subroutine for the
!> initialization of the solution for a Riemann problem.
!********************************************************************************

MODULE init_2d

   USE parameters_2d, ONLY : wp
   USE parameters_2d, ONLY : verbose_level
   USE parameters_2d, ONLY : n_solid , n_add_gas
   USE parameters_2d, ONLY : n_stoch_vars , n_pore_vars
   USE state_2d, ONLY : state_type

   IMPLICIT none

   REAL(wp), ALLOCATABLE :: q_init(:,:,:)

   REAL(wp), ALLOCATABLE :: thickness_init(:,:)

   !> Initial thickness of erodible layer (solid+voids)
   REAL(wp), ALLOCATABLE :: erodible_init(:,:)

CONTAINS

   SUBROUTINE init_empty(state)

      USE constitutive_parameters_2d, ONLY : T_ambient

      USE state_conversion_2d, ONLY : qp_to_qc

      USE geometry_2d, ONLY : comp_cells_x , comp_cells_y

      USE parameters_2d, ONLY : n_vars

      IMPLICIT NONE

      CLASS(state_type), INTENT(INOUT) :: state

      INTEGER :: j,k

      REAL(wp) :: qp_init(n_vars+2)

      WRITE(*,*) 'Initialization with zero thickness flow'

      qp_init(1:n_vars+2) = 0.0_wp
      qp_init(4) = T_ambient

      DO j = 1,comp_cells_x

         DO k = 1,comp_cells_y

            CALL qp_to_qc( qp_init(1:n_vars+2) , state%q(1:n_vars,j,k) )


         END DO

      END DO

      RETURN

   END SUBROUTINE init_empty

   !******************************************************************************
   !> \brief Collapsing volume initialization
   !
   !> This subroutine initialize the solution for a collpasing volume. Values for
   !>  the initial state (x, y, r, T, h, alphas) are read from the input file.
   !> \date 2019_12_11
   !
   !> @author
   !> Mattia de' Michieli Vitturi
   !
   !******************************************************************************

   SUBROUTINE collapsing_volume(state)

      USE state_conversion_2d, ONLY : qp_to_qc,                               &
           eval_mixture_properties_from_volume_fractions

      USE geometry_2d, ONLY : compute_cell_fract

      USE geometry_2d, ONLY : comp_cells_x , comp_cells_y

      USE parameters_2d, ONLY : n_vars

      USE parameters_2d, ONLY : x_collapse , y_collapse , r_collapse , T_collapse , &
         h_collapse , alphas_collapse , alphag_collapse

      IMPLICIT NONE

      CLASS(state_type), INTENT(INOUT) :: state

      INTEGER :: j,k

      REAL(wp) :: qp_init(n_vars+2) ,  qp0_init(n_vars+2)
      REAL(wp) :: xs_collapse(n_solid), xg_collapse(n_add_gas), xl_collapse
      REAL(wp) :: alphas_local(n_solid), alphal_local
      REAL(wp) :: rho_m, inv_rhom, rho_c, xc, sp_heat_c, sp_heat_mix

      REAL(wp), ALLOCATABLE :: cell_fract(:,:)

      ALLOCATE(cell_fract(comp_cells_x,comp_cells_y))

      CALL compute_cell_fract(x_collapse,y_collapse,r_collapse,r_collapse,0.0_wp,cell_fract)

      ! values outside the collapsing volume
      qp0_init = 0.0_wp
      qp0_init(1) = 0.0_wp                  ! h
      qp0_init(2) = 0.0_wp                  ! hu
      qp0_init(3) = 0.0_wp                  ! hv
      qp0_init(4) = T_collapse              ! T
      qp0_init(5:4+n_solid) = 0.0_wp        ! solid mass fractions
      qp0_init(4+n_solid+1:4+n_solid+n_add_gas) = 0.0_wp        ! gas mass fractions
      qp0_init(5+n_solid+n_add_gas:4+n_solid+n_add_gas+n_stoch_vars) = 0.0_wp
      qp0_init(5+n_solid+n_add_gas+n_stoch_vars:4+n_solid+n_add_gas+n_stoch_vars+ &
         n_pore_vars) =  0.0_wp
      qp0_init(n_vars+1:n_vars+2) = 0.0_wp  ! u,v

      ! values within the collapsing volume
      qp_init = 0.0_wp
      qp_init(2) = 0.0_wp
      qp_init(3) = 0.0_wp
      qp_init(4) = T_collapse

      qp_init(n_vars+1:n_vars+2) = 0.0_wp

      alphas_local = alphas_collapse(1:n_solid)
      alphal_local = 0.0_wp
      CALL eval_mixture_properties_from_volume_fractions(                     &
           T_collapse, alphag_collapse(1:n_add_gas), alphas_local,            &
           alphal_local, rho_m, inv_rhom, rho_c, xs_collapse, xg_collapse,    &
           xl_collapse, xc, sp_heat_c, sp_heat_mix)

      qp_init(5:4+n_solid) = xs_collapse
      qp_init(4+n_solid+1:4+n_solid+n_add_gas) = xg_collapse

      DO j = 1,comp_cells_x

         DO k = 1,comp_cells_y

            IF ( cell_fract(j,k) .GT. 0.0_wp ) THEN

               qp_init(1) = cell_fract(j,k) * h_collapse

               qp_init(5+n_solid+n_add_gas:4+n_solid+n_add_gas+n_stoch_vars) =    &
                  1.0_wp

               qp_init(5+n_solid+n_add_gas+n_stoch_vars:4+n_solid+n_add_gas       &
                  +n_stoch_vars+n_pore_vars) = 1.0_wp

               CALL qp_to_qc( qp_init(1:n_vars+2) , state%q(1:n_vars,j,k) )

            ELSE

               CALL qp_to_qc( qp0_init(1:n_vars+2) , state%q(1:n_vars,j,k) )

            END IF

         END DO

      END DO

      RETURN

   END SUBROUTINE collapsing_volume

END MODULE init_2d
