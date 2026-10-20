! Standalone Fortran mirror of kt_v2_original_1d.py for its Ritter test,
! using reconstructed u and h_face*u_face for the transported momentum.
program kt_ritter_1d
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none
  integer, parameter :: wp=real64, nx=120, nface=nx+1
  real(wp), parameter :: dx=0.25_wp, x0=0.0_wp, dam_x=10.0_wp
  real(wp), parameter :: h0=1.0_wp, grav=9.81_wp, theta=1.3_wp
  real(wp), parameter :: cfl=0.14_wp, dt0=1.0e-4_wp
  real(wp), parameter :: max_dt=2.0e-2_wp, dt_output=1.0e-1_wp
  real(wp), parameter :: dry_tol=1.0e-12_wp, t_end=1.5_wp
  real(wp) :: q(nx,2), q1(nx,2), q2(nx,2), qnext(nx,2)
  real(wp) :: l0(nx,2), l1(nx,2), l2(nx,2)
  real(wp) :: xedge(0:nx), h_exact(nx), h_num(nx)
  real(wp) :: t, dt, dt_old, dt_old_old, speed, min_h, mass, mass0
  real(wp) :: l1_error, linf_error, exact_mass, output_time
  integer :: i, nstep

  do i=0,nx
     xedge(i)=x0+real(i,wp)*dx
  end do

  ! Cell averages of the left wet reservoir/right dry bed initial state.
  do i=1,nx
     q(i,1)=h0*max(0.0_wp,min(1.0_wp,(dam_x-xedge(i-1))/dx))
     q(i,2)=0.0_wp
  end do
  mass0=sum(q(:,1))*dx
  t=0.0_wp
  dt_old=dt0
  dt_old_old=dt0
  nstep=0
  min_h=minval(q(:,1))

  do while (t < t_end-1.0e-14_wp)
     call spatial_rhs(q,l0,speed)
     dt=min(max_dt,cfl*dx/max(speed,1.0e-12_wp),t_end-t)
     output_time=(floor(t/dt_output+1.0e-12_wp)+1)*dt_output
     dt=min(dt,output_time-t)
     dt=min(dt,1.1_wp*0.5_wp*(dt_old+dt_old_old))

     q1=q+dt*l0
     call spatial_rhs(q1,l1,speed)
     q2=0.75_wp*q+0.25_wp*(q1+dt*l1)
     call spatial_rhs(q2,l2,speed)
     qnext=(1.0_wp/3.0_wp)*q+(2.0_wp/3.0_wp)*(q2+dt*l2)

     min_h=min(min_h,minval(qnext(:,1)))
     if (minval(qnext(:,1)) < -dry_tol) then
        i=minloc(qnext(:,1),dim=1)
        write(*,'(A,ES14.6,A,I0,A,ES14.6)') &
             'FAIL: negative h=',qnext(i,1),' in cell ',i, &
             ' at t=',t+dt
        error stop 1
     end if
     q=qnext
     t=t+dt
     nstep=nstep+1
     dt_old_old=dt_old
     dt_old=dt
  end do

  h_num=q(:,1)
  call ritter_cell_averages(xedge,t,h_exact)
  exact_mass=sum(h_exact)*dx
  mass=sum(h_num)*dx
  l1_error=sum(abs(h_num-h_exact))*dx
  linf_error=maxval(abs(h_num-h_exact))
  write(*,'(A)') 'Kurganov-Petrova (2007) 1-D Fortran diagnostic'
  write(*,'(A,ES14.6,A,I0,A,ES14.6)') &
       'case=ritter, momentum_reconstruction=u, t=',t, &
       ' s, steps=',nstep,', min(raw h)=',min_h
  write(*,'(A,ES14.6,A,ES14.6,A,ES14.6,A)') &
       'mass initial/final=',mass0,' / ',mass,' m^2'
  write(*,'(A,ES14.6,A,ES14.6,A,ES14.6,A)') &
       'mass error vs exact=',mass-exact_mass,' m^2, L1(cell-average h)=', &
       l1_error,' m^2, Linf=',linf_error,' m'

contains

  pure function mm3(a,b,c) result(v)
    real(wp), intent(in) :: a,b,c
    real(wp) :: v
    if ((a > 0.0_wp).and.(b > 0.0_wp).and.(c > 0.0_wp)) then
       v=min(a,b,c)
    else if ((a < 0.0_wp).and.(b < 0.0_wp).and.(c < 0.0_wp)) then
       v=max(a,b,c)
    else
       v=0.0_wp
    end if
  end function mm3

  subroutine limited_slope(v,slope)
    real(wp), intent(in) :: v(nx)
    real(wp), intent(out) :: slope(nx)
    real(wp) :: backward, forward, centered
    integer :: j
    slope=0.0_wp
    slope(1)=(v(2)-v(1))/dx
    slope(nx)=(v(nx)-v(nx-1))/dx
    do j=2,nx-1
       backward=(v(j)-v(j-1))/dx
       forward=(v(j+1)-v(j))/dx
       centered=(v(j+1)-v(j-1))/(2.0_wp*dx)
       slope(j)=mm3(theta*backward,centered,theta*forward)
    end do
  end subroutine limited_slope

  pure function desingularized_u(h,hu) result(u)
    real(wp), intent(in) :: h,hu
    real(wp) :: u, eps, denominator
    eps=dx**4
    denominator=sqrt(h**4+max(h**4,eps))
    u=0.0_wp
    if (denominator > 0.0_wp) u=sqrt(2.0_wp)*h*hu/denominator
  end function desingularized_u

  subroutine spatial_rhs(state,rhs,max_speed)
    real(wp), intent(in) :: state(nx,2)
    real(wp), intent(out) :: rhs(nx,2),max_speed
    real(wp) :: hslope(nx),uslope(nx),ucenter(nx)
    real(wp) :: hleft(nface),hright(nface),uleft(nface),uright(nface)
    real(wp) :: mleft(nface),mright(nface)
    real(wp) :: am(nface),ap(nface),denom
    real(wp) :: flux(nface,2),fl(2),fr(2),ul(2),ur(2)
    real(wp) :: hwest,heast
    integer :: j

    call limited_slope(state(:,1),hslope)
    ucenter=0.0_wp
    do j=1,nx
       if (state(j,1) > dry_tol) &
            ucenter(j)=desingularized_u(state(j,1),state(j,2))
    end do
    call limited_slope(ucenter,uslope)

    ! Face traces follow the Python K&P reconstruction and dry correction.
    hleft(1)=state(1,1)-0.5_wp*dx*hslope(1)
    hright(nface)=state(nx,1)+0.5_wp*dx*hslope(nx)
    do j=1,nx
       hwest=state(j,1)-0.5_wp*dx*hslope(j)
       heast=state(j,1)+0.5_wp*dx*hslope(j)
       if (heast < 0.0_wp) then
          heast=0.0_wp
          hwest=2.0_wp*state(j,1)-heast
       else if (hwest < 0.0_wp) then
          hwest=0.0_wp
          heast=2.0_wp*state(j,1)-hwest
       end if
       hright(j)=max(hwest,0.0_wp)
       hleft(j+1)=max(heast,0.0_wp)
    end do

    uleft(1)=ucenter(1)-0.5_wp*dx*uslope(1)
    uright(nface)=ucenter(nx)+0.5_wp*dx*uslope(nx)
    do j=1,nx
       uright(j)=ucenter(j)-0.5_wp*dx*uslope(j)
       uleft(j+1)=ucenter(j)+0.5_wp*dx*uslope(j)
    end do
    do j=1,nface
       if (hleft(j) <= dry_tol) uleft(j)=0.0_wp
       if (hright(j) <= dry_tol) uright(j)=0.0_wp
       mleft(j)=hleft(j)*uleft(j)
       mright(j)=hright(j)*uright(j)
       am(j)=min(uleft(j)-sqrt(grav*hleft(j)), &
            uright(j)-sqrt(grav*hright(j)),0.0_wp)
       ap(j)=max(uleft(j)+sqrt(grav*hleft(j)), &
            uright(j)+sqrt(grav*hright(j)),0.0_wp)
       denom=ap(j)-am(j)
       flux(j,:)=0.0_wp
       if (denom > 0.0_wp) then
          fl=[mleft(j),mleft(j)*uleft(j)+0.5_wp*grav*hleft(j)**2]
          fr=[mright(j),mright(j)*uright(j)+0.5_wp*grav*hright(j)**2]
          ul=[hleft(j),mleft(j)]
          ur=[hright(j),mright(j)]
          flux(j,:)=(ap(j)*fl-am(j)*fr+ap(j)*am(j)*(ur-ul))/denom
       end if
    end do
    do j=1,nx
       rhs(j,:)=-(flux(j+1,:)-flux(j,:))/dx
    end do
    max_speed=max(maxval(abs(am)),maxval(abs(ap)))
  end subroutine spatial_rhs

  subroutine ritter_cell_averages(edges,time,hbar)
    real(wp), intent(in) :: edges(0:nx),time
    real(wp), intent(out) :: hbar(nx)
    real(wp), parameter :: xg(12)=[ &
         -0.9815606342467192_wp,-0.9041172563704749_wp, &
         -0.7699026741943047_wp,-0.5873179542866175_wp, &
         -0.3678314989981802_wp,-0.1252334085114689_wp, &
          0.1252334085114689_wp, 0.3678314989981802_wp, &
          0.5873179542866175_wp, 0.7699026741943047_wp, &
          0.9041172563704749_wp, 0.9815606342467192_wp]
    real(wp), parameter :: wg(12)=[ &
         0.0471753363865118_wp,0.1069393259953184_wp, &
         0.1600783285433462_wp,0.2031674267230659_wp, &
         0.2334925365383548_wp,0.2491470458134028_wp, &
         0.2491470458134028_wp,0.2334925365383548_wp, &
         0.2031674267230659_wp,0.1600783285433462_wp, &
         0.1069393259953184_wp,0.0471753363865118_wp]
    real(wp) :: c0,x,xi,hpoint,total
    integer :: j,k
    c0=sqrt(grav*h0)
    do j=1,nx
       total=0.0_wp
       do k=1,12
          x=0.5_wp*(edges(j)-edges(j-1))*xg(k) &
               +0.5_wp*(edges(j)+edges(j-1))
          if (time <= 0.0_wp) then
             hpoint=h0*max(0.0_wp,min(1.0_wp,(dam_x-edges(j-1))/dx))
          else if (x <= dam_x-c0*time) then
             hpoint=h0
          else if (x < dam_x+2.0_wp*c0*time) then
             xi=(x-dam_x)/time
             hpoint=((2.0_wp*c0-xi)/3.0_wp)**2/grav
          else
             hpoint=0.0_wp
          end if
          total=total+wg(k)*hpoint
       end do
       hbar(j)=0.5_wp*total
    end do
  end subroutine ritter_cell_averages

end program kt_ritter_1d
