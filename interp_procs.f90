subroutine traject(u1,u2,v1,v2,w1,w2,el1,el2)
  !==============================================================================|
  !  Approximate particle postition from (x0,y0,z0) to (xn,yn,zn) using the      |
  !  4-stage Runge-Kutta itteritave method.                                      |
  !------------------------------------------------------------------------------|
  !  u1, v1, w1: Velocity field at time T0                                       |
  !  u2, v2, w2: Velocity field at time T0 + DTI                                 |
  !  el1, el2:   Free surface elevation                                          |
  !==============================================================================|
  use mod_config
  use mod_tracking
  use mod_flow_field
  implicit none
  !------------------------------------------------------------------------------|
  real(DP), dimension(ELEMENTS,SIGLAY), intent(in) :: u1, u2
  real(DP), dimension(ELEMENTS,SIGLAY), intent(in) :: v1, v2
  real(DP), dimension(ELEMENTS,SIGLAY), intent(in) :: w1, w2
  real(DP), dimension(NODES),           intent(in) :: el1, el2
  !------------------------------------------------------------------------------|
  integer                               :: i, nactive
  integer                               :: ns
  integer,  dimension(NDRFT)            :: ai
  real(DP), dimension(ELEMENTS,SIGLAY)  :: u, v, w
  real(DP), dimension(NODES)            :: el
  real(DP), dimension(:),   allocatable :: x0, y0, z0
  real(DP), dimension(:),   allocatable :: x, y, z
  real(DP), dimension(:),   allocatable :: hp, elp
  real(DP), dimension(:,:), allocatable :: vx, vy, vz
  real,     dimension(:),   allocatable :: depth
  integer,  dimension(:),   allocatable :: host
  logical,  dimension(:),   allocatable :: indomain
  !------------------------------------------------------------------------------|
  integer, parameter                    :: mstage = 4
  real,    parameter, dimension(mstage) :: a_rk = [0.0, 0.5, 0.5, 1.0]
  real,    parameter, dimension(mstage) :: b_rk = [1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0]
  real,    parameter, dimension(mstage) :: c_rk = [0.0, 0.5, 0.5, 1.0]
  !==============================================================================|

  !------------------------------------------------------------------------------|
  !  Select active particles                                                     |
  !------------------------------------------------------------------------------|
  nactive = 0
  do i = 1,NDRFT
    if (LAG_INDOMAIN(i) .and. LAG_START(i) <= SIM_TIME .and. LAG_STOP(i) > SIM_TIME) then
      nactive = nactive + 1
      ai(nactive) = i
    end if
  end do

  !------------------------------------------------------------------------------|
  !  Allocate working memory                                                     |
  !------------------------------------------------------------------------------|
  allocate(x0(nactive), y0(nactive), z0(nactive))
  allocate(x(nactive), y(nactive), z(nactive))
  allocate(hp(nactive), elp(nactive))
  allocate(host(nactive), indomain(nactive))
  allocate(vx(0:mstage,nactive))
  allocate(vy(0:mstage,nactive))
  allocate(vz(0:mstage,nactive))
  if (F_DEPTH) allocate(depth(nactive))

  !------------------------------------------------------------------------------|
  !  Initialize stage variables                                                  |
  !------------------------------------------------------------------------------|
  vx(0,:) = 0.0_dp
  vy(0,:) = 0.0_dp
  vz(0,:) = 0.0_dp

  do i = 1,nactive
    x0(i)   = LAG_XP(ai(i))
    y0(i)   = LAG_YP(ai(i))
    z0(i)   = LAG_ZP(ai(i))
    host(i) = LAG_HOST(ai(i))
    if (F_DEPTH) depth(i) = LAG_D(ai(i))
  end do
  indomain = .true.

  if (F_DEPTH) then
    call interp_elh(nactive,host,x0,y0,el1,hp,elp)
    z0 = 0.0_dp - depth/(hp + elp)
  end if

  !------------------------------------------------------------------------------|
  !  Loop over RK Stages                                                         |
  !------------------------------------------------------------------------------|
  do ns = 1,mstage

    !------------------------------------------------------------------------------|
    !  Particle position at stage N                                                |
    !------------------------------------------------------------------------------|
    x = x0 + a_rk(ns)*float(DTI)*vx(ns - 1,:)
    y = y0 + a_rk(ns)*float(DTI)*vy(ns - 1,:)
    z = a_rk(ns)*float(DTI)*vz(ns - 1,:) ! Displacement in meters
    do i = 1,nactive
      call fhe(x(i),y(i),host(i),indomain(i))
    end do

    !------------------------------------------------------------------------------|
    !  Calculate velocity field for stage N using c_rk coefficients                |
    !------------------------------------------------------------------------------|
    u  = (1.0 - c_rk(ns))*u1  + c_rk(ns)*u2
    v  = (1.0 - c_rk(ns))*v1  + c_rk(ns)*v2
    w  = (1.0 - c_rk(ns))*w1  + c_rk(ns)*w2
    el = (1.0 - c_rk(ns))*el1 + c_rk(ns)*el2

    !------------------------------------------------------------------------------|
    !  We need the particle's sigma position to determine the stage N velocity     |
    !------------------------------------------------------------------------------|
    call interp_elh(nactive,host,x,y,el,hp,elp)
    if (F_DEPTH) then
      z = 0.0_dp - depth/(hp + elp)
    else
      z = z0 + z/(hp + elp)
    end if

    if (F_DEPTH) where (z < -1.0) indomain = .false.
    where (z < -1.0) z = -2.0_dp - z ! Reflect off bottom
    where (z > 0.0) z = 0.0_dp ! Can't move above water surface

    !------------------------------------------------------------------------------|
    !  Evaluate velocity (u,v,w) at stage ns particle position                     |
    !------------------------------------------------------------------------------|
    call interp_v(nactive,host,x,y,z,u,v,w,vx(ns,:),vy(ns,:),vz(ns,:))
  end do

  !------------------------------------------------------------------------------|
  !  Sum stage contributions to get updated particle positions                   |
  !------------------------------------------------------------------------------|
  x = x0
  y = y0
  z = 0.0_dp
  do ns = 1,mstage
    x = x + float(DTI)*vx(ns,:)*b_rk(ns)
    y = y + float(DTI)*vy(ns,:)*b_rk(ns)
    z = z + float(DTI)*vz(ns,:)*b_rk(ns)
  end do

  do i = 1,nactive
    call fhe(x(i),y(i),host(i),indomain(i))
  end do

  !------------------------------------------------------------------------------|
  !  Convert the particle's new z position to sigma depth                        |
  !------------------------------------------------------------------------------|
  call interp_elh(nactive,host,x,y,el2,hp,elp)
  if (F_DEPTH) then
    z = 0.0_dp - depth/(hp + elp)
  else
    z = z0 + z/(hp + elp)
  end if

  if (F_DEPTH) where (z < -1.0) indomain = .false.
  where (z < -1.0) z = -2.0_dp - z ! Reflect off bottom
  where (z > 0.0) z = 0.0_dp ! Can't move above water surface

  !------------------------------------------------------------------------------|
  !  Update particle tracking data arrays                                        |
  !------------------------------------------------------------------------------|
  do i = 1,nactive
    LAG_XP(ai(i))       = x(i)
    LAG_YP(ai(i))       = y(i)
    LAG_ZP(ai(i))       = z(i)
    LAG_HOST(ai(i))     = host(i)
    LAG_INDOMAIN(ai(i)) = indomain(i)
  end do
end subroutine traject

!==============================================================================|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!==============================================================================|

subroutine random_walk(ell)
  !==============================================================================|
  !  Apply the random walk algoritham to the currently active particles          |
  !==============================================================================|
  use mod_tracking
  implicit none
  !------------------------------------------------------------------------------|
  real(DP), dimension(NODES), intent(in) :: ell
  !------------------------------------------------------------------------------|
  integer                               :: i, nactive
  integer,  dimension(NDRFT)            :: ai
  real(DP), dimension(:),   allocatable :: x, y, z
  real(DP), dimension(:),   allocatable :: dzm
  real(DP), dimension(:),   allocatable :: depth
  integer,  dimension(:),   allocatable :: host
  logical,  dimension(:),   allocatable :: indomain
  real(DP), dimension(:),   allocatable :: hp, elp
  real,     dimension(:,:), allocatable :: rnd
  real,     dimension(:),   allocatable :: rnd_vect
  !==============================================================================|

  !------------------------------------------------------------------------------|
  !  Select active particles                                                     |
  !------------------------------------------------------------------------------|
  nactive = 0
  do i = 1,NDRFT
    if (LAG_INDOMAIN(i) .and. LAG_START(i) <= SIM_TIME .and. LAG_STOP(i) > SIM_TIME) then
      nactive = nactive + 1
      ai(nactive) = i
    end if
  end do

  !------------------------------------------------------------------------------|
  !  Allocate working memory                                                     |
  !------------------------------------------------------------------------------|
  allocate(x(nactive), y(nactive), z(nactive))
  allocate(dzm(nactive))
  allocate(host(nactive), indomain(nactive))
  allocate(hp(nactive), elp(nactive))
  allocate(rnd(nactive,3))
  allocate(rnd_vect(nactive*3))
  if (F_DEPTH) allocate(depth(nactive))

  !------------------------------------------------------------------------------|
  !  Initalize working memory                                                    |
  !------------------------------------------------------------------------------|
  do i = 1,nactive
    x(i)   = LAG_XP(ai(i))
    y(i)   = LAG_YP(ai(i))
    z(i)   = LAG_ZP(ai(i))
    host(i) = LAG_HOST(ai(i))
    if (F_DEPTH) depth(i) = LAG_D(ai(i))
  end do
  indomain = .true.

  if (F_DEPTH) then
    call interp_elh(nactive,host,x,y,ell,hp,elp)
    z = 0.0_dp - depth/(hp + elp)
  end if

  !------------------------------------------------------------------------------|
  !  Add a random varriance to the particle positions                            |
  !------------------------------------------------------------------------------|
  rnd_vect = reshape(rnd, [size(rnd)])
  call gauss_random(size(rnd_vect), rnd_vect)
  rnd = reshape(rnd_vect, [size(rnd, 1), size(rnd, 2)])
  x = x + rnd(:,1)*((float(DTI)*2*K_XY)**0.5)
  y = y + rnd(:,2)*((float(DTI)*2*K_XY)**0.5)
  dzm =   rnd(:,3)*((float(DTI)*2*K_Z)**0.5)

  do i = 1,nactive
    call fhe(x(i),y(i),host(i),indomain(i))
  end do

  !------------------------------------------------------------------------------|
  !  Convert the particle's new z position to sigma depth                        |
  !------------------------------------------------------------------------------|
  call interp_elh(nactive,host,x,y,ell,hp,elp)
  if (F_DEPTH) then
    z = 0.0_dp - depth/(hp + elp)
  else
    z = z + dzm/(hp + elp)
  end if

  if (F_DEPTH) where (z < -1.0) indomain = .false.
  where (z < -1.0) z = -2.0_dp - z ! Reflect off bottom
  where (z > 0.0) z = 0.0_dp ! Can't move above water surface

  !------------------------------------------------------------------------------|
  !  Update particle tracking data arrays                                        |
  !------------------------------------------------------------------------------|
  do i = 1,nactive
    LAG_XP(ai(i))       = x(i)
    LAG_YP(ai(i))       = y(i)
    LAG_ZP(ai(i))       = z(i)
    LAG_HOST(ai(i))     = host(i)
    LAG_INDOMAIN(ai(i)) = indomain(i)
  end do
end subroutine random_walk

!==============================================================================|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!==============================================================================|

subroutine interp_v(np,host,x,y,z,uin,vin,win,u,v,w) 
  !==============================================================================|
  !  Give a point (x,y,z), obtain a linear interpolation                         |
  !  of the provided velocity field at the point                                 |
  !                                                                              |
  !  RETURN:                                                                     |
  !     u,v,w (velocity field at x,y,z)                                          |
  !==============================================================================|
  use mod_flow_field
  implicit none
  !------------------------------------------------------------------------------|
  integer,                              intent(in)  :: np
  integer,  dimension(np),              intent(in)  :: host
  real(DP), dimension(ELEMENTS,SIGLAY), intent(in)  :: uin, vin, win
  real(DP), dimension(np),              intent(in)  :: x, y, z
  real(DP), dimension(np),              intent(out) :: u, v, w
  !------------------------------------------------------------------------------|
  integer                     :: e1, e2, e3, k1, k2
  real(DP)                    :: x0c, y0c
  real(DP)                    :: dudx, dudy, dvdx, dvdy, dwdx, dwdy
  real(DP)                    :: ue01, ue02, ve01, ve02, we01, we02
  real(DP)                    :: zf1, zf2
  real(DP), dimension(SIGLAY) :: dzp
  integer                     :: npp, i, j, k
  !==============================================================================|

  u = 0.0_dp
  v = 0.0_dp
  w = 0.0_dp

  do j = 1,np
    i   = host(j)
    e1  = NBE(i,1)
    e2  = NBE(i,2)
    e3  = NBE(i,3)
    x0c = x(j) - XC(i)
    y0c = y(j) - YC(i)

    !------------------------------------------------------------------------------|
    !  Find sigma layer above and below particle                                   |
    !------------------------------------------------------------------------------|
    if(z(j) > Z_LEV(2)) then ! Particle near surface
      k1  = 1
      k2  = 1
      zf1 = 1.0_dp 
      zf2 = 0.0_dp
    else if(z(j) < Z_LEV(SIGLEV-1)) then ! Particle near bottom
      k1 = SIGLAY
      k2 = SIGLAY
      zf1 = 1.0_dp
      zf2 = 0.0_dp
    else
      dzp = Z_LAY - z(j)
      do npp = 1,SIGLAY-1
        if (dzp(npp) > 0 .and. dzp(npp + 1) < 0) then
          k1 = npp
        end if
      end do
      k2 = k1 + 1
      if (Z_LAY(k1) == Z_LAY(k2)) then
        zf1 = 0
        zf2 = 0
      else
        zf1 = (z(j) - Z_LAY(k2))/(Z_LAY(k1) - Z_LAY(k2))
        zf2 = (Z_LAY(k1) - z(j))/(Z_LAY(k1) - Z_LAY(k2))
      end if
    end if

    !------------------------------------------------------------------------------|
    !  Linear interpolation of particle velocity in sigma level above particle     |
    !------------------------------------------------------------------------------|
    k = k1
    dudx = A1U(i,1)*uin(i,k) + A1U(i,2)*uin(e1,k) + A1U(i,3)*uin(e2,k) + A1U(i,4)*uin(e3,k)
    dudy = A2U(i,1)*uin(i,k) + A2U(i,2)*uin(e1,k) + A2U(i,3)*uin(e2,k) + A2U(i,4)*uin(e3,k)
    dvdx = A1U(i,1)*vin(i,k) + A1U(i,2)*vin(e1,k) + A1U(i,3)*vin(e2,k) + A1U(i,4)*vin(e3,k)
    dvdy = A2U(i,1)*vin(i,k) + A2U(i,2)*vin(e1,k) + A2U(i,3)*vin(e2,k) + A2U(i,4)*vin(e3,k)
    dwdx = A1U(i,1)*win(i,k) + A1U(i,2)*win(e1,k) + A1U(i,3)*win(e2,k) + A1U(i,4)*win(e3,k)
    dwdy = A2U(i,1)*win(i,k) + A2U(i,2)*win(e1,k) + A2U(i,3)*win(e2,k) + A2U(i,4)*win(e3,k)
    ue01 = uin(i,k) + dudx*x0c + dudy*y0c
    ve01 = vin(i,k) + dvdx*x0c + dvdy*y0c
    we01 = win(i,k) + dwdx*x0c + dwdy*y0c

    !------------------------------------------------------------------------------|
    !  Linear interpolation of particle velocity in sigma level below particle     |
    !------------------------------------------------------------------------------|
    k = k2
    dudx = A1U(i,1)*uin(i,k) + A1U(i,2)*uin(e1,k) + A1U(i,3)*uin(e2,k) + A1U(i,4)*uin(e3,k)
    dudy = A2U(i,1)*uin(i,k) + A2U(i,2)*uin(e1,k) + A2U(i,3)*uin(e2,k) + A2U(i,4)*uin(e3,k)
    dvdx = A1U(i,1)*vin(i,k) + A1U(i,2)*vin(e1,k) + A1U(i,3)*vin(e2,k) + A1U(i,4)*vin(e3,k)
    dvdy = A2U(i,1)*vin(i,k) + A2U(i,2)*vin(e1,k) + A2U(i,3)*vin(e2,k) + A2U(i,4)*vin(e3,k)
    dwdx = A1U(i,1)*win(i,k) + A1U(i,2)*win(e1,k) + A1U(i,3)*win(e2,k) + A1U(i,4)*win(e3,k)
    dwdy = A2U(i,1)*win(i,k) + A2U(i,2)*win(e1,k) + A2U(i,3)*win(e2,k) + A2U(i,4)*win(e3,k)
    ue02 = uin(i,k) + dudx*x0c + dudy*y0c
    ve02 = vin(i,k) + dvdx*x0c + dvdy*y0c
    we02 = win(i,k) + dwdx*x0c + dwdy*y0c

    !------------------------------------------------------------------------------|
    !  Interpolate particle velocity between two sigma layers                      |
    !------------------------------------------------------------------------------|
    u(j) = ue01*zf1 + ue02*zf2
    v(j) = ve01*zf1 + ve02*zf2
    w(j) = we01*zf1 + we02*zf2
  end do
end subroutine interp_v

!==============================================================================|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!==============================================================================|

subroutine interp_elh(np,host,x,y,el,hp,elp)
  !==============================================================================|
  !  Given a point (x,y), obtain a linear interpolation                          |
  !  of the provided free surface/bathymetry (hin,ein) at the point              |
  !                                                                              |
  !  RETURN:                                                                     |
  !     hp(bathymetry at x,y) and elp (free surface elevation at x,y)            |
  !==============================================================================|
  use mod_flow_field
  implicit none
  !------------------------------------------------------------------------------|
  integer,                    intent(in)  :: np
  integer,  dimension(np),    intent(in)  :: host
  real(DP), dimension(np),    intent(in)  :: x, y
  real(DP), dimension(NODES), intent(in)  :: el
  real(DP), dimension(np),    intent(out) :: hp, elp
  !------------------------------------------------------------------------------|
  integer,  dimension(np) :: n1, n2, n3
  real(DP), dimension(np) :: h0, hx, hy
  real(DP), dimension(np) :: e0, ex, ey
  real(DP), dimension(np) :: x0c, y0c
  !==============================================================================|

  !------------------------------------------------------------------------------|
  !  Linearly interpolate free surface height and bathymetry                     |
  !------------------------------------------------------------------------------|
  n1  = NV(host,1)
  n2  = NV(host,2)
  n3  = NV(host,3)
  x0c = x - XC(host)
  y0c = y - YC(host)

  !------------------------------------------------------------------------------|
  !  Linear interpolation of bathymetry                                          |
  !------------------------------------------------------------------------------|
  h0 = AW0(host,1)*H(n1) + AW0(host,2)*H(n2) + AW0(host,3)*H(n3)
  hx = AWX(host,1)*H(n1) + AWX(host,2)*H(n2) + AWX(host,3)*H(n3)
  hy = AWY(host,1)*H(n1) + AWY(host,2)*H(n2) + AWY(host,3)*H(n3)
  hp = h0 + hx*x0c + hy*y0c

  !------------------------------------------------------------------------------|
  !  Linear interpolation of free surface height                                 |
  !------------------------------------------------------------------------------|
  e0 = AW0(host,1)*el(n1) + AW0(host,2)*el(n2) + AW0(host,3)*el(n3)
  ex = AWX(host,1)*el(n1) + AWX(host,2)*el(n2) + AWX(host,3)*el(n3)
  ey = AWY(host,1)*el(n1) + AWY(host,2)*el(n2) + AWY(host,3)*el(n3)
  elp = e0 + ex*x0c + ey*y0c
end subroutine interp_elh

!==============================================================================|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!==============================================================================|

subroutine gauss_random(n, rnd)
  !==============================================================================|
  !  Returns an array of random numbers, conforming to a standard gaussian       |
  !  normal distribution. That is, a normal distribution with zero mean and      |
  !  unit varriance: N(0,1)                                                      |
  !------------------------------------------------------------------------------|
  !  The provided random numbers are generated using George Marsaglia's KISS     |
  !  (Keep It Simple Stupid) pseudo random number generator, implemented by      |
  !  random_number(). The uniformly distributed numbers provided are then fit    |
  !  to a normal distribution via the Marsaglia polar method of the Box-Muller   |
  !  transform.                                                                  |
  !==============================================================================|
  implicit none
  !------------------------------------------------------------------------------|
  integer,            intent(in)  :: n
  real, dimension(n), intent(out) :: rnd
  !------------------------------------------------------------------------------|
  integer :: i
  real    :: s, x, y, temp
  !==============================================================================|

  do i = 1,n,2
    do
      call random_number(x)
      call random_number(y)
      x = (x - 0.5)*2
      y = (y - 0.5)*2
      s = x**2 + y**2
      if (s < 1) exit
    end do

    temp = sqrt(-2.0*log(s)/s)
    rnd(i) = x*temp
    if (i < n) rnd(i + 1) = y*temp
  end do
end subroutine gauss_random

