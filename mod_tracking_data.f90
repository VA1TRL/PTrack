module mod_tracking_data
  !==============================================================================|
  !  Lagrangian particle tracking information for the particles                  |
  !==============================================================================|
  use mod_io
  use mod_config
  implicit none
  save
  !------------------------------------------------------------------------------|
  integer :: NDRFT         ! Number of particles being tracked
  real    :: SIM_START_MJD ! Simulation real-time start date (mjd)

  integer,  allocatable, dimension(:) :: LAG_ID       ! User-supplied particle identifier
  integer,  allocatable, dimension(:) :: LAG_HOST     ! Element containing particle
  logical,  allocatable, dimension(:) :: LAG_INDOMAIN ! Particle is in the domain (Y/N)
  integer,  allocatable, dimension(:) :: LAG_STOP     ! End time for the particle track (s)
  integer,  allocatable, dimension(:) :: LAG_START    ! Release time for the particle (s)
  real,     allocatable, dimension(:) :: LAG_D        ! Initial depth of the particle (m)
  real(DP), allocatable, dimension(:) :: LAG_XP       ! X position of particle (m)
  real(DP), allocatable, dimension(:) :: LAG_YP       ! Y position of particle (m)
  real(DP), allocatable, dimension(:) :: LAG_ZP       ! Z position of particle (sigma)
  !==============================================================================|
contains

  subroutine read_seed
    !==============================================================================|
    !  Read the initial particle positions from the input file                     |
    !==============================================================================|
    !   VARIABLE       | TYPE | DESCRIPTION                                        |
    !  ----------------|------|--------------------------------                    |
    !   number         | INT  | Arbitrary identifier                               |
    !   x              | REAL | Domain co-ordinates (meters)                       |
    !   y              | REAL | Domain co-ordinates (meters)                       |
    !   z              | REAL | Particle depth (meters)                            |
    !   release        | REAL | Time to release particle into the simulation (mjd) |
    !   end            | REAL | Time to remove particle from the simulation (mjd)  |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    integer                         :: fileid
    integer                         :: i
    real, allocatable, dimension(:) :: start_in, stop_in, xp, yp
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Retrieve the number of particles in the tracking simulation                 |
    !------------------------------------------------------------------------------|
    call nc_open_file(SEEDFN, .false., fileid)
    call nc_dim(fileid, "number", NDRFT)

    !------------------------------------------------------------------------------|
    !  Initialize particle tracking arrays                                         |
    !------------------------------------------------------------------------------|
    allocate(LAG_ID(NDRFT), LAG_HOST(NDRFT))
    allocate(LAG_INDOMAIN(NDRFT), LAG_START(NDRFT), LAG_STOP(NDRFT))
    allocate(LAG_XP(NDRFT), LAG_YP(NDRFT))
    allocate(LAG_ZP(NDRFT))
    allocate(LAG_D(NDRFT))
    allocate(start_in(NDRFT), stop_in(NDRFT))
    allocate(xp(NDRFT), yp(NDRFT))

    LAG_HOST     = 0
    LAG_INDOMAIN = .true.

    !------------------------------------------------------------------------------|
    !  Read in the initial particle position file                                  |
    !------------------------------------------------------------------------------|
    call nc_read_var(fileid, "number",  LAG_ID)
    call nc_read_var(fileid, "x",       xp)
    call nc_read_var(fileid, "y",       yp)
    call nc_read_var(fileid, "z",       LAG_D)
    call nc_read_var(fileid, "release", start_in)
    call nc_read_var(fileid, "end",     stop_in)
    call nc_close_file(fileid)

    LAG_XP = xp
    LAG_YP = yp

    !------------------------------------------------------------------------------|
    !  Convert external mjd to internal seconds                                    |
    !------------------------------------------------------------------------------|
    SIM_START_MJD = minval(start_in)
    LAG_START     = int((start_in - SIM_START_MJD)*86400.0)
    LAG_STOP      = int((stop_in - SIM_START_MJD)*86400.0)

    !------------------------------------------------------------------------------|
    !  Locate the domain element the particles reside in                           |
    !------------------------------------------------------------------------------|
    do i = 1,NDRFT
      call fhe_robust(LAG_XP(i), LAG_YP(i), LAG_HOST(i), LAG_INDOMAIN(i))
    end do
  end subroutine read_seed

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine write_track(time, el)
    !==============================================================================|
    !  Write particle track to output file                                         |
    !==============================================================================|
    !    VARIABLE      | TYPE | DESCRIPTION                                        |
    !  ----------------|------|--------------------------------                    |
    !   time           | REAL | Time of particle record (mjd)                      |
    !   number         | INT  | Arbitrary identifier                               |
    !   x              | REAL | Domain co-ordinates (meters)                       |
    !   y              | REAL | Domain co-ordinates (meters)                       |
    !   z              | REAL | Particle depth (meters)                            |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    integer,                          intent(in) :: time
    real(DP), dimension(:), optional, intent(in) :: el
    !------------------------------------------------------------------------------|
    integer                             :: outid, i
    integer                             :: records
    real(DP), dimension(:), allocatable :: hp, elp
    real(DP), dimension(NDRFT)          :: z_pos
    real                                :: sim_time_mjd
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open the output file for writing                                            |
    !------------------------------------------------------------------------------|
    call nc_open_file(OUTFN, .true., outid)
    call nc_dim(outid, "time", records)

    !------------------------------------------------------------------------------|
    !  Shift z-coordinate to output domain                                         |
    !------------------------------------------------------------------------------|
    if (present(el)) then
      allocate(hp(NDRFT), elp(NDRFT))
      call interp_elh(NDRFT, LAG_HOST, LAG_XP, LAG_YP, el, hp, elp)
      if (F_DEPTH) then
        z_pos = 0.0_dp - LAG_D/(hp + elp)
      else
        z_pos = LAG_ZP
      end if
      if (P_REL_B) z_pos = -1.0_dp - z_pos
      z_pos = 0.0_dp - z_pos*(hp + elp)
      deallocate(hp, elp)
    else if (OUT_SIGMA) then
      z_pos = LAG_ZP
    else
      z_pos = 0.0_dp
    end if

    !------------------------------------------------------------------------------|
    !  Append particle records to the output file                                  |
    !------------------------------------------------------------------------------|
    records = records + 1
    sim_time_mjd = float(time)/86400.0 + SIM_START_MJD
    call nc_1d_write(outid, "time", records, sim_time_mjd)
    call nc_2d_write(outid, "x", records, NDRFT, real(LAG_XP))
    call nc_2d_write(outid, "y", records, NDRFT, real(LAG_YP))
    call nc_2d_write(outid, "z", records, NDRFT, real(z_pos))
    call nc_close_file(outid)
  end subroutine write_track

end module mod_tracking_data

