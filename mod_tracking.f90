module mod_tracking
  !==============================================================================|
  !  Lagrangian particle tracking information for the particles                  |
  !==============================================================================|
  use mod_prec
  use mod_config
  use mod_flow_field
  implicit none
  save
  !------------------------------------------------------------------------------|
  integer,  allocatable, dimension(:) :: LAG_ID       ! User-supplied particle identifier
  integer,  allocatable, dimension(:) :: LAG_HOST     ! Element containing particle
  logical,  allocatable, dimension(:) :: LAG_INDOMAIN ! Particle is in the domain (Y/N)
  integer,  allocatable, dimension(:) :: LAG_STOP     ! End time for the particle track (s)
  integer,  allocatable, dimension(:) :: LAG_START    ! Release time for the particle (s)
  real,     allocatable, dimension(:) :: LAG_D        ! Initial depth of the particle (m)
  real(DP), allocatable, dimension(:) :: LAG_XP       ! X position of particle (m)
  real(DP), allocatable, dimension(:) :: LAG_YP       ! Y position of particle (m)
  real(DP), allocatable, dimension(:) :: LAG_ZP       ! Z position of particle (sigma)

  integer :: NDRFT         ! Number of particles being tracked

  integer :: SIM_TIME      ! Curent simulation time (s)
  integer :: SIM_END       ! Durration of the simulation (s)
  real    :: SIM_START_MJD ! Simulation real-time start date (mjd)
  !==============================================================================|
contains

  subroutine init_tracking
    !==============================================================================|
    !  Read in lagrangian control parameters and initial lagrangian positions      |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    real :: end_time_mjd
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Initialize the random number generator                                      |
    !------------------------------------------------------------------------------|
    call random_seed()

    !------------------------------------------------------------------------------|
    !  Read initial partical positions from the seed file                          |
    !------------------------------------------------------------------------------|
    call read_seed()

    !------------------------------------------------------------------------------|
    !  Find particle tracking simulated end time                                   |
    !------------------------------------------------------------------------------|
    SIM_TIME = 0
    SIM_END  = maxval(LAG_STOP)
    end_time_mjd = float(SIM_END)/86400.0 + SIM_START_MJD
    if (.not.check_field_time(SIM_START_MJD, end_time_mjd)) then
      write(*,*) "ERROR: Provided times for the particle tracks are not covered by the NetCDF flow-field!"
      stop
    end if

    !------------------------------------------------------------------------------|
    !  Print statistics on lagrangian tracking to output                           |
    !------------------------------------------------------------------------------|
    write(*,*) "-- Lagrangian Tracking Informaion --"
    write(*,*) "# Particles to track   : ",NDRFT
    write(*,*) "# Sim. start time (mjd): ",SIM_START_MJD
    write(*,*) "# Sim. end time (mjd)  : ",end_time_mjd
    write(*,*) "# Time step (s)        : ",DTI
  end subroutine init_tracking

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine run_tracking
    !==============================================================================|
    !  Update particle positions, calculate scalar fields and particle velocities  |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    real(DP), dimension(ELEMENTS,SIGLAY) :: u, up
    real(DP), dimension(ELEMENTS,SIGLAY) :: v, vp
    real(DP), dimension(ELEMENTS,SIGLAY) :: w, wp
    real(DP), dimension(NODES)           :: el, elp
    real(DP), dimension(:), allocatable  :: h_part, el_part
    !==============================================================================|
    write(*,*) "== Running Tracking Simulation =="

    !------------------------------------------------------------------------------|
    !  Get velocity field at the beginning of the simulation (time 0)              |
    !------------------------------------------------------------------------------|
    call get_flow_at(up,vp,wp,elp,SIM_START_MJD)

    !------------------------------------------------------------------------------|
    !  Calculate the particle's initial sigma position                             |
    !------------------------------------------------------------------------------|
    allocate(h_part(NDRFT))
    allocate(el_part(NDRFT))
    call interp_elh(NDRFT,LAG_HOST,LAG_XP,LAG_YP,elp,h_part,el_part)
    LAG_ZP = 0.0_dp - LAG_D/(h_part + el_part)
    deallocate(h_part)
    deallocate(el_part)

    !------------------------------------------------------------------------------|
    !  Write initial particle positions to output file                             |
    !------------------------------------------------------------------------------|
    call write_track(elp)

    !------------------------------------------------------------------------------|
    !  Loop over the tracking period                                               |
    !------------------------------------------------------------------------------|
    do
      if (SIM_TIME > SIM_END) exit
      SIM_TIME = SIM_TIME + DTI

      !------------------------------------------------------------------------------|
      !  Get the velocity field for the current time step                            |
      !------------------------------------------------------------------------------|
      call get_flow_at(u,v,w,el,SIM_START_MJD + float(SIM_TIME)/86400.0)

      !------------------------------------------------------------------------------|
      !  Run the tracking simulation                                                 |
      !------------------------------------------------------------------------------|
      call traject(up,u,vp,v,wp,w,elp,el)
      if (P_RND_WALK) call random_walk(el)

      !------------------------------------------------------------------------------|
      !  Write particle records to file                                              |
      !------------------------------------------------------------------------------|
      if (mod(SIM_TIME,DTOUT) == 0) call write_track(el)

      !------------------------------------------------------------------------------|
      !  Time step update of velocity fields                                         |
      !------------------------------------------------------------------------------|
      up  = u
      vp  = v
      wp  = w
      elp = el
    end do

    write(*,*) "== Finished Particle Tracking =="
  end subroutine run_tracking

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine read_seed
    !==============================================================================|
    !  Read the initial partical positions from the input file                     |
    !==============================================================================|
    !  file format (By column, each row defines a diffrent particle)               |
    !                                                                              |
    !   PARAMETER      | TYPE | DESCRIPTION                                        |
    !  ----------------|------|--------------------------------                    |
    !   ID             | INT  | Arbitrary identifier                               |
    !   X              | REAL | Domain co-ordinates (meters)                       |
    !   Y              | REAL | Domain co-ordinates (meters)                       |
    !   Z              | REAL | Particle depth (meters)                            |
    !   Release Time   | REAL | Time to release particle into the simulation (mjd) |
    !   Track End Time | REAL | Time to remove particle from the simulation (mjd)  |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    integer                         :: i
    integer                         :: inlag
    integer                         :: stat
    logical                         :: fexist
    character(len=80)               :: junk
    real, allocatable, dimension(:) :: start_in, stop_in
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Count the number of particles being tracked                                 |
    !------------------------------------------------------------------------------|
    inquire(file=trim(STARTSEED),exist=fexist)
    if(.not.fexist) then
      write(*,*) "ERROR: Lagrangian particle initial position file: "
      write(*,*) trim(STARTSEED)," does not exist"
      stop
    end if

    NDRFT = 0
    open(unit=inlag,file=trim(STARTSEED),status='old')
    do
      read(inlag,*,iostat=stat) junk
      if (stat /= 0) exit
      if (len_trim(junk) > 0) NDRFT = NDRFT + 1
    end do
    rewind(inlag)

    !------------------------------------------------------------------------------|
    !  Initialize particle tracking arrays                                         |
    !------------------------------------------------------------------------------|
    allocate(LAG_ID(NDRFT),LAG_HOST(NDRFT))
    allocate(LAG_INDOMAIN(NDRFT),LAG_START(NDRFT),LAG_STOP(NDRFT))
    allocate(LAG_XP(NDRFT),LAG_YP(NDRFT))
    allocate(LAG_ZP(NDRFT))
    allocate(LAG_D(NDRFT))
    allocate(start_in(NDRFT), stop_in(NDRFT))

    LAG_HOST     = 0
    LAG_INDOMAIN = .true.

    !------------------------------------------------------------------------------|
    !  Scan in the initial particle position file                                  |
    !------------------------------------------------------------------------------|
    do i = 1,NDRFT
      read(inlag,*,iostat=stat) LAG_ID(i),LAG_XP(i),LAG_YP(i),LAG_D(i),start_in(i),stop_in(i)
      if (stat /= 0) then
        write(*,*) "ERROR: Problem reading the initial position file: ",trim(STARTSEED)
        stop
      end if
    end do
    close(inlag)

    !------------------------------------------------------------------------------|
    !  Convert external mjd to internal seconds                                    |
    !------------------------------------------------------------------------------|
    SIM_START_MJD = minval(start_in)
    LAG_START = int((start_in - SIM_START_MJD)*86400.0)
    LAG_STOP  = int((stop_in - SIM_START_MJD)*86400.0)

    !------------------------------------------------------------------------------|
    !  Locate the particle in the domain                                           |
    !------------------------------------------------------------------------------|
    do i = 1,NDRFT
      call fhe_robust(LAG_XP(i),LAG_YP(i),LAG_HOST(i),LAG_INDOMAIN(i))
    end do
  end subroutine read_seed

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine write_track(el)
    !==============================================================================|
    !  Write particle track to output file                                         |
    !==============================================================================|
    !  file format (By column, each row defines a diffrent particle)               |
    !                                                                              |
    !   PARAMETER      | TYPE | DESCRIPTION                                        |
    !  ----------------|------|--------------------------------                    |
    !   ID             | INT  | Arbitrary identifier                               |
    !   X              | REAL | Domain co-ordinates (meters)                       |
    !   Y              | REAL | Domain co-ordinates (meters)                       |
    !   Z              | REAL | Particle depth (meters)                            |
    !   TIME           | REAL | Time of particle record (mjd)                      |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    real(DP), dimension(NODES), intent(in) :: el
    !------------------------------------------------------------------------------|
    integer                    :: i
    integer                    :: outf
    logical                    :: fexist
    real(DP), dimension(NDRFT) :: hp, elp
    real(DP), dimension(NDRFT) :: z_pos
    real                       :: sim_time_mjd
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open the output file for writing                                            |
    !------------------------------------------------------------------------------|
    inquire(file=trim(OUTFN),exist=fexist)
    if (fexist) then
      open(outf,file=trim(OUTFN),status="old",position="append")
    else
      open(outf,file=trim(OUTFN),status="new")
    end if

    !------------------------------------------------------------------------------|
    !  Shift z-coordinate to output domain                                         |
    !------------------------------------------------------------------------------|
    call interp_elh(NDRFT,LAG_HOST,LAG_XP,LAG_YP,el,hp,elp)
    if (F_DEPTH) then
      z_pos = 0.0_dp - LAG_D/(hp + elp)
    else
      z_pos = LAG_ZP
    end if
    if (P_REL_B) z_pos = -1.0_dp - z_pos
    if (.not.OUT_SIGMA) z_pos = 0.0_dp - z_pos*(hp + elp)

    !------------------------------------------------------------------------------|
    !  Append particle records to the output file                                  |
    !------------------------------------------------------------------------------|
    sim_time_mjd = float(SIM_TIME)/86400.0 + SIM_START_MJD
    do i = 1,NDRFT
      write(outf,100) LAG_ID(i),LAG_XP(i),LAG_YP(i),z_pos(i),sim_time_mjd
    end do

    close(outf)
    return
100 format(i10,f14.2,f14.2,f9.2,f12.4)
  end subroutine write_track

end module mod_tracking

