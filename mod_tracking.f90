module mod_tracking
  !==============================================================================|
  !  Lagrangian particle tracking information for the particles                  |
  !==============================================================================|
  use mod_config
  use mod_flow_field
  use mod_tracking_data
  implicit none
  save
  !------------------------------------------------------------------------------|
  integer :: SIM_TIME      ! Current simulation time (s)
  integer :: SIM_END       ! Duration of the simulation (s)
  !==============================================================================|
contains

  subroutine init_tracking
    !==============================================================================|
    !  Read in Lagrangian control parameters and initial Lagrangian positions      |
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
    !  Read initial particle positions from the seed file                          |
    !------------------------------------------------------------------------------|
    call read_seed()
    call nc_new_outfile(OUTFN, NDRFT, LAG_ID, GRIDFN)

    !------------------------------------------------------------------------------|
    !  Find particle tracking simulated end time                                   |
    !------------------------------------------------------------------------------|
    SIM_TIME     = 0
    SIM_END      = maxval(LAG_STOP)
    end_time_mjd = float(SIM_END)/86400.0 + SIM_START_MJD
    if (.not.check_field_time(SIM_START_MJD, end_time_mjd)) then
      write(*,*) "ERROR: Provided times for the particle tracks are not covered by the NetCDF flow-field!"
      stop
    end if

    !------------------------------------------------------------------------------|
    !  Print statistics on lagrangian tracking to output                           |
    !------------------------------------------------------------------------------|
    write(*,*) "-- Lagrangian Tracking Information --"
    write(*,*) "Particles to track   : ", NDRFT
    write(*,*) "Sim. start time (mjd): ", SIM_START_MJD
    write(*,*) "Sim. end time (mjd)  : ", end_time_mjd
    write(*,*) "Time step (s)        : ", DTI
    write(*,*) "Output time step (s) : ", DTOUT
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
    write(*,*) "== Running 3D Tracking Simulation =="

    !------------------------------------------------------------------------------|
    !  Get velocity field at the beginning of the simulation (time 0)              |
    !------------------------------------------------------------------------------|
    call get_flow_at(SIM_START_MJD, up, vp, wp, elp)

    !------------------------------------------------------------------------------|
    !  Calculate the particle's initial sigma position                             |
    !------------------------------------------------------------------------------|
    allocate(h_part(NDRFT))
    allocate(el_part(NDRFT))
    call interp_elh(NDRFT, LAG_HOST, LAG_XP, LAG_YP, elp, h_part, el_part)
    LAG_ZP = 0.0_dp - LAG_D/(h_part + el_part)
    deallocate(h_part)
    deallocate(el_part)

    !------------------------------------------------------------------------------|
    !  Write initial particle positions to output file                             |
    !------------------------------------------------------------------------------|
    call write_track(SIM_TIME, elp)

    !------------------------------------------------------------------------------|
    !  Loop over the tracking period                                               |
    !------------------------------------------------------------------------------|
    do
      if (SIM_TIME >= SIM_END) exit
      SIM_TIME = SIM_TIME + DTI

      !------------------------------------------------------------------------------|
      !  Get the velocity field for the current time step                            |
      !------------------------------------------------------------------------------|
      call get_flow_at(SIM_START_MJD + float(SIM_TIME)/86400.0, u, v, w, el)

      !------------------------------------------------------------------------------|
      !  Run the tracking simulation                                                 |
      !------------------------------------------------------------------------------|
      call traject(SIM_TIME, up, u, vp, v, wp, w, elp, el)
      if (P_RND_WALK) call random_walk(SIM_TIME, el)

      !------------------------------------------------------------------------------|
      !  Write particle records to file                                              |
      !------------------------------------------------------------------------------|
      if (mod(SIM_TIME,DTOUT) == 0) call write_track(SIM_TIME, el)

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

  subroutine run_2d_tracking
    !==============================================================================|
    !  Update particle positions, calculate scalar fields and particle velocities  |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    real(DP), dimension(ELEMENTS,1) :: u, up
    real(DP), dimension(ELEMENTS,1) :: v, vp
    !==============================================================================|
    write(*,*) "== Running 2D Tracking Simulation =="

    !------------------------------------------------------------------------------|
    !  Get velocity field at the beginning of the simulation (time 0)              |
    !------------------------------------------------------------------------------|
    call get_flow_at(SIM_START_MJD, up, vp)

    !------------------------------------------------------------------------------|
    !  Write initial particle positions to output file                             |
    !------------------------------------------------------------------------------|
    call write_track(SIM_TIME)

    !------------------------------------------------------------------------------|
    !  Loop over the tracking period                                               |
    !------------------------------------------------------------------------------|
    do
      if (SIM_TIME >= SIM_END) exit
      SIM_TIME = SIM_TIME + DTI

      !------------------------------------------------------------------------------|
      !  Get the velocity field for the current time step                            |
      !------------------------------------------------------------------------------|
      call get_flow_at(SIM_START_MJD + float(SIM_TIME)/86400.0, u, v)

      !------------------------------------------------------------------------------|
      !  Run the tracking simulation                                                 |
      !------------------------------------------------------------------------------|
      call traject_2d(SIM_TIME, reshape(up,[ELEMENTS]), reshape(u,[ELEMENTS]), reshape(vp,[ELEMENTS]), reshape(v,[ELEMENTS]))
      if (P_RND_WALK) call random_walk(SIM_TIME)

      !------------------------------------------------------------------------------|
      !  Write particle records to file                                              |
      !------------------------------------------------------------------------------|
      if (mod(SIM_TIME,DTOUT) == 0) call write_track(SIM_TIME)

      !------------------------------------------------------------------------------|
      !  Time step update of velocity fields                                         |
      !------------------------------------------------------------------------------|
      up = u
      vp = v
    end do

    write(*,*) "== Finished Particle Tracking =="
  end subroutine run_2d_tracking

end module mod_tracking

