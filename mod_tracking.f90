module mod_tracking
  !==============================================================================|
  !  Lagrangian particle tracking information for the particles                  |
  !==============================================================================|
  use mod_prec
  use mod_config
  implicit none
  save
  !==============================================================================|
  integer,  allocatable, dimension(:) :: LAG_ID       ! Unique identifier
  integer,  allocatable, dimension(:) :: LAG_HOST     ! Element containing particle
  logical,  allocatable, dimension(:) :: LAG_INDOMAIN ! Particle is in the domain
  integer,  allocatable, dimension(:) :: LAG_RT       ! Release time for the particle (s)
  integer,  allocatable, dimension(:) :: LAG_AGE      ! Age of the particle (s)
  integer,  allocatable, dimension(:) :: LAG_TRACK    ! Lifetime of the particle (s)
  real(DP), allocatable, dimension(:) :: LAG_XP       ! X position of particle (m)
  real(DP), allocatable, dimension(:) :: LAG_YP       ! Y position of particle (m)
  real(DP), allocatable, dimension(:) :: LAG_ZP       ! Z position of particle (sigma)
  real(DP), allocatable, dimension(:) :: LAG_D        ! Initial depth of the particle (m)

  integer :: LAG_TIME ! Curent simulation time (s)
  integer :: LAG_DUR  ! Duration of simulation (s)

  integer :: ISLAG    ! Start record from NetCDF flow-field file
  integer :: IELAG    ! End record from NetCDF flow-field file
contains

  subroutine init_tracking
    !==============================================================================|
    !  Read in lagrangian control parameters and initial lagrangian positions      |
    !==============================================================================|
    use mod_flow_field
    implicit none
    !------------------------------------------------------------------------------|
    integer :: np_out, tdrift
    !==============================================================================|

    !------------------------------------------------------------------------------|
    ! Initialize particle tracking arrays                                          |
    !------------------------------------------------------------------------------|
    allocate(LAG_ID(NDRFT),LAG_HOST(NDRFT))
    allocate(LAG_INDOMAIN(NDRFT),LAG_RT(NDRFT),LAG_AGE(NDRFT),LAG_TRACK(NDRFT))
    allocate(LAG_XP(NDRFT),LAG_YP(NDRFT))
    allocate(LAG_ZP(NDRFT))
    allocate(LAG_D(NDRFT))

    LAG_HOST     = 0
    LAG_AGE      = 0
    LAG_INDOMAIN = .true.

    !------------------------------------------------------------------------------|
    !  Calculate particle tracking start iteration                                 |
    !------------------------------------------------------------------------------|
    ISLAG    = int(DAYST*86400.0_dp)/INSTP + 1

    !------------------------------------------------------------------------------|
    !  Read initial partical positions from seed file                              |
    !------------------------------------------------------------------------------|
    call read_seed

    !------------------------------------------------------------------------------|
    !  Calculate particle tracking end iteration                                   |
    !------------------------------------------------------------------------------|
    LAG_DUR  = maxval(LAG_RT + LAG_TRACK)
    IELAG    = LAG_DUR/INSTP + ISLAG

    !------------------------------------------------------------------------------|
    !  Print statistics on lagrangian tracking to output                           |
    !------------------------------------------------------------------------------|
    write(*,*) "-- Lagrangian Tracking Informaion --"
    write(*,*) "# Points               : ",NDRFT
    write(*,*) "# Start itteration     : ",ISLAG
    write(*,*) "# End itteration       : ",IELAG
    write(*,*) "# Forcing time step    : ",INSTP
    write(*,*) "# Internal time step   : ",DTI

    !------------------------------------------------------------------------------|
    !  Write initial particle positions to output file                             |
    !------------------------------------------------------------------------------|
    call write_track

    return
  end subroutine init_tracking

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine run_tracking
    !==============================================================================|
    !  Update particle positions, calculate scalar fields and particle velocities  |
    !==============================================================================|
    use mod_flow_field
    implicit none
    !------------------------------------------------------------------------------|
    real(DP), dimension(ELEMENTS,SIGLAY) :: u_start, v_start, w_start
    real(DP), dimension(ELEMENTS,SIGLAY) :: u_end, v_end, w_end
    real(DP), dimension(ELEMENTS,SIGLAY) :: u, v, w
    real(DP), dimension(NODES)           :: el, elp
    real(DP), dimension(ELEMENTS,SIGLAY) :: up, vp, wp
    real(DP), dimension(NODES)           :: el_start
    real(DP), dimension(NODES)           :: el_end
    !------------------------------------------------------------------------------|
    integer                   :: record
    integer                   :: np_out
    integer                   :: out_step
    integer                   :: i1, i2, it, p
    real(DP)                  :: tmp1, tmp2
    logical, dimension(NDRFT) :: active
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Get velocity field at the beginning of the simulation (time 0)              |
    !------------------------------------------------------------------------------|
    call read_flow(u_start,v_start,w_start,el_start,ISLAG)
    up = u_start
    vp = v_start
    wp = w_start
    elp = el_start

    !------------------------------------------------------------------------------|
    !  Loop over the tracking period                                               |
    !------------------------------------------------------------------------------|
    write(*,*) "== Running Tracking Simulation =="
    LAG_TIME = 0
    do record = (ISLAG+1),IELAG

      !------------------------------------------------------------------------------|
      !  Read velocity field from NetCDF file                                        |
      !------------------------------------------------------------------------------|
      call read_flow(u_end,v_end,w_end,el_end,record)

      !------------------------------------------------------------------------------|
      ! Loop within one forcing interval                                             |
      !------------------------------------------------------------------------------|
      i1 = 1
      i2 = INSTP/DTI ! caution of time step here
      do it = i1,i2
        LAG_TIME = LAG_TIME + DTI

        tmp2 = float(it - i1 + 1)/float(i2 - i1 + 1)
        tmp1 = float(it - i2)/float(i1 - 1 - i2)

        u  = tmp1*u_start  + tmp2*u_end
        v  = tmp1*v_start  + tmp2*v_end
        w  = tmp1*w_start  + tmp2*w_end
        el = tmp1*el_start + tmp2*el_end

        !------------------------------------------------------------------------------|
        ! Run the tracking simulation                                                  |
        !------------------------------------------------------------------------------|
        call traject(up,u,vp,v,wp,w,elp,el)

        !------------------------------------------------------------------------------|
        ! Write particle records to file                                               |
        !------------------------------------------------------------------------------|
        if (mod(LAG_TIME,DTOUT) == 0) call write_track

        !--Time step update of velocity fields
        up  = u
        vp  = v
        wp  = w
        elp = el
      end do

      !------------------------------------------------------------------------------|
      ! Hourly update of physical fields                                             |
      !------------------------------------------------------------------------------|
      u_start = u_end
      v_start = v_end
      w_start = w_end
      el_start = el_end

      write(*,*) "Finished: ",record,", sim time (s): ",LAG_TIME
    end do
    write(*,*) "== Finished Particle Tracking =="
    return
  end subroutine run_tracking

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine read_seed
    !==============================================================================|
    ! Read the initial partical positions from the input file                      |
    !==============================================================================|
    ! file format (By column, each row defines a diffrent particle)                |
    !                                                                              |
    !   PARAMETER      | TYPE | DESCRIPTION                                        |
    !  ----------------|------|--------------------------------                    |
    !   ID             | INT  | Arbitrary identifier                               |
    !   X              | REAL | Domain co-ordinates (meters)                       |
    !   Y              | REAL | Domain co-ordinates (meters)                       |
    !   Z              | REAL | Particle depth (meters)                            |
    !   Release Delay  | REAL | Delay until the particle is released (seconds)     |
    !   Track Duration | REAL | Time to track the particle for (seconds)           |
    !==============================================================================|
    use mod_flow_field
    implicit none
    !------------------------------------------------------------------------------|
    integer :: i, j, k, nhi
    integer :: inlag
    logical :: fexist
    integer :: stat
    real(DP), dimension(ELEMENTS,SIGLAY) :: u, v, w
    real(DP), dimension(NODES)           :: el
    real(DP), dimension(NDRFT)           :: hp, elp
    real(DP), dimension(NDRFT)           :: rt_in, track_in
    !==============================================================================|

    !------------------------------------------------------------------------------|
    ! Scan in the initial particle position file                                   |
    !------------------------------------------------------------------------------|
    inquire(file=trim(STARTSEED),exist=fexist)
    if(.not.fexist) then
      write(*,*) "ERROR: Lagrangian particle initial position file: "
      write(*,*) trim(STARTSEED)," does not exist"
      stop
    end if

    open(unit=inlag,file=trim(STARTSEED),status='old')
    do i = 1,NDRFT
      read(inlag,*,iostat=stat) LAG_ID(i),LAG_XP(i),LAG_YP(i),LAG_D(i),rt_in(i),track_in(i)
      if (stat /= 0) then
        write(*,*) "ERROR: Lagrangian particle initial position file: ",trim(STARTSEED)
        write(*,*) ,"does not contain ",NDRFT," records."
      stop
      end if
    end do
    close(inlag)

    !------------------------------------------------------------------------------|
    !  Convert external hours to internal seconds                                  |
    !------------------------------------------------------------------------------|
    LAG_RT = int(rt_in*3600.0_dp)
    LAG_TRACK = int(track_in*3600.0_dp)

    !------------------------------------------------------------------------------|
    !  Locate the particle in the domain                                           |
    !------------------------------------------------------------------------------|
    do i = 1,NDRFT
      call fhe_robust(LAG_XP(i),LAG_YP(i),LAG_HOST(i),LAG_INDOMAIN(i))
    end do

    !------------------------------------------------------------------------------|
    !  Shift z co-ordinate to model domain (sigma depth)                           |
    !------------------------------------------------------------------------------|
    call read_flow(u,v,w,el,ISLAG)
    call interp_elh(NDRFT,LAG_HOST,LAG_XP,LAG_YP,el,hp,elp)
    LAG_ZP = LAG_D/(hp + elp)
    if (P_REL_B) LAG_ZP = 1.0_dp - LAG_ZP
    return
  end subroutine read_seed

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine write_track
    !==============================================================================|
    ! Write particle track to output file                                          |
    !==============================================================================|
    ! file format (By column, each row defines a diffrent particle)                |
    !                                                                              |
    !   PARAMETER      | TYPE | DESCRIPTION                                        |
    !  ----------------|------|--------------------------------                    |
    !   ID             | INT  | Arbitrary identifier                               |
    !   X              | REAL | Domain co-ordinates (meters)                       |
    !   Y              | REAL | Domain co-ordinates (meters)                       |
    !   Z              | REAL | Domain co-ordinates (meters)                       |
    !   TIME           | REAL | Time position was recorded (hours)                 |
    !   AGE            | REAL | Duration of particle track (hours)                 |
    !==============================================================================|
    use mod_flow_field
    implicit none
    !------------------------------------------------------------------------------|
    integer :: i
    integer :: outf
    logical :: fexist
    real(DP), dimension(ELEMENTS,SIGLAY) :: u, v, w
    real(DP), dimension(NODES)           :: el
    real(DP), dimension(NDRFT)           :: hp, elp
    real(DP), dimension(NDRFT)           :: z_pos
    !==============================================================================|

    inquire(file=trim(OUTFN),exist=fexist)
    if (fexist) then
      open(outf,file=trim(OUTFN),status="old",position="append")
    else
      open(outf,file=trim(OUTFN),status="new")
    end if

    !------------------------------------------------------------------------------|
    !  Shift z-coordinate to output domain                                         |
    !------------------------------------------------------------------------------|
    z_pos = LAG_ZP
    if (P_REL_B) z_pos = 1.0_dp - z_pos
    if (.not.OUT_SIGMA) then
      i = LAG_TIME/INSTP + ISLAG
      call read_flow(u,v,w,el,i)
      call interp_elh(NDRFT,LAG_HOST,LAG_XP,LAG_YP,el,hp,elp)
      z_pos = z_pos*(hp + elp)
    end if
    
    !------------------------------------------------------------------------------|
    !  Append particle records to the output file                                  |
    !------------------------------------------------------------------------------|
    do i = 1,NDRFT
      write(outf,100) LAG_ID(i),LAG_XP(i),LAG_YP(i),z_pos(i),LAG_TIME,LAG_AGE(i)
    end do

    close(outf)
    return
100 format(i10,f14.2,f14.2,f9.2,i10,i10)
  end subroutine write_track

end module mod_tracking

