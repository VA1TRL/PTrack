module mod_config
  !==============================================================================|
  !  User parameters and control variables                                       |
  !==============================================================================|
  use mod_prec
  implicit none
  save
  !==============================================================================|
  integer :: DTI   ! Time step for particle track resolution (seconds)
  integer :: DTOUT ! Time step for output records (seconds)

  logical :: F_DEPTH
  logical :: P_REL_B
  logical :: OUT_SIGMA
  logical :: P_RND_WALK

  real    :: K_XY
  real    :: K_Z

  character(len=80) :: CASENAME   ! Name of curent run configuration file
  character(len=80) :: OUTFN      ! Name of output data file
  character(len=80) :: GRIDFN     ! Name of NetCDF flow-field/grid input file
  character(len=80) :: SEEDFN     ! Name of particle seed (initial location) file
contains

  subroutine init_model
    !==============================================================================|
    !  Initalize model configuration parameters from CASENAME                      |
    !==============================================================================|
    use mod_inp
    implicit none
    !------------------------------------------------------------------------------|
    integer :: iscan
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  DTI : Internal simulation time step (seconds)                               |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME, "DTI", iscal = DTI)
    if (iscan /= 0) then
      write(*,*) "ERROR reading DTI from: ", CASENAME
      call pscanmsg(iscan)
      stop
    end if

    !------------------------------------------------------------------------------|
    !  DTOUT : Output interval (seconds)                                           |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME, "DTOUT", iscal = DTOUT)
    if (iscan /= 0) then
      write(*,*) "ERROR reading DTOUT from: ", CASENAME
      call pscanmsg(iscan)
      stop
    end if

    !------------------------------------------------------------------------------|
    !  F_DEPTH : Run simulation holding particle depth constant                    |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME, "F_DEPTH", lval = F_DEPTH)
    if (iscan /= 0) then
      F_DEPTH = .false.
    end if

    !------------------------------------------------------------------------------|
    !  P_REL_B : Particle positions relitive to the bottom (instead of surface)    |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME, "P_REL_B", lval = P_REL_B)
    if (iscan /= 0) then
      P_REL_B = .false.
    end if

    !------------------------------------------------------------------------------|
    !  OUT_SIGMA : Output particle z position as sigma depth isnstead of meters    |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME, "OUT_SIGMA", lval = OUT_SIGMA)
    if (iscan /= 0) then
      OUT_SIGMA = .false.
    end if

    !------------------------------------------------------------------------------|
    !  P_RND_WALK : Apply a random walk behaviour to the active particles          |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME, "P_RND_WALK", lval = P_RND_WALK)
    if (iscan /= 0) then
      P_RND_WALK = .false.
    end if

    !------------------------------------------------------------------------------|
    !  GRIDFN : NetCDF input filename, containing grid and flow field data         |     
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME, "GRIDFN", cval = GRIDFN)
    if (iscan /= 0) then
      write(*,*) "ERROR reading GRIDFN from: ", CASENAME
      call pscanmsg(iscan)
      stop
    end if

    !------------------------------------------------------------------------------|
    !  OUTFN : Output file name                                                    |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME, "OUTFN", cval = OUTFN)
    if (iscan /= 0) then
      write(*,*) "ERROR reading OUTFN from: ", CASENAME
      call pscanmsg(iscan)
      stop
    end if

    !------------------------------------------------------------------------------|
    !  SEEDFN : Partical location input filename                                   |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME, "SEEDFN", cval = SEEDFN)
    if (iscan /= 0) then
      write(*,*) "ERROR reading SEEDFN from: ", CASENAME
      call pscanmsg(iscan)
      stop
    end if

    if (P_RND_WALK) then
      !------------------------------------------------------------------------------|
      !  K_XY : Horizontal particle diffusivity                                      |
      !------------------------------------------------------------------------------|
      iscan = scan_file(CASENAME, "K_XY", fscal = K_XY)
      if (iscan /= 0) then
        write(*,*) "ERROR reading K_XY: ", CASENAME
        call pscanmsg(iscan)
        stop
      end if

      !------------------------------------------------------------------------------|
      !  K_Z : Vertical particle diffusivity                                         |
      !------------------------------------------------------------------------------|
      iscan = scan_file(CASENAME, "K_Z", fscal = K_Z)
      if (iscan /= 0) then
        write(*,*) "ERROR reading K_Z from: ", CASENAME
        call pscanmsg(iscan)
        stop
      end if
    end if
  end subroutine init_model

end module mod_config

