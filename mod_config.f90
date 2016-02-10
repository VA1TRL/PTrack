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

  logical :: F_DEPTH    ! Run simulation holding particle depth constant
  logical :: P_REL_B    ! Particle positions relative to the bottom (instead of surface)
  logical :: OUT_SIGMA  ! Output particle z position as sigma depth instead of meters
  logical :: P_RND_WALK ! Apply a random walk behaviour to the active particles
  logical :: P_2D_MODEL ! Use a 2D, vertically averaged, flow field

  real    :: K_XY ! Horizontal particle diffusivity
  real    :: K_Z  ! Vertical particle diffusivity

  character(len=80) :: CASENAME   ! Name of current run configuration file
  character(len=80) :: OUTFN      ! Name of output data file
  character(len=80) :: GRIDFN     ! Name of NetCDF flow-field/grid input file
  character(len=80) :: SEEDFN     ! Name of particle seed (initial location) file
contains

  subroutine init_model
    !==============================================================================|
    !  Initialize model configuration parameters from CASENAME                     |
    !==============================================================================|
    use mod_inp
    implicit none
    !------------------------------------------------------------------------------|
    integer :: iscan
    !==============================================================================|

    iscan = scan_file(CASENAME, "DTI", iscal = DTI)
    if (iscan /= 0) then
      write(*,*) "ERROR reading DTI from: ", CASENAME
      call pscanmsg(iscan)
      stop
    end if

    iscan = scan_file(CASENAME, "DTOUT", iscal = DTOUT)
    if (iscan /= 0) then
      write(*,*) "ERROR reading DTOUT from: ", CASENAME
      call pscanmsg(iscan)
      stop
    end if

    iscan = scan_file(CASENAME, "GRIDFN", cval = GRIDFN)
    if (iscan /= 0) then
      write(*,*) "ERROR reading GRIDFN from: ", CASENAME
      call pscanmsg(iscan)
      stop
    end if

    iscan = scan_file(CASENAME, "OUTFN", cval = OUTFN)
    if (iscan /= 0) then
      write(*,*) "ERROR reading OUTFN from: ", CASENAME
      call pscanmsg(iscan)
      stop
    end if

    iscan = scan_file(CASENAME, "SEEDFN", cval = SEEDFN)
    if (iscan /= 0) then
      write(*,*) "ERROR reading SEEDFN from: ", CASENAME
      call pscanmsg(iscan)
      stop
    end if

    iscan = scan_file(CASENAME, "F_DEPTH", lval = F_DEPTH)
    if (iscan /= 0) F_DEPTH = .false.

    iscan = scan_file(CASENAME, "P_REL_B", lval = P_REL_B)
    if (iscan /= 0) P_REL_B = .false.

    iscan = scan_file(CASENAME, "OUT_SIGMA", lval = OUT_SIGMA)
    if (iscan /= 0) OUT_SIGMA = .false.

    iscan = scan_file(CASENAME, "P_RND_WALK", lval = P_RND_WALK)
    if (iscan /= 0) P_RND_WALK = .false.

    iscan = scan_file(CASENAME, "P_2D_MODEL", lval = P_2D_MODEL)
    if (iscan /= 0) P_2D_MODEL = .false.
    if (P_2D_MODEL) F_DEPTH = .false.

    if (P_RND_WALK) then
      iscan = scan_file(CASENAME, "K_XY", fscal = K_XY)
      if (iscan /= 0) then
        write(*,*) "ERROR reading K_XY: ", CASENAME
        call pscanmsg(iscan)
        stop
      end if

      iscan = scan_file(CASENAME, "K_Z", fscal = K_Z)
      if (iscan /= 0) then
        write(*,*) "ERROR reading K_Z from: ", CASENAME
        call pscanmsg(iscan)
        stop
      end if
    end if
  end subroutine init_model

end module mod_config

