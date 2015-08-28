module mod_config
  !==============================================================================|
  !  User parameters and control variables                                       |
  !==============================================================================|
  use mod_prec
  implicit none
  save
  !==============================================================================|
  integer :: DTI     ! Time step for particle track resolution (seconds)
  integer :: DTOUT   ! Time step for output records (seconds)

  logical :: F_DEPTH
  logical :: P_SIGMA
  logical :: P_REL_B
  logical :: OUT_SIGMA

  character(len=80) :: CASENAME   ! Name of curent run configuration file
  character(len=80) :: OUTFN      ! Name of output data file
  character(len=80) :: GRIDFN     ! Name of NetCDF flow-field/grid input file
  character(len=80) :: STARTSEED  ! Name of particle seed (initial location) file
contains

  subroutine init_model
    !==============================================================================|
    !  Initalize model configuration parameters from CASENAME                      |
    !==============================================================================|
    use mod_inp
    implicit none
    !------------------------------------------------------------------------------|
    integer            :: iscan
    integer            :: i
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  DTI : Internal simulation time step (seconds)                               |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME,"DTI",iscal = DTI)
    if(iscan /= 0) then
      write(*,*) "ERROR reading DTI from: ", CASENAME
      i = PScanMsg(iscan)
      stop 
    end if

    !------------------------------------------------------------------------------|
    !  DTOUT : Output interval (seconds)                                           |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME,"DTOUT",iscal = DTOUT)
    if(iscan /= 0) then
      write(*,*) "ERROR reading DTOUT from: ", CASENAME
      i = PScanMsg(iscan)
      stop 
    end if

    !------------------------------------------------------------------------------|
    !  P_SIGMA : Run vertical simulation over sigma levels, instead of meters      |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME,"P_SIGMA",lval = P_SIGMA)
    if(iscan /= 0) then
      write(*,*) "ERROR reading P_SIGMA from: ", CASENAME
      i = PScanMsg(iscan)
      stop 
    end if

    !------------------------------------------------------------------------------|
    !  F_DEPTH : Run simulation holding particle depth constant                    |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME,"F_DEPTH",lval = F_DEPTH)
    if(iscan /= 0) then
      write(*,*) "ERROR reading F_DEPTH from: ", CASENAME
      i = PScanMsg(iscan)
      stop 
    end if

    !------------------------------------------------------------------------------|
    !  P_REL_B : Particle positions relitive to the bottom (instead of surface)    |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME,"P_REL_B",lval = P_REL_B)
    if(iscan /= 0) then
      write(*,*) "ERROR reading P_REL_B from: ", CASENAME
      i = PScanMsg(iscan)
      stop 
    end if

    !------------------------------------------------------------------------------|
    !  OUT_SIGMA : Output particle z position as sigma depth isnstead of meters    |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME,"OUT_SIGMA",lval = OUT_SIGMA)
    if(iscan /= 0) then
      write(*,*) "ERROR reading OUT_SIGMA from: ", CASENAME
      i = PScanMsg(iscan)
      stop 
    end if

    !------------------------------------------------------------------------------|
    !  GRIDFN : NetCDF input filename, containing grid and flow field data         |     
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME,"GRIDFN",cval = GRIDFN)
    if(iscan /= 0) then
      write(*,*) "ERROR reading GRIDFN from: ", CASENAME
      i = PScanMsg(iscan)
      stop 
    end if

    !------------------------------------------------------------------------------|
    !  OUTFN : Output file name                                                    |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME,"OUTFN",cval = OUTFN)
    if(iscan /= 0) then
      write(*,*) "ERROR reading OUTFN from: ", CASENAME
      i = PScanMsg(iscan)
      stop 
    end if

    !------------------------------------------------------------------------------|
    !  STARTSEED : Partical location input filename                                |
    !------------------------------------------------------------------------------|
    iscan = scan_file(CASENAME,"STARTSEED",cval = STARTSEED)
    if(iscan /= 0) then
      write(*,*) "ERROR reading STARTSEED from: ", CASENAME
      i = PScanMsg(iscan)
      stop 
    end if
  end subroutine init_model

end module mod_config

