!==============================================================================|
!              Lagrangian particle tracking off-line program                   |
!==============================================================================|
!  Usage:                                                                      |
!    PTrack run.dat                                                            |
!------------------------------------------------------------------------------|
!  Where run.dat is the name of the run configuration file                     |
!==============================================================================|

program particle_traj
  !==============================================================================|
  use mod_config
  use mod_flow_field
  use mod_tracking
  implicit none
  !------------------------------------------------------------------------------|
  logical               :: fexist
  character(len=12)     :: dmy1, dmy2, dmy3
  integer, dimension(8) :: stc, ndc
  real                  :: duration
  !==============================================================================|

  !------------------------------------------------------------------------------|
  !  Get CASENAME from command line                                              |
  !------------------------------------------------------------------------------|
  call getarg(1,CASENAME)

  inquire(file=CASENAME,exist=fexist)
  if(.not.fexist)then
    write(*,*) 'ERROR: Parameter file: ',CASENAME,' does not exist'
    stop
  end if

  !------------------------------------------------------------------------------|
  !  Read configuration parameters controlling model run                         |
  !------------------------------------------------------------------------------|
  call init_model

  !------------------------------------------------------------------------------|
  !  Read domain information from NetCDF input file                              |
  !------------------------------------------------------------------------------|
  call ncd_read_grid

  write(*,*) "-- Domain Information --"
  write(*,*) "# Nodes               : ",NODES
  write(*,*) "# Elements            : ",ELEMENTS
  write(*,*) "# Sigma layers        : ",SIGLAY

  !------------------------------------------------------------------------------|
  !  Run the lagrangian particle tracking model                                  |
  !------------------------------------------------------------------------------|
  call init_tracking
  call run_tracking

  close(0)
end program particle_traj

