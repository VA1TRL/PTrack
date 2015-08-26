module mod_flow_field
  !==============================================================================|
  !  Grid/Flow-Field data storage and access functions                           |
  !==============================================================================|
  use mod_prec
  use mod_config
  use netcdf
  implicit none
  save
  !==============================================================================|

  !------------------------------------------------------------------------------|
  !  Grid metrics                                                                |
  !------------------------------------------------------------------------------|
  integer :: ELEMENTS  ! Number of elements (triangles) in the mesh
  integer :: NODES     ! Number of nodes (triangle corners) in the mesh
  integer :: SIGLEV    ! Number of sigma levels
  integer :: SIGLAY    ! Number of sigma layers
  integer :: N_TIME    ! Number of time records in NetCDF file

  !------------------------------------------------------------------------------|
  !  Grid co-ordinates                                                           |
  !------------------------------------------------------------------------------|
  real,    allocatable, dimension(:)   :: XP, YP       ! Node (x,y) co-ordinates (meters)
  real,    allocatable, dimension(:)   :: XC, YC       ! Element, triangle center, (x,y) co-ordinates (meters)
  real,    allocatable, dimension(:)   :: Z_LEV, Z_LAY ! Sigma co-ordinates
  real,    allocatable, dimension(:)   :: H            ! Bathymetric depth at node (meters)
  integer, allocatable, dimension(:)   :: NTVE         ! Number of elements sourounding each node
  integer, allocatable, dimension(:,:) :: NV           ! Corners nodes of the elements/grid triangles
  integer, allocatable, dimension(:,:) :: NBE          ! Elements sourounding each element
  integer, allocatable, dimension(:,:) :: NBVE         ! Elements sourounding each node

  !------------------------------------------------------------------------------|
  !  Shape coefficient arrays and control volume metrics                         |
  !------------------------------------------------------------------------------|
  real, allocatable, dimension(:,:) :: A1U
  real, allocatable, dimension(:,:) :: A2U
  real, allocatable, dimension(:,:) :: AWX
  real, allocatable, dimension(:,:) :: AWY
  real, allocatable, dimension(:,:) :: AW0
contains

  subroutine ncd_read_grid
    !==============================================================================|
    !  Read grid shape/dimmension information from NetCDF file                     |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    integer                           :: ierr
    integer                           :: dimid, varid, fid
    integer                           :: maxelem
    character(len=NF90_MAX_NAME)      :: temp
    real, allocatable, dimension(:,:) :: tmp
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    ierr = nf90_open(trim(GRIDFN),NF90_NOWRITE,fid)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Get model dimensions                                                        |
    !------------------------------------------------------------------------------|
    ierr = nf90_inq_dimid(fid,"nele",dimid)
    call handle_ncerror(ierr)
    ierr = nf90_inquire_dimension(fid,dimid,temp,ELEMENTS)
    call handle_ncerror(ierr)

    ierr = nf90_inq_dimid(fid,"node",dimid)
    call handle_ncerror(ierr)
    ierr = nf90_inquire_dimension(fid,dimid,temp,NODES)
    call handle_ncerror(ierr)

    ierr = nf90_inq_dimid(fid,"siglay",dimid)
    call handle_ncerror(ierr)
    ierr = nf90_inquire_dimension(fid,dimid,temp,SIGLAY)
    call handle_ncerror(ierr)

    ierr = nf90_inq_dimid(fid,"siglev",dimid)
    call handle_ncerror(ierr)
    ierr = nf90_inquire_dimension(fid,dimid,temp,SIGLEV)
    call handle_ncerror(ierr)

    ierr = nf90_inq_dimid(fid,"time",dimid)
    call handle_ncerror(ierr)
    ierr = nf90_inquire_dimension(fid,dimid,temp,N_TIME)
    call handle_ncerror(ierr)

    ierr = nf90_inq_dimid(fid,"maxelem",dimid)
    call handle_ncerror(ierr)
    ierr = nf90_inquire_dimension(fid,dimid,temp,maxelem)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Allocate space for grid information                                         |
    !------------------------------------------------------------------------------|
    allocate(XP(NODES),YP(NODES))
    allocate(XC(ELEMENTS),YC(ELEMENTS))
    allocate(H(NODES))
    allocate(Z_LEV(SIGLEV),Z_LAY(SIGLAY))
    allocate(A1U(ELEMENTS,4),A2U(ELEMENTS,4))
    allocate(AW0(ELEMENTS,3))
    allocate(AWX(ELEMENTS,3),AWY(ELEMENTS,3))
    allocate(NV(ELEMENTS,3),NBE(ELEMENTS,3))
    allocate(NTVE(NODES),NBVE(NODES,maxelem))

    !------------------------------------------------------------------------------|
    !  Get node co-ordinates                                                       |
    !------------------------------------------------------------------------------|
    ierr = nf90_inq_varid(fid,"x",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,XP)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"y",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,YP)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"xc",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,XC)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"yc",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,YC)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"h",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,H)
    call handle_ncerror(ierr)

    allocate(tmp(NODES,SIGLEV))
    ierr = nf90_inq_varid(fid,"siglev",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,tmp)
    call handle_ncerror(ierr)
    Z_LEV = tmp(1,:)
    deallocate(tmp)

    allocate(tmp(NODES,SIGLAY))
    ierr = nf90_inq_varid(fid,"siglay",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,tmp)
    call handle_ncerror(ierr)
    Z_LAY = tmp(1,:)
    deallocate(tmp)

    ierr = nf90_inq_varid(fid,"nv",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,NV)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"nbe",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,NBE)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"ntve",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,NTVE)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"nbve",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,NBVE)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Get interpolation parameters                                                |
    !------------------------------------------------------------------------------|
    ierr = nf90_inq_varid(fid,"a1u",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,A1U)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"a2u",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,A2U)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"aw0",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,AW0)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"awx",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,AWX)
    call handle_ncerror(ierr)

    ierr = nf90_inq_varid(fid,"awy",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,AWY)
    call handle_ncerror(ierr)

    ! Close file
    ierr = nf90_close(fid)
    call handle_ncerror(ierr)
  end subroutine ncd_read_grid

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine read_flow(u,v,w,el,time)
    !==============================================================================|
    !  Read flow field vectors from NetCDF file                                    |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    real, dimension(ELEMENTS,SIGLAY), intent(out) :: u, v, w
    real, dimension(NODES),           intent(out) :: el
    integer,                          intent(in)  :: time
    !------------------------------------------------------------------------------|
    integer :: ierr
    integer :: fid, varid
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    ierr = nf90_open(trim(GRIDFN),NF90_NOWRITE,fid)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Read flow field data from NetCDF file at specified time index               |
    !------------------------------------------------------------------------------|

    ! Free surface elevation
    ierr = nf90_inq_varid(fid,"zeta",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,el,[1,time],[NODES,1])
    call handle_ncerror(ierr)

    ! Eastward warter velocity
    ierr = nf90_inq_varid(fid,"u",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,u,[1,1,time],[ELEMENTS,SIGLAY,1])
    call handle_ncerror(ierr)

    ! Northward water velocity
    ierr = nf90_inq_varid(fid,"v",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,v,[1,1,time],[ELEMENTS,SIGLAY,1])
    call handle_ncerror(ierr)

    ! Upward water velocity
    ierr = nf90_inq_varid(fid,"ww",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,w,[1,1,time],[ELEMENTS,SIGLAY,1])
    call handle_ncerror(ierr)

    ! Close file
    ierr = nf90_close(fid)
    call handle_ncerror(ierr)
  end subroutine read_flow

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine handle_ncerror(errid)
    !==============================================================================|
    !  Check for and handle NetCDF file I/O errors                                 |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    integer, intent(in) :: errid
    !==============================================================================|

    if (errid /= NF90_NOERR) then
      write(*,*) "ERROR: Accessing NetCDF file: ",trim(GRIDFN)
      write(*,*) trim(nf90_strerror(errid))
      stop
    end if
  end subroutine handle_ncerror

end module mod_flow_field

