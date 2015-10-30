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

  !------------------------------------------------------------------------------|
  !  NetCDF flow-field records covering the current simulation time              |
  !------------------------------------------------------------------------------|
  real, allocatable, dimension(:,:) :: U_START, U_END
  real, allocatable, dimension(:,:) :: V_START, V_END
  real, allocatable, dimension(:,:) :: W_START, W_END
  real, allocatable, dimension(:)   :: EL_START, EL_END
  real                              :: T_START, T_END
  integer                           :: NC_RECORD
contains

  subroutine init_flow_field
    !==============================================================================|
    !  Read grid shape/dimmension information from NetCDF file and allocate space  |
    !  for the data arrays.                                                        |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    integer                      :: ierr
    integer                      :: maxelem
    integer                      :: dimid, fid
    character(len=NF90_MAX_NAME) :: temp
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

    ierr = nf90_inq_dimid(fid,"maxelem",dimid)
    call handle_ncerror(ierr)
    ierr = nf90_inquire_dimension(fid,dimid,temp,maxelem)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Close file                                                                  |
    !------------------------------------------------------------------------------|
    ierr = nf90_close(fid)
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
    allocate(U_START(ELEMENTS,SIGLAY),U_END(ELEMENTS,SIGLAY))
    allocate(V_START(ELEMENTS,SIGLAY),V_END(ELEMENTS,SIGLAY))
    allocate(W_START(ELEMENTS,SIGLAY),W_END(ELEMENTS,SIGLAY))
    allocate(EL_START(NODES),EL_END(NODES))

    !------------------------------------------------------------------------------|
    !  Initialize grid information arrays                                          |
    !------------------------------------------------------------------------------|
    call ncd_read_grid
    NC_RECORD = -1
  end subroutine init_flow_field

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine ncd_read_grid
    !==============================================================================|
    !  Read the finite element mesh from the NetCDF flow-field file                |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    integer                           :: ierr
    integer                           :: dimid, varid, fid
    real, allocatable, dimension(:,:) :: tmp
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    ierr = nf90_open(trim(GRIDFN),NF90_NOWRITE,fid)
    call handle_ncerror(ierr)

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

    !------------------------------------------------------------------------------|
    !  Close file                                                                  |
    !------------------------------------------------------------------------------|
    ierr = nf90_close(fid)
    call handle_ncerror(ierr)
  end subroutine ncd_read_grid

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine get_flow_at(u,v,w,el,time)
    !==============================================================================|
    !  Get the flow field corrosponding to the provided time (MJD format) by       |
    !  interpolating between the time-adjacent NetCDF records.                     |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    real(DP), dimension(ELEMENTS,SIGLAY), intent(out) :: u, v, w
    real(DP), dimension(NODES),           intent(out) :: el
    real,                                 intent(in)  :: time
    !------------------------------------------------------------------------------|
    real :: frac_start, frac_end
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Initialize flow field data on the first run of this subroutine              |
    !------------------------------------------------------------------------------|
    if (NC_RECORD < 0) then
      NC_RECORD = timeIndex(time)
      call read_flow(U_START,V_START,W_START,EL_START,T_START,NC_RECORD)
      NC_RECORD = NC_RECORD + 1
      call read_flow(U_END,V_END,W_END,EL_END,T_END,NC_RECORD)
    end if

    !------------------------------------------------------------------------------|
    !  Update forcing records if nessesairy                                        |
    !------------------------------------------------------------------------------|
    if (time > T_END) then
      write(*,*) "Done Processing: ",NC_RECORD,", simulation time (mjd): ",time

      U_START   = U_END
      V_START   = V_END
      W_START   = W_END
      EL_START  = EL_END
      T_START   = T_END
      NC_RECORD = NC_RECORD + 1

      call read_flow(U_END,V_END,W_END,EL_END,T_END,NC_RECORD)
    end if

    !------------------------------------------------------------------------------|
    !  Interpolate velocity field values to the current simulation time            |
    !------------------------------------------------------------------------------|
    frac_start = (T_END - time)/(T_END - T_START)
    frac_end   = (time - T_START)/(T_END - T_START)

    u  = frac_start*U_START  + frac_end*U_END
    v  = frac_start*V_START  + frac_end*V_END
    w  = frac_start*W_START  + frac_end*W_END
    el = frac_start*EL_START + frac_end*EL_END
  end subroutine get_flow_at

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine read_flow(u,v,w,el,time,n)
    !==============================================================================|
    !  Read flow field vectors from NetCDF file                                    |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    real, dimension(ELEMENTS,SIGLAY), intent(out) :: u, v, w
    real, dimension(NODES),           intent(out) :: el
    real,                             intent(out) :: time
    integer,                          intent(in)  :: n
    !------------------------------------------------------------------------------|
    integer            :: ierr
    integer            :: fid, varid
    real, dimension(1) :: tmp
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    ierr = nf90_open(trim(GRIDFN),NF90_NOWRITE,fid)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Read flow field data from NetCDF file at specified time index               |
    !------------------------------------------------------------------------------|

    ! Time step
    ierr = nf90_inq_varid(fid,"time",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,tmp,[n],[1])
    time = tmp(1)
    call handle_ncerror(ierr)

    ! Free surface elevation
    ierr = nf90_inq_varid(fid,"zeta",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,el,[1,n],[NODES,1])
    call handle_ncerror(ierr)

    ! Eastward warter velocity
    ierr = nf90_inq_varid(fid,"u",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,u,[1,1,n],[ELEMENTS,SIGLAY,1])
    call handle_ncerror(ierr)

    ! Northward water velocity
    ierr = nf90_inq_varid(fid,"v",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,v,[1,1,n],[ELEMENTS,SIGLAY,1])
    call handle_ncerror(ierr)

    ! Upward water velocity
    ierr = nf90_inq_varid(fid,"ww",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,w,[1,1,n],[ELEMENTS,SIGLAY,1])
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Close file                                                                  |
    !------------------------------------------------------------------------------|
    ierr = nf90_close(fid)
    call handle_ncerror(ierr)
  end subroutine read_flow

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  function timeIndex(time) result(n)
    !==============================================================================|
    !  Get the NetCDF record index corrosponding to the provided time              |
    !                                                                              |
    !  Note: Indecies start at 1, and the provided record will be the NetCDF       |
    !  record at or immidiatly before the provided time.                           |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    real, intent(in) :: time
    integer          :: n
    !------------------------------------------------------------------------------|
    real, allocatable, dimension(:) :: field_times
    character(len=NF90_MAX_NAME)    :: temp
    integer                         :: sizet
    integer                         :: ierr
    integer                         :: fid, dimid, varid
    integer                         :: i
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    ierr = nf90_open(trim(GRIDFN),NF90_NOWRITE,fid)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Get the time records                                                        |
    !------------------------------------------------------------------------------|
    ierr = nf90_inq_dimid(fid,"time",dimid)
    call handle_ncerror(ierr)
    ierr = nf90_inquire_dimension(fid,dimid,temp,sizet)
    call handle_ncerror(ierr)
    allocate(field_times(sizet))

    ierr = nf90_inq_varid(fid,"time",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,field_times)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Close file                                                                  |
    !------------------------------------------------------------------------------|
    ierr = nf90_close(fid)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Find the NetCDF record that covers the provided time                        |
    !------------------------------------------------------------------------------|
    if (field_times(1) <= time) then
      do i = 2,sizet
        if (field_times(i) > time) then
          n = i - 1
          return
        end if
      end do
    end if

    write(*,*) "ERROR: Could not find time ",time," in NetCDF flow-field file!"
    stop
  end function timeIndex

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  function check_field_time(time_start, time_end) result(isVallid)
    !==============================================================================|
    !  Check if the provided timeline is covered by the NetCDF dataset             |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    real, intent(in) :: time_start, time_end
    logical          :: isVallid
    !------------------------------------------------------------------------------|
    real, allocatable, dimension(:) :: field_times
    character(len=NF90_MAX_NAME)    :: temp
    integer                         :: sizet
    integer                         :: ierr
    integer                         :: fid, dimid, varid
    integer                         :: i
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    ierr = nf90_open(trim(GRIDFN),NF90_NOWRITE,fid)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Get the time records                                                        |
    !------------------------------------------------------------------------------|
    ierr = nf90_inq_dimid(fid,"time",dimid)
    call handle_ncerror(ierr)
    ierr = nf90_inquire_dimension(fid,dimid,temp,sizet)
    call handle_ncerror(ierr)
    allocate(field_times(sizet))

    ierr = nf90_inq_varid(fid,"time",varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fid,varid,field_times)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Close file                                                                  |
    !------------------------------------------------------------------------------|
    ierr = nf90_close(fid)
    call handle_ncerror(ierr)

    !------------------------------------------------------------------------------|
    !  Perform the check                                                           |
    !------------------------------------------------------------------------------|
    if (field_times(1) > time_start .or. field_times(sizet) <= time_end) then
      isVallid = .false.
      return 
    end if
    isVallid = .true.
    return
  end function check_field_time

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

