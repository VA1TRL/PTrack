module mod_flow_field
  !==============================================================================|
  !  Grid/Flow-Field data storage and access functions                           |
  !==============================================================================|
  use mod_prec
  use mod_config
  use mod_io
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
    integer :: maxelem
    integer :: fileid
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    call nc_open_file(GRIDFN, .false., fileid)

    !------------------------------------------------------------------------------|
    !  Get model dimensions                                                        |
    !------------------------------------------------------------------------------|
    call nc_dim(fileid, "nele",    ELEMENTS)
    call nc_dim(fileid, "node",    NODES)
    call nc_dim(fileid, "siglay",  SIGLAY)
    call nc_dim(fileid, "siglev",  SIGLEV)
    call nc_dim(fileid, "maxelem", maxelem)

    !------------------------------------------------------------------------------|
    !  Close file                                                                  |
    !------------------------------------------------------------------------------|
    call nc_close_file(fileid)

    !------------------------------------------------------------------------------|
    !  Allocate space for grid information                                         |
    !------------------------------------------------------------------------------|
    allocate(XP(NODES), YP(NODES))
    allocate(XC(ELEMENTS), YC(ELEMENTS))
    allocate(H(NODES))
    allocate(Z_LEV(SIGLEV), Z_LAY(SIGLAY))
    allocate(A1U(ELEMENTS,4), A2U(ELEMENTS,4))
    allocate(AW0(ELEMENTS,3))
    allocate(AWX(ELEMENTS,3), AWY(ELEMENTS,3))
    allocate(NV(ELEMENTS,3), NBE(ELEMENTS,3))
    allocate(NTVE(NODES), NBVE(NODES,maxelem))
    allocate(U_START(ELEMENTS,SIGLAY), U_END(ELEMENTS,SIGLAY))
    allocate(V_START(ELEMENTS,SIGLAY), V_END(ELEMENTS,SIGLAY))
    allocate(W_START(ELEMENTS,SIGLAY), W_END(ELEMENTS,SIGLAY))
    allocate(EL_START(NODES), EL_END(NODES))

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
    integer                           :: fileid
    real, allocatable, dimension(:,:) :: tmp
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    call nc_open_file(GRIDFN, .false., fileid)

    !------------------------------------------------------------------------------|
    !  Get node co-ordinates                                                       |
    !------------------------------------------------------------------------------|
    call nc_read_var(fileid, "x",    XP)
    call nc_read_var(fileid, "y",    YP)
    call nc_read_var(fileid, "xc",   XC)
    call nc_read_var(fileid, "yc",   YC)
    call nc_read_var(fileid, "h",    H)
    call nc_read_var(fileid, "nv",   NV)
    call nc_read_var(fileid, "nbe",  NBE)
    call nc_read_var(fileid, "ntve", NTVE)
    call nc_read_var(fileid, "nbve", NBVE)

    allocate(tmp(NODES,SIGLEV))
    call nc_read_var(fileid, "siglev", tmp)
    Z_LEV = tmp(1,:)
    deallocate(tmp)
    allocate(tmp(NODES,SIGLAY))
    call nc_read_var(fileid, "siglay", tmp)
    Z_LAY = tmp(1,:)
    deallocate(tmp)

    !------------------------------------------------------------------------------|
    !  Get interpolation parameters                                                |
    !------------------------------------------------------------------------------|
    call nc_read_var(fileid, "a1u", A1U)
    call nc_read_var(fileid, "a2u", A2U)
    call nc_read_var(fileid, "aw0", AW0)
    call nc_read_var(fileid, "awx", AWX)
    call nc_read_var(fileid, "awy", AWY)

    !------------------------------------------------------------------------------|
    !  Close file                                                                  |
    !------------------------------------------------------------------------------|
    call nc_close_file(fileid)
  end subroutine ncd_read_grid

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine get_flow_at(u, v, w, el, time)
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
      call read_flow(U_START, V_START, W_START, EL_START, T_START, NC_RECORD)
      NC_RECORD = NC_RECORD + 1
      call read_flow(U_END, V_END, W_END, EL_END, T_END, NC_RECORD)
    end if

    !------------------------------------------------------------------------------|
    !  Update forcing records if nessesairy                                        |
    !------------------------------------------------------------------------------|
    if (time > T_END) then
      write(*,*) "Done Processing: ", NC_RECORD, ", simulation time (mjd): ", time

      U_START   = U_END
      V_START   = V_END
      W_START   = W_END
      EL_START  = EL_END
      T_START   = T_END
      NC_RECORD = NC_RECORD + 1

      call read_flow(U_END, V_END, W_END, EL_END, T_END, NC_RECORD)
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

  subroutine read_flow(u, v, w, el, time, n)
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
    integer :: fileid
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    call nc_open_file(GRIDFN, .false., fileid)

    !------------------------------------------------------------------------------|
    !  Read flow field data from NetCDF file at specified time index               |
    !------------------------------------------------------------------------------|
    call nc_1d_read(fileid, "time", n, time)
    call nc_2d_read(fileid, "zeta", n, NODES,    el)
    call nc_3d_read(fileid, "u",    n, ELEMENTS, SIGLAY, u)
    call nc_3d_read(fileid, "v",    n, ELEMENTS, SIGLAY, v)
    call nc_3d_read(fileid, "ww",   n, ELEMENTS, SIGLAY, w)

    !------------------------------------------------------------------------------|
    !  Close file                                                                  |
    !------------------------------------------------------------------------------|
    call nc_close_file(fileid)
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
    integer                         :: sizet
    integer                         :: fileid
    integer                         :: i
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    call nc_open_file(GRIDFN, .false., fileid)

    !------------------------------------------------------------------------------|
    !  Get the time records                                                        |
    !------------------------------------------------------------------------------|
    call nc_dim(fileid, "time", sizet)
    allocate(field_times(sizet))
    call nc_read_var(fileid, "time", field_times)

    !------------------------------------------------------------------------------|
    !  Close file                                                                  |
    !------------------------------------------------------------------------------|
    call nc_close_file(fileid)

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

    write(*,*) "ERROR: Could not find time ", time, " in NetCDF flow-field file!"
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
    integer                         :: sizet
    integer                         :: fileid
    integer                         :: i
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open NetCDF data file                                                       |
    !------------------------------------------------------------------------------|
    call nc_open_file(GRIDFN, .false., fileid)

    !------------------------------------------------------------------------------|
    !  Get the time records                                                        |
    !------------------------------------------------------------------------------|
    call nc_dim(fileid, "time", sizet)
    allocate(field_times(sizet))
    call nc_read_var(fileid, "time", field_times)

    !------------------------------------------------------------------------------|
    !  Close file                                                                  |
    !------------------------------------------------------------------------------|
    call nc_close_file(fileid)

    !------------------------------------------------------------------------------|
    !  Perform the check                                                           |
    !------------------------------------------------------------------------------|
    if (field_times(1) > time_start .or. field_times(sizet) < time_end) then
      isVallid = .false.
      return 
    end if
    isVallid = .true.
    return
  end function check_field_time

end module mod_flow_field

