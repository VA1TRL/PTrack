module mod_io
  use netcdf

  interface nc_1d_var
    module procedure nc_1d_int_var, nc_1d_real_var
  end interface

contains

  subroutine nc_open_file(filepath, fileid)
    !==============================================================================|
    !  Prepare to read from the provided NetCDF file                               |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*), intent(in)  :: filepath
    integer,          intent(out) :: fileid
    !------------------------------------------------------------------------------|
    logical :: fexist
    integer :: ierr
    !==============================================================================|

    inquire(file=filepath, exist=fexist)
    if(.not.fexist) then
      write(*,*) "ERROR: NetCDF file ", filepath, " does not exist!"
      stop
    end if

    ierr = nf90_open(filepath, NF90_NOWRITE, fileid)
    call handle_ncerror(ierr)
  end subroutine nc_open_file

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine nc_close_file(fileid)
    !==============================================================================|
    !  Prepare to read from the provided NetCDF file                               |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    integer, intent(in) :: fileid
    !------------------------------------------------------------------------------|
    integer :: ierr
    !==============================================================================|

    ierr = nf90_close(fileid)
    call handle_ncerror(ierr)
  end subroutine nc_close_file

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine nc_dim(dimname, fileid, value)
    !==============================================================================|
    !  Read the value of the named dimmension from the provided NetCDF file        |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*), intent(in)  :: dimname
    integer,          intent(in)  :: fileid
    integer,          intent(out) :: value
    !------------------------------------------------------------------------------|
    integer                      :: ierr
    integer                      :: dimid
    character(len=NF90_MAX_NAME) :: temp
    !==============================================================================|

    ierr = nf90_inq_dimid(fileid, dimname, dimid)
    call handle_ncerror(ierr)
    ierr = nf90_inquire_dimension(fileid, dimid, temp, value)
    call handle_ncerror(ierr)
  end subroutine nc_dim

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine nc_1d_int_var(varname, fileid, nval, values)
    !==============================================================================|
    !  Read the value of the named variable from the provided NetCDF file          |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*),         intent(in)  :: varname
    integer,                  intent(in)  :: fileid
    integer,                  intent(in)  :: nval
    integer, dimension(nval), intent(out) :: values
    !------------------------------------------------------------------------------|
    integer                      :: ierr
    integer                      :: varid
    character(len=NF90_MAX_NAME) :: temp
    !==============================================================================|

    ierr = nf90_inq_varid(fileid, varname, varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fileid, varid, values, [1], [nval])
    call handle_ncerror(ierr)
  end subroutine nc_1d_int_var

  subroutine nc_1d_real_var(varname, fileid, nval, values)
    !==============================================================================|
    !  Read the value of the named variable from the provided NetCDF file          |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*),         intent(in)  :: varname
    integer,                  intent(in)  :: fileid
    integer,                  intent(in)  :: nval
    real,    dimension(nval), intent(out) :: values
    !------------------------------------------------------------------------------|
    integer                      :: ierr
    integer                      :: varid
    character(len=NF90_MAX_NAME) :: temp
    !==============================================================================|

    ierr = nf90_inq_varid(fileid, varname, varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fileid, varid, values, [1], [nval])
    call handle_ncerror(ierr)
  end subroutine nc_1d_real_var

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine nc_2d_var(varname, fileid, nval1, nval2, values)
    !==============================================================================|
    !  Read the value of the named variable from the provided NetCDF file          |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*),                intent(in)  :: varname
    integer,                         intent(in)  :: fileid
    integer,                         intent(in)  :: nval1, nval2
    integer, dimension(nval1,nval2), intent(out) :: values
    !------------------------------------------------------------------------------|
    integer                      :: ierr
    integer                      :: varid
    character(len=NF90_MAX_NAME) :: temp
    !==============================================================================|

    ierr = nf90_inq_varid(fileid, varname, varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fileid, varid, values, [1,1], [nval1,nval2])
    call handle_ncerror(ierr)
  end subroutine nc_2d_var

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine nc_3d_var(varname, fileid, nval1, nval2, nval3, values)
    !==============================================================================|
    !  Read the value of the named variable from the provided NetCDF file          |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*),                      intent(in)  :: varname
    integer,                               intent(in)  :: fileid
    integer,                               intent(in)  :: nval1, nval2, nval3
    integer, dimension(nval1,nval2,nval3), intent(out) :: values
    !------------------------------------------------------------------------------|
    integer                      :: ierr
    integer                      :: varid
    character(len=NF90_MAX_NAME) :: temp
    !==============================================================================|

    ierr = nf90_inq_varid(fileid, varname, varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fileid, varid, values, [1,1,1], [nval1,nval2,nval3])
    call handle_ncerror(ierr)
  end subroutine nc_3d_var

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
      write(*,*) "ERROR: Accessing NetCDF file: ", trim(nf90_strerror(errid))
      stop
    end if
  end subroutine handle_ncerror

end module mod_io

