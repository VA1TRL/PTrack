module mod_io
  use netcdf

  interface nc_1d_var
    module procedure nc_1d_int_var, nc_1d_real_var
  end interface

contains

  subroutine nc_open_file(filepath, writemode, fileid)
    !==============================================================================|
    !  Prepare to read from the provided NetCDF file                               |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*), intent(in)  :: filepath
    logical,          intent(in)  :: writemode
    integer,          intent(out) :: fileid
    !------------------------------------------------------------------------------|
    logical :: fexist
    integer :: ierr
    integer :: mode
    !==============================================================================|

    inquire(file=filepath, exist=fexist)
    if (.not.fexist) then
      write(*,*) "ERROR: NetCDF file ", filepath, " does not exist!"
      stop
    end if

    if (writemode) then
        mode = NF90_WRITE
    else
        mode = NF90_NOWRITE
    end if

    ierr = nf90_open(filepath, mode, fileid)
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

  subroutine nc_new_outfile(filepath, npoint, points)
    !==============================================================================|
    !  Create a new NetCDF file to recieve the simulation output                   |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*),           intent(in) :: filepath
    integer,                    intent(in) :: npoint
    integer, dimension(npoint), intent(in) :: points
    !------------------------------------------------------------------------------|
    integer :: ierr
    integer :: fid
    integer :: timeid, numid, numvarid, id
    !==============================================================================|

    ierr = nf90_create(filepath, NF90_CLOBBER, fid)
    call handle_ncerror(ierr)

    ierr = nf90_def_dim(fid, "time", NF90_UNLIMITED, timeid)
    call handle_ncerror(ierr)
    ierr = nf90_def_dim(fid, "number", npoint, numid)
    call handle_ncerror(ierr)

    ierr = nf90_def_var(fid, "time", NF90_FLOAT, [timeid], id)
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, id, "long_name", "Simulation Time")
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, id, "units", "days since 1858-11-17 00:00:00 +0:00")
    call handle_ncerror(ierr)

    ierr = nf90_def_var(fid, "number", NF90_INT, [numid], numvarid)
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, numvarid, "long_name", "Particle Identifier")
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, numvarid, "standard_name", "id")
    call handle_ncerror(ierr)

    ierr = nf90_def_var(fid, "x", NF90_FLOAT, [numid,timeid], id)
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, id, "long_name", "Domain X-Coordinate")
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, id, "units", "meters")
    call handle_ncerror(ierr)

    ierr = nf90_def_var(fid, "y", NF90_FLOAT, [numid,timeid], id)
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, id, "long_name", "Domain Y-Coordinate")
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, id, "units", "meters")
    call handle_ncerror(ierr)

    ierr = nf90_def_var(fid, "z", NF90_FLOAT, [numid,timeid], id)
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, id, "long_name", "Domain Z-Coordinate")
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, id, "standard_name", "depth")
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, id, "units", "meters")
    call handle_ncerror(ierr)
    ierr = nf90_put_att(fid, id, "positive", "down")
    call handle_ncerror(ierr)

    ierr = nf90_enddef(fid)
    call handle_ncerror(ierr)

    ierr = nf90_put_var(fid, numvarid, points, [1], [npoint])
    call handle_ncerror(ierr)

    ierr = nf90_close(fid)
    call handle_ncerror(ierr)
  end subroutine nc_new_outfile

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
    !==============================================================================|

    ierr = nf90_inq_varid(fileid, varname, varid)
    call handle_ncerror(ierr)
    ierr = nf90_get_var(fileid, varid, values, [1,1,1], [nval1,nval2,nval3])
    call handle_ncerror(ierr)
  end subroutine nc_3d_var

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine nc_1d_write(varname, fileid, location, value)
    !==============================================================================|
    !  Read the value of the named variable from the provided NetCDF file          |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*), intent(in) :: varname
    integer,          intent(in) :: fileid
    integer,          intent(in) :: location
    real,             intent(in) :: value
    !------------------------------------------------------------------------------|
    integer :: ierr
    integer :: varid
    !==============================================================================|

    ierr = nf90_inq_varid(fileid, varname, varid)
    call handle_ncerror(ierr)
    ierr = nf90_put_var(fileid, varid, [value], [location], [1])
    call handle_ncerror(ierr)
  end subroutine nc_1d_write

  subroutine nc_2d_write(varname, fileid, location, length, values)
    !==============================================================================|
    !  Read the value of the named variable from the provided NetCDF file          |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*),        intent(in) :: varname
    integer,                 intent(in) :: fileid
    integer,                 intent(in) :: location
    integer,                 intent(in) :: length
    real, dimension(length), intent(in) :: values
    !------------------------------------------------------------------------------|
    integer :: ierr
    integer :: varid
    !==============================================================================|

    ierr = nf90_inq_varid(fileid, varname, varid)
    call handle_ncerror(ierr)
    ierr = nf90_put_var(fileid, varid, values, [1,location], [length,1])
    call handle_ncerror(ierr)
  end subroutine nc_2d_write

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

