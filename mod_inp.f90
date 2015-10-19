module mod_inp
  use mod_prec
contains

  subroutine get_val(lnum, numchar, text_line, varname, vartype, logval, stringval, realval, intval, nval)
    !==============================================================================|
    !  Decompose input line into variable name and variable value(s)               |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    integer,                intent(in)    :: lnum, numchar
    character(len=numchar)                :: text_line
    character(len=20),      intent(out)   :: varname
    character(len=7),       intent(out)   :: vartype
    logical,                intent(out)   :: logval
    character(len=80),      intent(out)   :: stringval(150)
    real,                   intent(inout) :: realval(150)
    integer,                intent(inout) :: intval(150)
    integer,                intent(out)   :: nval
    !------------------------------------------------------------------------------|
    character(len=numchar) :: varval, temp, frag(200)
    character(len=80)      :: tstring
    character(len=6)       :: errstring
    character(len=16)      :: numchars 
    integer                :: length, eqloc, lvarval, dotloc
    integer                :: i, locex, np
    logical                :: onfrag
    !==============================================================================|
    frag     = " "
    numchars = "0123456789+-Ee. " 
    vartype  = "error"
    logval   = .false.
    length   = len_trim(text_line) 
    write(errstring, "(i6)") lnum
    locex = index(text_line, "!")

    !------------------------------------------------------------------------------|
    !  check for blank line of comment                                             |
    !------------------------------------------------------------------------------|
    if (length == 0 .or. locex == 1) then
      vartype = "no data"
      varname = "no data"
      return
    end if

    !------------------------------------------------------------------------------|
    !  Change commas to blanks                                                     |
    !------------------------------------------------------------------------------|
    do i = 1,length
      if (text_line(i:i) == ",") text_line(i:i) = " "
    end do

    !------------------------------------------------------------------------------|
    !  Removing trailing comments                                                  |
    !------------------------------------------------------------------------------|
    if (locex /= 0) then
      temp = text_line(1:locex - 1)
      text_line = temp
    end if

    !------------------------------------------------------------------------------|
    !  Ensure "=" exists and determine location                                    |
    !------------------------------------------------------------------------------|
    eqloc = index(text_line, "=")

    !------------------------------------------------------------------------------|
    !  Split off varname and varval strings                                        |
    !------------------------------------------------------------------------------|
    varname = text_line(1:eqloc - 1)
    varval  = adjustl(text_line(eqloc + 1:length))
    lvarval = len_trim(varval)

    !------------------------------------------------------------------------------|
    !  Determine type of varval                                                    |
    !------------------------------------------------------------------------------|

    ! check for logical
    if ((varval(1:1) == "T" .or. varval(1:1) == "F") .and. lvarval == 1) then
      vartype = "logical"
      if (varval(1:1) == "T") logval = .true.
      return
    end if

    ! check if it is a string (contains characters other than 0-9,+,-,e,E,.)
    do i = 1,lvarval
      if (index(numchars,varval(i:i)) == 0) vartype = "string" 
    end do

    ! Process string (may be multiple)
    if (vartype == "string") then
      tstring      = varval
      stringval(1) = tstring 
      nval         = 1
      onfrag       = .true.
      do i = 1,lvarval
        if (varval(i:i) /= " ") then
          frag(nval) = trim(frag(nval))//varval(i:i)
          onfrag = .true.
        else
          if (onfrag) nval = nval + 1
          onfrag = .false.
        end if
      end do
      do i = 1,nval
        stringval(i + 1) = trim(frag(i))
      end do
      return
    end if

    ! check if it is a float
    dotloc = index(varval, ".")
    if (dotloc /= 0) then
      vartype = "float"
    else
      vartype = "integer"
    end if

    !------------------------------------------------------------------------------|
    !  Fragment into strings for multiple values                                   |
    !------------------------------------------------------------------------------|
    np     = 1
    onfrag = .true.
    do i = 1,lvarval
      if(varval(i:i) /= " ")then 
        frag(np) = trim(frag(np))//varval(i:i)
        onfrag   = .true.
      else
        if(onfrag) np = np + 1
        onfrag = .false.
      end if
    end do

    !------------------------------------------------------------------------------|
    !  Extract number(s) from character strings                                    |
    !------------------------------------------------------------------------------|
    nval = np
    do i = 1,np
      temp = trim(frag(i))
      if (vartype == "float") then 
        read(temp,*) realval(i)
      else
        read(temp,*) intval(i)
      end if
    end do
  end subroutine get_val

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  function scan_file(fname, vname, iscal, fscal, ivec, fvec, cvec, nsze, cval, lval)
    !==============================================================================|
    !  Scan an input file for a variable                                           |
    !  Return value:                                                               |
    !       0 = File found, variable value found                                   |
    !      -1 = File does not exist or permissions are incorrect                   |
    !      -2 = Variable not found or improperly set                               |
    !      -3 = Variable is of different type, check input file                    |
    !      -4 = Vector provided bur data is scalar type                            |
    !      -5 = No datatype desired, exiting                                       |
    !                                                                              |
    !  Required input:                                                             |
    !      fname = file Name                                                       |
    !                                                                              |
    !  Optional (must provide one)                                                 |
    !      iscal = integer scalar                                                  |
    !      fscal = Float scalar                                                    |
    !      cval  = Character variable                                              |
    !      lval  = Logical variable                                                |
    !      ivec  = integer vector **                                               |
    !      fvec  = Float vector **                                                 |
    !      cvec  = String vector ** (strings of length 80)                         |
    !      nsze  = Array size ** (must be provided with ivec/fvec)                 |
    !==============================================================================|
    implicit none
    !------------------------------------------------------------------------------|
    character(len=*),  intent(in)              :: fname, vname
    integer,           intent(inout), optional :: iscal, ivec(*)
    real,              intent(inout), optional :: fscal, fvec(*)
    character(len=80), intent(inout), optional :: cval, cvec(*)
    logical,           intent(inout), optional :: lval
    integer,           intent(out),   optional :: nsze
    !------------------------------------------------------------------------------|
    integer             :: scan_file
    real                :: realval(150)
    integer             :: intval(150)
    character(len=20)   :: varname
    character(len=80)   :: stringval(150)
    character(len=100)  :: inpline
    character(len=2400) :: tline
    character(len=7)    :: vartype
    integer             :: i, nval, nset, nline, nrep
    logical             :: check, logval
    !==============================================================================|

    !------------------------------------------------------------------------------|
    !  Open the input file                                                         |
    !------------------------------------------------------------------------------|
    scan_file = 0
    inquire(file = trim(fname), exist = check)
    if (.not.check) then
      scan_file = -1
      return
    end if

    open(10, file = trim(fname))
    rewind(10)

    !------------------------------------------------------------------------------|
    !  Scan the file for the variable name                                         |
    !------------------------------------------------------------------------------|
    nset  = 0
    nline = 0
    do while(.true.)
      tline(1:len(tline)) = ' '
      nrep  = 0
      nline = nline + 1
      read(10, '(a)', end = 20) inpline
      tline(1:80) = inpline(1:80)

      !------------------------------------------------------------------------------|
      !  Process line continuations                                                  |
      !------------------------------------------------------------------------------|
110   continue
      i = len_trim(inpline)

      if (i /= 0) then
        if (inpline(i - 1:i) == '\\') then
          nrep = nrep + 1
          read(10, '(a)', end = 20) inpline
          nline = nline + 1
          tline(nrep*80 + 1:nrep*80 + 80) = inpline(1:80)
          goto 110
        end if
      end if

      !------------------------------------------------------------------------------|
      !  Remove line continuation character                                          |
      !------------------------------------------------------------------------------|
      if (nrep > 0) then
        do i = 2,len_trim(tline)
          if (tline(i - 1:i) == '\\') tline(i - 1:i) = '  '
        end do
      end if

      !------------------------------------------------------------------------------|
      !  Process the line                                                            |
      !------------------------------------------------------------------------------|
      call get_val(nline, len_trim(tline), adjustl(tline), varname, vartype, logval, stringval, realval, intval, nval)

      !------------------------------------------------------------------------------|
      !  If varname matches, process variable and error-check                        |
      !------------------------------------------------------------------------------|
      if (trim(to_upper(varname)) == trim(to_upper(vname))) then
        if (present(iscal)) then
          if (vartype == 'integer') then
            iscal = intval(1)
            close(10)
            return
          else
            scan_file = -3
          end if
        else if (present(fscal)) then
          if (vartype == 'float') then
            fscal = realval(1)
            close(10)
            return
          else
            scan_file = -3
          end if
        else if (present(cval)) then
          if (vartype == 'string') then
            cval = stringval(1)
            close(10)
            return
          else
            scan_file = -3
          end if
        else if (present(lval)) then
          if (vartype == 'logical') then
            lval = logval
            close(10)
            return
          else
            scan_file = -3
          end if
        else if (present(ivec)) then
          if (nval > 1) then
            if (vartype == 'integer') then
              ivec(1:nval) = intval(1:nval)
              nsze = nval
              close(10)
              return
            else
              scan_file = -3
            end if
          else
            scan_file = -4 
          end if
        else if (present(fvec)) then
          if (nval > 1) then
            if (vartype == 'float') then
              fvec(1:nval) = realval(1:nval)
              nsze = nval
              close(10)
              return
            else
              scan_file = -3
            end if
          else
            scan_file = -4 
          end if
        else if (present(cvec)) then
          if (nval > 0) then
            if (vartype == 'string') then
              cvec(1:nval) = stringval(2:nval+1)
              nsze = nval
              close(10)
              return
            else
              scan_file = -3
            end if
          else
            scan_file = -4
          end if
        else
          scan_file = -5
        end if
      end if

    end do
20  close(10) 
    scan_file = -2
    return 
  end function scan_file

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  function to_upper(strin) result(strout)
    implicit none
    character(len=*), intent(in) :: strin
    character(len=len(strin))    :: strout
    integer                      :: i, j

    do i = 1, len(strin)
      j = iachar(strin(i:i))
      if (j >= iachar("a") .and. j <= iachar("z")) then
        strout(i:i) = achar(iachar(strin(i:i)) - 32)
      else
        strout(i:i) = strin(i:i)
      end if
    end do
  end function to_upper

  !==============================================================================|
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
  !==============================================================================|

  subroutine pscanmsg(msgn)
    implicit none
    integer, intent(in) :: msgn
    character(len=50)   :: msg(5)

    data msg(1) /'File does not exist or permissions are incorrect'/
    data msg(2) /'Variable not found or improperly set'/
    data msg(3) /'Variable is of different type, check input file'/
    data msg(4) /'Vector provided bur data is scalar type'/
    data msg(5) /'No datatype desired, exiting'/

    print *, msg(-msgn)
  end subroutine pscanmsg

end module mod_inp

