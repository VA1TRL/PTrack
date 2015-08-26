MODULE MOD_INP

CONTAINS

  !==============================================================================!
  !  DECOMPOSE INPUT LINE INTO VARIABLE NAME AND VARIABLE VALUE(S)               !
  !==============================================================================!

  SUBROUTINE GET_VAL(LNUM,NUMCHAR,TEXT_LINE,VARNAME,VARTYPE,LOGVAL,STRINGVAL,&
       REALVAL,INTVAL,NVAL)

    !==============================================================================!
    USE MOD_PREC
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LNUM,NUMCHAR
    CHARACTER(LEN=NUMCHAR) :: TEXT_LINE
    CHARACTER(LEN=20), INTENT(OUT) :: VARNAME
    CHARACTER(LEN=7), INTENT(OUT) :: VARTYPE
    LOGICAL, INTENT(OUT) :: LOGVAL
    CHARACTER(LEN=80), INTENT(OUT) :: STRINGVAL(150)
    REAL, INTENT(INOUT) :: REALVAL(150)
    INTEGER, INTENT(INOUT) :: INTVAL(150)
    INTEGER, INTENT(OUT) :: NVAL
    !------------------------------------------------------------------------------!
    CHARACTER(LEN=NUMCHAR) :: VARVAL,TEMP,FRAG(200)
    CHARACTER(LEN=80) :: TSTRING
    CHARACTER(LEN=6) :: ERRSTRING
    CHARACTER(LEN=16) :: NUMCHARS 
    INTEGER LENGTH,EQLOC,LVARVAL,DOTLOC
    INTEGER I,LOCEX,NP
    LOGICAL ONFRAG

    !==============================================================================!
    FRAG = " "
    NUMCHARS = "0123456789+-Ee. " 
    VARTYPE = "error"
    LOGVAL = .FALSE.
    LENGTH = LEN_TRIM(TEXT_LINE) 
    WRITE(ERRSTRING,"(I6)") LNUM
    LOCEX = INDEX(TEXT_LINE,"!")

    !
    !-----------------------CHECK FOR BLANK LINE OR COMMENT------------------------!
    !
    IF(LENGTH == 0 .OR. LOCEX==1)THEN
       VARTYPE = "no data"
       VARNAME = "no data"
       RETURN
    END IF

    !
    !-----------------------CHANGE COMMAS TO BLANKS--------------------------------!
    !
    DO I=1,LENGTH
       IF(TEXT_LINE(I:I) == ",") TEXT_LINE(I:I) = " "
    END DO
    !
    !-----------------------REMOVING TRAILING COMMENTS-----------------------------!
    !
    IF(LOCEX /= 0)THEN
       TEMP = TEXT_LINE(1:LOCEX-1)
       TEXT_LINE = TEMP
    END IF
    !
    !--------------------ENSURE "=" EXISTS AND DETERMINE LOCATION------------------!
    !
    EQLOC = INDEX(TEXT_LINE,"=")
    !   IF(EQLOC == 0) CALL PERROR(6,'DATA LINE '//ERRSTRING//' MUST CONTAIN "=" ')

    !
    !--------------------SPLIT OFF VARNAME AND VARVAL STRINGS----------------------!
    !
    VARNAME = TEXT_LINE(1:EQLOC-1)
!    print *, VARNAME
    VARVAL  = ADJUSTL(TEXT_LINE(EQLOC+1:LENGTH))
!    print *, '&',trim(VARVAL),'&'
    LVARVAL = LEN_TRIM(VARVAL)
!    print *, LVARVAL
    !   IF(LVARVAL == 0) CALL PERROR(6,'IN DATA PARAMETER FILE', &
    !                'VARIABLE LINE'//ERRSTRING//' HAS NO ASSOCIATED VALUE')
    !
    !-----------------DETERMINE TYPE OF VARVAL-------------------------------------!
    !

    !
    !  CHECK FOR LOGICAL
    !
    IF((VARVAL(1:1) == "T" .OR. VARVAL(1:1) == "F") .AND. LVARVAL == 1)THEN 
       VARTYPE = "logical"
       IF(VARVAL(1:1) == "T") LOGVAL = .TRUE.
       RETURN
    END IF

    !
    !  CHECK IF IT IS A STRING  (CONTAINS CHARACTERS OTHER THAN 0-9,+,-,e,E,.)
    !
    DO I=1,LVARVAL
       IF(INDEX(NUMCHARS,VARVAL(I:I)) == 0) VARTYPE = "string" 
    END DO
!    print *, VARTYPE
    !
    !  PROCESS STRING (MAY BE MULTIPLE)
    !
    IF(VARTYPE == "string") THEN
       TSTRING = VARVAL
       STRINGVAL(1) = TSTRING 
       NVAL = 1
       ONFRAG = .TRUE.
       DO I=1,LVARVAL
          IF(VARVAL(I:I) /= " ")THEN
             FRAG(NVAL) = TRIM(FRAG(NVAL))//VARVAL(I:I)
             ONFRAG = .TRUE.
          ELSE
             IF(ONFRAG) NVAL = NVAL + 1
             ONFRAG = .FALSE.
          END IF
       END DO
       DO I=1,NVAL
          STRINGVAL(I+1) = TRIM(FRAG(I))
       END DO
       RETURN
    END IF

    !
    !  CHECK IF IT IS A FLOAT
    !

    DOTLOC = INDEX(VARVAL,".")
    IF(DOTLOC /= 0) THEN
       VARTYPE = "float"
    ELSE
       VARTYPE = "integer"
    END IF
!    print *, VARTYPE
    !
    !-----------------FRAGMENT INTO STRINGS FOR MULTIPLE VALUES---------------------!
    !
    NP = 1
    ONFRAG = .TRUE.
    DO I=1,LVARVAL
       IF(VARVAL(I:I) /= " ")THEN 
          FRAG(NP) = TRIM(FRAG(NP))//VARVAL(I:I)
          ONFRAG = .TRUE.
       ELSE
          IF(ONFRAG) NP = NP + 1
          ONFRAG = .FALSE.
       END IF
    END DO
    !
    !-----------------EXTRACT NUMBER(S) FROM CHARACTER STRINGS----------------------!
    !
!    print *, ONFRAG, trim(FRAG(1))
    NVAL = NP
!    print *, NP, FRAG
    DO I=1,NP
       TEMP = TRIM(FRAG(I))
!	   print *,'$', trim(TEMP),'$'
       IF(VARTYPE == "float") THEN 
          READ(TEMP,*)REALVAL(I)
          !print *, REALVAL
       ELSE
          READ(TEMP,*)INTVAL(I)
!	      print *,intval(i)
       END IF
    END DO
!  print *, REALVAL
  END SUBROUTINE GET_VAL


  !==============================================================================|

  FUNCTION SCAN_FILE(FNAME,VNAME,ISCAL,FSCAL,IVEC,FVEC,CVEC,NSZE,CVAL,LVAL)           

    !==============================================================================|
    !   Scan an Input File for a Variable                                          |
    !   RETURN VALUE:                                                              |
    !        0 = FILE FOUND, VARIABLE VALUE FOUND                                  |
    !       -1 = FILE DOES NOT EXIST OR PERMISSIONS ARE INCORRECT                  |
    !       -2 = VARIABLE NOT FOUND OR IMPROPERLY SET                              |
    !       -3 = VARIABLE IS OF DIFFERENT TYPE, CHECK INPUT FILE                   |
    !       -4 = VECTOR PROVIDED BUT DATA IS SCALAR TYPE                           |
    !       -5 = NO DATATYPE DESIRED, EXITING                                      |
    !							                                                   |
    !   REQUIRED INPUT:		        				                               |
    !        FNAME = File Name					                                   |
    !        FSIZE = Length of Filename					                           |
    !                                                                              | 
    !   OPTIONAL (MUST PROVIDE ONE)        					                       | 
    !        ISCAL = INTEGER SCALAR					                               |
    !        FSCAL = FLOAT SCALAR  						                           | 
    !        CVAL = CHARACTER VARIABLE                                             |
    !        LVAL = LOGICAL VARIABLE                                               |
    !        IVEC = INTEGER VECTOR **                                              |
    !        FVEC = FLOAT VECTOR **                                                |
    !        CVEC = STRING VECTOR ** (STRINGS OF LENGTH 80)                        |
    !      **NSZE = ARRAY SIZE (MUST BE PROVIDED WITH IVEC/FVEC)                   |
    !                                                                              | 
    !==============================================================================|

    USE MOD_PREC
    IMPLICIT NONE
    CHARACTER(LEN=*) :: FNAME,VNAME
    INTEGER, INTENT(INOUT), OPTIONAL :: ISCAL,IVEC(*)
    REAL, INTENT(INOUT), OPTIONAL    :: FSCAL,FVEC(*)
    CHARACTER(LEN=80), OPTIONAL      :: CVAL,CVEC(*)
    LOGICAL, INTENT(INOUT), OPTIONAL :: LVAL
    INTEGER, INTENT(out), OPTIONAL :: NSZE 

    !------------------------------------------------------------------------------|

    INTEGER :: SCAN_FILE
    REAL REALVAL(150)
    INTEGER  INTVAL(150)
    CHARACTER(LEN=20 ) :: VARNAME
    CHARACTER(LEN=80 ) :: STRINGVAL(150)
    CHARACTER(LEN=100 ) :: INPLINE
    CHARACTER(LEN=2400) :: TLINE
    CHARACTER(LEN=7  ) :: VARTYPE
    INTEGER I,NVAL,NSET,NLINE,NREP
    LOGICAL CHECK,LOGVAL


    SCAN_FILE = 0
    !==============================================================================|
    !            OPEN THE INPUT FILE                                               |
    !==============================================================================|
    INQUIRE(FILE=TRIM(FNAME),EXIST=CHECK)
    IF(.NOT.CHECK)THEN
       SCAN_FILE = -1
       RETURN
    END IF

    OPEN(10,FILE=TRIM(FNAME)) ; REWIND(10) 

    !==============================================================================|
    !            SCAN THE FILE FOR THE VARIABLE NAME                               |
    !==============================================================================|

    NSET = 0
    NLINE = 0
    DO WHILE(.TRUE.)
       TLINE(1:LEN(TLINE)) = ' ' 
       NREP  = 0
       NLINE = NLINE + 1
       READ(10,'(a)',END=20) INPLINE
!       print *, INPLINE
       TLINE(1:80) = INPLINE(1:80)

       !----PROCESS LINE CONTINUATIONS------------------------------------------------!
110    CONTINUE
       I = LEN_TRIM(INPLINE)
     
       IF(I /= 0)THEN


             IF( INPLINE(I-1:I) == '\\')THEN

                NREP = NREP + 1
                READ(10,'(a)',END=20) INPLINE
!		print *, '%',trim(INPLINE),'%'
                NLINE = NLINE + 1
                TLINE( NREP*80 + 1 : NREP*80 +80) = INPLINE(1:80)
!		print *, '%',trim(tline),'%'
                 GOTO 110
             END IF
          END IF
          !     IF(NREP > 4)CALL PERROR(6,"CANNOT HAVE > 4 LINE CONTINUATIONS")

          !----REMOVE LINE CONTINUATION CHARACTER \\-------------------------------------!
          IF(NREP > 0)THEN
             DO I=2,LEN_TRIM(TLINE)



                IF( TLINE(I-1:I) == '\\') TLINE(I-1:I) = '  '

             END DO
          END IF
!          print *, LEN_TRIM(TLINE)
          !----PROCESS THE LINE----------------------------------------------------------!
          CALL GET_VAL(NLINE,LEN_TRIM(TLINE),ADJUSTL(TLINE),VARNAME,VARTYPE,LOGVAL,&
               STRINGVAL,REALVAL,INTVAL,NVAL)

!          print *, I
          !----IF VARNAME MATCHES, PROCESS VARIABLE AND ERROR-CHECK----------------------!

          IF(TRIM(to_upper(VARNAME)) == TRIM(to_upper(VNAME)))THEN

             IF(PRESENT(ISCAL))THEN
                IF(VARTYPE == 'integer')THEN
                   ISCAL = INTVAL(1)
                   CLOSE(10)
                   RETURN
                ELSE
                   SCAN_FILE = -3
                END IF
             ELSE IF(PRESENT(FSCAL))THEN
                IF(VARTYPE == 'float')THEN
                   FSCAL = REALVAL(1)
                   CLOSE(10)
                   RETURN
                ELSE
                   SCAN_FILE = -3
                END IF
             ELSE IF(PRESENT(CVAL))THEN
                IF(VARTYPE == 'string')THEN
                   CVAL = STRINGVAL(1)
                   CLOSE(10)
                   RETURN
                ELSE
                   SCAN_FILE = -3
                END IF
             ELSE IF(PRESENT(LVAL))THEN
                IF(VARTYPE == 'logical')THEN
                   LVAL = LOGVAL
                   CLOSE(10)
                   RETURN
                ELSE
                   SCAN_FILE = -3
                END IF
             ELSE IF(PRESENT(IVEC))THEN
                IF(NVAL > 1)THEN
                   IF(VARTYPE == 'integer')THEN
                      IVEC(1:NVAL) = INTVAL(1:NVAL)
                      NSZE = NVAL
                      CLOSE(10)
                      RETURN
                   ELSE
                      SCAN_FILE = -3
                   END IF
                ELSE
                   SCAN_FILE = -4 
                END IF
             ELSE IF(PRESENT(FVEC))THEN
                IF(NVAL > 1)THEN
                   IF(VARTYPE == 'float')THEN
                      FVEC(1:NVAL) = REALVAL(1:NVAL)
                      NSZE = NVAL
                      CLOSE(10)
                      RETURN
                   ELSE
                      SCAN_FILE = -3
                   END IF
                ELSE
                   SCAN_FILE = -4 
                END IF
             ELSE IF(PRESENT(CVEC))THEN
                IF(NVAL > 0)THEN
                   IF(VARTYPE == 'string')THEN
                      CVEC(1:NVAL) = STRINGVAL(2:NVAL+1)
                      NSZE = NVAL
                      CLOSE(10)
                      RETURN
                   ELSE
                      SCAN_FILE = -3
                   END IF
                ELSE
                   SCAN_FILE = -4
                END IF
             ELSE
                SCAN_FILE = -5
             END IF
          END IF  !!VARIABLE IS CORRECT

       END DO !!LOOP OVER INPUT FILE
20     CLOSE(10) 
       SCAN_FILE = -2
       RETURN 
     END FUNCTION SCAN_FILE



function to_upper(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)

     implicit none

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

end function to_upper

function PScanMsg(msgn)
	!		print scan error message to stdout
    !       -1 = FILE DOES NOT EXIST OR PERMISSIONS ARE INCORRECT                  |
    !       -2 = VARIABLE NOT FOUND OR IMPROPERLY SET                              |
    !       -3 = VARIABLE IS OF DIFFERENT TYPE, CHECK INPUT FILE                   |
    !       -4 = VECTOR PROVIDED BUT DATA IS SCALAR TYPE                           |
    !       -5 = NO DATATYPE DESIRED, EXITING                                      |

	implicit none
	integer PScanMsg
	integer, intent(in) :: msgn
	character*50 msg(5)
	data msg(1)/'FILE DOES NOT EXIST OR PERMISSIONS ARE INCORRECT'/
	data msg(2)/'VARIABLE NOT FOUND OR IMPROPERLY SET'/
	data msg(3)/'VARIABLE IS OF DIFFERENT TYPE, CHECK INPUT FILE'/
	data msg(4)/'VECTOR PROVIDED BUT DATA IS SCALAR TYPE'/
	data msg(5)/'NO DATATYPE DESIRED, EXITING'/
	print *,msg(-msgn)
	PScanMsg=1
end function PScanMsg

     !==============================================================================|
   END MODULE MOD_INP
