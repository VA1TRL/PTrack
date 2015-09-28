subroutine fhe(x,y,host,indomain)
  !==============================================================================|
  !  Determine Which Element A Particle Resides in By Searching                  |
  !  Neighboring Elements.  Updates host component of Lagrangian Particle        |
  !  Type and update "found" flag to indicate whether the host                   |
  !  Has been found	                                                             |
  !==============================================================================|
  use mod_flow_field
  implicit none
  !------------------------------------------------------------------------------|
  real(DP), intent(in)    :: x, y
  integer,  intent(inout) :: host
  logical,  intent(out)   :: indomain
  logical                 :: isintriangle
  !------------------------------------------------------------------------------|
  integer                 :: i, j, iney, ncheck
  !==============================================================================|

  indomain = .true.

  ! Check if the particle has moved to a diffrent mesh element
  if (isintriangle(host,x,y)) then
    return
  end if

  ! Check neighbors
  do i = 1,3
    ncheck = NV(host,i)
    do j = 1,NTVE(ncheck)
      iney = NBVE(ncheck,j)
      if (isintriangle(iney,x,y)) then
        host = iney
        return
      end if
    end do
  end do

  ! Try harder
  call fhe_robust(x,y,host,indomain)
end subroutine fhe

!==============================================================================|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!==============================================================================|

subroutine fhe_robust(x,y,host,indomain)
  !==============================================================================|
  !  Find Home Element For Points (X,Y)                                          |
  !  Search Nearest Element to Progressively Further Elements. Updates Lagrangian|
  !  component "host" and marks Lagrangian component "ifound" with 1 if          |
  !  found.  returns logical variable "all_found" if all lagrangian variables    |
  !  have a known host element.  The host element may have been found prior to   |
  !  entry in this routine.                                                      |
  !==============================================================================|
  use mod_flow_field
  implicit none
  !------------------------------------------------------------------------------|
  real(DP), intent(in)    :: x, y
  integer,  intent(inout) :: host
  logical,  intent(inout) :: indomain
  logical                 :: isintriangle
  !------------------------------------------------------------------------------|
  integer                      :: min_loc
  real,    dimension(ELEMENTS) :: radlist
  real                         :: radlast
  integer, dimension(1)        :: locij
  integer                      :: peng
  !==============================================================================|

  if (.not.indomain) return
  radlist = (XC - x)**2 + (YC - y)**2 ! faster jm
  radlast = 0.0
  peng = 0

  do
    peng = peng + 1
    locij = minloc(radlist, radlist > radlast)
    min_loc = locij(1)

    if (min_loc == 0 .or. peng >= 32) then
      indomain = .false.
      return
    end if

    radlast = radlist(min_loc)

    if (isintriangle(min_loc,x,y)) then
      host = min_loc
      return
    end if
  end do
end subroutine fhe_robust

!==============================================================================|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!==============================================================================|

logical function isintriangle(ele,x,y) 
  !==============================================================================|
  !  Determine if Point (x0,y0) Is In triangle defined by nodes (xt(3),yt(3))    |
  !  Using Algorithm used for Scene Rendering in Computer Graphics               |
  !==============================================================================|
  use mod_flow_field
  implicit none
  !------------------------------------------------------------------------------|
  integer,  intent(in) :: ele
  real(DP), intent(in) :: x, y
  !------------------------------------------------------------------------------|
  real               :: f1, f2, f3
  real, dimension(3) :: xtri, ytri
  !==============================================================================|
 
  xtri = XP(NV(ele,1:3))
  ytri = YP(NV(ele,1:3))

  f1 = (y-ytri(1))*(xtri(2)-xtri(1)) - (x-xtri(1))*(ytri(2)-ytri(1))
  f2 = (y-ytri(3))*(xtri(1)-xtri(3)) - (x-xtri(3))*(ytri(1)-ytri(3))
  f3 = (y-ytri(2))*(xtri(3)-xtri(2)) - (x-xtri(2))*(ytri(3)-ytri(2))

  if (f1*f3 >= 0.0 .and. f3*f2 >= 0.0) then
    isintriangle = .true.
  else
    isintriangle = .false.
  end if
end function isintriangle

