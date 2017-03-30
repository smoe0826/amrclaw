!
!  :::::::::::::: BOUND :::::::::::::::::::::::::::::::::::::::::::
!     This routine sets the boundary values for a given grid 
!     at level level.
!     We are setting the values for a strip ng zones wide on
!     both borders.
!
!     Outputs from this routine:
!     The values around the border of the grid are inserted
!     directly into the enlarged valbig array.
!
!     This routine calls the routine filpatch
!     which for any block of mesh points on a given level,
!     intersects that block with all grids on that level and with
!     the physical boundaries, copies the values into the
!     appropriate intersecting regions, and interpolates the remaining
!     cells from coarser grids as required.
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine bound(time,nvar,ng,valbig,mitot,mptr,aux,naux)

  use amr_module, only: rnode, node, hxposs, cornxlo, cornxhi
  use amr_module, only: ndilo, ndihi
  use amr_module, only: nestlevel, xlower, xupper
  use amr_module, only: xperdom

  implicit none

  ! Input
  integer, intent(in) :: nvar, ng, mitot, mptr, naux
  real(kind=8), intent(in) :: time
  real(kind=8), intent(in out) :: valbig(nvar,mitot)
  real(kind=8), intent(in out) :: aux(naux,mitot)

  ! Locals
  integer :: ilo, ihi, level
  real(kind=8) :: xleft, xright, hx, xl, xr
  real(kind=8) :: xloWithGhost,  xhiWithGhost
  logical      :: patchOnly

  xleft  = rnode(cornxlo, mptr)
  xright = rnode(cornxhi, mptr)
  ilo    = node(ndilo, mptr)
  ihi    = node(ndihi, mptr)
  level  = node(nestlevel, mptr)
  hx     = hxposs(level)

  !print *,'here23=',shape(valbig)

  xloWithGhost = xleft  - ng*hx
  xhiWithGhost = xright + ng*hx
  ! used in filaptch for bc1amr: for patches it is called. for full grids called from bound below
  patchOnly = .false.

  !print *,'patch boundary=',xloWithGhost,xhiWithGhost,' level=',level
  !print *,'left'

  ! left boundary
  xl = xleft - ng*hx
  xr = xleft
  if ((xl < xlower) .and. xperdom) then
     !print *,'thing11'
     call  prefilrecur(level,nvar,valbig,aux,naux,time,mitot,1, &
          ilo-ng,ilo-1,ilo-ng,ihi+ng,patchOnly,0,1)
     !print *,'out11'
  else
     !print *,'thing1',level,'left=',xl,xl+hx,xl+ng*hx
     call filrecur(level,nvar,valbig,aux,naux,time,mitot,1,ilo-ng, &
          ilo-1,patchOnly,mptr,0,1)
     !print *,'out1',level
  endif

  !print *,'right'

  ! right boundary
  xl = xright
  xr = xright + ng*hx

  if ((xr .gt. xupper) .and. xperdom) then
     !print *,'thing21'
     call  prefilrecur(level,nvar,valbig,aux,naux,time,mitot, &
          mitot-ng+1,ihi+1,ihi+ng,ilo-ng,ihi+ng,patchOnly,1,1)
     !print *,'out21'
  else
     !print *,'thing2'
     call filrecur(level,nvar,valbig,aux,naux,time,mitot, &
          mitot-ng+1,ihi+1,ihi+ng,patchOnly,mptr,1,1)
     !print *,'out2'
  endif

  !print *,''

  ! set all exterior (physical)  boundary conditions for this grid at once
  ! used to be done from filpatch, but now only for recursive calls with new patch
  ! where the info matches. more efficient to do whole grid at once, and avoid copying
  call bc1amr(valbig,aux,mitot,nvar,naux,hx,level,time,    &
              xloWithGhost,xhiWithGHost)

end subroutine bound
