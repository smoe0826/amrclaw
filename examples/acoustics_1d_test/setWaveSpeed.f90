FUNCTION setWaveSpeed(meqn,npts,xivals,qvals)

  !This function returns the wave speed (eigenvalues)
  !of the PDE at every point in space. This specific
  !function is written for constant coefficient advection

  IMPLICIT NONE
  !Inputs
  INTEGER, INTENT(IN) :: meqn,npts
  DOUBLE PRECISION, DIMENSION(1:npts), INTENT(IN) :: xivals
  DOUBLE PRECISION, DIMENSION(1:npts,1:meqn), INTENT(IN) :: qvals

  ! Common block
  real(kind=8) :: rho, bulk, cc, zz

  common /cparam/ rho, bulk, cc, zz

  !Outputs
  DOUBLE PRECISION, DIMENSION(1:npts) :: setWaveSpeed

  Integer :: i

  Do i=1,npts
     setWaveSpeed(i)=1.d0
  EndDo

End FUNCTION setWaveSpeed
