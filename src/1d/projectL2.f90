SUBROUTINE projectL2(meqn,nQuad,maxPolyDegree,qvals,quadWeights,legendreQuad,coeffs)
  ! ===========================================================================
  ! Computes projection of function q(x) whos values are given by qvals onto Legendre basis polynomials
  ! for a single element. Note that qvals are assumed to be given at quadrature locations.
  ! INPUTS: meqn - number of equations solved
  !         maxPolyDegree - maximum Legendre polynomial degree used in basis
  !         nQuad - number of quadrature nodes per element
  !         quadWeights - Gaussian quadrature weights
  !         legendreQuad(k,i) - kth Legendre polynomial evaluated at ith quadrature node
  ! OUTPUTS: coeffs(k) - coefficient of kth degree Legendre polynomial in series expansion
  ! ===========================================================================
  use amr_module
  IMPLICIT NONE
  ! Inputs
  INTEGER :: meqn, nQuad, maxPolyDegree
  DOUBLE PRECISION, DIMENSION(0:nQuad,1:meqn), INTENT(IN) :: qvals
  DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: quadWeights
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad), INTENT(IN) :: legendreQuad
  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,1:meqn), INTENT(OUT) :: coeffs
  ! Local variables
  INTEGER :: k,m

  coeffs = 0D0
  DO m=1,meqn
    DO k=0,maxPolyDegree
      coeffs(k,m) = 0.5D0*(2D0*k+1D0)*SUM(quadWeights(:)*qvals(:,m)*legendreQuad(k,:))
      !prinT *,k,coeffs(k,m)
    ENDDO !k
  ENDDO !m

END SUBROUTINE projectL2
