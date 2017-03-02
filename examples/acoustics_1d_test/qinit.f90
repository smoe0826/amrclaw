subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.
    ! This default version prints an error message since it should
    ! not be used directly.  Copy this to an application directory and
    ! loop over all grid cells to set values of q(1:meqn, 1:mx).


    use amr_module, only: DGorder
    use modalDGmod
    implicit none


   
    integer, intent(in) :: meqn,mbc,mx,maux

    integer :: i,nQuad,k,maxPolyDegree,j
    integer :: num_eqn
    real(kind=8) :: beta, xcell
    common /cqinit/ beta

    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    double precision, dimension(:), allocatable :: quadNodes
    double precision, dimension(:), allocatable :: quadWeights
    double precision, dimension(:,:), allocatable :: legendreQuad
    DOUBLE PRECISION, DIMENSION(:), allocatable :: xvals
    DOUBLE PRECISION, DIMENSION(:,:), allocatable :: qvals
    DOUBLE PRECISION, DIMENSION(:,:),allocatable :: tmpCoeffs


    INTERFACE
      SUBROUTINE projectL2(meqn,nQuad,maxPolyDegree,qvals,quadWeights,legendreQuad,coeffs)
        INTEGER :: meqn, nQuad, maxPolyDegree
        DOUBLE PRECISION, DIMENSION(0:nQuad),INTENT(IN) :: quadWeights
        DOUBLE PRECISION, DIMENSION(0:nQuad,1:meqn),INTENT(IN) :: qvals
        DOUBLE PRECISION, DIMENSION(0:maxPolyDegree,0:nQuad),INTENT(IN) :: legendreQuad
        !Outputs
        DOUBLE PRECISION, DIMENSION(0:maxPolyDegree), INTENT(OUT) :: coeffs
      END SUBROUTINE projectL2
    END INTERFACE


    nQuad = DGorder
    maxPolyDegree = DGorder - 1

    num_eqn=meqn/DGorder

    allocate(quadNodes(0:nQuad),quadWeights(0:nQuad))
    allocate(legendreQuad(0:maxPolyDegree,0:nQuad))
    allocate(qvals(0:nQuad,1:meqn))
    allocate(xvals(1:nQuad+1))
    allocate(tmpCoeffs(0:maxPolyDegree,1:meqn))



      CALL gaussQuadNodes(nQuad+1,quadNodes)
      CALL gaussQuadWeights(nQuad+1,quadNodes,quadWeights)


      ! Precompute Legendre polynomial and its derivative values at quadrature locations
      DO i=0,nQuad
       DO k=0,maxPolyDegree
         legendreQuad(k,i) = legendre(quadNodes(i),k)
       ENDDO !k
      ENDDO !i


    q(:,:)=0.d0
      do i=1,mx
         do j=0,nQuad         
          xcell = xlower + (i-0.5d0)*dx+0.5d0*dx*quadNodes(j)
          !print *,'initialize=',i,xlower,xcell
          qvals(j,1) = dexp(-beta * (xcell-0.0d0)**2)
          qvals(j,2) = 0.d0
         end do
         CALL projectL2(meqn,nQuad,maxPolyDegree,qvals,quadWeights,legendreQuad,tmpCoeffs)
         !q(1,i) = dexp(-beta * (xcell-0.3d0)**2)  
         !q(1*DGorder+1,i) = 0.d0

          do k=1,num_eqn
               !print *,'after','k=',k,tmpCoeffs(:,k)
               q((k-1)*DGorder+1:k*DGorder,i) =tmpCoeffs(:,k)
          end do

      enddo

end subroutine qinit

