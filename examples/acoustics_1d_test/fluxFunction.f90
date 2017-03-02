SUBROUTINE fluxFunction(qvals,xpts,qpoint,flux,npts,meqn,dt,dx,translateMon,monQuad)
  ! ========================================================================
  ! User specified flux functions evaluated at xpts
  ! INPUTS: npts - number of points to evaluate at
  !         xpts(npts) - points to evaluate at
  !         qvals(npts,m) - mth equation solution values at xpts
  ! OUTPUT:
  !         fluxFunction(npts,m) - analytic flux function evaluated at xpts
  ! ========================================================================

  use amr_module
  IMPLICIT NONE
  ! Inputs


  INTEGER, INTENT(IN) :: npts
  INTEGER, INTENT(IN) :: meqn
  DOUBLE PRECISION, DIMENSION(1:npts), INTENT(IN) :: xpts
  DOUBLE PRECISION, DIMENSION(1:DGorder,1:meqn), INTENT(IN) :: qvals
  DOUBLE PRECISION, INTENT(IN) :: dx
  DOUBLE PRECISION, INTENT(IN) :: dt
  double precision, dimension(1:DGorder,1:DGorder), INTENT(IN) :: translateMon
  double precision, dimension(1:DGorder,1:npts), INTENT(IN) :: monQuad
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:npts,1:meqn) :: flux
  DOUBLE PRECISION, DIMENSION(1:npts,1:meqn) :: qpoint
  ! Local variables
  INTEGER :: i,j,k,l,m,h
  integer :: max_flux_derivs
  double precision tmp, test

  integer, dimension(1:2) :: b4

  Double precision, dimension(1:meqn,1:DGorder,1:DGorder) :: Q_mixed_derivs
  Double precision, dimension(1:meqn,1:DGorder,1:DGorder,1:npts) :: F_mixed_derivs

  double precision, dimension(1:DGorder,1:meqn):: testMat_eqn


  double precision, dimension(1:DGorder):: testing

  max_flux_derivs=DGorder-1

  Q_mixed_derivs = 0.d0 
  F_mixed_derivs = 0.d0


  testMat_eqn=matmul(translateMon,qvals)


  !print *,'eqn1=',xpts,qvals(:,1)
  !print *,'eqn2=',xpts,testMat_eqn(:,1)


  b4(1)=2
  b4(2)=1

  do m=1,meqn
     
     do h=1,DGorder
        
       Q_mixed_derivs(m,h,1)=testMat_eqn(h,m)
 
     enddo
  enddo

  do k=2,DGorder
    do m=1,meqn
      do h=0,DGorder-k
         !print *,'k=',k,'h=',h
         tmp=-(h+1.d0)/(k-1.d0)*Q_mixed_derivs(b4(m),h+2,k-1)
         Q_mixed_derivs(m,h+1,k)=2.d0/dx*tmp
         !Q_mixed_derivs(m,h+1,k)=tmp
      end do
    end do
  end do
  !if (abs(Q_mixed_derivs(1,1,1)) .ge. 1.d-5) then
    !print *,'dt/dx=',dt/dx,'dt=',dt,'dx=',dx
    !print *,'drivatives=','k=1,m=1',Q_mixed_derivs(2,1,1)
    !print *,'drivatives=','k=1,m=2',Q_mixed_derivs(2,2,1)
    !print *,'drivatives=','k=2,m=1',Q_mixed_derivs(2,1,2)
    !print *,'drivatives=','k=2,m=2',Q_mixed_derivs(2,2,2)
  !endif

  do i=1,npts
  do m=1,meqn
    tmp = 0.d0
    do k=0,max_flux_derivs
       do h=0,max_flux_derivs
          !print *,'eqn=',m,'fluxk=',k,'h=',h,dt**k,monQuad(1,i)**h,'point=',monQuad(1,i)
          tmp = tmp + ( dt**k/(1.d0+k)*(monQuad(1,i)**h)*Q_mixed_derivs(m,h+1,k+1))
       end do
    end do
    Flux(i,b4(m))=tmp


    !k=0
    !tmp = 0.d0
    !    do h=0,max_flux_derivs
    !      tmp = tmp + ( dt**k/(1.d0+k)*(monQuad(1,i)**h)*Q_mixed_derivs(m,h+1,k+1))
    !    end do
    qpoint(i,m)=tmp
  end do
  end do


    
END SUBROUTINE fluxFunction
