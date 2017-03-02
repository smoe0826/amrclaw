SUBROUTINE Riemann(meqn,qvalsl,qvalsr,Fl,Fr,xpt,RiemannSolve)
  !===============================!
  !This function evaluates an approximate
  !Local Lax Friedrichs flux. This is an
  !approximate Riemann solve 
  !===============================!
  IMPLICIT NONE
  ! Inputs
  Integer, Intent(IN) :: meqn
  DOUBLE PRECISION, DIMENSION(1), INTENT(IN) :: xpt
  Double Precision, Dimension(1:meqn), INTENT(IN) :: Fl,Fr
  Double Precision, Dimension(1:meqn), INTENT(IN) :: qvalsl,qvalsr
  DOUBLE PRECISION, DIMENSION(0:1,1:meqn), INTENT(INOUT) :: RiemannSolve
  ! Local Variables
  Double Precision, Dimension(1) :: s1,s2,smax

 
  !Interface blocks
  Interface

     FUNCTION setWaveSpeed(meqn,npts,xivals,qvals)
        IMPLICIT NONE
        !Inputs
        INTEGER, INTENT(IN) :: meqn,npts
        DOUBLE PRECISION, DIMENSION(1:npts), INTENT(IN) :: xivals
        DOUBLE PRECISION, DIMENSION(1:npts,1:meqn), INTENT(IN) :: qvals
        !Outputs
        Double Precision, Dimension(1:npts) :: setWaveSpeed
        Integer :: i
     End Function setWaveSpeed
  End Interface

  
  
  s1=setWaveSpeed(meqn,1,xpt,qvalsl)
  s2=setWaveSpeed(meqn,1,xpt,qvalsr)
  
  smax(1)=max(abs(s1(1)),abs(s2(1)))  
    
  RiemannSolve(0,:)=0.5d0*((Fl(:)+Fr(:))+smax(1)*(qvalsl(:)-qvalsr(:)))
  RiemannSolve(1,:)=0.5d0*((Fl(:)+Fr(:))+smax(1)*(qvalsl(:)-qvalsr(:)))
 

End SUBROUTINE Riemann
