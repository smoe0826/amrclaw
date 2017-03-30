c
c -------------------------------------------------------------
c
      subroutine stepgrid(q,fm,fp,mitot,mbc,dt,dtnew,dx,
     &                  nvar,xlow,time,mptr,maux,aux)
c
c          
c ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
c take a time step on a single grid. overwrite solution array q. 
c A modified version of the clawpack routine step is used.
c
c return fluxes in fm and fp.
c patch has room for ghost cells (mbc of them) around the grid.
c everything is the enlarged size (mitot by mjtot).
c
c mbc       = number of ghost cells  (= lwidth)
c mptr      = grid number  (for debugging)
c xlow      =  left corner of enlarged grid (including ghost cells).
c dt        = incoming time step
c dx        = mesh width for this grid
c dtnew     = return suggested new time step for this grid's soln.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use amr_module
      use modalDGmod
      implicit double precision (a-h,o-z)
      external rp1_ptwise

      common /comxyt/ dtcom,dxcom,tcom,icom

      dimension q(nvar,mitot)
      dimension fp(nvar,mitot)
      dimension fm(nvar,mitot)
      dimension aux(maux,mitot)

      double precision, dimension(:), allocatable :: quadNodes
      double precision, dimension(:), allocatable :: quadWeights
      double precision, dimension(:,:), allocatable :: legendreQuad
      double precision, dimension(:,:), allocatable :: translateMon
      double precision, dimension(:,:), allocatable :: monQuad
      double precision, dimension(:,:), allocatable :: dlegendreQuad
      double precision, dimension(:,:), allocatable :: localCoeffs
      double precision, dimension(:,:), allocatable :: numFlux
      double precision, dimension(:,:), allocatable :: numFlux1
      double precision, dimension(:,:), allocatable :: numFluxL
      double precision, dimension(:,:), allocatable :: numFluxR
      double precision, dimension(:), allocatable :: xivals
      double precision, dimension(:,:), allocatable :: qOut
      double precision, dimension(:,:), allocatable :: qvalsl
      double precision, dimension(:,:), allocatable :: qvalsr
      double precision, dimension(1) :: s1
      double precision, dimension(1,1) :: xpt
      double precision, dimension(1) :: xpt1
      double precision, dimension(:,:), allocatable :: coeffL
      double precision, dimension(:,:), allocatable :: coeffR
      double precision, dimension(:,:), allocatable :: fluxQuad
      double precision, dimension(:,:), allocatable :: qQuad
      double precision, dimension(:,:), allocatable :: fluxL
      double precision, dimension(:,:), allocatable :: fluxR

      double precision, dimension(:,:), allocatable :: qpointL
      double precision, dimension(:,:), allocatable :: qpointR
      double precision, dimension(:,:), allocatable :: fvalsL
      double precision, dimension(:,:), allocatable :: fvalsR
      double precision, dimension(:,:), allocatable :: qL
      double precision, dimension(:,:), allocatable :: qR
      double precision, dimension(:,:,:), allocatable :: fQuad


      double precision, dimension(:,:), allocatable :: testMat
      double precision, dimension(:,:), allocatable :: testMat_eqn

      double precision, dimension(:,:), allocatable :: k0
      double precision, dimension(:), allocatable :: k1

      logical    debug,  dump
      data       debug/.false./,  dump/.false./

c
c     # set tcom = time.  This is in the common block comxyt that could
c     # be included in the Riemann solver, for example, if t is explicitly
c     # needed there.

      nQuad = DGorder
      deg   = DGorder
      maxPolyDegree = DGorder-1
      meqn   = nvar
      num_eqn= meqn/DGorder

      allocate(k0(0:maxPolyDegree,1:meqn))
      allocate(k1(1:meqn))

      allocate(quadNodes(0:nQuad),quadWeights(0:nQuad))
      allocate(translateMon(0:maxPolyDegree,0:maxPolyDegree))
      allocate(legendreQuad(0:maxPolyDegree,0:nQuad))
      allocate(monQuad(0:maxPolyDegree,0:nQuad))
      allocate(dlegendreQuad(0:maxPolyDegree,0:nQuad))

      allocate(fluxQuad(1:nQuad+1,1:num_eqn))
      allocate(fluxL(1,1:num_eqn))
      allocate(fluxR(1,1:num_eqn))

      allocate(qQuad(1:nQuad+1,1:num_eqn))
      allocate(qpointL(1,1:num_eqn))
      allocate(qpointR(1,1:num_eqn))

      allocate(qL(1:mitot,1:num_eqn))
      allocate(qR(1:mitot,1:num_eqn))
      allocate(fQuad(1:mitot,1:nQuad+1,1:num_eqn))

      allocate(testMat(0:maxPolyDegree,0:maxPolyDegree))
      allocate(testMat_eqn(0:maxPolyDegree,1:num_eqn))

      CALL gaussQuadNodes(nQuad+1,quadNodes)
      CALL gaussQuadWeights(nQuad+1,quadNodes,quadWeights)



      ! Precompute Legendre polynomial and its derivative values at quadrature locations
      DO i=0,nQuad
       DO k=0,maxPolyDegree
         legendreQuad(k,i) = legendre(quadNodes(i),k)
         monQuad(k,i) = monomial(quadNodes(i),1)
         dlegendreQuad(k,i) = dlegendre(quadNodes(i),k)
       ENDDO !k
      ENDDO !i


      DO i=0,maxPolyDegree
       DO k=0,maxPolyDegree
         translateMon(k,i) = legendre2Mon(k,i)
       ENDDO !k
      ENDDO !i

    

      allocate(localCoeffs(0:maxPolyDegree,1:num_eqn))
      allocate(qOut(1:nQuad+1,1:meqn),xivals(1:nQuad+1))
      allocate(numFlux1(0:1,1:num_eqn),numFlux(0:1,1:num_eqn))
      allocate(numFluxL(1:mitot,1:num_eqn),numFluxR(1:mitot,1:num_eqn))
      allocate(qvalsl(1,1:num_eqn),qvalsr(1,1:num_eqn))
      allocate(coeffL(0:maxPolyDegree,1:num_eqn))
      allocate(coeffR(0:maxPolyDegree,1:num_eqn))
          


 
      tcom = time

      if (dump) then
         write(outunit,*) "dumping grid ",mptr," at time ",time
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i),ivar=1,nvar)
c    .                  ,(aux(ivar,i,j),ivar=1,maux)
 545        format(2i4,5e15.7)
         end do
         end do
      endif
c
      mx = mitot - 2*mbc
      mbig = mx       !# size for 1d scratch array
      xlowmbc = xlow + mbc*dx

c     # method(2:7) and mthlim
c     #    are set in the amr2ez file (read by amr)
c
      method(1) = 0
c
c
      call b4step1(mbc,mx,nvar,q,
     &             xlowmbc,dx,time,dt,maux,aux)
c
c
c     # take one step on the conservation law:
c
      call step1(mbig,nvar,maux,
     &           mbc,mx,
     &              q,aux,dx,dt,cflgrid,
     &              fm,fp,rp1_ptwise)
c
c
        write(outunit,1001) mptr, node(nestlevel,mptr),cflgrid
 1001   format(' Courant # of grid', i4,
     &        ' on level', i3, ' is  ', e10.3)
c

!$OMP  CRITICAL (cflm)

        cfl_level = dmax1(cfl_level,cflgrid)

!$OMP END CRITICAL (cflm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create the flux arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       i=1
       do k=1,num_eqn
               localCoeffs(:,k) = q((k-1)*DGorder+1:k*DGorder,i)
       end do

       xivals(:)=0.5D0*quadNodes(:)*dx+xlow+(i+0.5d0)*dx
       xpt(1,1)=1.d0
       call fluxFunction(localCoeffs,xivals,qpointR,fluxR,1,num_eqn
     &    ,dt,dx,translateMon,xpt)

       numFluxR(i,:)=fluxR(1,:)
       qR(i,:)=qpointR(1,:)

       xpt1(1)=1.d0
       CALL evaluateExpansion(localCoeffs,xpt1,qvalsr,1,DGorder-1
     & ,num_eqn) 
!!!!!!!!!!!!!!!!!!!!!!!Start Main Loop!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=2,mitot-1
          localCoeffs=0.d0
          do k=1,num_eqn
               localCoeffs(:,k) = q((k-1)*DGorder+1:k*DGorder,i)
          end do

          xivals(:)=0.5D0*quadNodes(:)*dx+xlow+(i+0.5d0)*dx

c          CALL evaluateExpansion(localCoeffs,quadNodes(:),qOut,nQuad+1
c     &      ,deg,num_eqn)

 
          call fluxFunction(localCoeffs,xivals,qQuad,fluxQuad,nQuad+1
     &    ,num_eqn,dt,dx,translateMon,monQuad)
 
          
          fQuad(i,:,:)=fluxQuad(:,:)


          xpt(1,1)=-1.d0
          call fluxFunction(localCoeffs,xivals,qpointL,fluxL,1,num_eqn
     &    ,dt,dx,translateMon,xpt)
  
          numFluxL(i,:)=fluxL(1,:)
          qL(i,:)=qpointL(1,:)

          xpt(1,1)=1.d0
          call fluxFunction(localCoeffs,xivals,qpointR,fluxR,1,num_eqn
     &    ,dt,dx,translateMon,xpt)

          numFluxR(i,:)=fluxR(1,:)
          qR(i,:)=qpointR(1,:)
         

c          testMat=matmul(legendreQuad,transpose(monQuad))
c          testMat_eqn=matmul(translateMon,localCoeffs)
      
c      do k=0,maxPolyDegree
c         print *,"k=",k,localCoeffs(k,:)!testMat_eqn(k,:)
c        print *,"k=",k,0.5D0*(2D0*k+1D0)*SUM(quadWeights(:)
c     & *testMat_eqn(k+1,1)*monQuad(k,:)*legendreQuad(2,:))
c      end do
       

c          xpt(1)=0.5d0*dx+xlow+(i+0.5d0)*dx
c          call Riemann(num_eqn,qvalsl,qvalsr,xpt,numFlux1(:,:))
c          numFluxR(i,:)=numFlux1(0,:)
c
c          coeffL=localCoeffs
        end do   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Main Time Stepping Loop 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c       # update q
        dtdx = dt/dx
        do 50 i=mbc+1,mitot-mbc



        xpt1(1)=-1.d0
        call Riemann(num_eqn,qR(i-1,:),qL(i,:)
     &    ,numFluxR(i-1,:),numFluxL(i,:),xpt1,numFlux1(:,:))
c        call AcousticsRiemann(num_eqn,qR(i-1,:),qL(i,:)
c     &    ,numFluxR(i-1,:),numFluxL(i,:),xpt1,numFlux1(:,:))

        numFlux(0,:)=numFlux1(1,:)


        xpt1(1)=1.d0
        call Riemann(num_eqn,qR(i,:),qL(i+1,:)
     &    ,numFluxR(i,:),numFluxL(i+1,:),xpt1,numFlux1(:,:))
c        call AcousticsRiemann(num_eqn,qR(i,:),qL(i+1,:)
c     &    ,numFluxR(i,:),numFluxL(i+1,:),xpt1,numFlux1(:,:))

        numFlux(1,:)=numFlux1(0,:)


          fluxQuad(:,:)=fQuad(i,:,:)

         if (mcapa.eq.0) then

          Do k=0,maxPolyDegree
            call forcingFunction(k,fluxQuad,quadWeights
     & ,dlegendreQuad(k,:),numFlux,meqn,nQuad,k1)
           k0(k,:)=k1(:)
          EndDo !k

          do k=1,num_eqn
             q((k-1)*DGorder+1:k*DGorder,i)=q((k-1)*DGorder+1:k*DGorder
     &       ,i)-dtdx*k0(:,k)
          end do

          do m=1,nvar
c
c            # no capa array.  Standard flux differencing:

c             q(m,i) = q(m,i) - dtdx
c            q(m,i) = q(m,i)
c     &           - dtdx * (fm(m,i+1) - fp(m,i))
          end do
         else
          do m=1,nvar
c            # with capa array.
c            q(m,i) = q(m,i)
c     &          - (dtdx * (fm(m,i+1) - fp(m,i))) / aux(mcapa,i)
          end do
         endif

 50      continue


c        print *,q

c
c
      if (method(4).eq.1) then
c        # with source term:   use Godunov splitting
         call src1(nvar,mbc,mx,xlowmbc,dx,
     &             q,maux,aux,time,dt)
         endif
c
c
c
c     # output fluxes for debugging purposes:
      if (debug) then
         write(dbugunit,*)" fluxes for grid ",mptr
c        do 830 j = mbc+1, mjtot-1
            do 830 i = mbc+1, mitot-1
               write(dbugunit,831) i,fm(1,i),fp(1,i)
               do 830 m = 2, meqn
                  write(dbugunit,832) fm(m,i),fp(m,i)
  831          format(2i4,4d16.6)
  832          format(8x,4d16.6)
  830    continue
      endif

c
c
c For variable time stepping, use max speed seen on this grid to 
c choose the allowable new time step dtnew.  This will later be 
c compared to values seen on other grids.
c

       if (cflgrid .gt. 0.d0) then
           dtnew = dt*cfl/cflgrid
         else
c          # velocities are all zero on this grid so there's no 
c          # time step restriction coming from this grid.
            dtnew = rinfinity
          endif

c     # give a warning if Courant number too large...
c
      if (cflgrid .gt. cflv1) then
            write(*,810) cflgrid
            write(outunit,810) cflgrid, cflv1
  810       format('*** WARNING *** Courant number  =', d12.4,
     &              '  is larger than input cfl_max = ', d12.4)
            endif
c
      if (dump) then
         write(outunit,*) "dumping grid ",mptr," after stepgrid"
         do i = 1, mitot
            write(outunit,545) i,(q(ivar,i),ivar=1,nvar)
         end do
      endif
      return
      end


