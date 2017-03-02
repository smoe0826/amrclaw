! ##############################################################
! Module containing Modal Discontinuous Galerkin basis functions,
! quadrature data for use in 1D simulations.
! ##############################################################

MODULE modalDGmod
        IMPLICIT NONE
        INTEGER, PARAMETER, PRIVATE :: DOUBLE=KIND(1D0)

        CONTAINS

    ! ########################################################################
    ! N-choose-k Function
    ! ########################################################################
            REAL(KIND=DOUBLE) FUNCTION choose(alpha,k)
                    IMPLICIT NONE
                    INTEGER, INTENT(IN) :: k
                    REAL(KIND=DOUBLE), INTENT(IN) :: alpha
                    INTEGER :: i
                    REAL(KIND=DOUBLE) :: HOLDER

                    HOLDER = 1D0

                    DO i = 1,k
                        HOLDER = HOLDER*((alpha-DBLE(k-i))/(DBLE(i)))
                    END DO
                    choose = HOLDER
                END FUNCTION choose

    ! ########################################################################
    ! Legendre Polynomial function of degree N
    ! ########################################################################
             REAL(KIND=DOUBLE) FUNCTION legendre(x,N)
                    IMPLICIT NONE
                    REAL(KIND=DOUBLE), INTENT(IN) :: x
                    REAL(KIND=DOUBLE) :: HOLDER
                    INTEGER, INTENT(IN) :: N
                    INTEGER :: k

                    HOLDER = 0.D0
                    DO k = 0,N
                        HOLDER = HOLDER + choose(DBLE(N),k)*choose((N+k-1)/2D0,N)*x**k
                    END DO

                    legendre = HOLDER*(2**N)

                END FUNCTION legendre

    ! ########################################################################
    ! Derivative of Legendre Polyomial of degree N
    ! ########################################################################
             REAL(KIND=DOUBLE) FUNCTION dlegendre(x,N)
                    IMPLICIT NONE
                    REAL(KIND=DOUBLE),INTENT(IN) :: x
                    REAL(KIND=DOUBLE) :: HOLDER
                    INTEGER, INTENT(IN) :: N ! Order of legendre polynomial
                    INTEGER :: k

                    HOLDER = 0.D0
                    DO k = 1,N
                        HOLDER = HOLDER + k*choose(DBLE(N),k)*choose((N+k-1)/2D0,N)*x**(k-1)
                    END DO

                    dlegendre = HOLDER*(2**N)

                END FUNCTION dlegendre

    ! ###########################################################################################################
    ! Subroutine for computing Gaussian quadrature nodes based on the derivative of Mth Order Legendre Polynomial
    ! For Modal DG, we require M=N+1 nodes, where N is the highest order of Legendre polynomial being used
    ! ###########################################################################################################
                SUBROUTINE gaussQuadNodes(M,nodes)
                    IMPLICIT NONE
                    INTEGER, INTENT(IN) :: M
                    REAL(KIND=DOUBLE), DIMENSION(0:M-1), INTENT(OUT) :: nodes
                    REAL(KIND=DOUBLE) :: xnew,xold,error,tol, PI
                    INTEGER :: k

                    PI = DACOS(-1D0)

                    tol  = 10D-12

                    DO k = 0,M-1
                        error = 1D0
                        xold = -1D0*DCOS(((2*k+1)/(2D0*M))*PI)

                        DO WHILE (error>tol)
                            xnew = xold - (legendre(xold,M))/(1D0*dlegendre(xold,M))
                            error = DABS(xnew-xold)
                            xold = xnew
                        END DO
                        nodes(k) = xold
                    END DO
                END SUBROUTINE gaussQuadNodes

  ! #######################################################################################################
  ! Computing weights associated with N+1 nodes for quadratures on [-1,1]
  ! For Modal DG, we require M=N+1 weights, where N is the highest order of Legendre polynomial being used
  ! #######################################################################################################
        SUBROUTINE gaussQuadWeights(M,nodes,wghts)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            REAL(KIND=DOUBLE), DIMENSION(0:M-1), INTENT(IN) :: nodes
            REAL(KIND=DOUBLE), DIMENSION(0:M-1), INTENT(OUT) :: wghts
            INTEGER :: k

            DO k = 0,M-1
              wghts(k) = 2D0/( (1-nodes(k)**2)*(dlegendre(nodes(k),M))**2 )
            END DO

        END SUBROUTINE gaussQuadWeights


    ! ########################################################################
    ! Legendre Polynomial of degree N coefficients
    ! ########################################################################
             REAL(KIND=DOUBLE) FUNCTION legendre2Mon(k,N)
                    IMPLICIT NONE
                    INTEGER, INTENT(IN) :: k
                    REAL(KIND=DOUBLE) :: HOLDER
                    INTEGER, INTENT(IN) :: N

                    HOLDER = choose(DBLE(N),k)*choose((N+k-1)/2D0,N)

                    legendre2Mon = HOLDER*(2**N)

                END FUNCTION legendre2Mon

    ! ########################################################################
    ! Monomial Evaluated at point
    ! ########################################################################
             REAL(KIND=DOUBLE) FUNCTION monomial(x,k)
                    IMPLICIT NONE
                    REAL(KIND=DOUBLE),INTENT(IN) :: x
                    REAL(KIND=DOUBLE) :: HOLDER
                    INTEGER, INTENT(IN) :: k ! Order of legendre polynomial

                    HOLDER = 0.D0
                    HOLDER = x**(k)

                    monomial = HOLDER

                END FUNCTION monomial
END MODULE modalDGmod
