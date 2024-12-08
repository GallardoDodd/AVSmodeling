PROGRAM euler_implicito
REAL, PARAMETER :: &
  Tc = 309.65, &        ! Temperatura del cos [K]
  k = 0.56, &           ! Conductivitat tèrmica [W/m·K]
  sigma = 0.472, &      ! Conductivitat elèctrica [S/m]
  D = 0.02, &           ! Amplitud del vas [m]
  V = 40.0, &           ! Potencial elèctric [V]
  c_v = 3686.0, &       ! Capacidad calorífica específica [J/kg·K]
  p = 1081.0            ! Densitat [kg/m^3]

INTEGER, PARAMETER :: N=101                          ! Nombre de punts de la idscretització espaial

REAL, PARAMETER :: &
    alpha=k/(c_v*p), &                               ! Coeficient de difusió [m^2/s]
    Tc_norm=1, &                                     ! Temperatura corporal normalitzada
    tin=0, &                                         ! Temps inicial (normalitzat)
    tfin=0.025, &                                    ! Temps final (normalitzat)
    xin=0, &                                         ! Coordenada espaial inicial (normalitzada)
    xfin=D*SQRT((0.5*sigma*V**2)/(Tc*k*D**2))        ! Coordenada espaial final (normalitzada)

REAL, PARAMETER :: &
    ax=(xfin-xin)/N, &          ! Increment de distància
    at1=(ax)**2, &              ! Increment temporal 1
    at2=0.5*(ax)**2             ! Increment temporal 2

INTEGER :: Nat1 = int((tfin-tin)/at1), Nat2 = int((tfin-tin)/at2), info     ! Nombre de punts per a cada cas i variable muda

INTEGER, ALLOCATABLE :: ipiv(:)

REAL, ALLOCATABLE :: &
    T(:,:), &           ! Matriu de temperatures
    xlin(:), &          ! Discretització espaial
    tlin(:), &          ! Discretització temporal
    M(:,:), &           ! Matriu euler_implícit
    M_inv(:,:), &       ! Inversa de la matriu euler_implícit
    b(:), &             ! Vector independent
    T_new(:), &         ! Matriu de temperatures temporal
    work(:)             ! Espai de treball per al mètode LU d'inversió de matrius

!----///---PRIMER CAS---///---:
!Es crea la malla de punts d'espai i temps
ALLOCATE(xlin(N), tlin(Nat1))
xlin = linspace(xin,xfin,N)
tlin = linspace(tin,tfin,Nat1)

!Calculem la matriu M i la seva inversa:
ALLOCATE(M(N-2,N-2),M_inv(N-2,N-2),ipiv(N-2),work(2*(N-2)))
M = 0.0
DO i = 1, N-2
    IF (i > 1) M(i, i-1) = -1.0     
    M(i, i) = 3.0                 
    IF (i < N-2) M(i, i+1) = -1.0   
END DO
CALL invertir(M, M_inv, N-2, ipiv, work, info)

!Condició que informa si la matriu construida és invertible o no
IF (info /= 0) THEN
    PRINT *, "La matriu no és invertible, info =", info
END IF
!Vector constant que conté les condicions inicials per a les iteracions
ALLOCATE(b(N-2))
b = at1
b(1)= Tc_norm + at1
b(N-2)= Tc_norm + at1

!Estructura del mapa de temperatures T(i,j): i és la component temporal, j la component espaial
ALLOCATE(T(Nat1,N))

!Condicions de contorn
T(:,1)=Tc_norm
T(:,N)=Tc_norm

!Condicions inicials
T(1,:)=Tc_norm
ALLOCATE(T_new(N-2))

!Iteracions aplicant l'equació matricial
DO i=2,Nat1
    T_new = matmul(M_inv,T(i-1,2:N-1)+b)
    T(i,2:N-1) = T_new
END DO
DEALLOCATE(T_new)

!Imprimim les dades com a llistes de dades que pugui entendre gnuplot o altre programa que les grafiqui:
PRINT*, "> Open dades_3D_imp_at1"
OPEN(unit=10,file='dades_3D_imp_at1.txt',status='replace',action='write')
PRINT*,"> Open dades_2D_imp_at1"
OPEN(unit=20,file='dades_2D_imp_at1.txt',status='replace',action='write')

!Dades gràfic 3D:
DO i=1,N
    DO j=1,Nat1
        WRITE(10, '(F10.6,F10.6,F10.3)') xlin(i), tlin(j), T(j,i)
    END DO
END DO

!Dades gràfic 2D a y =0.025
DO i=1,N
    WRITE(20,'(F10.6,F10.6)') xlin(i), T(Nat1,i)
END DO
CLOSE(10)
PRINT *, "> 'dades_3D_imp_at1' closed correctly"
CLOSE(20)
PRINT *, "> 'dades_2D_imp_at1' closed correctly"
DEALLOCATE(xlin,tlin,M,M_inv,ipiv,work,b,T)

!---///---SEGON CAS---///---:
!Es crea la malla de punts d'espai i temps
ALLOCATE(xlin(N), tlin(Nat2))
xlin = linspace(xin,xfin,N)
tlin = linspace(tin,tfin,Nat2)

!Calculem la matriu M i la seva inversa:
ALLOCATE(M(N-2,N-2),M_inv(N-2,N-2),ipiv(N-2),work(2*(N-2)))
M = 0.0
DO i = 1, N-2
    IF (i > 1) M(i, i-1) = -0.5     
    M(i, i) = 2.0                 
    IF (i < N-2) M(i, i+1) = -0.5   
END DO
CALL invertir(M, M_inv, N-2, ipiv, work, info)

!Condició que informa si la matriu construida és invertible o no
IF (info /= 0) THEN
    PRINT *, "La matriu no és invertible, info =", info
END IF
!Vector constant que conté les condicions inicials per a les iteracions
ALLOCATE(b(N-2))
b = at2
b(1)= 0.5*Tc_norm + at2
b(N-2)= 0.5*Tc_norm + at2

!Estructura del mapa de temperatures T(i,j): i és la component temporal, j la component espaial
ALLOCATE(T(Nat2,N))

!Condicions de contorn
T(:,1)=Tc_norm
T(:,N)=Tc_norm

!Condicions inicials
T(1,:)=Tc_norm
ALLOCATE(T_new(N-2))

!Iteracions aplicant l'equació matricial
DO i=2,Nat2
    T_new = matmul(M_inv,T(i-1,2:N-1)+b)
    T(i,2:N-1) = T_new
END DO
DEALLOCATE(T_new)

!Imprimim les dades com a llistes de dades que pugui entendre gnuplot o altre programa que les grafiqui:
PRINT*, "> Open dades_3D_imp_at2"
OPEN(unit=10,file='dades_3D_imp_at2.txt',status='replace',action='write')
PRINT*,"> Open dades_2D_imp_at2"
OPEN(unit=20,file='dades_2D_imp_at2.txt',status='replace',action='write')

!Dades gràfic 3D:
DO i=1,N
    DO j=1,Nat2
        WRITE(10, '(F10.6,F10.6,F10.3)') xlin(i), tlin(j), T(j,i)
    END DO
END DO

!Dades gràfic 2D a y =0.025
DO i=1,N
    WRITE(20,'(F10.6,F10.6)') xlin(i), T(Nat2,i)
END DO
CLOSE(10)
PRINT *, "> 'dades_3D_imp_at2' closed correctly"
CLOSE(20)
PRINT *, "> 'dades_2D_imp_at2' closed correctly"
DEALLOCATE(xlin,tlin,M,M_inv,ipiv,work,b,T)

!Subrutines que cridarem al executar el programa
CONTAINS
    FUNCTION linspace(valin, valfin, num) result(lista)
        REAL, INTENT(in) :: valin, valfin   
        INTEGER, INTENT(in) :: num 
        REAL, ALLOCATABLE :: lista(:)
        INTEGER :: k
        REAL :: paso
        
        paso = (valfin - valin) / (num - 1)
        ALLOCATE(lista(num))
        
        DO k = 1, num
            lista(k) = valin + paso*(k - 1)
        END DO
    END FUNCTION linspace
    SUBROUTINE invertir(A, A_inv, n, ipiv, work, info)
        REAL, INTENT(in) :: A(:,:)          
        REAL, INTENT(out) :: A_inv(:,:)     
        INTEGER, INTENT(in) :: n           
        INTEGER, INTENT(out) :: ipiv(:)     
        REAL, INTENT(out) :: work(:)       
        INTEGER, INTENT(out) :: info       

        A_inv = A
        CALL sgetrf(n, n, A_inv, n, ipiv, info)
        IF (info /= 0) RETURN 

        CALL sgetri(n, A_inv, n, ipiv, work, SIZE(work), info)
    END SUBROUTINE invertir
END PROGRAM euler_implicito