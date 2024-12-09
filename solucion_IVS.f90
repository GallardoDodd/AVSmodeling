module parametros

REAL, PARAMETER :: &
  Tc = 36.5, &        ! Temperatura del cos [K]
  k = 0.56, &           ! Conductivitat tèrmica [W/m·K]
  sigma = 0.472, &      ! Conductivitat elèctrica [S/m]
  D = 0.02, &           ! Amplitud del vas [m]
  V = 40.0, &           ! Potencial elèctric [V]
  c_v = 3686.0, &       ! Capacidad calorífica específica [J/kg·K]
  p = 1081.0            ! Densitat [kg/m^3]

INTEGER, ALLOCATABLE :: ipiv(:)
INTEGER, PARAMETER :: N=101                         ! Nombre de punts de la idscretització espaial
REAL :: tfin=0                        ! Temps final (normalitzat)
REAL, PARAMETER :: &
    alpha=k/(c_v*p), &                               ! Coeficient de difusió [m^2/s]
    Tc_norm=1, &                                     ! Temperatura corporal normalitzada
    tin=0, &                                         ! Temps inicial (normalitzat)
    xin=0, &                                         ! Coordenada espaial inicial (normalitzada)
    xfin=D*SQRT((0.5*sigma*V**2)/(Tc*k*D**2)), &     ! Coordenada espaial final (normalitzada)
    Tc_desnor=36.5
REAL :: ax, at2
INTEGER :: Nat2
REAL, ALLOCATABLE :: &
    T(:,:), &           ! Matriu de temperatures
    xlin(:), &          ! Discretització espaial
    tlin(:), &          ! Discretització temporal
    M1(:,:), &          ! Matriu 1 Crank_Nicholson
    M1_inv(:,:), &      ! Inversa de la matriu 1 Crank_Nicholson
    M2(:,:), &          ! Matriu 2 Crank_Nicholson
    b(:), &             ! Vector independent
    T_new(:), &         ! Matriu de temperatures temporal
    work(:), &           ! Espai de treball per al mètode LU d'inversió de matrius
    T_desnor(:,:), &
    xlin_desnor(:), &
    tlin_desnor(:)

end module parametros



PROGRAM solucion_IVS
use parametros
REAL, ALLOCATABLE :: &
    temps(:), &            ! Posició temporal
    Temp(:,:), &           ! Temperatura T(i,j): i és la component temporal, j la component espaial
    espai(:)               ! Posició espaial 
INTEGER, dimension(2) :: posmax
REAL :: maximo, T_max=80, aug_temp=0.01, maximo2
LOGICAL :: max_stop = .true., max_stop2=.true.


!Verifiquem que la temperatura no supera 80ºC en cap punt
DO WHILE (max_stop .and. max_stop2)
tfin = tfin + aug_temp
call cranck_nichols() !Obtenim la simulació de la temperatura fins un temps determinat

!Imposam que cap punt de la regió sana arribi a 50ºC
DO i=1,N
    if (xlin_desnor(i) < 0.0075 .or. xlin_desnor(i)>0.0125) THEN
    maximo2 = T_desnor(SIZE(tlin),i)
    if (maximo2 >= 50) THEN
        max_stop2 = (maximo2 < 50)
        !print *, tfin, maximo2 , xlin_desnor(i), max_stop2
    end if
    end if
END DO

!Imposam que la temperatura no superi els 90ºC
posmax = maxloc(T_desnor)
maximo = T_desnor(posmax(1),posmax(2))
max_stop = (maximo < T_max)


if (max_stop .and. max_stop2) THEN !Comprovem que no s'ha superat la temperatura màxima
    Temp = T_desnor
    espai = xlin_desnor
    temps = tlin_desnor
    DEALLOCATE(tlin_desnor)
    DEALLOCATE(xlin_desnor)
    DEALLOCATE(T_desnor)
    !print *, "vuelta con", maximo, "en tiempo", tfin
end if 
END DO
print *, "El temps més eficient trobat és", (tfin-aug_temp)*(2*D**2*c_v*p*Tc)/(sigma*v**2)
!----//IMPRIMIMOS LOS DATOS EN DOCMUENTO DE TEXTO PARA GRAFICAR//----
!Les dades com a llistes de dades que pugui entendre gnuplot o altre programa que les grafiqui:
OPEN(unit=10,file='resultat_problemaIVS-3D.txt',status='replace',action='write')
OPEN(unit=20,file='resultat_problemaIVS-2D.txt',status='replace',action='write')
!Dades gràfic 3D:
DO i=1,N
    DO j=1, SIZE(temps)
        WRITE(10, '(F10.6,F10.6,F10.3)') espai(i), temps(j), Temp(j,i)
    END DO
END DO
!Dades gràfic 2D a y =0.025
DO i=1,N
    WRITE(20,'(F10.6,F10.6)') espai(i), Temp(SIZE(temps),i)
END DO
CLOSE(10)
CLOSE(20)

!------------------------------------------------------------------//
!Subrutines que cridarem al executar el programa
CONTAINS
SUBROUTINE cranck_nichols()
use parametros
ax=(xfin-xin)/N
at2=0.5*(ax)**2  
Nat2 = int((tfin-tin)/at2)

!----///---PROBLEMA IVS---///----:
!Es crea la malla de punts d'espai i temps
ALLOCATE(xlin(N), tlin(Nat2))
xlin = linspace(xin,xfin,N)
tlin = linspace(tin,tfin,Nat2)

!Calculem la matriu M1 i la seva inversa:
ALLOCATE(M1(N-2,N-2),M1_inv(N-2,N-2),ipiv(N-2),work(2*(N-2)))
M1 = 0.0
DO i = 1, N-2
    IF (i > 1) M1(i, i-1) = -1.0     
    M1(i, i) = 4.0                 
    IF (i < N-2) M1(i, i+1) = -1.0   
END DO
CALL invertir(M1, M1_inv, N-2, ipiv, work, info)

!Condició que informa si la matriu construida és invertible o no
IF (info /= 0) THEN
    PRINT *, "La matriu M1 no és invertible, info =", info
END IF
DEALLOCATE(ipiv,work)

!Calculem la matriu M2
ALLOCATE(M2(N-2,N-2))
M2 = 0.0
DO i = 1, N-2
    IF (i > 1) M2(i, i-1) = 1.0                   
    IF (i < N-2) M2(i, i+1) = 1.0   
END DO

!Vector constant b
ALLOCATE(b(N-2))
b = 2*at2
b(1)= 2*(Tc_norm + at2)
b(N-2)= 2*(Tc_norm + at2)

!Estructura del mapa de temperatures T(i,j): i és la component temporal, j la component espaial
ALLOCATE(T(Nat2,N))

!Condicions de contorn
T(:,1)=Tc_norm
T(:,N)=Tc_norm

!Condicions inicials
T(1,:)=Tc_norm

!Iteracions aplicant l'equació matricial
DO i=2,Nat2
    T(i,2:N-1) = matmul(M1_inv,b + matmul(M2,T(i-1,2:N-1)))
END DO


!----//OBTENIM LES DADES REALS REVERTINT LA NORMALITZACIÓ//----
ALLOCATE(T_desnor(Nat2,N), xlin_desnor(N), tlin_desnor(Nat2))
DO i=1,N
    DO j=1,Nat2
        T_desnor(j,i) = T(j,i)*Tc_desnor
        xlin_desnor(i) = xlin(i)/SQRT((0.5*sigma*V**2)/(Tc*k*D**2))
        tlin_desnor(j) = tlin(j)*(2*D**2*c_v*p*Tc)/(sigma*V**2)
    END DO
END DO

DEALLOCATE(xlin,tlin,M1,M1_inv,M2,b,T)

END SUBROUTINE cranck_nichols


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
END PROGRAM solucion_IVS