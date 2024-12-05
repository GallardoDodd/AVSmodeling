program euler_implicito

! Consideramos un disco uniforme de altura h como el tejido enfermo
! c_v es calor específico, p la densidad, k la conductividad térmica
! i dif la difusión térmica
INTEGER, PARAMETER :: c_v=3686, p=1081, h=0.5, ampven=2
REAL, PARAMETER :: k=0.56, phi=0.472, dif=k/c_v/p

!Las condiciones iniciales i de frontera. El sistema se encuentra a la temperatura
!del cuerpo humano Tc inicialmente y en Tc siempre en la frontera.
!N es el nombre de puntos de espacio, ax la amplitud del espaciado en el espacio
!i at el espaciado en el tiempo.
INTEGER, PARAMETER :: N=101
CHARACTER(len=N) :: formato
REAL, PARAMETER :: Tc=36.5, tin=0, tfin=0.025,xin=0,xfin=ampven
REAL, PARAMETER :: ax=(xfin-xin)/N, at1=(ax)**2, at2=0.5*(ax)**2
REAL, PARAMETER :: gamma_1=1.0, gamma_2=0.5
REAL, PARAMETER :: beta_1=1/3, beta_2=0.5
INTEGER :: Nat1 = int((tfin-tin)/at1), Nat2=int((tfin-tin)/at2)
REAL, ALLOCATABLE :: T(:,:), xlin(:), tlin(:)

xlin = linspace(xin,xfin,N)
tlin = linspace(tin,tfin,Nat1)

!Estructura T(i,j): i és la component espaial, j la component temporal
ALLOCATE(T(N,Nat1))
!Condicions de contorn
T(1,:)=Tc
T(N,:)=Tc

!Condicions inicials
T(:,1)=Tc
!
DO i=2,N-1
    DO j=2,Nat1-1
        T(i,j)= beta_1*T(1,j-1) + beta_1*gamma_1*T(i-1,j) + beta_1*gamma_1*T(i+1,j) + beta_1*at1
    END DO
END DO

PRINT*, "> Open dades_3D_imp_at1"
OPEN(unit=10,file='dades_3D_imp_at1.txt',status='replace',action='write')
PRINT*,"> Open dades_2D_imp_at1"
OPEN(unit=20,file='dades_2D_imp_at1.txt',status='replace',action='write')
!Dades gràfic 3D:
DO i=1,N
    DO j=1,Nat1
        WRITE(10, '(F10.6,F10.6,F10.3)') xlin(i), tlin(j), T(i,j)
    END DO
END DO
!Dades gràfic 2D a y =0.025
DO i=1,Nat1
    WRITE(20,'(F10.6,F10.6)') xlin(i), T(i,Nat1)
END DO
CLOSE(10)
PRINT *, "> 'dades_3D_imp_at1' closed correctly"
CLOSE(20)
PRINT *, "> 'dades_2D_imp_at1' closed correctly"

contains
    function linspace(valin, valfin, num) result(lista)
        real, intent(in) :: valin, valfin   ! Valor final i inicial
        integer, intent(in) :: num ! Número de puntos
        real, allocatable :: lista(:)
        integer :: k
        real :: paso
        
        paso = (valfin - valin) / (num - 1)
        allocate(lista(num))
        
        do k = 1, num
            lista(k) = valin + paso*(k - 1)
        end do
    end function linspace
END program euler_implicito