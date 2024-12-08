program euler_explicito
! Consideramos un disco uniforme de altura h como el tejido enfermo
! c_v es calor específico, p la densidad, k la conductividad térmica,
! dif la difusión térmica i una amplitud del ventriculo ampven.
INTEGER, PARAMETER :: c_v=3686, p=1081, ampven=2
REAL, PARAMETER :: k=0.56, phi=0.472, dif=k/c_v/p,  h=0.5, pi=3.14159265

!Las condiciones iniciales i de frontera. El sistema se encuentra a la temperatura
!del cuerpo humano Tc inicialmente y en Tc siempre en la frontera.
!N es el nombre de puntos de espacio, ax la amplitud del espaciado en el espacio
!i at el espaciado en el tiempo.
INTEGER,  PARAMETER :: N=101
CHARACTER(len=N) :: formato
REAL, PARAMETER :: Tc=1, tin=0, tfin=0.025, xin=0, xfin=ampven !Condiciones de contorno
REAL, PARAMETER:: ax=(xfin-xin)/N , at1=0.51*(ax)**2, at2=0.49*(ax)**2, at3=0.25*(ax)**2 !Espaciado temporal segun el espacial en cada caso
INTEGER :: Nat1 =int((tfin-tin)/at1), Nat2 =int((tfin-tin)/at2), Nat3 =int((tfin-tin)/at3) !Número de puntos del espaciado temporal en cada caso

REAL, allocatable :: T(:,:), xlin(:), tlin(:), T_an(:,:), T_err(:,:)
xlin=linspace(xin,xfin,N)
tlin=linspace(tin,tfin,Nat1)
!----------------------------------------------------------------------------
!PRIMER CAS: at1
allocate(T(N,Nat1))
DO j=1,Nat1
    T(1,j)=Tc
    T(N,j)=Tc
end DO
DO i=1,N
    T(i,1)=Tc
end DO
DO j=1,Nat1-1
    DO i=2,N-1
        T(i,j+1)= T(i,j)+0.51*(T(i+1,j) - 2*T(i,j) + T(i-1,j))+at1
    end DO
end DO

!Guardamos los datos en archivos para luego graficar
PRINT *,"> Open temperatura_at2"
open(unit=10, file='temperatura_at1.txt', status='replace', action='write')
PRINT *,"> Open dades_2D_at1"
open(unit=20,file='dades_2D_at1.txt',status='replace',action='write')
!Datos gráfico 3D
DO i=1,N
DO j=1,Nat1
    write(10, '(F10.6, F10.6, F10.6)') xlin(i), tlin(j), T(i,j)
end DO
end Do

!Datos gráfico 2D a x = 0.025
DO i=1,N
    WRITE(20, '(F10.6,F10.6)') xlin(i), T(i,Nat1)
END DO

!Error numèrico de la solución en x=0.025
allocate(T_an(N,Nat1))
allocate(T_err(N,Nat1))
open(unit=30,file='Grafics_error/error_2D_at1.txt',status='replace',action='write')
DO i=1,N
    T_an(i,Nat1)= Tc + (1-exp(-pi*tfin**2))*sin(pi*xlin(i)) *  4/pi**3 +  (1-exp(-pi*9*tfin**2))*sin(9*pi*xlin(i)) *  4/(27*pi**3) !Solución analítica para la temperatura
end DO
DO j=1,Nat1
    DO i=1,N
        T_err(i,j) = -T_an(i,j) + T(i,j)  !Error numèrico: diferencia entre T() numericamente calculado i T_an analiticamente calculado
    end DO
end DO
DO i=1,N
    WRITE(30, '(F10.6,F10.6)') xlin(i), T_err(i,Nat1) !Datos gráfico error numèrico en y=0.025
END DO
print *, T_an
close(30)
close(10)
PRINT *,"> 'temperatura_at1' closed correctly"
close(20)
PRINT *,"> 'dades_2D_at1' closed correctly"


!----------------------------------------------------------------------------
!SEGON CAS: at2
tlin=linspace(tin,tfin,Nat2)
DEALLOCATE(T)
allocate(T(N,Nat2))
DO j=1,Nat2
    T(1,j)=Tc
    T(N,j)=Tc
end DO
DO i=1,N
    T(i,1)=Tc
end DO

DO j=1,Nat2-1
    DO i=2,N-1
        T(i,j+1)= T(i,j)+0.49*(T(i+1,j) - 2*T(i,j) + T(i-1,j))+at2
    end DO
end DO

!Unidad 10 para guardar los datos 3D, unidad 20 para el plot T vs z a t = 0.025
PRINT *,"> Open temperatura_at2"
open(unit=10, file='temperatura_at2.txt', status='replace', action='write')
PRINT *,"> Open dades_2D_at2"
open(unit=20,file='dades_2D_at2.txt',status='replace',action='write')

!Datos gráfico 3D
DO i=1,N
DO j=1,Nat2
    write(10, '(F10.6, F10.6, F10.6)') xlin(i), tlin(j), T(i,j)
end DO
end DO
!Datos gráfico 2D a y = 0.025
DO i=1,N
    WRITE(20, '(F10.6,F10.6)') xlin(i), T(i,Nat2)
END DO

close(10)
PRINT *,"> 'temperatura_at2' closed correctly"
close(20)
PRINT *,"> 'dades_2D_at2' closed correctly"



!----------------------------------------------------------------------------
!TERCER CAS: at3
tlin=linspace(tin,tfin,Nat3)
DEALLOCATE(T)
allocate(T(N,Nat3))
DO j=1,Nat3
    T(1,j)=Tc
    T(N,j)=Tc
end DO
DO i=1,N
    T(i,1)=Tc
end DO

DO j=1,Nat3-1
    DO i=2,N-1
        T(i,j+1)= T(i,j)+0.25*(T(i+1,j) - 2*T(i,j) + T(i-1,j))+at3
    end DO
end DO

!Guardamos los datos en archivos para luego graficar
PRINT *,"> Open temperatura_at3"
open(unit=10, file='temperatura_at3.txt', status='replace', action='write')
PRINT *,"> Open dades_2D_at3"
open(unit=20,file='dades_2D_at3.txt',status='replace',action='write')
!Datos gráfico 3D
DO i=1,N
DO j=1,Nat3
    write(10, '(F10.6, F10.6, F10.6)') xlin(i), tlin(j), T(i,j)
end DO
end Do

!Datos gráfico 2D a y = 0.025
DO i=1,N
    WRITE(20, '(F10.6,F10.6)') xlin(i), T(i,Nat3)
END DO

close(10)
PRINT *,"> 'temperatura_at3' closed correctly"
close(20)
PRINT *,"> 'dades_2D_at3' closed correctly"




!---- FUNCIONES ----
contains
function linspace(valin, valfin, num) result(lista)
   real, intent(in) :: valin, valfin   ! Valor final i inicial
    integer, intent(in) :: num ! Número de puntos
    real, allocatable :: lista(:)
    integer :: i
    real :: paso
    
    paso = (valfin - valin) / (num - 1)
    allocate(lista(num))
      
    do i = 1, num
       lista(i) = valin + paso*(i - 1)
    end do
end function linspace


    
END program