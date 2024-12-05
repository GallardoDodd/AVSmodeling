program euler_explicito
! Consideramos un disco uniforme de altura h como el tejido enfermo
! c_v es calor específico, p la densidad, k la conductividad térmica,
! dif la difusión térmica i una amplitud del ventriculo ampven.
INTEGER, PARAMETER :: c_v=3686, p=1081, ampven=2
REAL, PARAMETER :: k=0.56, phi=0.472, dif=k/c_v/p,  h=0.5

!Las condiciones iniciales i de frontera. El sistema se encuentra a la temperatura
!del cuerpo humano Tc inicialmente y en Tc siempre en la frontera.
!N es el nombre de puntos de espacio, ax la amplitud del espaciado en el espacio
!i at el espaciado en el tiempo.
INTEGER,  PARAMETER :: N=101
CHARACTER(len=N) :: formato
REAL, PARAMETER :: Tc=36.5, tin=0, tfin=0.025, xin=0, xfin=ampven !Condiciones de contorno
REAL, PARAMETER:: ax=(xfin-xin)/N , at1=0.51*(ax)**2, at2=0.49*(ax)**2, at3=0.25*(ax)**2 !Espaciado temporal segun el espacial en cada caso
INTEGER :: Nat1 =int((tfin-tin)/at1), Nat2 =int((tfin-tin)/at2), Nat3 =int((tfin-tin)/at3) !Número de puntos del espaciado temporal en cada caso
!PREGUNTA!!!!: Aquí calculo el número de puntos que tiene la malla temporal en cada caso, pero me sale un numero real no entero, i lo estoy truncando
! a causa de esto realmente nunca llegamos a 0,025 (tfin), pero si llegamos cerca, lo damos por valido? o añadimos al final de la malla otro valor extra que sea tfin
REAL, allocatable :: T(:,:), xlin(:), tlin(:)
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

!Datos gráfico 2D a y = 0.025
DO i=1,N
    WRITE(20, '(F10.6,F10.6)') xlin(i), T(i,Nat1)
END DO

close(10)
PRINT *,"> 'temperatura_at1' closed correctly"
close(20)
PRINT *,"> 'dades_2D_at1' closed correctly"

tlin=linspace(tin,tfin,Nat2)
!----------------------------------------------------------------------------

!SEGON CAS: at2
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

tlin=linspace(tin,tfin,Nat3)

!----------------------------------------------------------------------------
!TERCER CAS: at3
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

!funciones que luego ya moveremos a un archivo separado
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