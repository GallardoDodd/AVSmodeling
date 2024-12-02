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

!PRIMER CAS: at1
allocate(T(N,Nat1))
DO j=1,Nat1
    T(1,j)=Tc
    T(N,j)=Tc
end DO
DO i=2,N
    T(i,1)=Tc
end DO

DO j=1,Nat1-1
    DO i=2,N-1
        T(i,j+1)= T(i,j)+0.51*(T(i+1,j) - 2*T(i,j) + T(i-1,j))+at1
    end DO
end DO

!Guardamos los datos en archivos para luego graficar
open(unit=10, file='temperatura_at1.txt', status='replace', action='write')
!write(formato, '(A, I0, A)') '(', N, 'F10.6)'
!DO j=1,Nat1
!    write(10, formato) T(:,j)
!end DO

DO i=1,N
DO j=1,Nat1
    write(10, '(F10.6, F10.6, F10.6)') xlin(i), tlin(j), T(i,j)
end DO
end Do

close(10)

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