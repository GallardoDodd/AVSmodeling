program euler_explicito

! Consideramos un disco uniforme de altura h como el tejido enfermo
! c_v es calor específico, p la densidad, k la conductividad térmica
! i dif la difusión térmica
INTEGER, PARAMETER :: c_v=3686, p=1081, h=0.5
REAL, PARAMETER :: k=0.56, phi=0.472, dif=k/c_v/p

!Las condiciones iniciales i de frontera. El sistema se encuentra a la temperatura
!del cuerpo humano Tc inicialmente y en Tc siempre en la frontera.
!N es el nombre de puntos de espacio, ax la amplitud del espaciado en el espacio
!i at el espaciado en el tiempo.
REAL, PARAMETER :: Tc=36.5, N=101, ax=0
REAL :: at1=0.51*(ax)**2, at2=0.49*(ax)**2, at3=0.25*(ax)**2

!modificaciones
END program
