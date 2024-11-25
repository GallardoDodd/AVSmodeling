!Los modulos por ahora los dejo aquí, luego los pasamos a otro archivo i organizamos. att: Daniel
module linspace_modulo
    implicit none
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
  end module linspace_modulo !!esto no se usa ahora mismo, se usa la funcion linspace de abajo
  