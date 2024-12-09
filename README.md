Modelización del tratamiento de ablación cardíaca IVS

Este proyecto consistirá en desarrollar un modelo sencillo para describir de manera aproximada el proceso físico que ocurre durante un tratamiento de ablación cardíaca IVS. 
Utilizaremos los métodos numéricos EULER EXPLÍCITO, EULER IMPLÍCITO y el esquema CRANK-NICOLSON. Los compararemos entre si y respecto a la solución analítica del sistema planteado, 
para elegir el método con el que se obtiene mejor solución para nuestro modelo.

Para compilar los códigos de Fortran de cada método se debe usar:

- gfortran euler_explicito.f90 -o euler_explicito
Para Euler_Explicito
- gfortran euler_implicito.f90 -o euler_implicito
Para Euler_Implicito
- gfortran Crank_Nicholson.f90 -o crank_nicholson -llapack -lblas
Para Crank_Nicholson
- gfortran solucion_IVS.f90 -o solucion_IVS-compilado -llapack -lblas
Para Solucion_IVS
