
set xlabel "Pos. X"
set ylabel "Pos. t"

# Dibujar líneas verticales en x=2 y x=4
set arrow from 0.0075, 36 to 0.0075,55 nohead linecolor rgb "red"  # Línea vertical en x=2
set arrow from 0.0125, 36 to 0.0125,55 nohead linecolor rgb "red"  # Línea vertical en x=4

# Graficar una función (por ejemplo, y = sin(x))
plot "resultat_problemaIVS-2D.txt" using 1:2

# Guardar la salida en un archivo PNG
set terminal pngcairo
set output 'solucion_IVS-2D.png'
replot
set terminal x11  # Vuelve a la salida estándar
