
set xlabel "Pos. X"
set ylabel "Temperatura"

set arrow from 0.0075, 36 to 0.0075,51 nohead linecolor rgb "red"  # Línea vertical en x=2
set arrow from 0.0125, 36 to 0.0125,51 nohead linecolor rgb "red"  # Línea vertical en x=4

plot "resultat_problemaIVS-2D.txt" using 1:2

# Guardar la salida en un archivo PNG
set terminal pngcairo
set output 'solucion_IVS-2D.png'
replot
set terminal x11  # Vuelve a la salida estándar
