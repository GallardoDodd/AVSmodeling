set terminal pngcairo enhanced size 800,600
set output 'T_vs_z_at2.png'
set title "Temperatura vs posició a t=0.0025"
set xlabel "z"
set ylabel "T"
set style line 1 lc rgb "blue" lw 2
plot "dades_2D_at2.txt" using 1:2 with linespoints ls 1
set output