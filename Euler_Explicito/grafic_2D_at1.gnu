set terminal pngcairo enhanced size 800,600
set output 'T_vs_z_at1.png'
set title "Temperatura vs posició a t=0.0025" font "Times,20"
set xlabel "z"
set ylabel "T"
set style line 1 lc rgb "blue" lw 2
plot "dades_2D_at1.txt" using 1:2 with linespoints ls 1
set output