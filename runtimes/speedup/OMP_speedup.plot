set terminal png
set ylabel "speedup"
set xlabel "number of threads"
set xtics 1
set title "Speedup for OpenMP execution"
set output "OpenMP_STATIC.png"

plot x title 'Speedup max' with lines, \
     'omp_nbody_brute_force.txt' using 1:($2/$3) title 'OpenMP' with linespoints;