set terminal png
set ylabel "speedup"
set xlabel "number of processes"
set xtics 1
set title "Speedup for MPI execution"
set output "MPI.png"

plot x title 'Speedup max' with lines, \
     'mpi_nbody_brute_force.txt' using 1:($2/$3) title 'MPI' with linespoints;