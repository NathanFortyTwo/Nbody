total_time=2
nbparticule=2500

rm *.txt &>/dev/null
rm *.png &>/dev/null
rm particles.log &>/dev/null

../../hosts.sh

mpi_init=$(mpirun -np 1 --hostfile hosts ../../src/mpi_nbody_brute_force $nbparticule $total_time | grep "Simulation took" | awk '{print $3}') 
echo "1 $mpi_init $mpi_init" >> "mpi_nbody_brute_force.txt"

export OMP_NUM_THREADS=1
omp_init=$(../../src/omp_nbody_brute_force $nbparticule $total_time| grep "Simulation took" | awk '{print $3}') 
echo "1 $omp_init $omp_init" >> "omp_nbody_brute_force.txt"


for nb_thread in 2 3 4 5 6; do
    mpi_res=$(mpirun -np $nb_thread --hostfile hosts ../../src/mpi_nbody_brute_force $nbparticule $total_time | grep "Simulation took" | awk '{print $3}') 
    echo "$nb_thread $mpi_init $mpi_res" >> "mpi_nbody_brute_force.txt"

    export OMP_NUM_THREADS=$nb_thread
    omp_res=$(../../src/omp_nbody_brute_force $nbparticule $total_time| grep "Simulation took" | awk '{print $3}') 
    echo "$nb_thread $omp_init $omp_res" >> "omp_nbody_brute_force.txt"

done

gnuplot OMP_speedup.plot &>/dev/null
gnuplot MPI_speedup.plot &>/dev/null
