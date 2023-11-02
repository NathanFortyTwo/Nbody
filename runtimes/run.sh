total_time=2
rm *.txt
rm particles.log

for nbparticule in 500 1000 1500 2000 2500;  do

    ../parallel/sequential_nbody_brute_force $nbparticule $total_time | grep "Simulation took" | awk '{print $3}' >> "sequential_nbody_brute_force.txt"
    ../parallel/sequential_nbody_barnes_hut $nbparticule $total_time | grep "Simulation took" | awk '{print $3}' >> "sequential_nbody_barnes_hut.txt"

    for nb_thread in 2 4 6; do
        mpirun -np $nb_thread --hostfile hosts ../parallel/mpi_nbody_brute_force $nbparticule $total_time | grep "Simulation took" | awk '{print $3}' >> "mpi_nbody_brute_force_$nb_thread.txt"
        export OMP_NUM_THREADS=$nb_thread
        ../parallel/omp_nbody_brute_force $nbparticule $total_time| grep "Simulation took" | awk '{print $3}' >> "omp_nbody_brute_force_$nb_thread.txt"
    done

done