total_time=2
rm *.txt &>/dev/null
rm particles.log &>/dev/null

for nbparticule in 500 1000 1500 2000 2500;  do
    for nb_thread in 1 2 3 4 5 6; do
        export OMP_NUM_THREADS=$nb_thread
        ../../src/omp_nbody_brute_force $nbparticule $total_time| grep "Simulation took" | awk '{print $3}' >> "omp_nbody_brute_force_$nb_thread.txt"
    done

done