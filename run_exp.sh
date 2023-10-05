nbparticules=500
simul_time=5
#logdir=../results/$(date)

logdir=./results
algodir=./parallel
for algo in brute_force barnes_hut
do

    filename=$algo.txt

    rm "$logdir/$filename"
    touch "$logdir/$filename"
    
    make &> "/dev/null"

    for nbproc in 1 2 3 4 5 6 7 8
    do
        export OMP_NUM_THREADS=$nbproc
        timeval=$(./$algodir/nbody_$algo $nbparticules $simul_time | tail -n 1 | awk '{print $3}')
        echo "$nbproc:$timeval" >> "$logdir/$filename"
    done

done
exit 0

