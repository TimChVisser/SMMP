echo "min_recusrsion multicore_enabled cores length order time"

# for length in 1024 4096 9000 32768 70000 1572864 3500000

for length in 10000 100000 1000000 10000000 100000000 1000000000
do
    for order in -a -d -r
    do
        ./mergesort -l $length -i 10 -s 999 $order
     done
done


for min_l in 64 256 1024 16384 262144
do
    for length in 10000 100000 1000000 10000000 100000000 1000000000
	do
        for order in -a -d -r
        do
            for cores in -c1 -c2 -c3 -c4 -c5 -c6 -c7 -c8
            do
                 ./mergesort -s 999 -L $min_l -l $length -p $cores -i 10 $order
            done

        done
    done
done


