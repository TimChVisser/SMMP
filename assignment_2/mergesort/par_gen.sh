for length in 10000 100000 1000000 10000000 100000000 1000000000
do
    for order in -a -d -r
    do
        for cores in -c1 -c2 -c3 -c4 -c5 -c6 -c7 -c8
        do
                 ./mergesort -s 999 -L $1 -l $length -p $cores -i 10 $order
        done

    done
done