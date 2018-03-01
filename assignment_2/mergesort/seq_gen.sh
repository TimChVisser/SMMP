for length in 10000 100000 1000000 10000000 100000000 1000000000
do
    for order in -a -d -r
    do
        ./mergesort -l $length -i 10 -s 999 $order
     done
done