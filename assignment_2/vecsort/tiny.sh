echo "outer_cpu iner_cpu outer_size iner_max_size time"

 #long and tinny
for cpu in 1 2 3 4
do
    ./vecsort -s 999 -L 1024 -l 1000000 -m 1000  -c $cpu -v `expr "8 - $cpu"` -i 10 -d
done