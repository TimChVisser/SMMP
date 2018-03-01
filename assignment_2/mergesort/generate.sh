echo "min_recusrsion multicore_enabled cores length order time"

sh seq_gen.sh

for min_l in 64 256 1024 16384 262144
do
    sh par_gen.sh $min_l
done


