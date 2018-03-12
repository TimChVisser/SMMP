#!/bin/sh
# 3 4 5 6 7 8 9 10 11 12 13 14 15 16
run () {
   for cores in 1 2
	do
		for data in n100m100.pgm n1000m1000.pgm n5000m5000.pgm n20000m100.pgm n100m20000.pgm
		do
			m=`echo $data | grep -o [0-9]* | tail -n 1`
			n=`echo $data | grep -o [0-9]* | head -n 1`
			./heat_pth -n $n -m $m -i 100 -r 100 -c ../../$data -t ../../$data -p $cores | tail -n 1
		done
	done
}

make -q
run

echo "openMP\n"
cat main.c | sed -n 's/\/\/#define OPENMP/#define OPENMP/' > main.c
make -q
run

echo "PTHREAD\n"
cat main.c | sed -n 's/#define OPENMP/\/\/#define OPENMP/' > main.c
cat main.c | sed -n 's/\/\/#define PTHREAD/#define PTHREAD/' > main.c
make -q
run

cat main.c | sed -n 's/#define PTHREAD/\/\/#define PTHREAD/' > main.c
