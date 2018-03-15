#!/bin/sh

use_cores=1

run () {
   for data in n100m100.pgm n1000m1000.pgm n5000m5000.pgm n20000m100.pgm n100m20000.pgm
	do
		for cores in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
		do
			m=`echo $data | grep -o [0-9]* | tail -n 1`
			n=`echo $data | grep -o [0-9]* | head -n 1`
			./heat_pth -n $n -m $m -i 100 -r 100 -c ../../$data -t ../../$data -p $cores | tail -n 1
			if [ $use_cores = "0" ]; then
				break
			fi
		done

	done
}

cp ref2.c ref2.c.backup
make clean
make
use_cores=0
 run

use_cores=1
echo "openMP\n"
cat ref2.c.backup | sed 's/\/\/#define OPENMP/#define OPENMP/' > ref2.c
make clean
make
run

echo "PTHREAD\n"
cat ref2.c.backup | sed 's/\/\/#define PTHREAD/#define PTHREAD/' > ref2.c
make clean
make
run

mv ref2.c.backup ref2.c