#!/bin/bash -l

make

# for n in {13..16}
for n in {3..12}
do
    mkdir -p mult$n
    cp run.sh mult$n
    cp wolff.o mult$n
    cd mult$n
    qsub -N hex-$n -v mult=$n -j y run.sh
    cd ..
done
