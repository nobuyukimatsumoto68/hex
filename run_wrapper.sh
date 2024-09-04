#!/bin/bash -l

module load gcc/12.2.0
make

# for n in {13..16}
ns=(8 12)
# starts=(535 177)

len=${#ns[@]}

# for n in {4..8}

for (( i=0; i<$len; i++ ));
do
    n=${ns[$i]}
    dir=mult_n_$n
    mkdir -p ${dir}
    cp run.sh ${dir}
    cp wolff.o ${dir}
    cd ${dir}
    qsub -N hex_old-$n -v mult=$n -j y run.sh
    cd ..
done
