#!/bin/bash -l

module load gcc/12.2.0
make

ns=(8 12 16)
# starts=(535 177)

len=${#ns[@]}
echo $len

for (( i=0; i<$len; i++ ));
do
    n=${ns[$i]}
    # for j in {0..4..1}
    # do
    # tt=$(echo "$j*0.05" | bc -l)
    dir=mult${n}
    # start=${starts[$i]}
    mkdir -p ${dir}
    cp run.sh ${dir}
    cp wolff.o ${dir}
    cd ${dir}
    # qsub -N hex-${n} -v mult=$n -v start=$start -j y run.sh
    qsub -N hex-${n} -v mult=$n -j y run.sh
    echo $n $start
    cd ..
    # done
done
