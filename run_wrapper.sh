#!/bin/bash -l

module load gcc/13.2.0

# app='wolff.o'
app='wolff2.o'
make ${app}

# # ns=(8 12)
# ns=(4 6 8 10 12)
# # starts=(535 177)
# len=${#ns[@]}

# for (( i=0; i<$len; i++ ));
# do
#     n=${ns[$i]}
#     # for j in {0..4..1}
#     # do
#     # tt=$(echo "$j*0.05" | bc -l)
#     dir=mult${n}
#     mkdir -p ${dir}

#     # start=$(ls ${dir}/config*tautil0.208378_1.181769/ | cut -f2 -d'/' | sed 's/ckpoint//g' | sed 's/.dat//g')
#     # start=$(ls ${dir}/config*tautil0.500000_0.866025/ | cut -f2 -d'/' | sed 's/ckpoint//g' | sed 's/.dat//g')
#     cp run.sh ${dir}
#     cp ${app} ${dir}
#     cd ${dir}
#     # qsub -N hex-${n} -v mult=$n -v start=$start -j y run.sh
#     qsub -N hex-${n} -v mult=$n -v app=$app -j y run.sh
#     echo $n $start
#     cd ..
#     # done
# done

n=4
dir=mult${n}
mkdir -p ${dir}

cp run.sh ${dir}
cp ${app} ${dir}
cd ${dir}
# qsub -N hex-${n} -v mult=$n -v start=$start -j y run.sh
qsub -N hex-${n} -v mult=$n -v app=$app -j y run.sh
echo $n $start
