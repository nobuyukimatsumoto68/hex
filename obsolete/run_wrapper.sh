#!/bin/bash -l

source /projectnb/qfe/nmatsumo/scripts/load_modules_grid.sh

OBJ="./bin/ising_t2"
Lx_list=( 32 48 64 96 128 160 192 )

make
export OBJ

for Lx in "${Lx_list[@]}"
do
    export Lx
    # bash run.sh
    # dir="${Lx}_$((${Lx}*8))_1.000_1.000_0.000"
    # mkdir -p ${dir}
    # cd ${dir}
    qsub -N corr_${Lx} -o "log_${Lx}" ./run.sh
done

