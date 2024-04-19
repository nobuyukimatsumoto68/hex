#!/bin/bash -l

# https://www.bu.edu/tech/support/research/system-usage/running-jobs/submitting-jobs/
# to submit: qsub script.sh

#------ qsub options --------#
#$ -P qfe
#$ -M mtsmtnbyk@gmail.com
##### export environmental variables
#$ -V
##### run time limit. format: hh:mm:ss; default 12 hrs
#$ -l h_rt=24:00:00
##### merge error and output
#$ -j y
##### email options; begins (b), ends (e), is aborted (a), is suspended (s), or never (n) - default
#$ -m beas


# --------- job info -----------#

echo "start"
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "NSlots : $NSLOTS"
echo "Host name : $HOSTNAME"
echo "working directory : $TMPDIR"
echo "=========================================================="

#-------------#
# echo "load modules"
# source /projectnb/qfe/nmatsumo/scripts/load_modules_grid.sh
echo "Lx=${Lx}"
# export ...
#-------------#

#------- Program execution -------#

echo "running program"
pwd
./$OBJ ${Lx}
echo "finished"
