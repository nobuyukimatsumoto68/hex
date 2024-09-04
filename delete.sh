#!/bin/bash

jobs=$(qstat -u nmatsum | grep nmatsum | cut -d ' ' -f 1)

echo $jobs

for elem in "${jobs[@]}"
do
    qdel $elem
done
