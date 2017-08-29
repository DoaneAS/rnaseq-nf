#!/bin/bash -l

spack load pigz

fr1=$1
fr2=$2
pair_id=$3



salmon quant -l A -i index --threads ${NSLOTS} -o $pair_id  -1 <(gzip -cd ${fr1}) -2 <(gzip -cd ${fr2}) --numBootstraps 30

