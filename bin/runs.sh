#!/bin/bash

#config=swan-queue
#bench=ferret

bench=$1
config=$2
f=$3
t=$4

#./parsecmgmt -a build -p $bench -c $config

#0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11,15,18,22,26,30,19,23,27,31

#CPUS[0]=
CPUS[1]='0'
CPUS[2]='0,4'
CPUS[3]='0,4,8'
CPUS[4]='0,4,8,12'
CPUS[5]='0,4,8,12,1'
CPUS[6]='0,4,8,12,1,5'
CPUS[7]='0,4,8,12,1,5,9'
CPUS[8]='0,4,8,12,1,5,9,13'
CPUS[9]='0,4,8,12,1,5,9,13,16'
CPUS[10]='0,4,8,12,1,5,9,13,16,20'
CPUS[11]='0,4,8,12,1,5,9,13,16,20,24'
CPUS[12]='0,4,8,12,1,5,9,13,16,20,24,28'
CPUS[13]='0,4,8,12,1,5,9,13,16,20,24,28,17'
CPUS[14]='0,4,8,12,1,5,9,13,16,20,24,28,17,21'
CPUS[15]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25'
CPUS[16]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29'
CPUS[17]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2'
CPUS[18]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6'
CPUS[19]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10'
CPUS[20]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14'
CPUS[21]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3'
CPUS[22]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7'
CPUS[23]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11'
CPUS[24]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11,15'
CPUS[25]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11,15,18'
CPUS[26]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11,15,18,22'
CPUS[27]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11,15,18,22,26'
CPUS[28]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11,15,18,22,26,30'
CPUS[29]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11,15,18,22,26,30,19'
CPUS[30]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11,15,18,22,26,30,19,23'
CPUS[31]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11,15,18,22,26,30,19,23,27'
CPUS[32]='0,4,8,12,1,5,9,13,16,20,24,28,17,21,25,29,2,6,10,14,3,7,11,15,18,22,26,30,19,23,27,31'

mkdir -p output_${bench}_${config}

for cores in 1
do
for i in `seq $f $t` ; do
export NUM_THREADS=$cores
echo $NUM_THREADS $bench $i
./parsecmgmt -a run -p $bench -c $config -i native -s "time numactl --physcpubind=+${CPUS[$cores]} --membind=0-$[(NUM_THREADS+7)/8-1]" > output_${bench}_${config}/${bench}_${config}_${cores}_${i}.txt
done
done
