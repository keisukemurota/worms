#!/usr/bin/env bash

# Set the temperatures
Bs=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.3 1.5 1.8 2.0 2.5 3.0 3.5 4.0 4.6 5.0 6.0 7.0)
N=$(echo "10^5" | bc -l)
K=$(echo "10^4" | bc -l)
L=50
p=20

echo "N = $N, K = $K and L = $L"

if [ ! -d "build" ]; then
    echo "Error: build directory does not exist."
    exit 1
else
    cd build
fi

source /opt/materiapps-gcc/env.sh
source /opt/materiapps-gcc/alpscore/alpscorevars.sh

u_path="../python/rmsKit/array/torch/BLBQ1D_loc/J0_1_J1_1_hx_1_hz_0/-1_mel/Adam/lr_0.001_epoch_10000/loss_0.0000002/u"

if [ ! -d "$u_path" ]; then
    echo "Error: u_path directory does not exist."
    exit 1
fi

rm -f output.log

for B in ${Bs[@]}; do
    T=$(echo "1/$B" | bc -l)
    out=$(mpirun -n $p ./main_MPI -m BLBQ1D -L1 $L -T $T -N $N -K $K 2>&1)
    echo "T = $T done."
    c=$(echo "$out" | grep "Specific heat")
    e=$(echo "$out" | grep "Energy per site")
    as=$(echo "$out" | grep "Average sign")
    n=$(echo "$out" | grep "sweeps(in total)")
    echo "T = $T,$n,$c,$e,$as" >> output.log
done
