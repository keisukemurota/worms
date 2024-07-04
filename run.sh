#!/usr/bin/env bash

# Set the temperatures
Bs=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.5 2.0 2.5 3.0 3.5 4.0 5.0 6.0)
N=$(echo "10^6" | bc -l)
K=$(echo "10^4" | bc -l)
L=30

echo "N = $N, K = $K and L = $L"

if [ ! -d "build" ]; then
    echo "Error: build directory does not exist."
    exit 1
else
    cd build
fi

rm -f output.log

for B in ${Bs[@]}; do
    T=$(echo "1/$B" | bc -l)
    out=$(mpirun -n 6 ./main_MPI -m BLBQ1D -L1 $L -T $T -N $N -K $K 2>&1)
    echo "T = $T done."
    c=$(echo "$out" | grep "Specific heat")
    e=$(echo "$out" | grep "Energy per site")
    as=$(echo "$out" | grep "Average sign")
    echo "T = $T, $c, $e, $as" >> output.log
done
