#!/bin/bash

widths=(100 1000 5000 10000)
threads=(4 6 12 20)



m1_files=()
m2_files=()
m2t_files=()
matmul_files=()
matmult_files=()

for file in m_1/*; do
    m1_files+=("test/m_1/$file")
done

for file in m_2/*0.dat; do
    m2_files+=("test/m_2/$file")
done

for file in m_2/*t.dat; do
    m2t_files+=("test/m_2/$file")
done

for file in matmul/*0.dat; do
    matmul_files+=("test/matmul/$file")
done

for file in matmul/*t.dat; do
    matmult_files+=("test/matmul/$file")
done

omp_size="omp_size.dat"

#tesing for different input sizes
base_num_threads=12
num_impl=3
timefile="timing.dat"
impl=1
if [ ! -e output/ ]; then
    mkdir output
fi

# timing study for different input sizes for 12 threads
echo $impl
echo $num_impl
while [ $impl -le $num_impl ]; do
    for ((i=0;i<${#m1_files[@]};i++));do
        m1="${m1_files[$i]}"
        m2="${m2_files[$i]}"
        width="${widths[$i]}"         
        "../bin/omp.out -l $m1 -r $m2 -i $impl -t $base_num_threads -w $width > output/output.dat 2>> $timefile"
    done
    echo -e "\n" >> $timefile
    impl=$(expr $impl + 1)
done


