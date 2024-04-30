#! /bin/bash

output_dir=/home/joono/CGVM_FDoG_opencv/results/240430
input_dir=/home/joono/CGVM_FDoG_opencv/inputs/simple_data/inputs
infodraw_dir=/home/joono/CGVM_FDoG_opencv/inputs/simple_data/infodraws

for ((i=1; i<=10; i++))
do
    input_file_name="color_$i.png"
    infodraw_file_name="color_$i""_out.png"

    input_path="$input_dir/$input_file_name"
    infodraw_path="$infodraw_dir/$infodraw_file_name"

    ../build/FDoG_app $input_path $infodraw_path $output_dir

    if [ $(($i % 5)) -eq 0 ]; then
        echo "compress .npy files to .npz"
        python npy_to_npz.py --input_dir=../results/fpath_npzs --output_dir=../results/fpath_npzs
    fi

done

