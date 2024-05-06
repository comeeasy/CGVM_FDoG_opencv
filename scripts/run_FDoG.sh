#! /bin/bash

output_dir=/mnt/d/Datasets/simple_fpath_data
input_dir=/home/joono/CGVM_FDoG_opencv/inputs/simple_data/inputs
infodraw_dir=/home/joono/CGVM_FDoG_opencv/inputs/simple_data/infodraws

npy_dir=/mnt/d/Datasets/simple_fpath_data/fpath_npzs
npz_dir=/mnt/d/Datasets/simple_fpath_data/fpath_npzs

cd /home/joono/CGVM_FDoG_opencv/scripts

for ((i=1; i<=1000; i++))
do
    input_file_name="color_$i.png"
    infodraw_file_name="color_$i""_out.png"

    echo ""
    echo "Processing $input_file_name..."

    input_path="$input_dir/$input_file_name"
    infodraw_path="$infodraw_dir/$infodraw_file_name"

    ../build/FDoG_app $input_path $infodraw_path $output_dir

    if [ $(($i % 5)) -eq 0 ]; then
        python npy_to_npz.py --input_dir=$npy_dir --output_dir=$npz_dir
    fi

done

