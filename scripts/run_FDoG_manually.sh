#! /bin/bash

output_dir=/mnt/d/Datasets/hard_fpath_data
npy_dir=/mnt/d/Datasets/hard_fpath_data/fpath_npzs
npz_dir=/mnt/d/Datasets/hard_fpath_data/fpath_npzs

input_dir="/home/joono/CGVM_FDoG_opencv/inputs/origin_res/imgs/"
input_imgs=(
    "beach.bmp"
    "bently.bmp"
    "boat1.bmp"
    "boat2.bmp"
    "boy.bmp"
    "bridge.bmp"
    "cat.bmp"
    "girl_gen.bmp"
    "girl.bmp"
    "lizard.bmp"
    "milk.bmp"
    "simple01.bmp"
    "slrclub.bmp"
    "testpattern.bmp"    
)
infodraw_dir="/home/joono/CGVM_FDoG_opencv/inputs/origin_res/infodraw/"
infodraw_imgs=(
    "beach.bmp"
    "bently_out.bmp"
    "boat1.bmp"
    "boat2.bmp"
    "boy.bmp"
    "bridge.bmp"
    "cat.bmp"
    "girl_gen.bmp"
    "girl.bmp"
    "lizard.bmp"
    "milk.bmp"
    "simple01.bmp"
    "slrclub.bmp"
    "testpattern.bmp"   
)

cd /home/joono/CGVM_FDoG_opencv/scripts

# Ensure that both arrays have the same length
if [ ${#input_imgs[@]} -ne ${#infodraw_imgs[@]} ]; then
    echo "Error: Input and Infodraw arrays must have the same length."
    exit 1
fi

for ((i=0; i<${#input_imgs[@]}; i++))
do
    input_file_name="${input_imgs[i]}"
    infodraw_file_name="${infodraw_imgs[i]}"

    input_path="$input_dir$input_file_name"
    infodraw_path="$infodraw_dir$infodraw_file_name"

    echo ""
    echo "Processing $input_file_name..."

    ../build/FDoG_app "$input_path" "$infodraw_path" "$output_dir"

    python npy_to_npz.py --input_dir="$npy_dir" --output_dir="$npz_dir"
done

# for ((i=1; i<=1000; i++))
# do
#     input_path="$input_dir$input_img"
#     infodraw_path="$infodraw_dir$infodraw_img"

#     echo ""
#     echo "Processing $input_file_name..."

#     ../build/FDoG_app $input_path $infodraw_path $output_dir

#     if [ $(($i % 5)) -eq 0 ]; then
#         python npy_to_npz.py --input_dir=$npy_dir --output_dir=$npz_dir
#     fi

# done

