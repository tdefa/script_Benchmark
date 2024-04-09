#!/bin/bash

#SBATCH -N 1

#nombre threads max sur GPU 48

#SBATCH -n 2

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=40000

#SBATCH -J bt_25_inv4

#SBATCH --output="bt_25_inv4.out"
#SBATCH --partition=cbio-cpu

cd /cluster/CBIO/data1/data3/tdefard/st_seg_code/st_seg/scs

python  main_cluster.py \
--do_prerpo_csv_to_tsv 1 \
--do_prerpo 1 \
--do_train 1 \
--do_postprocess 1 \
--main_folder "/cluster/CBIO/data1/st_segmentation/open_vizgen/breast_comseg/scs/${1}" \
--regex_save $2 \
--image_file_name $3 \
--rna_filenames $4 \
--segmentation_labels "cellpose_label" \
--simulation_seg_path $5 \
--cellpose_diameter 90 \
--bin_size  15 \
--epoch 300 \
--background_threshold 10 \
--pixel_size 0.108 \

echo "Done from slurm"


