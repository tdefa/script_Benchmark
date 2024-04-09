#!/bin/bash

#SBATCH -N 1

#nombre threads max sur GPU 48

#SBATCH -n 2

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=80000

#SBATCH -J entire_image

#SBATCH --output="entire_image_gpu.out"
#  SBATCH --partition=cbio-cpu

#SBATCH --partition=cbio-gpu

#SBATCH --gres=gpu:1





module load cuda/11.3

cd /cluster/CBIO/data1/data3/tdefard/st_seg_code/st_seg/scs

python main_cluster.py \
--do_prerpo_csv_to_tsv 1 \
--do_prerpo 1 \
--do_train 1 \
--do_postprocess 1 \
--main_folder "/cluster/CBIO/data1/st_segmentation/baysor_data/scs/entire_image/" \
--regex_save $1 \
--image_file_name $2 \
--rna_filenames $3 \
--segmentation_labels "cellpose_labels" \
--cellpose_pretrained_model  "/cluster/CBIO/data1/st_segmentation/pretrained_cellpose/CP_baysor" \
--cellpose_diameter 90 \
--bin_size  15 \
--epoch 300 \
--background_threshold 10 \

echo "Done from slurm"





