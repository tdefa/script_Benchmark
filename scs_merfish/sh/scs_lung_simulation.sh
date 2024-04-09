#!/bin/bash

#SBATCH -N 1

#nombre threads max sur GPU 48

#SBATCH -n 2

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40000

#SBATCH -J bt_25_inv4

#SBATCH --output="bt_25_inv4.out"
#SBATCH --partition=cbio-cpu

## restric to node 7
#SBATCH --nodelist=node007



module load cuda/11.3

cd /cluster/CBIO/data1/data3/tdefard/st_seg_code/st_seg/scs

python main_cluster.py \
--do_prerpo_csv_to_tsv 1 \
--do_prerpo 1 \
--do_train 1 \
--do_postprocess 1 \
--main_folder "/cluster/CBIO/data1/data3/tdefard/T7/simulation/scs_input_v2/${1}" \
--regex_save $2 \
--image_file_name $3 \
--rna_filenames $4 \
--segmentation_labels "simulation_labels" \
--simulation_seg_path $5 \
--cellpose_pretrained_model  "/cluster/CBIO/data1/st_segmentation/pretrained_cellpose/CP_baysor" \
--cellpose_diameter 90 \
--bin_size  $6 \
--new_rna_csv_filename "${2}_${6}.tsv" \
--epoch 100 \
--background_threshold 15 \

echo "Done from slurm"






