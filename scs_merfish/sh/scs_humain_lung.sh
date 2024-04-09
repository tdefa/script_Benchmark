#!/bin/bash


#SBATCH -N 1

#nombre threads max sur GPU 48

#SBATCH -n 1

#SBATCH --cpus-per-task=1

#SBATCH -J cp_15

#SBATCH --output="cp_15.out"
#SBATCH --partition=cbio-cpu

#SBATCH --mem-per-cpu=50000

# ## SBATCH --gres=gpu:1

## specify node name





module load cuda/11.3


cd /cluster/CBIO/data1/data3/tdefard/st_seg_code/st_seg/scs

echo "Running scs on $1"
echo "Image file name: $3"
echo "RNA file names: $4"

python main_cluster.py \
--do_prerpo_csv_to_tsv 1 \
--do_prerpo 1 \
--do_train 1 \
--do_postprocess 1 \
--cellpose_model nuclei \
--cellpose_diameter 15 \
--segmentation_labels "cellpose_labels" \
--bin_size  15 \
--main_folder "/cluster/CBIO/data1/data3/tdefard/T7/sp_data/In_situ_Sequencing_16/scs_input_10_10/${1}" \
--regex_save $2 \
--image_file_name $3 \
--rna_filenames $4 \
--epoch 100 \
--background_threshold 10 \




echo "Done from slurm"




