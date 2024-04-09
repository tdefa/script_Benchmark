



#%%
import numpy as np
import pandas as pd
import tifffile
from tqdm import tqdm
import scipy
import argparse
import datetime
import random
from pathlib import Path



def get_aggrement_nuclei_cell(masks_nuclei, masks_cell):
    unique_nuclei = np.unique(masks_nuclei)
    list_nuc = []
    list_cell = []
    dict_nuc2cell = {}
    dict_cell2nuc = {}
    for nuc in tqdm(unique_nuclei):
        if nuc == 0:
            continue
        unique_cell_nuc = np.unique(masks_cell[masks_nuclei == nuc])
        if 0 in unique_cell_nuc:
            unique_cell_nuc = unique_cell_nuc[1:] ## if the nuclei must be 100% inside the cell
        print(f'unique cell in the nuc {unique_cell_nuc}')
        if len(unique_cell_nuc) == 1 and unique_cell_nuc[0] != 0:
            unique_nuc_cell = np.unique(masks_nuclei[masks_cell == unique_cell_nuc[0]])
            print(f'unique nuc in the cell {unique_nuc_cell}')
            print()
            if 0 in unique_nuc_cell:
                unique_nuc_cell = unique_nuc_cell[1:]
            if len(unique_nuc_cell) == 1:
                # new_mask[masks_cell == unique_cell_nuc[0]] = nuc
                list_nuc.append(nuc)
                list_cell.append(unique_cell_nuc[0])
                print(f'nuc added {len(list_cell)}')
                dict_nuc2cell[nuc] = unique_cell_nuc[0]
                dict_cell2nuc[unique_cell_nuc[0]] = nuc
    return list_nuc, list_cell, dict_nuc2cell, dict_cell2nuc


#%%

if __name__ == "__main__":


    #size_min_xy
    #size_min_z
    #path_nuclei_mask
    #path_cell_mask


    ############## check consistency of the cellpose segmentation with the prior nuclei

    e = datetime.datetime.now()
    print(e)
    date_str = f"{e.month}_d{e.day}_h{e.hour}_min{e.minute}_s{e.second}_r" + str(random.randint(0, 5000))
    parser = argparse.ArgumentParser(description='test')
    ###################
    ###### hyper-parameter
    ###################

    ## task to do
    parser.add_argument("--compute_filter", type=int,
                        default=0)
    parser.add_argument("--compute_list_cell", type=int,
                        default=0)


    parser.add_argument("--size_min_z_nuclei", type=float,
                        default=0)
    parser.add_argument("--size_min_xy_nuclei", type=int,
                        default=10)

    parser.add_argument("--size_min_z_cellbound", type=float,
                        default=0)
    parser.add_argument("--size_min_xy_cellbound", type=int,
                        default=10)

    parser.add_argument("--z_stack", type=int,
                        default=None)


    parser.add_argument("--path_nuclei_mask", type=str,
                        default="/home/tom/share/baysor_data/input_comseg_3Dseg_tom/prior/prior_stitched.tif")  # TODO
    parser.add_argument("--path_cell_mask", type=str,
                        default="/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/data_release_baysor_merfish_gut/data_analysis/cellpose/cell_boundaries/results/cellpose_membrane.tif")  # TODO
    parser.add_argument("--path_save", type=str,
                        default="/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/data_release_baysor_merfish_gut/prior_stitchedz4_filtered/")

    parser.add_argument("--port", default=3950)
    parser.add_argument("--mode", default='client')
    parser.add_argument("--host", default='127.0.0.2')

    args = parser.parse_args()

    args.path_save = Path(args.path_save)
    args.path_nuclei_mask = Path(args.path_nuclei_mask)
    args.path_cell_mask = Path(args.path_cell_mask)


    Path(args.path_save).mkdir(parents=True, exist_ok=True)
    #### writie the argument to a file

    args.image_name_save = f"{args.size_min_z_nuclei}_{args.size_min_xy_nuclei}_{args.size_min_z_cellbound}_{args.size_min_xy_cellbound}_{args.path_nuclei_mask.stem}"

    with open(Path(args.path_save) / "script_parameter_.txt", "w") as f:
        for k, v in args.__dict__.items():
            f.write(f"{k} : {v}\n")
            print(f"{k} : {v}\n")

    masks_nuclei_o = tifffile.imread(args.path_nuclei_mask)
    masks_cell_o = tifffile.imread(args.path_cell_mask)
    masks_cell = masks_cell_o.copy()
    masks_nuclei = masks_nuclei_o.copy()
    masks_cell = masks_cell.astype(float)


    if args.z_stack is None:
        masks_nuclei = masks_nuclei[args.z_stack]
        masks_cell = masks_cell[args.z_stack]
    new_masks_nuclei = np.zeros_like(masks_nuclei)
    new_masks_cell = np.zeros_like(masks_cell)
    unique_cell = np.unique(masks_cell)
    unique_nuclei = np.unique(masks_nuclei)

    if args.compute_filter:

        for cell in tqdm(unique_cell):
            if cell == 0:
                continue
            mask = masks_cell == cell
            mask = scipy.ndimage.minimum_filter(
                mask,
                size=10
            )
            new_masks_cell[mask] = cell
        for nuclei in tqdm(unique_nuclei):
            if nuclei == 0:
                continue
            mask = masks_nuclei == nuclei
            mask = scipy.ndimage.minimum_filter(
                mask,
                size=args.size_min_xy_nuclei
            )
            new_masks_nuclei[mask] = nuclei

        tifffile.imwrite(args.path_save / f"nuclei_min{args.size_min_xy_nuclei}.tif",
                         new_masks_nuclei)
        tifffile.imwrite(args.path_save / f"cell_min{args.size_min_xy_cellbound}.tif",
                         new_masks_cell)

    if args.compute_list_cell:
        new_masks_nuclei =  tifffile.imread(args.path_save / f"nuclei_min{args.size_min_xy_nuclei}.tif")
        new_masks_cell = tifffile.imread(args.path_save / f"cell_min{args.size_min_xy_cellbound}.tif")



        list_nuc, list_cell, dict_nuc2cell, dict_cell2nuc = get_aggrement_nuclei_cell(masks_nuclei=new_masks_nuclei
                                                                                     ,masks_cell=new_masks_cell
                                                                                      )

        np.save(args.path_save / f"dict_nuc2cell_min{args.size_min_xy_nuclei}_min{args.size_min_xy_cellbound}.npy",
                dict_nuc2cell)
        np.save(args.path_save / f"dict_cell2nuc_min{args.size_min_xy_nuclei}_min{args.size_min_xy_cellbound}.npy",
                dict_cell2nuc)
        np.save(args.path_save / f"list_nuc_min{args.size_min_xy_nuclei}_min{args.size_min_xy_cellbound}.npy",
                list_nuc)
        np.save(args.path_save / f"list_cell_min{args.size_min_xy_nuclei}_min{args.size_min_xy_cellbound}.npy",
                list_cell)


        new_masks_cell[~np.isin(new_masks_cell, list_cell)] = 0
        tifffile.imwrite(args.path_save  / f"cell_min{args.size_min_xy_nuclei}_min{args.size_min_xy_cellbound}.tif",
                         new_masks_cell)







