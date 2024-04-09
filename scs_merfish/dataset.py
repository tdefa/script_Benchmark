import os
import math
import numpy as np
import spateo as st
from pathlib import Path
import matplotlib.pyplot as plt
import tifffile
from tqdm import tqdm

from segmentation import compute_cellpose_segmentation
from dataset_utils import find_all_genes, compute_segmentation2center, compute_idx2exp, compute_all_exp_merged_bins, compute_offsets, get_adatasub
from plots import plot_alignment, plot_segmentation
import tifffile

def create_train_test_data(
        bin_file, 
        image_file, 
        save_folder,
        segmentation_labels,
        prealigned = True,
        align = None,
        patchsize = 0,
        n_neighbor = 50,
        bin_size=3 ,
        startx = 0,
        starty = 0,
        background_threshold = 10,
        cellpose_params = {},
        n_top_genes=2000,
        simulation_seg_path = None,
    ):

    fig_folder = Path(save_folder) / "fig"
    fig_folder.mkdir(parents=True, exist_ok=True)
    data_folder = Path(save_folder) / "data"
    data_folder.mkdir(parents=True, exist_ok=True)
    print(f"bin_file: {bin_file}, image_file: {image_file},")

    adatasub = get_adatasub(bin_file, image_file, prealigned,
                            startx=startx,
                            starty=starty,
                            patchsize=patchsize
                            )

    patch_ids = f'{startx}:{starty}:{patchsize}:{patchsize}'
    patchsizex = int(adatasub.X.shape[0])
    patchsizey = int(adatasub.X.shape[1])

    # Align staining image with bins
    before = adatasub.layers['stain'].copy()
    if align:
        st.cs.refine_alignment(adatasub, mode=align, transform_layers=['stain'])
    plot_alignment(fig_folder, adatasub, before, patch_ids)

    # Create segmentation

    ## cellpose parameters
    if segmentation_labels == 'cellpose_labels':
        gpu = cellpose_params["gpu"] if "gpu" in cellpose_params else True
        model_type = cellpose_params["model_type"] if "model_type" in cellpose_params else 'cyto'
        pretrained_model = cellpose_params["pretrained_model"] if "pretrained_model" in cellpose_params else None
        diameter = cellpose_params["diameter"] if "diameter" in cellpose_params else 90
        flow_threshold = cellpose_params["flow_threshold"] if "flow_threshold" in cellpose_params else 0.9
        adatasub.layers[segmentation_labels] = compute_cellpose_segmentation(adatasub, gpu=gpu, model_type=model_type, pretrained_model=pretrained_model,
                                      diameter=diameter,  flow_threshold=flow_threshold)

    elif segmentation_labels == 'watershed_labels':
        st.cs.mask_nuclei_from_stain(adatasub, otsu_classes=4, otsu_index=1)
        st.cs.find_peaks_from_mask(adatasub, 'stain', 7)
        st.cs.watershed(adatasub, 'stain', 5, out_layer=segmentation_labels)
    elif segmentation_labels == 'simulation_labels':
        mask_seg = tifffile.imread(simulation_seg_path)
        if mask_seg.ndim == 3:
            mask_seg = np.amax(mask_seg, 0)
        adatasub.layers[segmentation_labels] = mask_seg


    else:
        raise ValueError(f"segmentation_labels {segmentation_labels} not recognized")
    plot_segmentation(fig_folder, adatasub, segmentation_labels, patch_ids)

    adatasub.write(data_folder / f'spots{patch_ids}.h5ad')

    # Prepare data for network
    print('Prepare data for neural network...')
    segmentation2center = compute_segmentation2center(adatasub.layers[segmentation_labels])
    geneid, id2gene, genecnt = find_all_genes(bin_file)
    idx2exp = compute_idx2exp(bin_file, geneid, bin_size, startx, starty, patchsizex, patchsizey)
    all_exp_merged_bins, selected_index = compute_all_exp_merged_bins(
        genecnt, idx2exp, bin_size, patchsizex, patchsizey,
        n_top_genes=n_top_genes,
    )
    with open(data_folder / f'variable_genes_{patch_ids}.txt', 'w') as fw:
        for id in selected_index:
            fw.write(id2gene[id] + '\n')

    x_train_tmp = []
    x_train = []
    x_train_pos = []
    y_train = []
    y_binary_train = []
    x_train_bg_tmp = []
    x_train_bg = []
    x_train_pos_bg = []
    y_train_bg = []
    y_binary_train_bg = []
    x_test_tmp = []
    x_test= []
    x_test_pos = []

    offsets = compute_offsets(bin_size)
    for i in tqdm(range(adatasub.layers[segmentation_labels].shape[0])):
        for j in range(adatasub.layers[segmentation_labels].shape[1]):

            if (not i % bin_size == 0) or (not j % bin_size == 0):
                continue
            idx = int(math.floor(i / bin_size) * math.ceil(patchsizey / bin_size) + math.floor(j / bin_size))
            if adatasub.layers[segmentation_labels][i, j] > 0:
                if idx >= 0 and idx < all_exp_merged_bins.shape[0] and np.sum(all_exp_merged_bins[idx, :]) > 0:
                    x_train_sample = [all_exp_merged_bins[idx, :]]
                    x_train_pos_sample = [[i, j]]
                    y_train_sample = [segmentation2center[adatasub.layers[segmentation_labels][i, j]]] # cell id
                    for dx, dy in offsets:
                        if len(x_train_sample) == n_neighbor:
                            break
                        x = i + dx
                        y = j + dy
                        if x < 0 or x >= adatasub.layers[segmentation_labels].shape[0] or y < 0 or y >= adatasub.layers[segmentation_labels].shape[1]:
                            continue
                        idx_nb = int(math.floor(x / bin_size) * math.ceil(patchsizey / bin_size) + math.floor(y / bin_size))
                        if idx_nb >= 0 and idx_nb < all_exp_merged_bins.shape[0] and np.sum(all_exp_merged_bins[idx_nb, :]) > 0:
                            x_train_sample.append(all_exp_merged_bins[idx_nb, :])
                            x_train_pos_sample.append([x, y])
                    if len(x_train_sample) < n_neighbor:
                        continue
                    x_train_tmp.append(x_train_sample)
                    if len(x_train_tmp) > 500:
                        x_train.extend(x_train_tmp)
                        x_train_tmp = []

                    x_train_pos.append(x_train_pos_sample)
                    y_train.append(y_train_sample)
                    y_binary_train.append(1)
            else:
                if idx >= 0 and idx < all_exp_merged_bins.shape[0] and np.sum(all_exp_merged_bins[idx, :]) > 0:
                    backgroud = True
                    for nucleus in segmentation2center:
                        if (i - segmentation2center[nucleus][0]) ** 2 + (j - segmentation2center[nucleus][1]) ** 2 <= 900 or adatasub.layers['stain'][i, j] > background_threshold:
                            backgroud = False
                            break
                    if backgroud:
                        if len(x_train_bg) + len(x_train_bg_tmp) >= len(x_train) + len(x_train_tmp):
                            continue
                        x_train_sample = [all_exp_merged_bins[idx, :]]
                        x_train_pos_sample = [[i, j]]
                        y_train_sample = [[-1, -1]]
                        for dx, dy in offsets:
                            if len(x_train_sample) == n_neighbor:
                                break
                            x = i + dx
                            y = j + dy
                            if x < 0 or x >= adatasub.layers[segmentation_labels].shape[0] or y < 0 or y >= adatasub.layers[segmentation_labels].shape[1]:
                                continue
                            idx_nb = int(math.floor(x / bin_size) * math.ceil(patchsizey / bin_size) + math.floor(y / bin_size))
                            if idx_nb >= 0 and idx_nb < all_exp_merged_bins.shape[0] and np.sum(all_exp_merged_bins[idx_nb, :]) > 0:
                                x_train_sample.append(all_exp_merged_bins[idx_nb, :])
                                x_train_pos_sample.append([x, y])
                        if len(x_train_sample) < n_neighbor:
                            continue
                        x_train_bg_tmp.append(x_train_sample)
                        if len(x_train_bg_tmp) > 500:
                            x_train_bg.extend(x_train_bg_tmp)
                            x_train_bg_tmp = []
                        x_train_pos_bg.append(x_train_pos_sample)
                        y_train_bg.append(y_train_sample)
                        y_binary_train_bg.append(0)
                    else:
                        x_test_sample = [all_exp_merged_bins[idx, :]]
                        x_test_pos_sample = [[i, j]]
                        for dx, dy in offsets:
                            if len(x_test_sample) == n_neighbor:
                                break
                            x = i + dx
                            y = j + dy
                            exp_merge = np.zeros(len(selected_index))
                            if x < 0 or x >= adatasub.layers[segmentation_labels].shape[0] or y < 0 or y >= adatasub.layers[segmentation_labels].shape[1]:
                                continue
                            idx_nb = int(math.floor(x / bin_size) * math.ceil(patchsizey / bin_size) + math.floor(y / bin_size))
                            if idx_nb >= 0 and idx_nb < all_exp_merged_bins.shape[0] and np.sum(all_exp_merged_bins[idx_nb, :]) > 0:
                                x_test_sample.append(all_exp_merged_bins[idx_nb, :])
                                x_test_pos_sample.append([x, y])
                        if len(x_test_sample) < n_neighbor:
                            continue
                        x_test_tmp.append(x_test_sample)
                        if len(x_test_tmp) > 500:
                            x_test.extend(x_test_tmp)
                            x_test_tmp = []
                        x_test_pos.append(x_test_pos_sample)#


    x_train.extend(x_train_tmp)
    x_train_bg.extend(x_train_bg_tmp)
    x_test.extend(x_test_tmp)

    x_train = np.array(x_train)
    x_train_pos = np.array(x_train_pos)
    y_train = np.vstack(y_train)
    y_binary_train = np.array(y_binary_train)
    x_train_bg = np.array(x_train_bg)
    x_train_pos_bg = np.array(x_train_pos_bg)
    y_train_bg = np.vstack(y_train_bg)
    y_binary_train_bg = np.array(y_binary_train_bg)

    bg_index = np.arange(len(x_train_bg))
    np.random.shuffle(bg_index)
    x_train = np.vstack((x_train, x_train_bg[bg_index[:len(x_train)]]))
    x_train_pos = np.vstack((x_train_pos, x_train_pos_bg[bg_index[:len(x_train_pos)]]))
    y_train = np.vstack((y_train, y_train_bg[bg_index[:len(y_train)]]))
    y_binary_train = np.hstack((y_binary_train, y_binary_train_bg[bg_index[:len(y_binary_train)]]))
    x_test= np.array(x_test)
    x_test_pos = np.array(x_test_pos)

    print(x_train.shape, x_train_pos.shape, y_train.shape, y_binary_train.shape, x_train_bg.shape, x_train_pos_bg.shape, y_train_bg.shape, y_binary_train_bg.shape)
    print(x_test.shape, x_test_pos.shape)

    np.savez_compressed(data_folder / f'x_train_{patch_ids}.npz', x_train=x_train)
    np.savez_compressed(data_folder / f'x_train_pos_{patch_ids}.npz', x_train_pos=x_train_pos)
    np.savez_compressed(data_folder / f'y_train_{patch_ids}.npz', y_train=y_train)
    np.savez_compressed(data_folder / f'y_binary_train_{patch_ids}.npz', y_binary_train=y_binary_train)
    np.savez_compressed(data_folder / f'x_test_{patch_ids}.npz', x_test=x_test)
    np.savez_compressed(data_folder / f'x_test_pos_{patch_ids}.npz', x_test_pos=x_test_pos)

if __name__=="__main__":

    from plots import plot_rna_dist
    from dataset_utils import create_csv_file, merge_rna_count_dfs

    rna_filenames = ["all_z.csv"]
    image_file_name = "mosaic_DAPI_z3.tif"
    new_rna_csv_filename =  "preprocessed_rna_no_spot.tsv"
    main_folder = Path("/cluster/CBIO/data1/ablondel1/vallot/202304141840_1184_VallotAB5301andAG5988_VMSC09402/region_0/crop10000_2_3")
    spot_size = None
    pixel_size = 0.108
    image_size=(1000, 1000)

    image_file = main_folder / "images" / image_file_name
    save_folder = main_folder / "scs_pred"
    segmentation_labels = 'cellpose_labels'
    prealigned = True
    align = None
    patchsize = 0
    n_neighbor = 50
    bin_size=3 
    startx = 0
    starty = 0
    background_threshold = 10
    n_top_genes = 2000

    rna_folder = main_folder / "csv"
    bin_file = rna_folder / new_rna_csv_filename

    rna_dfs = []
    for rna_filename in rna_filenames:
        rna_filepath = os.path.join(rna_folder, rna_filename)
        rna_count_df = create_csv_file(
            rna_filepath, 
            spot_size=spot_size, 
            pixel_size=pixel_size, 
            image_size=image_size,
        )
        rna_dfs.append(rna_count_df)
        
    result_df = merge_rna_count_dfs(rna_dfs)
    result_df.to_csv(bin_file, sep='\t', index=False)
    plot_rna_dist(main_folder, result_df)

    create_train_test_data(
        bin_file, 
        image_file, 
        save_folder,
        segmentation_labels,
        prealigned = prealigned,
        align = align,
        patchsize = patchsize,
        n_neighbor = n_neighbor,
        bin_size = bin_size,
        startx = startx,
        starty = startx,
        background_threshold = background_threshold,
        n_top_genes = n_top_genes,
        )
