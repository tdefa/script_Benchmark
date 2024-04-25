


#%%

import numpy as np
import pandas as pd
import tifffile
from tqdm import tqdm
import scipy
from pathlib import Path

def compute_AJI(image_gr_cm,
                image_pred_cm,
                max_index = None,
                cell_index_backgroud =  0 ): ## to move to an other file

    from tqdm import tqdm

    """
    https://ieeexplore.ieee.org/document/7872382
    Parameters
    ----------
    dico_pred_cm :  { cell_id :list of index of mol predicted as belonging to the cell}}
    image_pred_cm :  { cell_id :list of index of mol belonging to the cell}}
    max_index : do not take molecule after a certain indice, to exclude centroid for example
    Returns
    -------    
    """

    dico_cell_pred_rna = {}
    dico_m2b_index_match = {}
    list_used_pred = []
    max_iou_list = []
    C_image = 0
    U_image = 0

    for cell_gr in tqdm(image_gr_cm):  # forEach ground truth nucleus  Gi do :j←arg maxk(|Gi∩Sk|/|Gi∪Sk|)
        mol_gr = image_gr_cm[cell_gr]
        if max_index is not None: ## remove centroid
            mol_gr = np.array(mol_gr)
            mol_gr = list(mol_gr[mol_gr < max_index])
        iou_list = []
        C_list = []
        U_list = []
        list_cell_pred = list(image_pred_cm.keys())
        for cell_pred_index in range(len(list_cell_pred)):
            if cell_index_backgroud == list_cell_pred[cell_pred_index]:
                iou_list.append(0)
                C_list.append(0)
                U_list.append(0)
                continue
            mol_pred = list(image_pred_cm[list_cell_pred[cell_pred_index]])
            if max_index is not None:  ## remove centroid
                mol_pred = np.array(mol_pred)
                mol_pred = list(mol_pred[mol_pred < max_index])
            C_list.append(len(set(list(mol_pred)).intersection(set(list(mol_gr)))))
            U_list.append(len(set(list(mol_pred)).union(set(list(mol_gr))))
                          )
            if len(set(list(mol_pred)).union(set(list(mol_gr)))) > 0 and len(mol_gr) > 0:

                iou = len(set(list(mol_pred)).intersection(set(list(mol_gr)))) / len(
                    set(list(mol_pred)).union(set(list(mol_gr))))
            else:
                iou = -1

            # iou not defined
            # print(iou)

            iou_list.append(iou)
        if np.max(iou_list) == 0:
            index_cell_max_overlap = np.argmax(iou_list)
            max_iou_list.append(0)
            dico_cell_pred_rna[cell_gr] = []

        elif np.max(iou_list) == -1:
            dico_cell_pred_rna[cell_gr] = []

        else:
            index_cell_max_overlap = np.argmax(iou_list)
            max_iou_list.append(iou_list[index_cell_max_overlap])
            C_image += C_list[index_cell_max_overlap]
            U_image += U_list[index_cell_max_overlap]
            list_used_pred.append(list_cell_pred[index_cell_max_overlap])
            dico_cell_pred_rna[cell_gr] = image_pred_cm[list_cell_pred[index_cell_max_overlap]]
            dico_m2b_index_match[cell_gr] = list_cell_pred[index_cell_max_overlap]

    list_unused_list = list(set(list_cell_pred) - set(list_used_pred))

    for cell_unused in list_unused_list:
        U_image += len(image_pred_cm[cell_unused])

    #list_aji_global.append(C_image / U_image)
    return C_image / U_image, C_image, U_image, dico_cell_pred_rna, max_iou_list, dico_m2b_index_match




def compute_perf(pred_df,
                 dict_gr_cm,
        cell_column_name = "cell",
        max_index = None,
                 ):

    if max_index is None:
        max_index = len(pred_df)

    dict_pred_cm = {}
    for i, row in tqdm(pred_df.iterrows()):
        cell = row[cell_column_name]
        if cell in dict_pred_cm:
            dict_pred_cm[cell].append(i)
        else:
            dict_pred_cm[cell] = [i]
        if i > max_index:
            break




    aji, C_image, U_image, dico_cell_pred_rna, max_iou_list, dico_m2b_index_match = compute_AJI(
                image_gr_cm= dict_gr_cm,
                image_pred_cm= dict_pred_cm,
                max_index = None,
                cell_index_backgroud =  0
    )
    return aji, C_image, U_image, dico_cell_pred_rna, max_iou_list, dico_m2b_index_match



#%%




if __name__ == "__local__" :

    str_save = "190324"

    ##########

    cell_seg  = tifffile.imread("/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/data_release_baysor_merfish_gut/mask_for_metrics/cellpose_membrane_min.tif")
    list_unique_cell = list(np.unique(cell_seg))
    df_input = pd.read_csv("/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/input_comseg/df/mouse1.csv")

    cell_from_seg = []
    for row in tqdm(df_input.iterrows()):
        if cell_seg.ndim == 3:
            x = row[1].x
            y = row[1].y
            z = row[1].z
            cell_from_seg.append(cell_seg[z, y, x])
        else:
            raise KeyError("do not support 2D image for validation, should we ?")

    dict_gr_cm = {}
    for i in tqdm(list(range(len(cell_from_seg)))):
        if cell_from_seg[i] in dict_gr_cm:
            dict_gr_cm[cell_from_seg[i]].append(i)
        else:
            dict_gr_cm[cell_from_seg[i]] = [i]
    del dict_gr_cm[0]

    dict_res = {}



    #### comseg test
    df_input = pd.read_csv("/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/input_comseg/df/mouse1.csv")

    pred_df = pd.read_csv("/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/input_comseg/df/results/2_d18_h19_min41_s28_r3259/cell0_r10_rmax8_small_p.csv")
    aji, C_image, U_image, dico_cell_pred_rna, max_iou_list, dico_m2b_index_match = compute_perf(
        pred_df = pred_df,
        dict_gr_cm=dict_gr_cm,
                 cell_column_name="cell",
        max_index=len(df_input),

    )

    list_cell_cyto = list(dico_m2b_index_match.keys())
    list_cell_comseg  = list(dico_m2b_index_match.values())
    list_unique_cell_comseg = list(pred_df.cell.unique())
    list_unique_cyto = list(np.unique(cell_seg))
    assert np.isin(list_cell_cyto, list_unique_cyto).all()
    assert np.isin(list_cell_comseg, list_unique_cell_comseg).all()

    inverse_dico_m2b_index_match = {}
    for k, v in dico_m2b_index_match.items():
        inverse_dico_m2b_index_match[v] = k


    list_unique_cell = list(pred_df.cell.unique())
    list_unique_cyto = list(np.unique(cell_seg))
    np.isin(list_cell_cyto, list_unique_cyto).all()

    model_name = "comseg"

    dict_res[model_name] = {}
    dict_res[model_name]['aji'] = aji
    dict_res[model_name]['C_image'] = C_image
    dict_res[model_name]['U_image'] = U_image
    dict_res[model_name]['dico_cell_pred_rna'] = dico_cell_pred_rna
    dict_res[model_name]['max_iou_list'] = max_iou_list
    dict_res[model_name]['dico_m2b_index_match'] = dico_m2b_index_match
    print(np.mean(max_iou_list))
    print(len(max_iou_list))
    print()
    np.save(f"/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/data_release_baysor_merfish_gut/mask_for_metrics/dict_res_{str_save}.npy", dict_res)


    ### pciseq
    pred_df = pd.read_csv("/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/pciseq/results/geneData.csv")
    model_name = "pciseq"

    aji, C_image, U_image, dico_cell_pred_rna, max_iou_list, dico_m2b_index_match = compute_perf(
                                                                                                pred_df,
                                                                                                dict_gr_cm=dict_gr_cm,
                                                                                                cell_column_name="neighbour"
                                                                                                 )
    dict_res[model_name] = {}
    dict_res[model_name]['aji'] = aji
    dict_res[model_name]['C_image'] = C_image
    dict_res[model_name]['U_image'] = U_image
    dict_res[model_name]['dico_cell_pred_rna'] = dico_cell_pred_rna
    dict_res[model_name]['max_iou_list'] = max_iou_list
    dict_res[model_name]['dico_m2b_index_match'] = dico_m2b_index_match
    print(np.mean(max_iou_list))
    print(len(max_iou_list))
    print()


    ### baysor + dapi


    pred_df = pd.read_csv(
        "/home/tom/Bureau/phd/simulation/baysor/baysor_ubuntu-latest_x64_build/Baysor/bin/mouse_iliuem_dapi/segmentation.csv")
    model_name = "baysor_dapi"

    aji, C_image, U_image, dico_cell_pred_rna, max_iou_list, dico_m2b_index_match = compute_perf(pred_df,
                                                                                                 dict_gr_cm=dict_gr_cm,

                                                                                                 cell_column_name="cell")
    dict_res[model_name] = {}
    dict_res[model_name]['aji'] = aji
    dict_res[model_name]['C_image'] = C_image
    dict_res[model_name]['U_image'] = U_image
    dict_res[model_name]['dico_cell_pred_rna'] = dico_cell_pred_rna
    dict_res[model_name]['max_iou_list'] = max_iou_list
    dict_res[model_name]['max_iou_list'] = max_iou_list
    dict_res[model_name]['dico_m2b_index_match'] = dico_m2b_index_match
    print(np.mean(max_iou_list))
    print(len(max_iou_list))
    print()
    np.save(f"/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/data_release_baysor_merfish_gut/mask_for_metrics/dict_res_{str_save}.npy", dict_res)




    ######### watershed
    pred_df = pd.read_csv(
        "/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/watershed/df.csv")
    model_name = "watershed"

    aji, C_image, U_image, dico_cell_pred_rna, max_iou_list, dico_m2b_index_match = compute_perf(pred_df,
                                                                                                 dict_gr_cm=dict_gr_cm,
                                                                                                 cell_column_name="cell_id")
    dict_res[model_name] = {}
    dict_res[model_name]['aji'] = aji
    dict_res[model_name]['C_image'] = C_image
    dict_res[model_name]['U_image'] = U_image
    dict_res[model_name]['dico_cell_pred_rna'] = dico_cell_pred_rna
    dict_res[model_name]['max_iou_list'] = max_iou_list
    dict_res[model_name]['max_iou_list'] = max_iou_list
    dict_res[model_name]['dico_m2b_index_match'] = dico_m2b_index_match


    print(np.mean(max_iou_list))
    print(len(max_iou_list))
    print()

    np.save(f"/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/data_release_baysor_merfish_gut/mask_for_metrics/dict_res_{str_save}.npy", dict_res)



    ### scs


    pred_df = pd.read_csv(
        "/home/tom/share/baysor_data/scs/entire_image/all1/results/df_res.csv")
    model_name = "scs_e300_all"
    aji, C_image, U_image, dico_cell_pred_rna, max_iou_list, dico_m2b_index_match = compute_perf(pred_df,
                                                                                                 dict_gr_cm=dict_gr_cm,
                                                                                                 cell_column_name="cell_pred")
    dict_res[model_name] = {}
    dict_res[model_name]['aji'] = aji
    dict_res[model_name]['C_image'] = C_image
    dict_res[model_name]['U_image'] = U_image
    dict_res[model_name]['dico_cell_pred_rna'] = dico_cell_pred_rna
    dict_res[model_name]['max_iou_list'] = max_iou_list
    dict_res[model_name]['max_iou_list'] = max_iou_list
    dict_res[model_name]['dico_m2b_index_match'] = dico_m2b_index_match


    print(np.mean(max_iou_list))
    print(len(max_iou_list))
    print()

    np.save(f"/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/data_release_baysor_merfish_gut/mask_for_metrics/dict_res_{str_save}.npy", dict_res)


    ### baysor
    pred_df = pd.read_csv("/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/data_release_baysor_merfish_gut/data_analysis/baysor/segmentation/segmentation.csv")
    model_name = "baysor"
    aji, C_image, U_image, dico_cell_pred_rna, max_iou_list, dico_m2b_index_match = compute_perf(pred_df,
                                                                                                 dict_gr_cm=dict_gr_cm,
                                                                                                 cell_column_name="cell")
    dict_res[model_name] = {}
    dict_res[model_name]['aji'] = aji
    dict_res[model_name]['C_image'] = C_image
    dict_res[model_name]['U_image'] = U_image
    dict_res[model_name]['dico_cell_pred_rna'] = dico_cell_pred_rna
    dict_res[model_name]['max_iou_list'] = max_iou_list
    print(np.mean(max_iou_list))
    print(len(max_iou_list))
    print()

    np.save(f"/media/tom/Transcend1/doi_10_5061_dryad_jm63xsjb2__v20210916/data_release_baysor_merfish_gut/mask_for_metrics/dict_res_{str_save}.npy", dict_res)


    for model in dict_res.keys():
        print(model)
        list_iou = dict_res[model]['max_iou_list']
        print(f'mean iou {np.mean(list_iou)}')
        print()












