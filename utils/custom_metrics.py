


import numpy as np

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
    dico_gr_cm :  { cell_id :list of index of mol belonging to the cell}}
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
    for cell_gr in image_gr_cm:  # forEach ground truth nucleus  Gi do :j←arg maxk(|Gi∩Sk|/|Gi∪Sk|)
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