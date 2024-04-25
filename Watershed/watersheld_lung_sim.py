





#############################
# Example of how we benchmarked watershed
##################################


import numpy as np
import tifffile
import anndata as ad
from skimage.segmentation import watershed
from scipy import ndimage as ndi
import list_markers_hardcoded

from utils.count_matrix_metrics import sctransform_from_parameters
from tqdm import tqdm
from sklearn import metrics
from utils.custom_metrics import  compute_AJI
from pathlib import Path
import pandas as pd


def get_watersheld(ground_truth_mask, nuclei,
                   pixel_size=np.array([0.3, 0.103, 0.103])):
    if ground_truth_mask is not None:
        empty_space = (ground_truth_mask > 0).astype(np.uint8)
    else:
        empty_space = None
    inverted_mask = np.ones(nuclei.shape, dtype = np.uint8) - (nuclei != 0).astype(np.uint8)
    distance = ndi.distance_transform_edt(inverted_mask,
                                          sampling=pixel_size)
    labels_with_empty = watershed(image=distance,
                                  markers=nuclei,
                                  mask=empty_space)
    return labels_with_empty






if __name__ == '__main__':

    ## load df

    path_folder = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0"
    path_nuclei = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/remove20/nuclei/"
    path_ind_cyto  = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/individual_cytoplasms_rd/"
    path_save_watershed = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/remove20/saved_watershed_seg_wit_gr_mask/"


    path_df_folder = "/media/tom/T7/simulation/exp_same_cyto/mouse_lung_simulation/df_spots/"
    dico_dico_commu = {}
    for path_df in Path(path_df_folder).glob("*.csv"):
        dico_dico_commu[path_df.stem] = {}
        dico_dico_commu[path_df.stem]["df_spots_label"] = pd.read_csv(path_df)


    path_anndata = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/NON_CONV_PAPER26/anndata/"

    scale = [0.3, 0.103, 0.103]
    list_pred = []
    list_gr = []
    list_iou = []
    last_letter = None
    list_aji = []

    #####################
    # Compute Voronoid tesselation
    ######################

    last_letter = -4
    for image_name in tqdm(list(dico_dico_commu.keys())[:30]):
        print(image_name)
        nuclei = tifffile.imread(path_nuclei + image_name+ ".tiff" )#+ ".tiff")
        ground_truth_mask =None
        unique_nuclei_raw = np.unique(nuclei)
        labels_with_empty = get_watersheld(ground_truth_mask, nuclei,
                       pixel_size=np.array(scale))
        np.save(path_save_watershed+image_name, labels_with_empty)




    ###############################
    ### Compute expression matrix
    ##############################

    selected_genes = list_markers_hardcoded.nsforest0_thalassa
    selected_genes = np.unique(selected_genes)
    path_dico_simulation = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/dico_simulation/"

    scrna_centroids = np.load(
        "/centroid/scrna_centroidsnsforest0_thalassadico_simulationnornmalized_True.npy",
        allow_pickle=True)
    scrna_unique_clusters = np.load(
        "/centroid/scrna_unique_clustersnsforest0_thalassadico_simulation_nsforest0_thalassa.npy",
        allow_pickle=True)
    param_sctransform_anndata = np.load(
        "/centroid/param_sctransform_anndatansforest0_thalassa.npy",
        allow_pickle=True)

    gene_index_dico = {}
    for gene_id in range(len(selected_genes)):
        gene_index_dico[selected_genes[gene_id]] = gene_id
    list_count_matrix = []
    dico_dico_spots_pos = {}
    for image_name in tqdm(list(dico_dico_commu.keys())[0:30]):
        dico_dico_spots_pos[image_name] = {}
        try:
            labels_with_empty = np.load(path_save_watershed + image_name + ".npy")
        except:
            labels_with_empty = np.load(path_save_watershed + image_name + ".tiff.npy")
        ground_truth_mask = np.load(path_ind_cyto+image_name  + ".tiff.npy",allow_pickle=True)
        nuclei  = tifffile.imread(path_nuclei + image_name+ ".tiff")
        unique_resize_nuc = list(np.unique(nuclei))
        if 0 in unique_resize_nuc:
            unique_resize_nuc = unique_resize_nuc[1:]
        dico_ground_truth = {}
        dico_simulation = np.load(Path(path_dico_simulation) /  (image_name+ ".tiff.npy"), allow_pickle=True ).item()
        for nuc in unique_resize_nuc:
            dico_ground_truth[nuc] = dico_simulation['dico_cell_index'][nuc]["type"]
        df_spots_label = dico_dico_commu[image_name]["df_spots_label"]
        ### compute the cunt matrix (count dico)
        dico_count_matrix = {}
        list_x = list(df_spots_label.x)
        list_y = list(df_spots_label.y)
        list_z = list(df_spots_label.z)
        list_df_gene = list(df_spots_label.gene)
        list_cell_pred = ["-1"] * len(list_df_gene)

        image_gr_cm = {} # { cell_id :list of index of mol belonging to the cell}}
        image_pred_cm = {} # { cell_id :list of index of mol belonging to the cell}}
        count_matrix  = []
        list_nuc_pred = []
        dico_spots_pos = {}
        for nucleus in tqdm(unique_resize_nuc):
            dico_count_matrix[nucleus] = {}
            dico_count_matrix[nucleus]["list_gene"] = []
            dico_count_matrix[nucleus]["list_gene_index_pred"] = []
            dico_count_matrix[nucleus]["list_gene_gr"] = []
            dico_count_matrix[nucleus]["list_gene_index_gr"] = []
            dico_spots_pos[nucleus] = []

            gr_type = dico_ground_truth[nucleus]
            #assert dico_dico_commu[image_name]['dico_ground_truth'][nucleus] == gr_type
            dico_count_matrix[nucleus]['type'] = gr_type
        for gene_index in range(len(list_df_gene)):
            nuc_pred = labels_with_empty[list_z[gene_index], list_y[gene_index], list_x[gene_index]]
            if nuc_pred != 0 and nuc_pred in unique_resize_nuc:
                dico_count_matrix[nuc_pred]["list_gene"].append(list_df_gene[gene_index])
                dico_count_matrix[nuc_pred]["list_gene_index_pred"].append(gene_index)
                dico_spots_pos[nuc_pred].append([list_z[gene_index], list_y[gene_index], list_x[gene_index]])
                #list_cell_pred[gene_index] = nuc_pred
                #list_nuc_pred.append(nuc_pred)
            nuc_gr = ground_truth_mask[list_z[gene_index], list_y[gene_index], list_x[gene_index]]
            if nuc_gr != 0 and nuc_gr in unique_resize_nuc:
                dico_count_matrix[nuc_gr]["list_gene_gr"].append(list_df_gene[gene_index])
                dico_count_matrix[nuc_gr]["list_gene_index_gr"].append(gene_index)
                #list_nuc_pred.append(nuc_pred)

        for nucleus in tqdm(unique_resize_nuc):
            list_gene = dico_count_matrix[nucleus]["list_gene"]
            expression_vector = np.zeros(len(selected_genes))
            for str_gene in list_gene:
                expression_vector[gene_index_dico[str_gene]] += 1
            dico_count_matrix[nucleus]["expression_vector"] = expression_vector
            count_matrix.append(expression_vector)

        dico_dico_commu[image_name]["dico_count_matrix"] = dico_count_matrix
        #df_spots_label["cell_pred"] = list_cell_pred
        dico_dico_commu[image_name]['df_spots_label'] = df_spots_label
        list_count_matrix += count_matrix
        dico_dico_spots_pos[image_name] = dico_spots_pos





    param_sctransform = param_sctransform_anndata

    # prediction
    list_aji = []
    list_total_iou = []
    list_pred = []
    list_gr = []
    mol_list_cell_true = []
    mol_list_cell_pred = []

    mol_list_cell_type_true = []
    mol_list_cell_type_pred = []

    nb_total_mol = 0
    nb_well_classif = 0
    error_list_molecule = []
    not_catch_list_molecule = []
    all_molecule = []
    total_molecule = []

    list_total_prc_not_catch = []
    list_total_prc_error = []


    for image_name in tqdm(list(dico_dico_commu.keys())[:]):
        print()
        dico_count_matrix = dico_dico_commu[image_name]["dico_count_matrix"]
        df_spots_label = dico_dico_commu[image_name]["df_spots_label"]

        dico_nucleus_celltype = {}
        dico_nucleus_celltype["-1"] = "-1"
        for nucleus in dico_count_matrix:
            if np.sum(dico_count_matrix[nucleus]["expression_vector"]) == 0:
                gr_type = dico_count_matrix[nucleus]['type']
                list_pred.append("np pred")
                list_gr.append(gr_type)
                continue
            expression_vector = dico_count_matrix[nucleus]["expression_vector"]
            array_of_vect = expression_vector.reshape(1, len(expression_vector))

            if  len(param_sctransform) > 1 is not None or param_sctransform.item() :
                norm_expression_vectors = sctransform_from_parameters(
                    np.array(param_sctransform),
                    array_of_vect)

            else:
                norm_expression_vectors = array_of_vect
            dico_count_matrix[nucleus]["norm_expression_vectors"] = norm_expression_vectors
            correlation_array =  metrics.pairwise.cosine_similarity(norm_expression_vectors,
                                                                   np.array(scrna_centroids))[0]
            index_cluster_max = np.argmax(correlation_array)
            pred_rna_seq = scrna_unique_clusters[index_cluster_max]
            dico_nucleus_celltype[nucleus] = pred_rna_seq
            dico_count_matrix[nucleus]['pred'] = pred_rna_seq
            gr_type = dico_count_matrix[nucleus]['type']
            nb_total_mol += np.sum(array_of_vect)
            if pred_rna_seq == gr_type:
                nb_well_classif += np.sum(array_of_vect)
            mol_list_cell_type_true += [gr_type]* int(np.sum(array_of_vect))
            mol_list_cell_type_pred += [pred_rna_seq] * int(np.sum(array_of_vect))
            list_pred.append(pred_rna_seq)
            list_gr.append(gr_type)
        image_gr_cm = {}
        image_pred_cm =  {}
        for nuc in dico_count_matrix:
            image_gr_cm[nuc] = dico_count_matrix[nuc]["list_gene_index_gr"]
            image_pred_cm[nuc] = dico_count_matrix[nuc]["list_gene_index_pred"]
        aji, C_image, U_image, dico_cell_pred_rna, list_iou, _ = compute_AJI(image_gr_cm=image_gr_cm,
                                                                image_pred_cm= image_pred_cm,
                                                                max_index = None,
                                                                cell_index_backgroud =  0)
        error_list =  []
        prc_not_catch_list = []
        list_iou2  = []
        prc_not_catch_list_molecule = []
        for nuc_index in dico_count_matrix:
            try:
                iou = len(
                    set(image_pred_cm[nuc_index]).intersection(image_gr_cm[nuc_index])) / len(
                    set(image_pred_cm[nuc_index]).union(image_gr_cm[nuc_index]))
                list_iou2.append(iou)
                total_molecule += list(image_gr_cm[nuc_index])
            except Exception as e:
                    print(e)
            try:
                error_list_molecule += list(set(image_pred_cm[nuc_index]) - set(image_gr_cm[nuc_index]))
                prc_false = len(
                    set(image_pred_cm[nuc_index]) - set(image_gr_cm[nuc_index])) / len(
                    set(image_pred_cm[nuc_index]))
                error_list.append(prc_false)
            except Exception as e:
                print(e)
                error_list.append(0)
            try:
                #prc_not_catch_list_molecule += list(set(image_gr_cm[nuc_index]) - set(image_pred_cm[nuc_index]))
                all_molecule += list(image_gr_cm[nuc_index])
                prc_not_catch = len(
                    set(image_gr_cm[nuc_index]) - set(image_pred_cm[nuc_index])) / len(
                    image_gr_cm[nuc_index])
                prc_not_catch_list.append(prc_not_catch)
            except Exception as e:
                print(e)
        list_total_iou.append(list_iou2)
        list_total_prc_not_catch += prc_not_catch_list
        list_total_prc_error += error_list
        list_aji.append(aji)
        acc = metrics.accuracy_score(list_pred, list_gr)
        print(f'acc cell type nucleus {acc}"')
    print(f'mean aji {round(np.mean(list_aji), 4)}')
    print(f'mean list_total_iou {round(np.mean(np.concatenate(list_total_iou)), 4)}')
    print(f'medain list_total_iou {round(np.median(np.concatenate(list_total_iou)), 4)}')
    print(f'median  total prc_not_catch_list {round(np.median(list_total_prc_not_catch), 4)}')
    print(f'median total error_list {round(np.median(list_total_prc_error), 4)}')
    print(f'mean total prc_not_catch_list {round(np.mean(list_total_prc_not_catch), 4)}')
    print(f'mean total error_list {round(np.mean(list_total_prc_error), 4)}')
    print(f'acc cell type nucleus {acc}"')