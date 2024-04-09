








##############################
# Example of how we benchmarked Baysor
#############################
import sys
import os
import numpy as np
import pandas as pd
import tifffile
from sklearn import metrics
from pathlib import Path
sys.path.insert(1, os.getcwd() )

sys.path.insert(1, os.getcwd() + "/code/")
import list_markers_hardcoded

from utils.data_processing import sctransform_from_parameters
from utils.custom_metrics import  compute_AJI
from tqdm import tqdm
import anndata  as ad

#%% extrat information from cvs result baysor




#%%


if __name__ == "__main__":


    ##################""
    ## Preprocessing
    ##################

    ### cell type centroid
    scrna_centroids = np.load(
        "/centroid/scrna_centroidsnsforest0_thalassadico_simulationnornmalized_True.npy",
        allow_pickle=True)

    ### cell type name
    scrna_unique_clusters = np.load(
        "/centroid/scrna_unique_clustersnsforest0_thalassadico_simulation_nsforest0_thalassa.npy",
        allow_pickle=True)

    ### normalization parameter
    param_sctransform_anndata = np.load(
        "/centroid/param_sctransform_anndatansforest0_thalassa.npy",
        allow_pickle=True)
    scale = "set_with_mask'"
    nb_cluster = 5
    cf_prior = 1
    script_name = f"max_{scale}_nb_cluster{nb_cluster}"
    path_save_csv = "/home/tom/Bureau/phd/simulation/baysor/baysor_csv_input/" + script_name + '/'
    path_save_result = "/home/tom/Bureau/phd/simulation/baysor/baysor_csv_input/" + script_name + '/'
    save_res_paper = f"/home/tom/Bureau/phd/simulation/paper_result/baysor/" + script_name + '/'
    Path(save_res_paper).mkdir(exist_ok=True)
    path_nuclei_mask = "/media/tom/T7/regular_grid/simu1912/cube2D_step100/remove20/nuclei_tif/"
    path_nuclei = "/media/tom/T7/regular_grid/simu1912/cube2D_step100/remove20/nuclei/"
    gene_list  = np.unique(list_markers_hardcoded.nsforest0_thalassa)
    REAL_DATA = False


    last_letter = None
    dico_scale = {"x": 1, 'y': 1, "z":1}
    path_to_csv_folder = "/media/tom/T7/regular_grid/simu1912/cube2D_step100/dataframe_folder/ns0_talassa_max3/"
    use_um_scale = True
    Path(path_save_result).mkdir(exist_ok=True)
    twoD = False
    dico_param = {'scale' : scale, 'cf_prior' :cf_prior, 'nb_cluster' : nb_cluster, 'use_um_scale' : use_um_scale, 'twoD' : twoD}
    np.save(Path(save_res_paper) / "dico_param.npy", dico_param)




    ##################
    ## Generate Baysor input dataframe
    ##################
    for path_to_csv in Path(path_to_csv_folder).glob('*csv'):
        assert last_letter is None or last_letter < 0
        print(path_to_csv.name)
        df_spots_label = pd.read_csv(path_to_csv)
        df_spots_label[':in_nucleus'] = df_spots_label['in_nucleus']
        if use_um_scale:
            df_spots_label = df_spots_label.rename(columns={"cell": "cell_gr"})
            df_spots_label.x = df_spots_label.x * dico_scale["x"]
            df_spots_label.y = df_spots_label.y * dico_scale["y"]
            if not twoD:
                df_spots_label.z = df_spots_label.z * dico_scale["z"] + 0.01
        df_spots_label.to_csv(path_save_csv + path_to_csv.name ) #, sep='\t', encoding='utf-8')


    ##################
    # generate .sh script to apply baysor
    ##############"

    for path_to_csv in Path(path_to_csv_folder).glob('*csv'):
        image_name = path_to_csv.stem
        print(image_name)
        with open(f"/home/tom/Bureau/phd/simulation/baysor/baysor_ubuntu-latest_x64_build/Baysor/bin/{script_name}.sh", "a+") as file:
            file.write(f"\n mkdir {path_save_result}{image_name} \n"
                       f"./Baysor run  -x x -y y -z z --gene gene -c config.toml -o {path_save_result}{image_name} "
                         f"-p {path_save_csv}{image_name}.csv "
                         f"{path_nuclei_mask}{image_name}.tif "
                         f"--n-clusters {nb_cluster} --save-polygons GeoJSON -p --prior-segmentation-confidence {cf_prior} --no-ncv-estimation \n"
                       )








    REAL_DATA = False
    dico_gr_cm = {}
    dico_pred_cm = {}
    list_aji_global = []
    dico_cell_vec = {}
    error_list = []
    prc_not_catch_list = []
    list_iou2 = []
    list_list_iou = []
    list_max_iou = []
    gene_index_dico = {}
    dico_cell_pred_rna_all = {}
    count_matrix_pred_all = []



    #################
    # COMPUTE IOU
    #################

    #### dictionary with the ground_truth
    path_simu_only = "/media/tom/T7/regular_grid/simu1912/cube2D_step100/dico_simulation_nsforest0_thalassa_v0/"


    list_images_path = [img_name.stem for img_name in Path(path_to_csv_folder).glob('*csv')]
    for gene_id in range(len(gene_list)):
        gene_index_dico[gene_list[gene_id]] = gene_id
    for image_path in tqdm(list_images_path):
        image_name = str(image_path).split("/")[-1]
        print(image_name)
        segmentation = pd.read_csv(path_save_result + f"{image_name[:last_letter]}/segmentation.csv")
        list_cell_pred = np.array(segmentation.cell)#### generate count matrix pred ####
        dico_pred_cm[image_name] = {}
        list_cell_pred = np.array(segmentation.cell)
        list_gene = np.array(segmentation.gene)
        index_mol_pred = np.array(range(len(list_cell_pred)))
        list_unique_cell_pred = np.unique(list_cell_pred)
        if 0 in list_unique_cell_pred:
            list_unique_cell_pred = list_unique_cell_pred[1:]
        for cell_id in list_unique_cell_pred:
            dico_pred_cm[image_name][cell_id] = index_mol_pred[list_cell_pred == cell_id]
        # print(len(dico_pred_cm[image_name]))
        image_pred_cm = dico_pred_cm[image_name]
        count_matrix_pred = []
        dico_cell_vec_local = {}
        for cell in image_pred_cm:
            dico_cell_vec_local[cell] = {}
            dico_cell_vec_local[cell]["genes"] = list_gene[image_pred_cm[cell]]

            expression_vector = np.bincount(
                [gene_index_dico[gene] for gene in dico_cell_vec_local[cell]["genes"]],
                minlength=len(gene_index_dico))
            count_matrix_pred.append(expression_vector)
        if not REAL_DATA:
            dico_simulation = np.load(path_simu_only + image_name + ".npy",
                                      allow_pickle=True).item()
            dico_gr_cm[image_name] = {}
            list_cell_gr = np.array(segmentation.cell_gr)
            index_mol_gr = np.array(range(len(list_cell_gr)))
            list_unique_cell = np.unique(np.load(path_nuclei + image_name + ".npy"))
            if 0 in list_unique_cell:
                assert list_unique_cell[0] == 0
                list_unique_cell = list_unique_cell[1:]
            for cell_id in list_unique_cell:
                dico_gr_cm[image_name][cell_id] = index_mol_gr[list_cell_gr == cell_id]
            image_gr_cm = dico_gr_cm[image_name]
            aji, C_image, U_image, dico_cell_pred_rna, max_iou_list, dico_b2m_index_match =  compute_AJI(
                                                image_gr_cm=image_gr_cm,
                                                image_pred_cm = image_pred_cm,
                                                max_index = None)
            list_max_iou += max_iou_list
            list_aji_global.append(aji)
            dico_cell_pred_rna_all[image_name] = dico_cell_pred_rna


            for nuc_index in image_gr_cm:
                try:
                    iou = len(
                        set(dico_cell_pred_rna[nuc_index]).intersection(image_gr_cm[nuc_index])) / len(
                        set(dico_cell_pred_rna[nuc_index]).union(image_gr_cm[nuc_index]))
                    list_iou2.append(iou)
                except Exception as e:
                    print("union is empty")
                    print(e)
                try:
                    prc_false = len(
                        set(dico_cell_pred_rna[nuc_index]) - set(image_gr_cm[nuc_index])) / len(
                        set(dico_cell_pred_rna[nuc_index]))
                    error_list.append(prc_false)
                except Exception as e:
                    print(e)

                try:
                    prc_not_catch = len(
                        set(image_gr_cm[nuc_index]) - set(dico_cell_pred_rna[nuc_index])) / len(
                        image_gr_cm[nuc_index])
                    prc_not_catch_list.append(prc_not_catch)
                except Exception as e:
                    print(e)





    #################
    # COMPUTE acc
    #################

    list_cell_type_pred = []
    list_cell_type_gr = []
    count_matrix_pred = []
    for image_path in tqdm(list_images_path):
        image_name = str(image_path).split("/")[-1]
        print(image_name)
        dico_cell_vec[image_name] =  {}
        segmentation = pd.read_csv(path_save_result + f"{image_name[:last_letter]}/segmentation.csv")
        dico_simulation = np.load(path_simu_only + image_name +".npy",
                                  allow_pickle=True).item()
        list_gene = np.array(segmentation.gene)
        dico_cell_pred_rna = dico_cell_pred_rna_all[image_name]
        for cell in dico_cell_pred_rna:
            dico_cell_vec[image_name][cell] = {}
            dico_cell_vec[image_name][cell]["genes"] = list_gene[dico_cell_pred_rna[cell]]
            expression_vector = np.bincount(
                [gene_index_dico[gene] for gene in dico_cell_vec[image_name][cell]["genes"]],
                minlength=len(gene_index_dico))
            dico_cell_vec[image_name][cell]["vec"] = expression_vector
            if expression_vector.ndim != 2:
                expression_vector = expression_vector.reshape(1, len(expression_vector))
            if param_sctransform_anndata is not None and  not param_sctransform_anndata.shape == ():
                norm_expression_vectors = sctransform_from_parameters(
                    np.array(param_sctransform_anndata),
                    expression_vector)
            else:
                norm_expression_vectors = expression_vector
            if np.sum(expression_vector) == 0:
                pred_rna_seq = "None"
            else:
                correlation_array = metrics.pairwise.cosine_similarity(norm_expression_vectors,
                                                                       np.array(scrna_centroids))[0]
                index_cluster_max = np.argmax(correlation_array)
                pred_rna_seq = scrna_unique_clusters[index_cluster_max]
                dico_cell_vec[image_name][cell]["pred_cell_type"] = pred_rna_seq
                dico_cell_vec[image_name][cell]["gr"] = dico_simulation['dico_cell_index'][cell]['type']
                list_cell_type_pred.append(pred_rna_seq)
                list_cell_type_gr.append(dico_simulation['dico_cell_index'][cell]['type'])
        acc = metrics.accuracy_score(y_true=list_cell_type_gr,
                                     y_pred=list_cell_type_pred)
        print(f'acc {acc}')
    from sklearn import metrics
    acc = metrics.accuracy_score(y_true=list_cell_type_gr,
                                 y_pred=list_cell_type_pred)
    dico_f1_cell_type = {k: v for k, v in zip(scrna_unique_clusters, metrics.f1_score(list_cell_type_gr,
                                                                                          list_cell_type_pred,
                                                                                          labels=scrna_unique_clusters,
                                                                                          average=None))}
    print(script_name)
    print(f'mean aji {np.mean(list_aji_global)}')
    print(f' acc from aji {acc}')
    print(f'median prc_not_catch_list {round(np.median(prc_not_catch_list), 4)}')
    print(f'median error_list {round(np.median(error_list), 4)}')
    print(f'mean prc_not_catch_list {round(np.mean(prc_not_catch_list), 4)}')
    print(f'mean error_list {round(np.mean(error_list), 4)}')
    print(dico_f1_cell_type)
    print(f'mean iou {np.mean(list_iou2)}')
    print(f'median iou {np.median(list_iou2)}')
    print(f'mean prc_not_catch_list {round(np.mean(prc_not_catch_list), 4)}')
    print(f'mean error_list {round(np.mean(error_list), 4)}')
