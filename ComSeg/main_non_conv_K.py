


import matplotlib
#matplotlib.use('Qt5Agg')
import sys
sys.path += ["/home/tom/Bureau/phd/simulation/ComSeg_pkg/src"]



import numpy as np
import scanpy as sc
import random
import tifffile
from comseg import dataset
import importlib
import comseg
import comseg
from comseg import model
from comseg import dictionary
import importlib

from pathlib import Path
from tqdm import tqdm
from comseg.utils.preprocessing import sctransform_from_parameters
importlib.reload(dataset)

import argparse
import datetime
#importlib.reload(dataset)

if __name__ == '__main__':

    e = datetime.datetime.now()
    print(e)
    date_str = f"{e.month}_d{e.day}_h{e.hour}_min{e.minute}_s{e.second}_r" + str(random.randint(0, 5000))
    parser = argparse.ArgumentParser(description='test')
    ###################
    ###### hyper-parameter
    ###################
    parser.add_argument("--max_cell_radius", type=float,
                        default=15)
    parser.add_argument("--mean_cell_diameter", type=float,
                        default=15)

    parser.add_argument("--path_dataset_folder", type=str,
                        default="/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/test_dataframe/")  # TODO
    parser.add_argument("--path_to_mask_prior", type=str,
                        default="/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/remove20/nuclei/")  # TODO
    parser.add_argument("--path_dict_cell_centroid", type=str,
                        default="/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/remove20/dico_centroid/")
    parser.add_argument("--path_simulation_gr", type=str,
                        default='/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/dico_simulation/')  # TODO

    parser.add_argument("--path_centroid", type=str,
                        default='/home/tom/Bureau/phd/simulation/mycode/centroid/')  # TODO



    parser.add_argument("--port", default=3950)
    parser.add_argument("--mode", default='client')
    parser.add_argument("--host", default='127.0.0.2')

    args = parser.parse_args()

    #### HYPERPARAMETER ####
    max_cell_radius = args.max_cell_radius
    mean_cell_diameter = args.mean_cell_diameter



    path_dataset_folder = args.path_dataset_folder
    ##path to your prior segmentation mask
    path_to_mask_prior = args.path_to_mask_prior

    path_dict_cell_centroid = args.path_dict_cell_centroid
    path_simulation_gr = args.path_simulation_gr


    args.path_save = args.path_dataset_folder + "results/" + date_str + "/"
    Path(args.path_save).mkdir(parents=True, exist_ok=True)


    ### define the gene you want to study, you can restrict it to few genes you want to study.


    ## scale/ pixel size in um
    dict_scale = {"x": 0.103, 'y': 0.103, "z": 0.3}

    ### create the dataset object
    dataset_non_conv = comseg.dataset.ComSegDataset(
        path_dataset_folder=path_dataset_folder,
        path_to_mask_prior=path_to_mask_prior,
        dict_scale=dict_scale,
        mask_file_extension=".npy",
        mean_cell_diameter=mean_cell_diameter
    )

    ### add prior knowledge, here using nucleus segmentation mask
    dataset_non_conv.add_prior_from_mask(prior_keys_name='in_nucleus', overwrite=True)

    ### compute the co-expression correlation at the dataset scale
    if "dict_co_expression" not in dataset_non_conv.__dict__:
        dico_proba_edge, count_matrix = dataset_non_conv.compute_edge_weight(  # in micrometer
            images_subset=None,
            distance="pearson",
        )

    import seaborn as sns
    from matplotlib import pyplot as plt

    corr_matrix = []

    for gene0 in dataset_non_conv.dict_co_expression:
        list_corr_gene0 = []
        for gene1 in dataset_non_conv.dict_co_expression:
            list_corr_gene0.append(dataset_non_conv.dict_co_expression[gene0][gene1]  )
        corr_matrix.append(list_corr_gene0)
    list_gene = list(dataset_non_conv.dict_co_expression.keys())
    # plotting the heatmap for correlation
    ax = sns.heatmap(corr_matrix, xticklabels=list_gene, yticklabels=list_gene, )
    plt.title("Correlation matrix of the dataset")
    plt.show()
    # corr_matrix_false = np.array(corr_matrix).copy()
    import comseg
    from comseg import model
    from comseg import dictionary
    import importlib

    import comseg
    from comseg import model
    from comseg import dictionary
    import importlib
    from comseg import clustering

    Comsegdict = dictionary.ComSegDict(
        dataset=dataset_non_conv,
        mean_cell_diameter=mean_cell_diameter,
        community_detection="with_prior",
        # weights_name="weight",
        prior_name="in_nucleus",
        seed=None,
    )
    Comsegdict.compute_community_vector()

    Comsegdict.compute_insitu_clustering(
        size_commu_min=3,
        norm_vector= True,
        ### parameter clustering
        n_pcs=30,
        n_comps=30,
        clustering_method="leiden",
        n_neighbors=20,
        resolution=1,
        palette=None,
        min_merge_correlation=0.9,
        nb_min_cluster=15,
    )

    import scanpy as sc
    import random

    palette = {}
    for i in range(-1, 500):
        palette[str(i)] = "#" + "%06x" % random.randint(0, 0xFFFFFF)
    adata = Comsegdict.in_situ_clustering.anndata_cluster
    adata.obs["leiden_merged"] = adata.obs["leiden_merged"].astype(str)
    sc.tl.umap(adata)
    fig_ledien = sc.pl.umap(adata, color=["leiden_merged"], palette=palette, legend_loc='on data',
                            )

    Comsegdict.add_cluster_id_to_graph(clustering_method="leiden_merged")

    Comsegdict.classify_centroid(
        path_dict_cell_centroid=path_dict_cell_centroid,
        file_extension=".npy")

    Comsegdict.associate_rna2landmark(
        key_pred="leiden_merged",
        prior_name='in_nucleus',
        distance='distance',
        max_cell_radius=max_cell_radius)

    Comsegdict.anndata_from_comseg_result(
        key_cell_pred='cell_index_pred',
    )

    ############# compute_result metrics #############

    dico_dico_commu = {}
    selected_genes = dataset_non_conv.selected_genes
    gene_index_dico = {}
    for gene_id in range(len(selected_genes)):
        gene_index_dico[selected_genes[gene_id]] = gene_id

    list_image = [img for img in Comsegdict if  str(type(Comsegdict[img]))]

    for img_name in tqdm(list_image):
        dico_dico_commu[img_name] = {}
        cell_unique = np.unique(tifffile.imread(path_to_mask_prior + img_name + ".tiff"))
        G = Comsegdict[img_name].G
        ###################
        dico_expression_m_nuc = {}
        for cell in cell_unique:
            dico_expression_m_nuc[cell] = []
        for node_id, node_data in G.nodes(data=True):
            if node_data['cell_index_pred'] != 0:
                try:
                    dico_expression_m_nuc[node_data['cell_index_pred']].append(node_id)
                except:
                    print('cell with no rna')
        dico_expression_m_nuc_vec = {}
        for nuc_index in dico_expression_m_nuc:
            expression_vector = np.bincount(
                [gene_index_dico[G.nodes[ind_node]["gene"]] for ind_node in dico_expression_m_nuc[nuc_index]
                 if G.nodes[ind_node]["gene"] != 'centroid'], minlength=len(gene_index_dico))
            if expression_vector.ndim != 2:
                expression_vector = expression_vector.reshape(1, len(expression_vector))
            dico_expression_m_nuc_vec[nuc_index] = expression_vector
        dico_dico_commu[img_name]["dico_expression_m_nuc"] = dico_expression_m_nuc
        dico_dico_commu[img_name]["dico_expression_m_nuc_vec"] = dico_expression_m_nuc_vec
        dico_dico_commu[img_name]["G"] = G
        dico_gr_cm = {}
        for cell in cell_unique:
            dico_gr_cm[cell] = [node_index for node_index, n_d in G.nodes(data=True) if
                                "cell" in n_d and n_d["cell"] == cell]
        dico_dico_commu[img_name]["dico_gr_cm"] = dico_gr_cm
        print(img_name + ' erere')







    key_nodes_cell_id_pred = "cell_index_pred"
    key_dico_particule_index = "dico_expression_m_nuc"  # "dico_expression_m_prior_commu"
    key_dico_expression = "dico_expression_m_nuc_vec"
    molecule_label_pred = "leiden2"

    #################
    ### iou / aji ###
    #################

    total_list_iou = []
    total_list_iou_from_aji = []
    total_prc_not_catch_list = []
    total_error_list = []
    total_aji = []

    total_list_iou_classify_only = []
    total_list_iou_from_aji_classify_only = []
    total_prc_not_catch_list_classify_only = []
    total_error_list_classify_only = []

    for image_name in list(dico_dico_commu.keys()):
        list_iou = []
        prc_not_catch_list = []
        error_list = []
        G = dico_dico_commu[image_name]["G"]
        dico_particule_index = dico_dico_commu[image_name][key_dico_particule_index]
        dico_gr_cm = dico_dico_commu[image_name]["dico_gr_cm"]
        list_iou = []
        prc_not_catch_list = []
        error_list = []
        list_iou_classify_only = []
        prc_not_catch_list_classify_only = []
        error_list_classify_only = []
        for nuc_index in dico_particule_index:
            try:
                iou = len(
                    set(dico_particule_index[nuc_index]).intersection(dico_gr_cm[nuc_index])) / len(
                    set(dico_particule_index[nuc_index]).union(dico_gr_cm[nuc_index]))
                list_iou.append(iou)
            except ZeroDivisionError:
                print(ZeroDivisionError)
            except KeyError:
                print(nuc_index)
                print("inconsistency between  df_spots_label and dico_coord_map to check")
            try:
                prc_false = len(
                    set(dico_particule_index[nuc_index]) - set(dico_gr_cm[nuc_index])) / len(
                    set(dico_particule_index[nuc_index]))
                error_list.append(prc_false)

            except ZeroDivisionError:
                print(ZeroDivisionError)
            except KeyError:
                print(nuc_index)
                print("inconsistency between  df_spots_label and dico_coord_map to check")
            try:
                prc_not_catch = len(
                    set(dico_gr_cm[nuc_index]) - set(dico_particule_index[nuc_index])) / len(
                    dico_gr_cm[nuc_index])
                prc_not_catch_list.append(prc_not_catch)
            except KeyError:
                print(nuc_index)
                print("inconsistency between  df_spots_label and dico_coord_map to check")
            except ZeroDivisionError:
                print(ZeroDivisionError)

        total_list_iou += list_iou
        total_prc_not_catch_list += prc_not_catch_list
        total_error_list += error_list

        total_list_iou_classify_only += list_iou_classify_only
        total_prc_not_catch_list_classify_only += prc_not_catch_list_classify_only
        total_error_list_classify_only += error_list_classify_only
        print(f"mean iou {np.mean(list_iou)}")
        print(f"mean prc_not_catch_list {np.mean(prc_not_catch_list)}")
        print(f"mean error_list {np.mean(error_list)}")
        print()
        print()
        print()




    print(f"mean total iou {round(np.mean(total_list_iou), 4)}")


    print(f"mean total prc_not_catch_list {round(np.mean(total_prc_not_catch_list), 4)}")
    print(f"mean total error_list {np.mean(total_error_list)}")
    print()
    print(f"median iou {round(np.median(total_list_iou), 4)}")
    print(f"median prc_not_catch_list {round(np.median(total_prc_not_catch_list), 4)}")
    print(f"median error_list {round(np.median(total_error_list), 4)}")
    print()

    ###### cell id mol matchinf
    from sklearn.metrics.cluster import adjusted_rand_score
    from sklearn import metrics


    for image_name in list(dico_dico_commu.keys()):
        print(image_name)
        cell_unique = np.unique(tifffile.imread(path_to_mask_prior + image_name + ".tiff"))
        dico_simulation = np.load( path_simulation_gr + image_name + '.tiff.npy', allow_pickle = True).item()


        dico_ground_truth = {}
        for nuc in cell_unique[1:]:
            dico_ground_truth[nuc] = dico_simulation["dico_cell_index"][nuc]['type']
            dico_dico_commu[image_name]["dico_ground_truth"] = dico_ground_truth

    seq_scrna_centroids = np.load(Path(args.path_centroid)  / "scrna_centroidsnsforest0_thalassadico_simulationnornmalized_True.npy",
        allow_pickle=True)
    seq_scrna_unique_clusters =  np.load(Path(args.path_centroid)  / "scrna_unique_clustersnsforest0_thalassadico_simulation_nsforest0_thalassa.npy",
        allow_pickle=True)
    seq_param_sctransform_anndata =  np.load(Path(args.path_centroid)  / "param_sctransform_anndatansforest0_thalassa.npy",
        allow_pickle=True)

    y_true_cell_type_classif_only = []
    y_pred_cell_type_classif_only = []
    y_true_cell_type_total = []
    y_pred_cell_type_total = []
    y_true_cell_type_with_no_pred = []
    y_pred_cell_type_with_no_pred = []
    total_cell = 0
    min_rna_to_be_pred = 1
    from sklearn import metrics
    for image_name in list(dico_dico_commu.keys()):
        print(image_name)
        dico_expression = dico_dico_commu[image_name][key_dico_expression]
        G = dico_dico_commu[image_name]["G"]

        dico_ground_truth = dico_dico_commu[image_name]["dico_ground_truth"]
        total_cell += len(dico_ground_truth)
        y_true_cell_type  = []
        y_pred_cell_type = []
        dico_cell_index_pred = {}
        for nuc_index in dico_expression:  ## take into account other non pred list
            if True:
                if nuc_index == 1:
                    if 1 not in dico_ground_truth:
                        print("problem index start to 2")
                        continue
                if nuc_index == 0:
                    print('problem to chck 0 nuclei index')
                    continue
                ground_truth_type = dico_ground_truth[nuc_index]
            expression_vector = dico_expression[nuc_index]
            if np.sum(expression_vector) < min_rna_to_be_pred:
                if True:
                    y_true_cell_type_with_no_pred.append(ground_truth_type)
                y_pred_cell_type_with_no_pred.append("not predicted")
                dico_cell_index_pred[nuc_index] = {'type': "not predicted"}
                continue
            if expression_vector.ndim != 2:
                expression_vector = expression_vector.reshape(1, len(expression_vector))
            if seq_param_sctransform_anndata is not None and not seq_param_sctransform_anndata.shape == ():
                norm_expression_vectors = sctransform_from_parameters(
                    np.array(seq_param_sctransform_anndata),
                    expression_vector)
            else:
                norm_expression_vectors = expression_vector

            correlation_array = metrics.pairwise.cosine_similarity(norm_expression_vectors,
                                                                       np.array(seq_scrna_centroids))[0]
            index_cluster_max = np.argmax(correlation_array)
            pred_rna_seq = seq_scrna_unique_clusters[index_cluster_max]
            dico_cell_index_pred[nuc_index] = {'type': pred_rna_seq}
            #np.save("/media/tom/T7/Stitch/dico_cell_index_pred", dico_cell_index_pred)
            if True:
                y_true_cell_type.append(ground_truth_type)
                y_pred_cell_type.append(pred_rna_seq)
                y_true_cell_type_with_no_pred.append(ground_truth_type)
                y_pred_cell_type_with_no_pred.append(pred_rna_seq)
        acc = metrics.accuracy_score(y_true=y_true_cell_type,
                                     y_pred=y_pred_cell_type)
        print(f"{image_name}  {acc}")
        y_true_cell_type_total += y_true_cell_type
        y_pred_cell_type_total += y_pred_cell_type
        acc = metrics.accuracy_score(y_true=y_true_cell_type_total,
                                     y_pred=y_pred_cell_type_total)
        print(f"acc without taking into account no pred cell {round(acc, 4)}")
        acc_all_cell = metrics.accuracy_score(y_true=y_true_cell_type_with_no_pred,
                                              y_pred=y_pred_cell_type_with_no_pred)
        print(f"acc with taking into account no pred cell {round(acc_all_cell, 4)}")
        print(f"% of predicted {len(y_true_cell_type_total) / len(y_true_cell_type_with_no_pred)}")
        dico_f1_cell_type = {k: v for k, v in
                             zip(seq_scrna_unique_clusters, metrics.f1_score(y_true_cell_type_total,
                                                                             y_pred_cell_type_total,
                                                                             labels=seq_scrna_unique_clusters,
                                                                             average=None))}
        print(f"acc without taking into account no pred cell {round(acc, 4)}")
        acc_all_cell = metrics.accuracy_score(y_true=y_true_cell_type_with_no_pred,
                                              y_pred=y_pred_cell_type_with_no_pred)
        print(f"acc with taking into account no pred cell {round(acc_all_cell, 4)}")
        print(f"% of predicted {len(y_true_cell_type_total) / len(y_true_cell_type_with_no_pred)}")
        print(f"mean iou {round(np.mean(total_list_iou), 4)}")
        print(f"mean prc_not_catch_list {round(np.mean(total_prc_not_catch_list), 4)}")
        print(f"mean error_list {np.mean(total_error_list)}")
        print()
        print(f"median iou {round(np.median(total_list_iou), 4)}")
        print(f"median prc_not_catch_list {round(np.median(total_prc_not_catch_list), 4)}")
        print(f"median error_list {round(np.median(total_error_list), 4)}")
        print()
        print(f"acc without taking into account no pred cell {round(acc, 4)}")
        print()
        print(args)
        Comsegdict.save(Path(args.path_save) / ("Comsegdict." + date_str + "pickle.txt"))
        print("model saved")

