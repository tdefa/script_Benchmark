



import matplotlib
#matplotlib.use('Qt5Agg')
import sys

import sys

sys.path += ['/home/tom/Bureau/phd/simulation/ComSeg_pkg/src']


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
importlib.reload(model)
importlib.reload(dictionary)

import datetime
import argparse

#importlib.reload(dataset)

if __name__ == '__main__':

    e = datetime.datetime.now()
    date_str = f"{e.month}_d{e.day}_h{e.hour}_min{e.minute}_s{e.second}_r" + str(random.randint(0, 5000))
    parser = argparse.ArgumentParser(description='test')


    ###################
    ###### hyper-parameter
    ###################
    parser.add_argument("--max_cell_radius", type=float,
                        default=15)  # TODO
    parser.add_argument("--mean_cell_diameter", type=float,
                        default=10)  # TODO

    parser.add_argument("--path_dataset_folder", type=str,
                        default="/media/tom/T7/regular_grid/simu1912/elbow_cube/dataframe_max3/")  # TODO
    parser.add_argument("--path_to_mask_prior", type=str,
                        default="/media/tom/T7/regular_grid/simu1912/elbow_cube/remove20/nuclei/")  # TODO
    parser.add_argument("--path_dict_cell_centroid", type=str,
                        default="/media/tom/T7/regular_grid/simu1912/elbow_cube/remove20/dico_centroid/")
    parser.add_argument("--path_simulation_gr", type=str,
                        default='/media/tom/T7/regular_grid/simu1912/elbow_cube/dico_simulation_nsforest0_thalassa/')  # TODO
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

    args.k_nearest_neighbors = 10
    args.n_neighbors_centroid = args.k_nearest_neighbors
    args.max_dist_centroid = mean_cell_diameter / 4



    args.path_save = args.path_dataset_folder + "results/" + date_str + "/"
    Path(args.path_save).mkdir(parents=True, exist_ok=True)

    ## scale/ pixel size in um
    dict_scale = {"x": 0.150, 'y': 0.150, "z": 0.3}

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
            list_corr_gene0.append(dataset_non_conv.dict_co_expression[gene0][gene1])
        corr_matrix.append(list_corr_gene0)
    list_gene = list(dataset_non_conv.dict_co_expression.keys())
    # plotting the heatmap for correlation
    ax = sns.heatmap(corr_matrix, xticklabels=list_gene, yticklabels=list_gene, )
    plt.show()

    import comseg
    from comseg import model
    from comseg import dictionary
    import importlib

    importlib.reload(comseg)
    importlib.reload(model)
    importlib.reload(dictionary)

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



    ###################""



    dico_dico_commu = {}
    selected_genes = dataset_non_conv.selected_genes
    gene_index_dico = {}
    for gene_id in range(len(selected_genes)):
        gene_index_dico[selected_genes[gene_id]] = gene_id

    list_image = [img for img in Comsegdict if  str(type(Comsegdict[img]))]

    for img_name in tqdm(list_image):
        dico_dico_commu[img_name] = {}
        cell_unique = np.unique(np.load(path_to_mask_prior + img_name + ".npy"))
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
