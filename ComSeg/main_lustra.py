


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
                        default=40)  # TODO
    parser.add_argument("--mean_cell_diameter", type=float,
                        default=20)  # TODO

    parser.add_argument("--path_dataset_folder", type=str,
                        default="/media/tom/T7/lustra_all/2023-09-06_LUSTRA-14rounds/image_per_pos/dataframe_ComSeg/")  # TODO
    parser.add_argument("--path_to_mask_prior", type=str,
                        default="")  # TODO
    parser.add_argument("--path_dict_cell_centroid", type=str,
                        default="/media/tom/T7/lustra_all/2023-09-06_LUSTRA-14rounds/image_per_pos/dico_centroid_ComSeg/")
    parser.add_argument("--path_simulation_gr", type=str,
                        default='/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/dico_simulation/')  # TODO

    parser.add_argument("--path_centroid", type=str,
                        default='/home/tom/Bureau/phd/simulation/mycode/centroid/')  # TODO

    parser.add_argument("--weight_mode", type=str,
                        default='original')  # TO
    parser.add_argument("--eps_min_weight", type=float, default=0)

    parser.add_argument("--port", default=3950)
    parser.add_argument("--mode", default='client')
    parser.add_argument("--host", default='127.0.0.2')

    args = parser.parse_args()
    args.path_save = args.path_dataset_folder + "results/" + date_str + "/"
    Path(args.path_save).mkdir(parents=True, exist_ok=True)


    print(args)

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
    args.max_dist_centroid = mean_cell_diameter / 2




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


    ####### compute the co-expression
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
    from comseg import clustering
    import importlib

    importlib.reload(comseg)
    importlib.reload(model)
    importlib.reload(dictionary)
    importlib.reload(clustering)

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
        n_pcs=0,
        n_comps=0,
        clustering_method="leiden",
        n_neighbors=10,
        resolution=0.1,
        palette=None,
        min_merge_correlation=0.8,
        nb_min_cluster=4,
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

    adata = Comsegdict.final_anndata
    print(adata[np.sum(adata.X, axis=1) > 5, :])
    #sc.tl.pca(adata, svd_solver='arpack', n_comps = 0)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=0)
    sc.tl.leiden(adata,  resolution=0.5)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=["leiden"], palette=palette, legend_loc='on data')
    Comsegdict.save(Path(args.path_save) / ("Comsegdict." + date_str + "pickle.txt"))
    ## adata object save in comseg
    import pickle
    for image in Comsegdict:
        G = Comsegdict[image].G
        with open(Path(args.path_save) / (image + "." + date_str + "pickle.txt"), 'wb') as handle:
            pickle.dump(G, handle, -1)
    adata.write_h5ad(Path(args.path_save) / ("adata." + date_str + "h5ad"))



