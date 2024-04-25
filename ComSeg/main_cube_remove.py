



import random
from comseg import dataset
import comseg
from pathlib import Path
from tqdm import tqdm
import argparse
import datetime
#importlib.reload(dataset)
import comseg
from comseg import model
from comseg import dictionary
import scanpy as sc
import random

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
                        default="/media/tom/T7/regular_grid/simu1912/cube2D_step100/dataframe_folder/ns0_talassa_max3/")  # TODO
    parser.add_argument("--path_to_mask_prior", type=str,
                        default="/media/tom/T7/regular_grid/simu1912/cube2D_step100/remove20/nuclei/")  # TODO
    parser.add_argument("--path_dict_cell_centroid", type=str,
                        default="/media/tom/T7/regular_grid/simu1912/cube2D_step100/remove20/dico_centroid/")
    parser.add_argument("--path_simulation_gr", type=str,
                        default='/media/tom/T7/regular_grid/simu1912/cube2D_step100/dico_simulation_nsforest0_thalassa_v0/')  # TODO
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
    ## scale/ pixel size in um
    dict_scale = {"x":  0.150, 'y':  0.150, "z": 0.3}

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
    ax = sns.heatmap(corr_matrix, xticklabels=list_gene, yticklabels=list_gene, )
    plt.show()



    Comsegdict = dictionary.ComSegDict(
        dataset=dataset_non_conv,
        mean_cell_diameter=mean_cell_diameter,
        community_detection="with_prior",
        prior_name="in_nucleus",
        seed=None,
    )
    Comsegdict.compute_community_vector()

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

