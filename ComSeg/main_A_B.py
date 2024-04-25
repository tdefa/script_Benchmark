



import numpy as np
import random
from comseg import dataset

import comseg


from pathlib import Path
import datetime
import argparse

if __name__ == '__main__':
#%%


    e = datetime.datetime.now()
    date_str = f"{e.month}_d{e.day}_h{e.hour}_min{e.minute}_s{e.second}_r" + str(random.randint(0, 5000))
    for nb_B in  [10, 100]:

        parser = argparse.ArgumentParser(description='test')
        ###################
        ###### hyper-parameter
        ###################
        parser.add_argument("--max_cell_radius", type=float,
                            default=15)  # TODO
        parser.add_argument("--mean_cell_diameter", type=float,
                            default=10)  # TODO

        parser.add_argument("--path_dataset_folder", type=str,
                            default=f"/media/tom/T7/regular_grid/simu1912/cube2D_step100/dataframe_folder/A_B/A100_B{nb_B}/")  # TODO
        parser.add_argument("--path_to_mask_prior", type=str,
                            default="/media/tom/T7/regular_grid/simu1912/cube2D_step100/nuclei/")  # TODO
        parser.add_argument("--path_dict_cell_centroid", type=str,
                            default= "/media/tom/T7/regular_grid/simu1912/cube2D_step100/dico_centroid/")
        parser.add_argument("--path_simulation_gr", type=str,
                            default=f'/media/tom/T7/regular_grid/simu1912/cube2D_step100/dico_simulation_A100_B{nb_B}/')  # TODO


        parser.add_argument("--port", default=3950)
        parser.add_argument("--mode", default='client')
        parser.add_argument("--host", default='127.0.0.2')
        args = parser.parse_args()

        #### HYPERPARAMETER ####
        max_cell_radius = args.max_cell_radius
        mean_cell_diameter = args.mean_cell_diameter



        path_dataset_folder = args.path_dataset_folder
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
                                                mean_cell_diameter = mean_cell_diameter
                                    )

        ### add prior knowledge, here using nucleus segmentation mask
        dataset_non_conv.add_prior_from_mask(prior_keys_name='in_nucleus', overwrite=True)


        ### compute the co-expression correlation at the dataset scale
        if "dict_co_expression" not in dataset_non_conv.__dict__:
            dico_proba_edge, count_matrix = dataset_non_conv.compute_edge_weight(  # in micrometer
                images_subset=None,
                distance="pearson",
            )

        import comseg
        from comseg import model
        from comseg import dictionary
        import importlib
        from comseg import clustering


        Comsegdict = dictionary.ComSegDict(
            dataset=dataset_non_conv,
                     mean_cell_diameter= mean_cell_diameter,
                     community_detection="with_prior",
                    # weights_name="weight",
                     prior_name="in_nucleus",
                     seed=None,
        )
        Comsegdict.compute_community_vector()




        #############""


        Comsegdict.compute_insitu_clustering(
            size_commu_min=3,
            norm_vector=False,
            ### parameter clustering
            n_pcs=0,
            n_comps=0,
            clustering_method="leiden",
            n_neighbors=30,
            resolution=1,
            palette=None,
            nb_min_cluster=1,
            min_merge_correlation=0.8,
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
        from matplotlib import pyplot as plt
        image_name = 'mask_cyto0'
        nuclei = np.load(path_to_mask_prior + image_name   + ".npy")
        from comseg.utils.plot import plot_result

        G = Comsegdict[image_name].G
        plot_result(nuclei,
                    G,
                    key_node='index_commu',
                    title=None,
                    dico_cell_color=None,
                    figsize=(15, 15),
                    spots_size=5,
                    plot_outlier=True)
        plt.show()

        Comsegdict.classify_centroid(
                        path_dict_cell_centroid = path_dict_cell_centroid,
                        file_extension = ".npy")

        Comsegdict.associate_rna2landmark(
            key_pred="leiden_merged",
            prior_name='in_nucleus',
            distance='distance',
            max_cell_radius=max_cell_radius)


        Comsegdict.anndata_from_comseg_result(
                                       key_cell_pred='cell_index_pred',
                                       )


        image_name = 'mask_cyto0'
        nuclei = np.load(path_to_mask_prior + image_name   + ".npy")
        from comseg.utils.plot import plot_result

        G = Comsegdict[image_name].G
        plot_result(nuclei,
                    G,
                    key_node='cell_index_pred',
                    title=None,
                    dico_cell_color=None,
                    figsize=(15, 15),
                    spots_size=5,
                    plot_outlier=True)
        plt.show()

















