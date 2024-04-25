

import pickle
from pathlib import Path
import scanpy as sc
import random
import comseg
from comseg import model
from comseg import dictionary
import importlib

import matplotlib
# matplotlib.use('Qt5Agg')
import sys
import os
import psutil
import pandas as pd

import numpy as np
import random
import tifffile
from comseg.dataset import ComSegDataset
import comseg

import importlib

from pathlib import Path
from tqdm import tqdm
from comseg.utils.preprocessing import sctransform_from_parameters

import argparse
import datetime
import scanpy as sc
# importlib.reload(ComSegDataset)

from sklearn.utils.random import sample_without_replacement
from comseg.utils.preprocessing import run_sctransform
import seaborn as sns

import seaborn as sns
from matplotlib import pyplot as plt
import random
from memory_profiler import profile
import time





if __name__ == '__main__':

    time_start = time.time()

    e = datetime.datetime.now()

    date_str = f"{e.month}_d{e.day}_h{e.hour}_min{e.minute}_s{e.second}_r" + str(random.randint(0, 5000))
    parser = argparse.ArgumentParser(description='test')
    ###################
    ###### hyper-parameter
    ###################
    parser.add_argument("--max_cell_radius", type=float,
                        default=15)
    parser.add_argument("--mean_cell_diameter", type=float,
                        default=10)

    parser.add_argument("--path_dataset_folder", type=str,
                        default="/media/tom/Transcend/doi_10_5061_dryad_jm63xsjb2__v20210916/input_comseg_3Dseg_tom/df/")  # TODO
    parser.add_argument("--path_to_mask_prior", type=str,
                        default="/media/tom/Transcend/doi_10_5061_dryad_jm63xsjb2__v20210916/input_comseg_3Dseg_tom/prior/")  # TODO
    parser.add_argument("--path_dict_cell_centroid", type=str,
                        default="")




    parser.add_argument("--port", default=3950)
    parser.add_argument("--mode", default='client')
    parser.add_argument("--host", default='127.0.0.2')

    pid = os.getpid()
    # Create a Process object
    process = psutil.Process(pid)

    args = parser.parse_args()

    #### HYPERPARAMETER ####
    max_cell_radius = args.max_cell_radius
    mean_cell_diameter = args.mean_cell_diameter



    path_dataset_folder = args.path_dataset_folder
    ##path to your prior segmentation mask
    path_to_mask_prior = args.path_to_mask_prior

    path_dict_cell_centroid = args.path_dict_cell_centroid


    args.path_save = str(Path(args.path_dataset_folder) / ("results/" + date_str + "/"))
    Path(args.path_save).mkdir(parents=True, exist_ok=True)
    #### writie the argument to a file


    with open(Path(args.path_save) / "script_parameter.txt", "w") as f:
        for k, v in args.__dict__.items():
            f.write(f"{k} : {v}\n")


    ### define the gene you want to study, you can restrict it to few genes you want to study.


    ## scale/ pixel size in um
    dict_scale = {"x": 0.17, 'y': 0.17, "z": 1.5}

    ### create the dataset object
    dataset = comseg.dataset.ComSegDataset(
        path_dataset_folder=path_dataset_folder,
        path_to_mask_prior=path_to_mask_prior,
        dict_scale=dict_scale,
        mask_file_extension=".tif",
        mean_cell_diameter=mean_cell_diameter
    )
    cpu_usage = process.cpu_percent(interval=None)
    print(f"CPU usage: {cpu_usage} percent")
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss} bytes")
    ### add prior knowledge, here using nucleus segmentation mask
    dataset.add_prior_from_mask(prior_keys_name='in_nucleus',
                                         overwrite=True,
                                         compute_centroid = True) # compute cell centroid

    cpu_usage = process.cpu_percent(interval=None)
    print(f"CPU usage: {cpu_usage} percent")
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss} bytes")
    ### compute the co-expression correlation at the dataset scale
    try :
        assert False
        dict_co_expression = np.load(Path(args.path_to_mask_prior) /'dict_co_expression_n40_50000.npy', allow_pickle=True).item()
        dataset.__dict__["dict_co_expression"] = dict_co_expression
    except Exception as e:
        print(e)
        print("in co_expression")
        print()



    if "dict_co_expression" not in dataset.__dict__:

        dico_proba_edge, count_matrix = dataset.compute_edge_weight(  # in micrometer
            images_subset=None,
            n_neighbors=40,
            sampling=True,
            sampling_size=10000
        )
        ########## plot corr
        cpu_usage = process.cpu_percent(interval=None)
        print(f"CPU usage: {cpu_usage} percent")
        mem_info = process.memory_info()
        print(f"Memory usage: {mem_info.rss} bytes")

    print(" coexpression added")
    corr_matrix = []

    np.save(Path(args.path_to_mask_prior) /'dict_co_expression_n40_50000.npy', dataset.dict_co_expression)

    Comsegdict = dictionary.ComSegDict(
        dataset=dataset,
        mean_cell_diameter=args.mean_cell_diameter,
        community_detection="with_prior",
        prior_name="in_nucleus",
    )

    cpu_usage = process.cpu_percent(interval=None)
    print(f"CPU usage: {cpu_usage} percent")
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss} bytes")

    Comsegdict.compute_community_vector()
    cpu_usage = process.cpu_percent(interval=None)
    print(f"CPU usage: {cpu_usage} percent")
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss} bytes")

    Comsegdict.compute_insitu_clustering(
        size_commu_min=3,
        norm_vector=True,
        ### parameter clustering
        n_pcs=3,
        n_comps=3,
        clustering_method="leiden",
        n_neighbors=20,
        resolution=1,
        n_clusters_kmeans=4,
        palette=None,
        nb_min_cluster=0,
        min_merge_correlation=0.8,
    )
    cpu_usage = process.cpu_percent(interval=None)
    print(f"CPU usage: {cpu_usage} percent")
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss} bytes")




    palette = {}
    for i in range(-1, 500):
        palette[str(i)] = "#" + "%06x" % random.randint(0, 0xFFFFFF)
    adata = Comsegdict.in_situ_clustering.anndata_cluster
    adata.obs["leiden_merged"] = adata.obs["leiden_merged"].astype(int)
    #sc.tl.umap(adata)
    #fig_ledien = sc.pl.umap(adata, color=["leiden_merged"], palette=palette, legend_loc='on data',
     #                       )

    cpu_usage = process.cpu_percent(interval=None)
    print(f"CPU usage: {cpu_usage} percent")
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss} bytes")

    Comsegdict.add_cluster_id_to_graph(clustering_method="leiden_merged")

    ### get a csv spot/cluster

    gene_list = []
    x_list = []
    y_list = []
    z_list = []
    leiden = []
    cell_id = []

    img_name = list(Comsegdict.keys())[0]
    for node in Comsegdict[img_name].G.nodes:
        gene_list.append(Comsegdict[img_name].G.nodes[node]["gene"])
        x_list.append(Comsegdict[img_name].G.nodes[node]["x"])
        y_list.append(Comsegdict[img_name].G.nodes[node]["y"])
        z_list.append(Comsegdict[img_name].G.nodes[node]["z"])
        leiden.append(Comsegdict[img_name].G.nodes[node]["leiden_merged"])

    dictio = {'gene': gene_list, 'x': x_list, 'y': y_list,  'z': z_list,
              "leiden": leiden}
    df = pd.DataFrame(dictio)

    df.to_csv(
        Path(args.path_save) / "leiden0.csv")

    Comsegdict.classify_centroid(
        path_dict_cell_centroid=None,
        n_neighbors=15,
        dict_in_pixel=True,
        max_dist_centroid=None,
        key_pred="leiden_merged",
        distance="ngb_distance_weights",
        convex_hull_centroid=True,
        file_extension=".tiff.npy"
    )

    cpu_usage = process.cpu_percent(interval=None)
    print(f"CPU usage: {cpu_usage} percent")
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss} bytes")

    Comsegdict.associate_rna2landmark(
        key_pred="leiden_merged",
        prior_name='in_nucleus',
        distance='distance',
        max_cell_radius=8)

    time_end = time.time() - time_start
    print(f"Time elapsed: {time_end} seconds")

    gene_list = []
    x_list = []
    y_list = []
    z_list = []
    leiden = []
    cell_index_pred_list = []

    img_name = list(Comsegdict.keys())[0]
    for node in Comsegdict[img_name].G.nodes:
        gene_list.append(Comsegdict[img_name].G.nodes[node]["gene"])
        x_list.append(Comsegdict[img_name].G.nodes[node]["x"])
        y_list.append(Comsegdict[img_name].G.nodes[node]["y"])
        z_list.append(Comsegdict[img_name].G.nodes[node]["z"])
        leiden.append(Comsegdict[img_name].G.nodes[node]["leiden_merged"])
        cell_index_pred_list.append(Comsegdict[img_name].G.nodes[node]["cell_index_pred"])

    dictio = {'gene': gene_list, 'x': x_list, 'y': y_list, 'z': z_list,
              "leiden": leiden, "cell": cell_index_pred_list}
    df = pd.DataFrame(dictio)
    df.to_csv(
        Path(args.path_save) / "cell0_r10_rmax8_small_p.csv")

    adata = Comsegdict.in_situ_clustering.anndata_cluster
    adata.obs["leiden_merged"] = adata.obs["leiden_merged"].astype(int)
    #sc.tl.umap(adata)
    #fig_ledien = sc.pl.umap(adata, color=["leiden_merged"], palette=palette, legend_loc='on data',
     #                       )
    ### vizulaize  point cloud with napari
    Comsegdict.anndata_from_comseg_result(
                               key_cell_pred='cell_index_pred',
                               )
    cpu_usage = process.cpu_percent(interval=None)
    print(f"CPU usage: {cpu_usage} percent")
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss} bytes")

    filename = Path(args.path_save) / "result.h5ad"
    with open(filename, 'wb') as handle:
        pickle.dump(Comsegdict.final_anndata, handle, -1)

    print('done')
