


#%%


import numpy as np
import scanpy as sc
from scipy.spatial.distance import cdist
import ssam


def select_genes_for_sct(vec = None,
                 genes = None,
                 min_expr = 0.01,
                 min_cell = 5):
    """
    Select the gene where it is possible to apply sctransform
    default value from original vst :https://github.com/satijalab/sctransform/blob/master/R/vst.R
    :param vec:
    :param analysis: need only if vec is None
    :param genes:
    :param min_expr: default 0.01
    :param min_cell: default 5
    :return:
    """
    if type(vec) != type(np.array([])):
        vec = vec.toarray() #sparse object for ex
    bool_index  = np.sum(vec > min_expr, axis=0) >= min_cell
    if bool_index.ndim == 2:
        bool_index = bool_index[0]
    return bool_index, np.array(genes)[bool_index]

def get_cluster_centroid(anndata, sct_normalize = True, n_comps=20,
                         _neighbors=15,  n_pcs=20,
                         resolution=0.5, plot_umap = True,
                         min_rna = 10):


    selected_genes = list(anndata.var['features'])
    anndata = anndata[list(np.array(np.sum(anndata.X, axis=1)).reshape(-1) > min_rna), :]
    if sct_normalize:
        bool_index, selected_genes = select_genes_for_sct(
            vec=anndata.X,
            genes=selected_genes,
            min_expr=0.01,
            min_cell=5)
        anndata = anndata[:, bool_index]
        norm_expression_vectors, param_sctransform = ssam.run_sctransform(anndata.X,
                                                                          debug_path=None)
        anndata.X = norm_expression_vectors

    sc.tl.pca(anndata, svd_solver='arpack', n_comps=n_comps)
    sc.pp.neighbors(anndata, n_neighbors=_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(anndata, resolution=resolution)
    if plot_umap:
        sc.tl.umap(anndata)
        sc.pl.umap(anndata, color=['leiden'], palette=None, legend_loc='on data')
    list_centroids = []
    unique_cluster = np.unique(anndata.obs["leiden"])
    for cluster in unique_cluster:
        list_centroids.append(np.mean(anndata[anndata.obs["leiden"] == cluster, :].X, axis=0))
        print(cluster, np.sum(anndata.obs["leiden"] == cluster))

    return list_centroids, bool_index, selected_genes


def dist_to_centroid(anndata,
                     list_centroids,
                     selected_genes,
                     sct_normalize = True,
                     sampling  = True,
                     sampling_size = 100000,
                     min_rna = 3):


    selected_genes = list(selected_genes)
    anndata = anndata[:, selected_genes]
    anndata = anndata[list(np.sum(anndata.X, axis=1) > min_rna), :]
    if sct_normalize:
        if sampling: # sampling
            count_matrix = anndata.X
            count_matrix = count_matrix[np.random.choice(len(count_matrix), sampling_size, replace=False)]
            bool_index, new_selected_genes = select_genes_for_sct(
                vec=count_matrix,
                genes=selected_genes,
                min_expr=0.01,
                min_cell=5)
            count_matrix = count_matrix[:, bool_index]
            anndata = anndata[:, bool_index]
            _, param_sctransform = ssam.run_sctransform(count_matrix,
                                                          debug_path=None)

            norm_expression_vectors = sctransform_from_parameters(
                np.array(param_sctransform),
                anndata.X)
            anndata.X = norm_expression_vectors
        else:
            bool_index, new_selected_genes = select_genes_for_sct(
                vec=anndata.X,
                genes=selected_genes,
                min_expr=0.01,
                min_cell=5)
            anndata = anndata[:, bool_index]
            norm_expression_vectors, param_sctransform = ssam.run_sctransform(anndata.X,
                                                                              debug_path=None)
            anndata.X = norm_expression_vectors
    else:
        bool_index = [True] * len(selected_genes)


####################""

    list_centroids = np.array(list_centroids)[~np.isnan(np.sum(list_centroids, axis=1))]
    list_centroids = list_centroids[:, bool_index]
    matrix_dist = cdist(anndata.X, list_centroids, metric="cosine")
    list_dict_to_centroid = np.min(matrix_dist, axis=1)
    anndata.obs['dist_to_centroid'] = list_dict_to_centroid
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4))
    # plot the cumulative histogram
    values, bins, patches = ax.hist(list_dict_to_centroid, range=(0, 2), bins=np.arange(0, 2, 0.01), density=True, histtype='step',
                                    cumulative=True, label='Empirical')
    plt.show()
    area = sum(np.diff(bins) * values)
    print("area ", sum(np.diff(bins) * values))
    return area,  values, bins, list_dict_to_centroid





