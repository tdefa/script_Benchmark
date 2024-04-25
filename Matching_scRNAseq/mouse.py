




import list_markers_hardcoded
import scanpy as sc
import ssam
from utils.count_matrix_metrics import dist_to_centroid
from comseg.dictionary import ComSegDict
from pathlib import Path
import numpy as np

from utils.data_processing import sctransform_from_parameters


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



if __name__ == '__main__':



    selected_genes = list_markers_hardcoded.oligo_pool2_6g
    anndata_seq = sc.read(
        "/home/tom/Bureau/phd/markers_selection/data/20220321_lung_merge_27samples_raw_selected_with_subtype.h5ad")
    anndata_seq = anndata_seq[anndata_seq.obs.condition == 'IR_17Gy_5M', selected_genes]
    anndata_seq.X = anndata_seq.X.toarray()
    anndata_seq = anndata_seq[np.sum(anndata_seq.X, axis=1) > 5, :]
    count_matrix = anndata_seq.X
    bool_index, new_selected_genes = select_genes_for_sct(
        vec=count_matrix,
        genes=list(anndata_seq.var_names),
        min_expr=0.01,
        min_cell=5)

    _, param_sctransform = ssam.run_sctransform(count_matrix.toarray(),
                                                debug_path=None)


    anndata_seq.X = anndata_seq.X.toarray()
    norm_expression_vectors = sctransform_from_parameters(
        np.array(param_sctransform),
        anndata_seq.X.toarray())
    anndata_seq.X = norm_expression_vectors
    anndata = anndata_seq.copy()
    anndata.X = anndata.X.toarray()

    plot_umap = True
    n_comps = 0
    n_pcs = 0
    n_comps = 0
    resolution = 0.1
    _neighbors = 30
    # sc.tl.pca(anndata, svd_solver='arpack', n_comps=n_comps)
    sc.pp.neighbors(anndata, n_neighbors=_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(anndata, resolution=resolution, random_state=5)
    list_centroids = []
    unique_cluster = np.unique(anndata.obs["leiden"])
    for cluster in unique_cluster:
        list_centroids.append(np.median(anndata[anndata.obs["leiden"] == cluster, :].X, axis=0))
        print(cluster, np.sum(anndata.obs["leiden"] == cluster))

    selected_genes = list_markers_hardcoded.oligo_pool2_6g


####################### comseg

    model = ComSegDict()

    path_parameter = '/media/tom/T7/lustra_all/2023-09-06_LUSTRA-14rounds/image_per_pos/dataframe_ComSeg/results/11_d10_h16_min15_s31_r3763/Comsegdict.11_d10_h16_min15_s31_r3763pickle.txt'
    model.load(path_parameter)
    model.anndata_from_comseg_result()
    adata_comseg = model.final_anndata
    print(adata_comseg[np.sum(adata_comseg.X, axis=1) > 5, :])


    area_comseg =  dist_to_centroid(anndata = adata_comseg,
                             list_centroids = list_centroids,
                             selected_genes  = selected_genes,
                             sct_normalize = True,
                            sampling=False,
                            sampling_size = 50000,
                            min_rna = 5)


    ######## baysor


    adata_baysor10 = sc.read_h5ad("/home/tom/Bureau/phd/simulation/baysor/baysor_csv_input/max3_50_nb_cluster5_lustra_with_mask/count_matrix_pred_all.h5ad")
    print(adata_baysor10[np.sum(adata_baysor10.X, axis=1) > 5, :])

    area_baysor10 = dist_to_centroid(anndata=adata_baysor10,
                                     list_centroids=list_centroids,
                                     selected_genes=selected_genes,
                                     sct_normalize=True,
                                     sampling=False,
                                     sampling_size=50000, min_rna=5)





    ####### watershed


    adata_watershed = sc.read_h5ad("/media/tom/T7/lustra_all/2023-09-06_LUSTRA-14rounds/image_per_pos/watershed/anndata/stich0_max6_20_20p.h5ad")
    print(adata_watershed[np.sum(adata_watershed.X, axis = 1) > 5, :])
    area_watershed  =  dist_to_centroid(anndata = adata_watershed,
                             list_centroids = list_centroids,
                             selected_genes  = selected_genes,
                             sct_normalize = True,
                                      sampling=False,
                            sampling_size = 50000, min_rna = 5)

    adata_pciSeq = sc.read_h5ad("/media/tom/T7/lustra_all/2023-09-06_LUSTRA-14rounds/image_per_pos/pciseq_res/14juin_pciseq_result_231023.h5ad")
    print(adata_pciSeq[np.sum(adata_pciSeq.X, axis = 1) > 5, :])
    adata_pciSeq  =  dist_to_centroid(anndata = adata_pciSeq,
                             list_centroids = list_centroids,
                             selected_genes  = selected_genes,
                             sct_normalize = True,
                                      sampling=False,
                            sampling_size = 50000, min_rna = 5)


    from matplotlib import pyplot
    bins = area_comseg[-1]
    values_comseg, bins, patches = pyplot.hist(area_comseg[-1], area_comseg[-2], alpha=0.25, label='comseg', cumulative=True, histtype='step', color = "blue")
    values_baysor, bins, patches =pyplot.hist(area_baysor10[-1], area_baysor10[-2], alpha=0.25, label='baysor', cumulative=True, histtype='step', color = "red")
    values_watershed, bins, patches =pyplot.hist(area_watershed[-1], area_watershed[-2], alpha=0.25, label='watershed', cumulative=True, histtype='step', color = "red")
    values_pciSeq, bins, patches =pyplot.hist(adata_pciSeq[-1], adata_pciSeq[-2], alpha=0.25, label='watershed', cumulative=True, histtype='step', color = "green")
    pyplot.legend(loc='upper right')
    pyplot.show()


    import numpy as np
    import matplotlib.pyplot as plt

    # set width of bar

    # set width of bar
    barWidth = 0.20

    # set height of bar




    comseg_bar = [values_comseg[25], values_comseg[50], values_comseg[75]]
    baysor_bar =  [values_baysor[25], values_baysor[50], values_baysor[75]]
    watershed_bar =  [values_watershed[25], values_watershed[50], values_watershed[75]]
    pciSeq_bar =  [values_pciSeq[25], values_pciSeq[50], values_pciSeq[75]]

    folder_save = '/home/tom/Bureau/phd/simulation/paper_result/figure/'
    plt.rcParams['xtick.labelsize'] = 55
    plt.rcParams['ytick.labelsize'] = 55
    plt.rc_context({"axes.labelsize": 55, })

    barWidth = 0.20
    #fig, ax = plt.subplots(figsize=(18, 15))
    plt.rcParams['xtick.labelsize'] = 55
    plt.rcParams['ytick.labelsize'] = 55
    plt.rc_context({"axes.labelsize": 55, })

    barWidth = 0.20
    fig, ax = plt.subplots(figsize=(25, 15))
    # Set position of bar on X axis
    br1 = np.arange(len(comseg_bar))
    br2 = [x + barWidth for x in br1]
    br3 = [x + barWidth for x in br2]
    br4 = [x + barWidth for x in br3]

    # Make the plot
    plt.bar(br1, comseg_bar, color="#386cb0", width=barWidth,
            edgecolor='grey', label='comseg')
    plt.bar(br2, watershed_bar, color='#fdc086', width=barWidth,
            edgecolor='grey', label='watershed')
    plt.bar(br3, pciSeq_bar, color='#7fc97f', width=barWidth,
            edgecolor='grey', label='pciSeq')
    plt.bar(br4, baysor_bar, color='#beaed4', width=barWidth,
            edgecolor='grey', label='baysor')

    plt.xticks([r + barWidth for r in range(len(comseg_bar))],
               ['0.25', '0.50', '0.75'], fontsize=55)
    for spline_v in plt.gca().spines.values():
        spline_v.set_visible(False)
    plt.show()

    image_format = 'svg'  # e.g .png, .svg, etc.
    file_name = Path(folder_save) / 'matching_lung_mice1011_r20_bis.svg'
    fig.savefig(file_name, format=image_format, dpi=400, bbox_inches='tight')
    image_format = 'png'  # e.g .png, .svg, etc.
    file_name = Path(folder_save) / 'matching_lung_mice1011_r20_bis.png'
    fig.savefig(file_name, format=image_format, dpi=400, bbox_inches='tight')
