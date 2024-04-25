


#%%


import anndata as ad

import ssam
from utils.count_matrix_metrics import dist_to_centroid
from pathlib import Path
import numpy as np
import comseg
import scanpy


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




    ref_pwc_13 = scanpy.read("/media/tom/T7/sp_data/In_Situ_Sequencing/single_cell/andata_pwc13.h5ad")

    anndata = ref_pwc_13
    sct_normalize = True
    n_comps=20
    _neighbors=20
    n_pcs=20
    resolution=0.5
    plot_umap = True


    bool_index, new_selected_genes = select_genes_for_sct(
        vec=anndata.X,
        genes=list(anndata.var['features']),
        min_expr=0.01,
        min_cell=5)
    selected_genes_centroid = new_selected_genes
    anndata = anndata[:, bool_index]

    norm_expression_vectors, param_sctransform = ssam.run_sctransform(anndata.X,
                                                                          debug_path=None)



    ## take only the genes that were use in the paper
    anndata_pciseq = scanpy.read("/media/tom/T7/sp_data/In_situ_Sequencing_16/single_cell/pciseq_from_paper.h5ad")
    selected_genes = list(anndata_pciseq.var["features"])
    selected_genes2 = list(ref_pwc_13.var["features"])
    selected_genes = list(set(selected_genes).intersection(set(selected_genes2)))
    selected_genes_centroid = selected_genes






    list_centroids = []
    unique_cluster = np.unique(anndata.obs["seurat_clusters"])
    for cluster in unique_cluster:
        list_centroids.append(np.median(anndata[anndata.obs["seurat_clusters"] == cluster, selected_genes_centroid].X, axis=0))
        print(cluster, np.sum(anndata.obs["seurat_clusters"] == cluster))


    ########### pciseq

    anndata_pciseq = scanpy.read("/media/tom/T7/sp_data/In_situ_Sequencing_16/single_cell/pciseq_from_paper.h5ad")
    print(anndata_pciseq[np.sum(anndata_pciseq.X, axis = 1) > 5, selected_genes_centroid])
    area_pci =  dist_to_centroid(anndata = anndata_pciseq,
                             list_centroids = list_centroids,
                             selected_genes  = selected_genes_centroid,
                             sct_normalize = True,
                             sampling_size = 50000, min_rna = 5)
    set(anndata_pciseq.var["features"])



    ########### baysor
    #adata_baysor = scanpy.read("/home/tom/Bureau/phd/simulation/baysor/baysor_csv_input/baysor_iss_scale_15_pwc13_138_priordefault_csv_cluster/all_adata.h5ad")
    adata_baysor = scanpy.read("/home/tom/Bureau/phd/simulation/baysor/baysor_csv_input/cluster_lung_tissue_mask2D/count_matrix_pred_all.h5ad")
    print(adata_baysor[np.sum(adata_baysor.X, axis = 1) > 5, selected_genes_centroid])
    area_baysor  =  dist_to_centroid(anndata = adata_baysor,
                             list_centroids = list_centroids,
                             selected_genes  = selected_genes_centroid,
                             sct_normalize = True,
                             sampling_size = 50000, min_rna = 5)

    ###### watershed
    anndata_watershed = scanpy.read("/media/tom/T7/sp_data/In_situ_Sequencing_16/anndata/watershed20_size16_10_10_v1.h5ad")
    print(anndata_watershed[np.sum(anndata_watershed.X, axis = 1) > 5, selected_genes_centroid])
    area_w = dist_to_centroid(anndata = anndata_watershed,
                     list_centroids = list_centroids,
                     selected_genes = selected_genes_centroid,
                     sct_normalize = True,
                         sampling  = False,
                     sampling_size = 50000,
                              min_rna=5)







    anndata_comseg0  = scanpy.read("/media/tom/T7/sp_data/In_situ_Sequencing_16/comsegres/resolution05/r45/1/adata.11_d9_h19_min45_s38_r2984.h5ad")
    anndata_comseg1  = scanpy.read("/media/tom/T7/sp_data/In_situ_Sequencing_16/comsegres/resolution05/r45/2/adata.11_d9_h19_min45_s38_r2419.h5ad")#r35 R05
    anndata_comseg2  = scanpy.read("/media/tom/T7/sp_data/In_situ_Sequencing_16/comsegres/resolution05/r45/3/adata.11_d9_h19_min45_s38_r1946.h5ad") #r35 R05
    anndata_comseg3  = scanpy.read("/media/tom/T7/sp_data/In_situ_Sequencing_16/comsegres/resolution05/r45/4/adata.11_d9_h19_min45_s38_r3796.h5ad")  #r35 R05
    anndata_comseg4  = scanpy.read("/media/tom/T7/sp_data/In_situ_Sequencing_16/comsegres/resolution05/r45/5/adata.11_d9_h19_min45_s38_r2532.h5ad")


    X  = np.concatenate([anndata_comseg0.X, anndata_comseg1.X, anndata_comseg2.X, anndata_comseg3.X, anndata_comseg4.X])
    adata_full_comseg = ad.AnnData(X)
    adata_full_comseg.var["features"] = anndata_comseg4.var["features"]
    adata_full_comseg.var_names = anndata_comseg4.var_names
    print(adata_full_comseg)
    print(adata_full_comseg[np.sum(adata_full_comseg.X, axis = 1) > 5,   selected_genes_centroid])
    selected_genes = list(selected_genes)
    min_rna = 5
    adata_full_comseg = adata_full_comseg[:, selected_genes]
    adata_full_comseg = adata_full_comseg[list(np.sum(adata_full_comseg.X, axis=1) > min_rna), :]
    print(adata_full_comseg)
    area_comseg = dist_to_centroid(anndata =adata_full_comseg,
                     list_centroids = list_centroids,
                     selected_genes = selected_genes_centroid,
                         sct_normalize = True,
                         sampling  = False,
                         sampling_size = 50000,
                         min_rna = 5)





    from matplotlib import pyplot
    values_comseg, bins, patches  = pyplot.hist(area_comseg[-1], area_comseg[-2], alpha=0.25, label='comseg', cumulative=True, color = "blue")
    values_baysor, bins, patches  = pyplot.hist(area_baysor[-1], area_baysor[-2], alpha=0.25, label='baysor', cumulative=True, color = "red")
    values_watershed, bins, patches  = pyplot.hist(area_w[-1], area_w[-2], alpha=0.25, label='baysor', cumulative=True, color = "red")
    values_pciSeq, bins, patches  = pyplot.hist(area_pci[-1], area_pci[-2], alpha=0.25, label='baysor', cumulative=True, color = "red")
    pyplot.hist(area_pci[-1], area_pci[-2], alpha=0.25, label='pciSeq',  cumulative=True, color = "green")
    pyplot.legend(loc='upper right')
    pyplot.show()


    import numpy as np
    import matplotlib.pyplot as plt

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

    # Adding Xticks
    #plt.xlabel('cosine distance threshold', fontweight='bold', fontsize=25)
    #plt.ylabel('number of cell matched with scRNA-seq', fontweight='bold', fontsize=25)
    plt.xticks([r + barWidth for r in range(len(comseg_bar))],
               ['0.25', '0.50', '0.75'], fontsize=55)
    #plt.yticks(fontsize=20)
    #plt.legend(fontsize=20)
    #plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off')
    for spline_v in plt.gca().spines.values():
        spline_v.set_visible(False)
    plt.show()
    image_format = 'svg'  # e.g .png, .svg, etc.
    file_name = Path(folder_save) / 'matching_lung_humain0911.svg'
    fig.savefig(file_name, format=image_format, dpi=400, bbox_inches='tight')
    image_format = 'png'  # e.g .png, .svg, etc.
    file_name = Path(folder_save) / 'matching_lung_humain0911.png'
    fig.savefig(file_name, format=image_format, dpi=400, bbox_inches='tight')

