


####################"
# Here is an example of how we benchmarked pciSeq
###################

from pathlib import Path
import scipy
from scipy import sparse
import anndata as ad
import pciSeq
import pandas as pd
import numpy as np
import scanpy

if __name__ == '__main__':

    path_save = "/home/tom/Bureau/phd/simulation/paper_result/pcsiseq/IRM17_MOUSE/"
    ### for 2D version
    mask_nuclei = np.load('/media/tom/T7/lustra_all/2023-09-06_LUSTRA-14rounds/image_per_pos/from_cluster/final_mask/14juin_final_masks.npy')
    mip_maks_nuclei = np.max(mask_nuclei, axis=0)
    unique_cell_id = np.unique(mip_maks_nuclei)[1:]
    mip_maks_nuclei = scipy.sparse.coo_matrix(mip_maks_nuclei)

    dico_dico_pred2gr = {}
    for i in range(len(unique_cell_id)):
        dico_dico_pred2gr[i+1] = unique_cell_id[i]
    np.save(path_save + "dico_dico_pred2gr.npy", dico_dico_pred2gr)



    dico_dico_commu = np.load(
        "/media/tom/T7/lustra_all/2023-09-06_LUSTRA-14rounds/image_per_pos/comseg_imput/14juin_dico_dico_commu_mask_removal_23j_6round.npy",
        allow_pickle=True).item()
    df_spots_label = dico_dico_commu['stich0']['df_spots_label']
    df_spots_label = df_spots_label.rename(columns={"gene": "Gene", "z": 'z_stack' })
    df_spots_label = df_spots_label[['y', 'x', 'z_stack', 'Gene']]

    df_spots_label.x = df_spots_label.x
    df_spots_label.y = df_spots_label.y
    df_spots_label.z_stack = df_spots_label.z_stack
    selected_genes = np.unique(df_spots_label.Gene)
    gene_index_dico = {}
    for gene_id in range(len(selected_genes)):
        gene_index_dico[selected_genes[gene_id]] = gene_id


    #### load scRNAseq
    #from pciSeq.src.diagnostics.utils import is_redis_running, validate_redis

    anndata = scanpy.read(
        "/home/tom/Bureau/phd/markers_selection/data/20220321_lung_merge_27samples_raw_selected_with_subtype.h5ad")[:, selected_genes]



    column_name = list(anndata.obs.cell_ID)
    last_letter = 4
    rows = anndata.X.T
    row_name = list(anndata.var['features'])
    df_scRNAseq = pd.DataFrame(data=rows.toarray(),
                               index=row_name,
                               columns=column_name)
    cellData, geneData = pciSeq.fit(df_spots_label, mip_maks_nuclei, df_scRNAseq)
    list_gene_name = list(cellData.Genenames)
    list_gene_count = list(cellData.CellGeneCount)
    geneData.to_csv(Path(path_save) / "geneData.csv")
    cellData.to_csv(Path(path_save) / "cellData.csv")
    list_pred_cell = list(geneData.neighbour)
    list_pred_cell_prob = np.array([i[0] for i in geneData.neighbour_prob])
    dico_cell_nb_rna = {}
    for cell_index in np.unique(list_pred_cell):
        dico_cell_nb_rna[cell_index] = 0

    for index_cell in range(len(list_pred_cell)):
        if  list_pred_cell_prob[index_cell] > 0.5:
            dico_cell_nb_rna[list_pred_cell[index_cell]] += 1
    selected_genes = np.unique(np.concatenate(list_gene_name))
    gene_index_dico = {}
    for gene_id in range(len(selected_genes)):
        gene_index_dico[selected_genes[gene_id]] = gene_id
    pci_count_matrix = np.zeros([len(cellData), len(gene_index_dico)])
    for indl in range(len(list_gene_name)):
        if indl not in dico_cell_nb_rna.keys():
            continue
        if dico_cell_nb_rna[indl] <= 5:
            continue
        expression_vector = np.zeros(len(gene_index_dico))
        for gene_index in range(len(list_gene_name[indl])):
            gene = list_gene_name[indl][gene_index]
            expression_vector[gene_index_dico[gene]] = list_gene_count[indl][gene_index]
        pci_count_matrix[indl] = expression_vector
    adata = ad.AnnData(pci_count_matrix)
    adata.var["features"] = selected_genes
    adata.var_names = selected_genes
    print(adata[np.sum(adata.X, axis=1) > 5, :])
