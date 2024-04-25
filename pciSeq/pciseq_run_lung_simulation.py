








####################
## Here we displaya example of how we benchmarked pciSeq
####################


#%%
from pathlib import Path
import os
import numpy as np
import pandas as pd
import skimage.color
import matplotlib.pyplot as plt
from scipy.sparse import load_npz, coo_matrix
import pciSeq
import scipy
from pathlib import Path
from sklearn import metrics
import tifffile
from utils.custom_metrics import  compute_AJI
from tqdm import tqdm
import list_markers_hardcoded
from utils.data_processing import sctransform_from_parameters

####
#%%



path_dico_simulation = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/dico_simulation/"
path_nuclei = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/remove20/nuclei/"
path_cyto = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/individual_cytoplasms_rd/"
path_save_res = "/home/tom/Bureau/phd/simulation/paper_result/pcsiseq/lung_simu_remove20/max_filter3/res"
Path(path_save_res).mkdir(exist_ok=True, parents= True)

path_df_folder = "/media/tom/T7/simulation/exp_same_cyto/same_param1_4_0/dataframe_cluster_filter_max3/"
dico_dico_commu = {}
for path_df in Path(path_df_folder).glob("*.csv"):
    dico_dico_commu[path_df.stem + ".tiff.npy"] = {}
    dico_dico_commu[path_df.stem + ".tiff.npy" ]["df_spots_label"] = pd.read_csv(path_df)

# single cell data
import scanpy

anndata = scanpy.read("/home/tom/Bureau/phd/markers_selection/data/20220321_lung_merge_27samples_raw_selected_with_subtype.h5ad")



do_pred_cosine = True
scrna_centroids = np.load(
    "/home/tom/Bureau/phd/simulation/mycode/centroid/scrna_centroidsnsforest0_thalassadico_simulationnornmalized_True.npy",
    allow_pickle=True)
scrna_unique_clusters = np.load(
    "/home/tom/Bureau/phd/simulation/mycode/centroid/scrna_unique_clustersnsforest0_thalassadico_simulation_nsforest0_thalassa.npy",
    allow_pickle=True)
param_sctransform_anndata = np.load(
    "/home/tom/Bureau/phd/simulation/mycode/centroid/param_sctransform_anndatansforest0_thalassa.npy",
    allow_pickle=True)

gene_list  = np.unique(list_markers_hardcoded.nsforest0_thalassa)
gene_index_dico = {}
for gene_id in range(len(gene_list)):
    gene_index_dico[gene_list[gene_id]] = gene_id


###ns forest
column_name = list(anndata.obs.cell_ID)


last_letter = 4


rows = anndata.X.T
row_name = list(anndata.var['features'])
df_scRNAseq = pd.DataFrame(data=rows.toarray(),
                           index=row_name,
                           columns=column_name)


#### load nuclei segmentation  array under the form of scipy.sparse.coo.coo_matrix


total_y_true = []
total_y_pred = []
total_y_pred_cosine = []
total_y_pred_cosine_int = []

total_aji = []
show_plot = False
total_error_list = []
total_prc_not_catch_list = []
list_iou_from_aji = []
total_list_iou = []

dico_dico_pred2gr = {}



for path_to_simu in tqdm(list(Path(path_dico_simulation).glob("*.npy"))):
    image_name = str(path_to_simu).split("/")[-1]
    print(image_name)
    if image_name not in dico_dico_commu:
        continue
    try:
        #mask_nuclei = np.load(path_nuclei+image_name)
        mask_nuclei = tifffile.imread(path_nuclei  + image_name[:-last_letter])
    except FileNotFoundError as e:
        print(e)
        continue
    #mask_nuclei = tifffile.imread(path_nuclei + "dapi_maskdapi_" +image_name[:-last_letter])
    cyto =    np.load(path_cyto+image_name)



    mask_nuclei = np.amax(mask_nuclei, 0)
    unique_nuclei_raw = np.unique(mask_nuclei)
    unique_resize_nuc = np.unique(mask_nuclei)
    for nuc in unique_nuclei_raw:
        if nuc not in unique_resize_nuc:
            mask_nuclei[mask_nuclei == nuc] = 0
    #mask_nuclei = np.amax(mask_nuclei, 0)
    unique_cell_id = np.unique(mask_nuclei)[1:]
    dico_psc_id_mask_id = {}
    for i in range(len(unique_cell_id)):
        dico_psc_id_mask_id[i+1] = unique_cell_id[i]
    dico_dico_pred2gr[image_name] = dico_psc_id_mask_id
    mask_nuclei = scipy.sparse.coo_matrix(mask_nuclei)
    #mask_nuclei = [coo_matrix(d) for d in mask_nuclei]
    if show_plot:
        rgb_label_image = skimage.color.label2rgb(mask_nuclei.toarray(), bg_label=0)
        _dpi = 72
        plt.figure(figsize=(1500/_dpi, 1500/_dpi), dpi=_dpi)
        imgplot = plt.imshow(rgb_label_image)
        plt.show()

    #### load spot dataframe

    df_spots_label = dico_dico_commu[image_name]['df_spots_label']
    df_spots_label = df_spots_label.rename(columns={"gene": "Gene", 'z': 'z_stack'})

    #### try with scaling

    if show_plot:
        my_dpi = 72
        plt.figure(figsize=(1500/my_dpi, 1500/my_dpi), dpi=my_dpi)
        imgplot = plt.imshow(rgb_label_image)
        plt.scatter(df_spots_label.x, df_spots_label.y, s=1)
        plt.show()
    pciSeq.attach_to_log()

    #cellData, geneData = pciSeq.fit(df_spots_label, mask_nuclei, df_scRNAseq)

    df_spots_label_input = df_spots_label[['y', 'x', 'z_stack', 'Gene']].copy()

    cellData, geneData = pciSeq.fit(df_spots_label, mask_nuclei, df_scRNAseq)

    ### compute cell type calling accuracy
    dico_simulation = np.load(path_to_simu, allow_pickle = True).item()
    dico_cell_index = dico_simulation['dico_cell_index']
    list_cell_id = list(cellData.Cell_Num)
    list_cell_type = list(cellData.ClassName)
    list_cell_type_max_probe = [np.argmax(i) for i in  list(cellData.Prob)]
    y_true = []
    y_pred = []
    y_pred_cosine = []
    y_pred_cosine_int = []

    for cell_id_index in dico_psc_id_mask_id:
        if unique_cell_id[cell_id_index-1] not in dico_cell_index:
            continue
        if list_cell_id[cell_id_index-1] == 1:
            if 1 not in dico_cell_index:
                print("problem index start to 2")
                continue
        y_true.append(dico_cell_index[dico_psc_id_mask_id[cell_id_index]]['type'])
        y_pred.append(list_cell_type[cell_id_index-1][list_cell_type_max_probe[cell_id_index-1]])



        if do_pred_cosine:
            norm_expression_vectors = np.zeros((len(gene_list)))
            for gene_name_index in range(len(list(cellData['Genenames'])[cell_id_index-1])):
                gene_name = list(cellData['Genenames'])[cell_id_index-1][gene_name_index]
                norm_expression_vectors[gene_index_dico[gene_name]] = list(cellData['CellGeneCount'])[cell_id_index-1][gene_name_index]
            norm_expression_vectors = norm_expression_vectors.reshape(1, -1)

            norm_expression_vectors = sctransform_from_parameters(
                np.array(param_sctransform_anndata),
                norm_expression_vectors)

            correlation_array = metrics.pairwise.cosine_similarity(norm_expression_vectors,
                                                                   np.array(scrna_centroids))[0]

            index_cluster_max = np.argmax(correlation_array)
            pred_rna_seq = scrna_unique_clusters[index_cluster_max]
            y_pred_cosine.append(pred_rna_seq)
        if  True:   #do_use_count_int
            expression_vectors_gene = geneData.Gene[geneData.neighbour == cell_id_index]

            expression_vector = np.bincount([gene_index_dico[gene] for gene in expression_vectors_gene],
                        minlength=len(gene_index_dico))

            expression_vector = expression_vector.reshape(1, -1)
            if np.sum(expression_vector) == 0:
                y_pred_cosine_int.append('Zero')
                continue
            norm_expression_vectors = sctransform_from_parameters(
                np.array(param_sctransform_anndata),
                expression_vector)
            correlation_array = metrics.pairwise.cosine_similarity(norm_expression_vectors,
                                                                   np.array(scrna_centroids))[0]
            index_cluster_max = np.argmax(correlation_array)
            pred_rna_seq = scrna_unique_clusters[index_cluster_max]
            y_pred_cosine_int.append(pred_rna_seq)

    acc = metrics.accuracy_score(y_true = y_true, y_pred=y_pred)
    acc_cosine = metrics.accuracy_score(y_true = y_true, y_pred=y_pred_cosine)
    acc_cosine_int = metrics.accuracy_score(y_true = y_true, y_pred=y_pred_cosine_int)

    print(f' {image_name} acc {acc}')
    print(f' {image_name} acc_cosine {acc_cosine}')
    print(f' {image_name} acc_cosine {acc_cosine_int}')


    total_y_true += y_true
    total_y_pred += y_pred
    total_y_pred_cosine += y_pred_cosine
    total_y_pred_cosine_int += y_pred_cosine_int

    ### gene count matrix
    assert list(df_spots_label.x) == list(geneData.x)
    list_pred_cell = list(geneData.neighbour)
    list_cell_gr = list(df_spots_label.cell)
    df_spots_label['pciseq_pred'] = list_pred_cell
    ari = metrics.adjusted_rand_score(labels_true=list_cell_gr,
                                      labels_pred=list_pred_cell)
    ami = metrics.adjusted_mutual_info_score(labels_true=list_cell_gr,
                                             labels_pred=list_pred_cell)
    vm = metrics.v_measure_score(labels_true=list_cell_gr,
                                 labels_pred=list_pred_cell)
    print(f"molecule  cell type ari {round(ari, 4)} ami {round(ami, 4)}  vm {round(vm, 4)}")
    image_pred_cm = {}
    image_gr_cm = {}

    for cell in list(np.unique(list_cell_gr)) + list(np.unique(list_pred_cell)):
        image_gr_cm[cell] = []
    for cell in  list(np.unique(list_pred_cell)):
        image_pred_cm[cell] = []

    for mol_ind in range(len(list_pred_cell)):
        image_gr_cm[list_cell_gr[mol_ind]].append(mol_ind)
        image_pred_cm[list_pred_cell[mol_ind]].append(mol_ind)
    aji,  C_image, U_image, dico_cell_pred_rna, list_iou, dico_matching = compute_AJI(image_pred_cm, image_gr_cm,
                                                                       max_index=None, cell_index_backgroud =  0)
    total_aji.append(aji)
    list_iou_from_aji.append(list_iou)
    print(f' aji {aji}')


    error_list = []
    prc_not_catch_list = []
    list_iou2 = []
    for nuc_index in image_pred_cm: ## this one is
        nuc_index
        try:
            iou = len(
                set(list(image_pred_cm[nuc_index])).intersection(set(list(image_gr_cm[dico_psc_id_mask_id[nuc_index]])))) / len(
                set(list(image_pred_cm[nuc_index])).union(set(list(image_gr_cm[dico_psc_id_mask_id[nuc_index]]))))
            list_iou2.append(iou)
        except Exception as e:
            print(e)

        try:
            prc_false = len(
                set(image_pred_cm[nuc_index]) - set(image_gr_cm[dico_psc_id_mask_id[nuc_index]])) / len(
                set(image_pred_cm[nuc_index]))
            error_list.append(prc_false)
        except Exception as e:
            print(e)
        try:
            prc_not_catch = len(
                set(image_gr_cm[dico_psc_id_mask_id[nuc_index]]) - set(image_pred_cm[nuc_index])) / len(
                image_gr_cm[dico_psc_id_mask_id[nuc_index]])
            prc_not_catch_list.append(prc_not_catch)
        except Exception as e:
            print(e)
    total_list_iou += list_iou2
    total_error_list += error_list
    total_prc_not_catch_list += prc_not_catch_list

    ### save df for ploting

    df_spots_label.to_csv(path_save_res + 'df_spots_label_' + image_name +".csv" )
    geneData.to_csv(path_save_res + 'geneData' + image_name )
    geneData.to_csv(path_save_res + 'geneData' + image_name )

    acc = metrics.accuracy_score(y_true=total_y_true, y_pred=total_y_pred)
    print(f'final acc {round(acc, 4)}')
    #print(f'mean list_iou_from_aji {round(np.mean(list_iou_from_aji), 4)}')
    print(f'mean list_iou2 {round(np.mean(total_list_iou), 4)}')







acc = metrics.accuracy_score(y_true = total_y_true, y_pred=total_y_pred)
print(f'final acc {round(acc, 4)}')
bool_zero = np.array(total_y_pred) != 'Zero'
acc = metrics.accuracy_score(y_true = np.array(total_y_true)[bool_zero], y_pred=np.array(total_y_pred)[bool_zero])
print(f'final acc without no predicted {round(acc, 4)}')
acc = metrics.accuracy_score(y_true = total_y_true, y_pred=total_y_pred_cosine)
print(f'final acc cosine {round(acc, 4)}')

acc = metrics.accuracy_score(y_true = total_y_true, y_pred=total_y_pred_cosine_int)
print(f'final acc cosine int {round(acc, 4)}')

print(f'mean list_iou2 {round(np.mean(total_list_iou), 4)}')
print(f'median list_iou2 {round(np.median(total_list_iou), 4)}')

print((f'mean aji {round(np.mean(total_aji), 4)}'))


print(f'mean total_error_list {round(np.mean(total_error_list), 4)}')
print((f'mean total_prc_not_catch_list  {round(np.mean(total_prc_not_catch_list), 4)}'))
print(f'median total_error_list {round(np.median(total_error_list), 4)}')
print((f'median total_prc_not_catch_list  {round(np.median(total_prc_not_catch_list), 4)}'))
