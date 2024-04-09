

import math
import spateo as st
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import lil_matrix

def create_csv_file(rna_file, pixel_size=0.108, spot_size = None, image_size=(10000, 10000)):
    
    # Preprocess df
    #df = pd.read_csv(rna_file, index_col=0)
    df = pd.read_csv(rna_file)
    print(f" df shape is {df.shape}")
    if not np.isin( list(df.columns), ['x', 'y', 'geneID']).all():
        try:
            columns_to_drop = ['x', 'y', 'fov', 'transcript_id', 'cell_id', 'global_z', 'barcode_id']
            if columns_to_drop in df.columns:
                df = df.drop(columns=columns_to_drop)
            df = df.rename(columns={"gene":"geneID", "global_x":"x", "global_y": "y"})
        except Exception as e:
            print(f"exception is {e}")
    assert np.isin(['x', 'y', 'geneID'], df.columns).all(), "The columns of the dataframe should be x, y, geneID"
    df['x'] = df['x'] / pixel_size
    df['y'] = df['y'] / pixel_size
    df['MIDCounts'] = 1

    # Create spots
    if spot_size is not None:
        x_max, y_max = image_size
        df['x_bin'] = pd.cut(df['x'], bins=np.arange(0, x_max, spot_size))
        df['y_bin'] = pd.cut(df['y'], bins=np.arange(0, y_max, spot_size))
        grouped_df = df.groupby(['x_bin', 'y_bin', 'geneID']).size().reset_index(name='MIDCounts')
        grouped_df = grouped_df[grouped_df['MIDCounts'] > 0]
        grouped_df['x'] = grouped_df['x_bin'].apply(lambda x: x.mid).astype(int)
        grouped_df['y'] = grouped_df['y_bin'].apply(lambda x: x.mid).astype(int)
        df = grouped_df.drop(columns = ["x_bin", "y_bin"])
    else:
        df['MIDCounts'] = 1

    df['x'] = df['x'].astype(int)
    df['y'] = df['y'].astype(int)

    return df[['geneID', 'x', 'y', 'MIDCounts']].reset_index(drop=True)

def merge_rna_count_dfs(rna_dfs):
    merged_df = pd.concat(rna_dfs)
    grouped_df = merged_df.groupby(['x', 'y', 'geneID'])['MIDCounts'].sum().reset_index()
    grouped_df = grouped_df[grouped_df['MIDCounts'] > 0]
    return grouped_df[['geneID', 'x', 'y', 'MIDCounts']].reset_index(drop=True)


def get_adatasub(bin_file, image_file, prealigned, startx=0, starty=0, patchsize=1000):
    if prealigned:
        adatasub = st.io.read_bgi_agg(bin_file, image_file, prealigned=True)
    else:
        adatasub = st.io.read_bgi_agg(bin_file, image_file)
    if int(patchsize) > 0:
        adatasub = adatasub[startx:startx+patchsize,starty:starty+patchsize]
    min_val, max_val = adatasub.layers["stain"].min(), adatasub.layers["stain"].max()
    adatasub.layers["stain"] = ((adatasub.layers["stain"]-min_val)/(max_val-min_val)*255).astype(int)
    adatasub.layers['unspliced'] = adatasub.X
    return adatasub


def compute_segmentation2center(segmentation_array):
    segmentation2x = {}
    segmentation2y = {}
    for i in range(segmentation_array.shape[0]):
        for j in range(segmentation_array.shape[1]):
            if segmentation_array[i, j] == 0:
                continue
            if segmentation_array[i, j] in segmentation2x:
                segmentation2x[segmentation_array[i, j]].append(i)
                segmentation2y[segmentation_array[i, j]].append(j)
            else:
                segmentation2x[segmentation_array[i, j]] = [i]
                segmentation2y[segmentation_array[i, j]] = [j]

    segmentation2center = {}
    for nucleus in segmentation2x:
        segmentation2center[nucleus] = [np.mean(segmentation2x[nucleus]), np.mean(segmentation2y[nucleus])]
    
    return segmentation2center

def compute_xmin_ymin(bin_file):
    xall = []
    yall = []
    with open(bin_file) as fr:
        header = fr.readline()
        for line in fr:
            gene, x, y, count = line.split()
            xall.append(int(x))
            yall.append(int(y))
    xmin = np.min(xall)
    ymin = np.min(yall)
    return xmin, ymin

def find_all_genes(bin_file):
    geneid = {}
    genecnt = 0
    id2gene = {}
    with open(bin_file) as fr:
        header = fr.readline()
        for line in fr:
            gene, x, y, count = line.split()
            if gene not in geneid:
                geneid[gene] = genecnt
                id2gene[genecnt] = gene
                genecnt += 1
    return geneid, id2gene, genecnt


def compute_idx2exp(bin_file, geneid, downrs, startx, starty, patchsizex, patchsizey):
    xmin, ymin = compute_xmin_ymin(bin_file)
    print(f'xmin: {xmin}, ymin: {ymin} from merge bin file {bin_file} are replace by  1 1')
    xmin = 1
    ymin = 1
    idx2exp = {}
    with open(bin_file) as fr:
        header = fr.readline()
        for line in fr:
            gene, x, y, count = line.split()
            x = int(x) - xmin
            y = int(y) - ymin
            if gene not in geneid:
                continue
            if int(x) < startx or int(x) >= startx + patchsizex or int(y) < starty or int(y) >= starty + patchsizey:
                continue
            idx = int(math.floor((int(x) - startx) / downrs) * math.ceil(patchsizey / downrs) + math.floor((int(y) - starty) / downrs))
            if idx not in idx2exp:
                idx2exp[idx] = {}
                idx2exp[idx][geneid[gene]] = int(count)
            elif geneid[gene] not in idx2exp[idx]:
                idx2exp[idx][geneid[gene]] = int(count)
            else:
                idx2exp[idx][geneid[gene]] += int(count)
    return idx2exp

def compute_all_exp_merged_bins(genecnt, idx2exp, downrs, patchsizex, patchsizey, n_top_genes=2000):

    all_exp_merged_bins = lil_matrix((int(math.ceil(patchsizex / downrs) * math.ceil(patchsizey / downrs)), genecnt), dtype=np.int8)
    for idx in idx2exp:
        for gid in idx2exp[idx]:
            all_exp_merged_bins[idx, gid] = idx2exp[idx][gid]
        
    all_exp_merged_bins = all_exp_merged_bins.tocsr()
    all_exp_merged_bins_ad = ad.AnnData(
        all_exp_merged_bins,
        obs=pd.DataFrame(index=[i for i in range(all_exp_merged_bins.shape[0])]),
        var=pd.DataFrame(index=[i for i in range(all_exp_merged_bins.shape[1])]),
    )
    sc.pp.highly_variable_genes(all_exp_merged_bins_ad, n_top_genes=n_top_genes, flavor='seurat_v3', span=1.0)
    selected_index = all_exp_merged_bins_ad.var[all_exp_merged_bins_ad.var.highly_variable].index
    selected_index = list(selected_index)
    selected_index = [int(i) for i in selected_index]
   
    #check total gene counts
    all_exp_merged_bins = all_exp_merged_bins.toarray()[:, selected_index]
    return all_exp_merged_bins, selected_index

def compute_offsets(bin_size, max_dis = 10):
    offsets = []
    for dis in range(1, max_dis + 1):
        for dy in range(-dis, dis + 1):
            offsets.append([-dis * bin_size, dy * bin_size])
        for dy in range(-dis, dis + 1):
            offsets.append([dis * bin_size, dy * bin_size])
        for dx in range(-dis + 1, dis):
            offsets.append([dx * bin_size, -dis * bin_size])
        for dx in range(-dis + 1, dis):
            offsets.append([dx * bin_size, dis * bin_size])
    return offsets