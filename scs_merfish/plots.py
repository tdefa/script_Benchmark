import spateo as st
import matplotlib.pyplot as plt
import numpy as np
import random

def plot_alignment(fig_folder, adatasub, before, patch_ids):
    fig, axes = plt.subplots(ncols=2, figsize=(16, 8), tight_layout=True)
    axes[0].imshow(before)
    st.pl.imshow(adatasub, 'unspliced', ax=axes[0], alpha=0.6, cmap='Reds', vmax=10, use_scale=False, save_show_or_return='return')
    axes[0].set_title('before alignment')
    st.pl.imshow(adatasub, 'stain', ax=axes[1], use_scale=False, save_show_or_return='return')
    st.pl.imshow(adatasub, 'unspliced', ax=axes[1], alpha=0.6, cmap='Reds', vmax=10, use_scale=False, save_show_or_return='return')
    axes[1].set_title('after alignment')
    plt.savefig(fig_folder / f'alignment_{patch_ids}.png')
    plt.show()

def plot_segmentation(fig_folder, adatasub, segmentation_labels, patch_ids):
    fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)
    st.pl.imshow(adatasub, 'stain', save_show_or_return='return', ax=ax)
    st.pl.imshow(adatasub, segmentation_labels, labels=True, alpha=0.5, ax=ax)
    plt.savefig(fig_folder / f'{segmentation_labels}_{patch_ids}.png')
    plt.show()

def plot_rna_dist(fig_folder, result_df):
    total_rna_df = result_df.groupby(['x', 'y'])['MIDCounts'].sum().reset_index()
    total_rna_df = total_rna_df[total_rna_df['MIDCounts'] > 0]
    plt.hist(total_rna_df["MIDCounts"])
    plt.title("Distribution of MID Counts per spot in RNA Data")
    plt.savefig(fig_folder / f'rna_distribution.png')
    plt.show()


def plot_segmentation_plus_rna(fig_folder, adatasub, result_df, patch_ids):
    fig, ax = plt.subplots(figsize=(8, 8), tight_layout=True)
    st.pl.imshow(adatasub, 'stain', save_show_or_return='return', ax=ax)
    st.pl.scatter(result_df, color='red', ax=ax, save_show_or_return='return')
    plt.savefig(fig_folder / f'rna_{patch_ids}.png')
    plt.show()


def plot_rna_on_image(fig_folder, df_spot,
                      images, patch_ids , figsize=(16, 16),
                      spots_size=1, pixel_size=0.108,gene_column="geneID",
                           x_key="x", y_key="y"):
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16, 16), tight_layout=True, )
    ax[0].imshow(images)
    ax[1].imshow(images)
    for gene in np.unique(df_spot[gene_column]):
        df_gene = df_spot[df_spot[gene_column] == gene]
        x_list = (df_gene[x_key] * (1 / pixel_size)).tolist()
        y_list = (df_gene[y_key] * (1 / pixel_size)).tolist()
        ## add random color in hexa
        color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
        ## add color in the scatter plot for each gene
        ax[1].scatter(y_list, x_list, s=spots_size, c=color)
    plt.savefig(fig_folder / f'rna_on_image_{patch_ids}.png')
    plt.show()
