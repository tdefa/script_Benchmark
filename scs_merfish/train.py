import os
import gc
import numpy as np
import anndata as ad
from pathlib import Path

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import tensorflow_addons as tfa
from tensorflow.keras.layers.experimental import preprocessing

from scs.transformer import dir_to_class, create_transformer_classifier, run_experiment


def train(startx, starty, patchsize, epochs, val_ratio, save_folder, class_num = 16):
 
    results_folder = Path(save_folder) / "results"
    results_folder.mkdir(parents=True, exist_ok=True)
    data_folder = Path(save_folder) / "data"
    data_folder.mkdir(parents=True, exist_ok=True)
    checkpoint_folder = Path(save_folder) / "ckpt" / f'model_{startx}_{starty}_{patchsize}'
    checkpoint_folder.mkdir(parents=True, exist_ok=True)

    patch_id = str(startx) + ':' + str(starty) + ':' + str(patchsize) + ':' + str(patchsize)

    x_train_ = np.load(data_folder / f'x_train_{patch_id}.npz')
    x_train_ = x_train_['x_train'].astype(np.float32)
    x_train_pos_ = np.load(data_folder /f'x_train_pos_{patch_id}.npz')
    x_train_pos_ = x_train_pos_['x_train_pos'].astype(np.int32)
    y_train_ = np.load(data_folder /f'y_train_{patch_id}.npz')
    y_train_ = y_train_['y_train']
    y_binary_train_ = np.load(data_folder / f'y_binary_train_{patch_id}.npz')
    y_binary_train_ = y_binary_train_['y_binary_train'].astype(np.int32)
    x_test = np.load(data_folder / f'x_test_{patch_id}.npz')
    x_test = x_test['x_test'].astype(np.float32)
    x_test_pos = np.load(data_folder / f'x_test_pos_{patch_id}.npz')
    x_test_pos = x_test_pos['x_test_pos'].astype(np.int32)

    x_train_select = []
    x_validation_select = []
    adata = ad.read_h5ad(data_folder / f'spots{patch_id}.h5ad')
    for i in range(len(x_train_pos_)):
        if x_train_pos_[i][0][0] > int(adata.X.shape[0] * (1 - np.sqrt(val_ratio))) and x_train_pos_[i][0][1] > int(adata.X.shape[1] * (1 - np.sqrt(val_ratio))):
            x_validation_select.append(i)
        else:
            x_train_select.append(i)

    for i in range(len(y_train_)):
        if y_train_[i][0] != -1:
            y_train_[i] = y_train_[i] - x_train_pos_[i][0]
        else:
            y_train_[i][0] = -9999
            y_train_[i][1] = -9999
    y_train_ = dir_to_class(y_train_, class_num)
    for i in range(len(x_train_pos_)):
        for j in range(1, len(x_train_pos_[i])):
            x_train_pos_[i][j] = x_train_pos_[i][j] - x_train_pos_[i][0]
        x_train_pos_[i][0] = x_train_pos_[i][0] - x_train_pos_[i][0]
    for i in range(len(x_test_pos)):
        for j in range(1, len(x_test_pos[i])):
            x_test_pos[i][j] = x_test_pos[i][j] - x_test_pos[i][0]
        x_test_pos[i][0] = x_test_pos[i][0] - x_test_pos[i][0]

    x_train = x_train_[x_train_select]
    x_train_pos = x_train_pos_[x_train_select]
    y_train = y_train_[x_train_select]
    y_binary_train = y_binary_train_[x_train_select]
    x_validation = x_train_[x_validation_select]
    x_validation_pos = x_train_pos_[x_validation_select]
    y_validation = y_train_[x_validation_select]
    y_binary_validation = y_binary_train_[x_validation_select]

    input_shape = (x_train.shape[1], x_train.shape[2])
    input_position_shape = (x_train_pos.shape[1], x_train_pos.shape[2])

    learning_rate = 0.001
    weight_decay = 0.0001
    batch_size = 10
    num_epochs = epochs
    num_patches = x_train.shape[1]
    projection_dim = 64
    num_heads = 1
    transformer_units = [
        projection_dim * 2,
        projection_dim,
    ]  # Size of the transformer layers
    transformer_layers = 8
    mlp_head_units = [1024, 256]  # Size of the dense layers of the final classifier

    transformer_classifier = create_transformer_classifier(class_num, input_shape, input_position_shape, num_patches, projection_dim, num_heads, transformer_units, transformer_layers, mlp_head_units)
    run_experiment(
        str(startx), str(starty), str(patchsize), 
        transformer_classifier, 
        x_train, x_train_pos, x_train_, x_train_pos_, y_train, y_train_, y_binary_train, 
        x_test, x_test_pos, 
        x_validation, x_validation_pos, y_validation, y_binary_validation, 
        learning_rate, weight_decay, batch_size, num_epochs,
        checkpoint_folder = checkpoint_folder / "ckpt",
        data_folder = data_folder, 
        results_folder = results_folder,
    )


if __name__=="__main__":
    patch_size = 0
    epochs = 100
    val_ratio = 0.0625
    save_folder = "/cluster/CBIO/data1/ablondel1/vallot/202304141840_1184_VallotAB5301andAG5988_VMSC09402/region_0/crop10000_2_3/scs_pred/"
    train(0, 0, patch_size, epochs, val_ratio, save_folder = save_folder)