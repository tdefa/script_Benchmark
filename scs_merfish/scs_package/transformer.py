#%%
import gc
import os
from pathlib import Path
import anndata as ad
from pathlib import Path
import numpy as np
import tensorflow as tf
import tensorflow_addons as tfa
from tensorflow.keras import layers
from tensorflow import keras
from tensorflow.keras.layers.experimental import preprocessing


def dir_to_class(y_dir, class_num):
    y_dir_class = []
    for i in range(len(y_dir)):
        x, y = y_dir[i]
        if x == -9999:
            y_vec = np.zeros(class_num)
            y_dir_class.append(y_vec)
        else:
            if y == 0 and x > 0:
                deg = np.arctan(float('inf'))
            elif y == 0 and x < 0:
                deg = np.arctan(-float('inf'))
            elif y == 0 and x == 0:
                deg = np.arctan(0)
            else:
                deg = np.arctan((x/y))
            if (x > 0 and y < 0) or (x <= 0 and y < 0):
                deg += np.pi
            elif x < 0 and y >= 0:
                deg += 2 * np.pi
            cla = int(deg / (2 * np.pi / class_num))
            y_vec = np.zeros(class_num)
            y_vec[cla] = 1
            y_dir_class.append(y_vec)
    return np.array(y_dir_class)

def mlp(x, hidden_units, dropout_rate):
    for units in hidden_units:
        x = layers.Dense(units, activation=tf.nn.gelu)(x)
        x = layers.Dropout(dropout_rate)(x)
    return x

def masked_categorical_cross_entropy(y_true, y_pred):
    cce = keras.losses.CategoricalCrossentropy(from_logits=True)
    return cce(y_true, y_pred, sample_weight=tf.math.reduce_sum(y_true, axis=-1)) #* weights

class PatchEncoder(layers.Layer):
    def __init__(self, num_patches, projection_dim):
        super(PatchEncoder, self).__init__()
        self.num_patches = num_patches
        self.projection = layers.Dense(units=projection_dim)
        self.position_embedding = layers.Dense(units=projection_dim)

    def call(self, patch, position):
        encoded = self.projection(patch) + self.position_embedding(position)
        return encoded

def create_transformer_classifier(class_num, input_shape, input_position_shape, num_patches, projection_dim, num_heads, transformer_units, transformer_layers, mlp_head_units):
    inputs = layers.Input(shape=input_shape)
    inputs_positions = layers.Input(shape=input_position_shape)
    encoded_patches = PatchEncoder(num_patches, projection_dim)(inputs, inputs_positions)

    for _ in range(transformer_layers):
        # Layer normalization 1.
        x1 = layers.LayerNormalization(epsilon=1e-6)(encoded_patches)
        # Create a multi-head attention layer.
        attention_output = layers.MultiHeadAttention(
            num_heads=num_heads, key_dim=projection_dim, dropout=0.1
        )(x1, x1)
        # Skip connection 1.
        x2 = layers.Add()([attention_output, encoded_patches])
        # Layer normalization 2.
        x3 = layers.LayerNormalization(epsilon=1e-6)(x2)
        # MLP.
        x3 = mlp(x3, hidden_units=transformer_units, dropout_rate=0.1)
        # Skip connection 2.
        encoded_patches = layers.Add()([x3, x2])

    # Create a [batch_size, projection_dim] tensor.
    representation = layers.LayerNormalization(epsilon=1e-6)(encoded_patches)
    # Add MLP.
    features = mlp(representation[:, 0, :], hidden_units=mlp_head_units, dropout_rate=0.5)
    # Classify outputs.
    pos = layers.Dense(class_num, name='pos_out')(features)
    binary = layers.Dense(1, activation='sigmoid', name='cat_out')(features)

    model = keras.Model(inputs=[inputs, inputs_positions], outputs=[pos, binary])
    return model

def run_experiment(
        startx, starty, patchsize, model, 
        x_train, x_train_pos, x_train_, x_train_pos_, y_train, y_train_, y_binary_train, 
        x_test, x_test_pos, x_validation, x_validation_pos, y_validation, y_binary_validation, 
        learning_rate, weight_decay, batch_size, num_epochs, 
        checkpoint_folder, data_folder, results_folder,
    save_freq, save_best_only=False):
    optimizer = tfa.optimizers.AdamW(
        learning_rate=learning_rate, weight_decay=weight_decay
    )
    model.compile(
        optimizer=optimizer,
        loss={
            'pos_out': masked_categorical_cross_entropy,
            'cat_out': keras.losses.BinaryCrossentropy(from_logits=False),
        },
        metrics={
            'pos_out': keras.metrics.CategoricalAccuracy(name="accuracy"),
        },
    )

    if len(x_validation[np.where(y_binary_validation == 1)]) < 100:
        checkpoint_callback = keras.callbacks.ModelCheckpoint(
            checkpoint_folder / '{epoch:03d}.keras',
            monitor="pos_out_accuracy",
            save_best_only=save_best_only,
            save_weights_only=True,
        )
    else:
        checkpoint_callback = keras.callbacks.ModelCheckpoint(
            checkpoint_folder / '{epoch:03d}.keras',
            monitor="pos_out_accuracy",
            save_best_only=save_best_only,
            save_weights_only=True,
        )

    #print(x_train)
    #print(x_train_pos)
    #print(y_train)
    #print(y_binary_train)
    #print(len(x_validation[np.where(y_binary_validation == 1)]))
    model.fit(
        x=[x_train, x_train_pos],
        y=[y_train, y_binary_train],
        batch_size=batch_size,
        epochs=num_epochs,
        validation_split=0.0,
        validation_data=([x_validation[np.where(y_binary_validation == 1)], x_validation_pos[np.where(y_binary_validation == 1)]], [y_validation[np.where(y_binary_validation == 1)], y_binary_validation[np.where(y_binary_validation == 1)]]),
        callbacks=[checkpoint_callback],
    )

    print('Inference on all the spots...')
    try:
        checkpoint_folder = Path("/home/tom/share/T7/2023-09-06_LUSTRA-14rounds/scs_input/15_5/ckpt/model_0_0_0")
        list_checkpoint = list(checkpoint_folder.glob("*.keras"))
        ## sort them by name
        list_checkpoint = sorted(list_checkpoint, key=lambda x: int(x.stem))
        print(f'Load the best model...{list_checkpoint[-1]}')

        model.load_weights(list_checkpoint[-1])
    except Exception as e:
        print(e)
    pred_centers_test_all = []
    pred_binary_test_all = []
    for i in range(int(len(x_test) / 10000) + 1):
        pred_centers_test_, pred_binary_test_ = model.predict(x = [x_test[i*10000: (i+1)*10000], x_test_pos[i*10000: (i+1)*10000]], batch_size=batch_size)
        pred_centers_test_all.append(pred_centers_test_)
        pred_binary_test_all.append(pred_binary_test_)
        gc.collect()
    pred_centers_test = np.vstack(pred_centers_test_all)
    pred_binary_test = np.vstack(pred_binary_test_all)
    #for i in range(len(x_test_pos)):
    #    print(x_test_pos[i][0], pred_binary_test[i], pred_centers_test[i])

    pred_centers_train_all = []
    pred_binary_train_all = []
    for i in range(int(len(x_train_) / 10000) + 1):
        pred_centers_train_, pred_binary_train_ = model.predict(x = [x_train_[i*10000: (i+1)*10000], x_train_pos_[i*10000: (i+1)*10000]], batch_size=batch_size)
        pred_centers_train_all.append(pred_centers_train_)
        pred_binary_train_all.append(pred_binary_train_)
        gc.collect()
    pred_centers_train = np.vstack(pred_centers_train_all)
    pred_binary_train = np.vstack(pred_binary_train_all)
    #for i in range(len(y_train_)):
    #    print(y_train_[i], pred_binary_train[i], pred_centers_train[i])
    x_train_pos__ = np.load(data_folder / f'x_train_pos_{startx}:{starty}:{patchsize}:{patchsize}.npz')
    x_train_pos__ = x_train_pos__['x_train_pos']
    x_test_pos_ = np.load(data_folder / f'x_test_pos_{startx}:{starty}:{patchsize}:{patchsize}.npz')
    x_test_pos_ = x_test_pos_['x_test_pos']
    #print(y_train_.shape, x_train_pos__.shape)

    print('Write prediction results...')
    with open(results_folder /f'spot_prediction_{startx}:{starty}:{patchsize}:{patchsize}.txt', 'w') as fw:
        for i in range(len(y_train_)):
            fw.write(str(x_train_pos__[i][0][0]) + '\t' + str(x_train_pos__[i][0][1]) + '\t' + str(pred_binary_train[i][0]) + '\t' + ':'.join([str(c) for c in pred_centers_train[i]]) + '\n')
        for i in range(len(x_test_pos_)):
            fw.write(str(x_test_pos_[i][0][0]) + '\t' + str(x_test_pos_[i][0][1]) + '\t' + str(pred_binary_test[i][0]) + '\t' + ':'.join([str(c) for c in pred_centers_test[i]]) + '\n')

    return

def get_train_val_test_dataset(patch_id, data_folder, val_ratio, class_num=16):
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

    return (
        x_train, x_train_pos, y_train, y_binary_train, 
        x_validation, x_validation_pos, y_validation, y_binary_validation,
        x_train_, x_train_pos_, y_train_, y_binary_train_, 
        x_test, x_test_pos,
    )

def train(startx, starty, patchsize, epochs, val_ratio, save_folder, class_num = 16, save_freq = None, save_best_only = True):
 
    data_folder = Path(save_folder) / "data"
    checkpoint_folder = Path(save_folder) / "ckpt" / f'model_{startx}_{starty}_{patchsize}'
    checkpoint_folder.mkdir(parents=True, exist_ok=True)

    results_folder = Path(save_folder) / "results"
    results_folder.mkdir(parents=True, exist_ok=True)

    patch_id = str(startx) + ':' + str(starty) + ':' + str(patchsize) + ':' + str(patchsize)

    (
        x_train, x_train_pos, y_train, y_binary_train, 
        x_validation, x_validation_pos, y_validation, y_binary_validation,
        x_train_, x_train_pos_, y_train_, y_binary_train_, 
        x_test, x_test_pos,
    ) = get_train_val_test_dataset(patch_id, data_folder, val_ratio=val_ratio)

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
        checkpoint_folder = checkpoint_folder  ,#/ "ckpt",
        data_folder = data_folder, 
        results_folder = results_folder,save_freq = save_freq,
        save_best_only = save_best_only
    )

#%%
if __name__=="__main__":
    patch_size = 0
    epochs = 100
    val_ratio = 0.0625
    save_folder = "/cluster/CBIO/data1/ablondel1/vallot/202304141840_1184_VallotAB5301andAG5988_VMSC09402/region_0/crop10000_2_2/crop1000_6_4/scs_pred/"
    train(0, 0, patch_size, epochs, val_ratio, save_folder = save_folder)