import tifffile
import spateo as st
from cellpose import models 

def compute_cellpose_segmentation(adata, gpu=True, model_type='cyto',pretrained_model=None,
                                  diameter=90, channels=[0, 0], flow_threshold=0.4,
                                  ):
    """

    Parameters
    ----------
    adata
    gpu
    model_type : choose in ['cyto', 'nuclei', TN1', 'TN2' ... check availalbe model here https://cellpose.readthedocs.io/en/latest/models.html#model-zoo
    pretrained_model : set model to None to use the pretrained_model. pretrained_model = str(path_to_model) to use a custom model
    diameter :  parameter that plays with the ration pixel size / cell size
    channels
    flow_threshold

    Returns
    -------

    """
    model = models.CellposeModel(gpu=gpu, model_type=model_type, pretrained_model=pretrained_model)
    masks, flows, styles = model.eval(
        adata.layers['stain'],
        diameter=diameter,
        channels=channels,
        flow_threshold=flow_threshold,
        do_3D=False,
        stitch_threshold=0
    )
    return masks

def compute_watershed_segmentation(adata, otsu_classes = 4, otsu_index=1):
    st.cs.mask_nuclei_from_stain(adata, otsu_classes = otsu_classes, otsu_index=otsu_index)
    st.cs.find_peaks_from_mask(adata, 'stain', 7)
    st.cs.watershed(adata, 'stain', 5, out_layer='watershed_labels')
    return adata.layers["watershed_labels"]
