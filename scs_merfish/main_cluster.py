
import os
from pathlib import Path
import pandas as pd
from scs_package.transformer import train
from scs_package.postprocessing import postprocess
from dataset import create_train_test_data
from dataset_utils import create_csv_file, merge_rna_count_dfs
import argparse
import random
import datetime
e = datetime.datetime.now()
date_str = f"{e.month}_d{e.day}_h{e.hour}_min{e.minute}_s{e.second}_r" + str(random.randint(0 ,5000))
import numpy as np
import tifffile
from plots import plot_rna_on_image
import datetime
print("not main_cluster.py")

if __name__ == '__main__':


    print("in main_cluster.py")
    parser = argparse.ArgumentParser(description='test')
    ## task
    parser.add_argument("--do_prerpo_csv_to_tsv", type=int, default=1)
    parser.add_argument("--do_prerpo", type=int, default=1)
    parser.add_argument("--do_train", type=int, default=1)
    parser.add_argument("--do_postprocess", type=int, default=1)
    ## prerpo parameters
    parser.add_argument("--rna_filenames", type=str,
                        default="08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_009.csv")
    parser.add_argument("--image_file_name", type=str, default="08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_009.tiff")
    parser.add_argument("--new_rna_csv_filename", type=str, default="preprocessed_rna_no_spot.tsv")
    parser.add_argument("--main_folder", type=str,
                        default="/cluster/CBIO/data1/data3/tdefard/T7/simulation/scs_input/")
    parser.add_argument("--regex_save", type=str, default="test_bis")
    parser.add_argument("--pixel_size", type=float, default=1)
    parser.add_argument("--image_size_0y", type=int, default=0)
    parser.add_argument("--image_size_1x", type=int, default=0)
    parser.add_argument("--segmentation_labels", type=str,
                        default='simulation_labels', help="cellpose_labels or watershed_labels or simulation_labels")
    parser.add_argument("--simulation_seg_path", type=str,
                        default="/cluster/CBIO/data1/data3/tdefard/T7/simulation/scs_input/images/08_IR5M_Pdgfra-Cy3_Serpine1-Cy5_009.tiff")
    ## add cellpose parameters
    parser.add_argument("--cellpose_model", type=str, default="nuclei")
    parser.add_argument("--cellpose_diameter", type=int, default=90)
    parser.add_argument("--cellpose_flow_threshold", type=float, default=0.9)
    parser.add_argument("--cellpose_pretrained_model", type=str, default=None)
    #### SCS parameters
    parser.add_argument("--prealigned", type=int, default=1)
    parser.add_argument("--align", type=int, default=0)
    parser.add_argument("--patchsize", type=int, default=0)
    parser.add_argument("--n_neighbor", type=int, default=50)
    parser.add_argument("--bin_size", type=int, default=8)
    parser.add_argument("--background_threshold", type=int, default=15)
    parser.add_argument("--n_top_genes", type=int, default=2000)
    parser.add_argument("--epoch", type=int, default=100)
    parser.add_argument("--val_ratio", type=float, default=0.0625)
    parser.add_argument("--dia_estimate", type=int, default=15)
    parser.add_argument("--cell_size_threshold", type=int, default=200)
    parser.add_argument("--port", default=3950)
    parser.add_argument("--mode", default='client')
    parser.add_argument("--host", default='127.0.0.2')
    args = parser.parse_args()
    args.prealigned = bool(args.prealigned)
    args.align = bool(args.align)

    cellpose_params = {
        "model_type": args.cellpose_model,
        "diameter": args.cellpose_diameter,
        "flow_threshold": args.cellpose_flow_threshold,
        "pretrained_model": args.cellpose_pretrained_model,
    }




    rna_filenames = [args.rna_filenames]
    #args.rna_filenames = ['z_0.0.csv', 'z_1.0.csv', 'z_2.0.csv', 'z_3.0.csv', 'z_4.0.csv', 'z_5.0.csv', 'z_6.0.csv']
    #rna_filenames = args.rna_filenames
    print("args.rna_filenames hardcoded ", args.rna_filenames)


    main_folder = Path(args.main_folder)
    image_file_name = args.image_file_name
    new_rna_csv_filename = args.new_rna_csv_filename
    pixel_size = args.pixel_size
    image_size = (args.image_size_0y, args.image_size_1x)

    image_file = main_folder / "images" / image_file_name

    if args.regex_save is not None:
        args.save_folder = main_folder / (args.regex_save)
    else:
        args.save_folder = main_folder / ("scs_pred" + date_str)


    Path(args.save_folder).mkdir(parents=True, exist_ok=True)
    rna_folder = main_folder / "csv"
    Path(rna_folder).mkdir(parents=True, exist_ok=True)

    bin_file = rna_folder / new_rna_csv_filename
    fig_folder = args.save_folder / "fig"
    Path(fig_folder).mkdir(parents=True, exist_ok=True)

    # Dataset creation parameters
    segmentation_labels = args.segmentation_labels
    prealigned = args.prealigned
    align = args.align
    patchsize = args.patchsize
    n_neighbor = args.n_neighbor
    bin_size = args.bin_size
    startx = 0
    starty = 0
    background_threshold = args.background_threshold
    n_top_genes = args.n_top_genes

    np.save(args.save_folder / "script_parameter", args)
    with open(args.save_folder / "script_parameter.txt", "w") as f:
        for k, v in args.__dict__.items():
            f.write(f"{k} : {v}\n")
            print(f"{k} : {v}")

    # Training parameters
    epochs = args.epoch
    val_ratio = args.val_ratio
    print("in code")
    if args.do_prerpo_csv_to_tsv:
        rna_dfs = []
        for rna_filename in rna_filenames:
            print(rna_filename)
            rna_filepath = os.path.join(rna_folder, rna_filename)
            rna_count_df = create_csv_file(
                rna_filepath,
                pixel_size=pixel_size,
                image_size=image_size, ### image size for what ?
            )
            rna_dfs.append(rna_count_df)
            #### plot RNA on image
            df_spot = pd.read_csv(rna_filepath)
            #list_columns_to_drop = ['global_z']
            #df_spot = df_spot.drop(columns = list_columns_to_drop)
            #df_spot = df_spot.rename(columns={"gene":"geneID", "global_x":"x", "global_y": "y"})
            #df_spot.to_csv(rna_filepath)
            images = tifffile.imread(image_file)
            plot_rna_on_image(fig_folder, df_spot,
                                images, patch_ids = rna_filename,
                              figsize=(16, 16), spots_size=1, pixel_size=args.pixel_size)
        result_df = merge_rna_count_dfs(rna_dfs)
        ## check if perent fodelr folder
        result_df.to_csv(bin_file, sep='\t', index=False)




    if args.do_prerpo:
        print("create_train_test_data")

        create_train_test_data(
            bin_file = bin_file,
            image_file = image_file,
            save_folder = args.save_folder,
            segmentation_labels = segmentation_labels,
            prealigned = prealigned,
            align = align,
            patchsize = patchsize,
            n_neighbor = n_neighbor,
            bin_size = bin_size,
            startx = startx,
            starty = startx,
            background_threshold = background_threshold,
            n_top_genes = n_top_genes,
            cellpose_params = cellpose_params,
            simulation_seg_path = args.simulation_seg_path,
        )

    if args.do_train:
        train(startx , starty, patchsize, epochs, val_ratio, save_folder = args.save_folder)
        print("done training")

    if args.do_postprocess:
        print("post processing")
        dia_estimate = args.dia_estimate
        bin_size_postprocess = args.bin_size
        postprocess(
                    startx = 0,
                    starty = 0,
                    patchsize = patchsize,
                    bin_size = bin_size_postprocess,
                    dia_estimate = dia_estimate,
                    save_folder = args.save_folder,
                    segmentation_labels = args.segmentation_labels,
                    cell_size_threshold = args.cell_size_threshold,
                    )
    print("done")
