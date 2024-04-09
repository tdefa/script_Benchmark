#!/bin/bash


cd /home/tom/Bureau/phd/simulation/baysor/baysor_ubuntu-latest_x64_build/Baysor/bin/

./Baysor run  -x x -y y -z z --gene gene -c config.toml -o /home/tom/Bureau/phd/simulation/baysor/baysor_ubuntu-latest_x64_build/Baysor/bin/mouse_ilieum_3D/ -m 30 -i 400 -p --n-clusters=4 --scale-std=50% --nuclei-genes=Neat1 --cyto-genes=Apob,Net1,Slc5a1,Mptx2 --exclude-genes=Blank* --prior-segmentation-confidence=0.95 "/media/tom/Transcend/doi_10_5061_dryad_jm63xsjb2__v20210916/input_comseg_3Dseg_tom/df/prior_stitched.csv" :in_nucleus