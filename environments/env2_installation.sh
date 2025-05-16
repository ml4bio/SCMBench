#!/bin/bash

conda create -n env2 python=3.11 -y
conda activate env2
pip install torch==2.6.0 torchvision==0.21.0 torchaudio==2.6.0 --index-url https://download.pytorch.org/whl/cu118
pip install scanpy==1.9.8
pip install pybedtools
pip install scvi-tools
pip install --quiet scvi-colab
pip install muon
pip install mudata