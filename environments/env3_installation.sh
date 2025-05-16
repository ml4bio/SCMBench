#!/bin/bash

conda create -n env3 python=3.10 -y
conda activate env3

pip install torch==2.1.0 torchvision==0.16.0 torchaudio==2.1.0 --index-url https://download.pytorch.org/whl/cu118
pip install torchtext
pip install scgpt "flash-attn<1.0.5"  # optional, recommended
pip install wandb