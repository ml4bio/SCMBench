#!/bin/bash

conda create -n env1 python=3.8 -y
conda activate env1
pip install torch==2.0.0 torchvision==0.15.1 torchaudio==2.0.1
pip install scanpy==1.9.8
pip3 install pamona
pip install scglue
pip install git+https://github.com/epurdom/cobolt.git#egg=cobolt
pip install mofapy2
git clone https://github.com/xianglin226/scMDC.git
pip install scmomat
pip3 install unioncom
pip install harmonypy
pip3 install unioncom