r"""
Performance evaluation metrics
"""

from typing import Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.spatial
import sklearn.metrics
import sklearn.neighbors
from anndata import AnnData
from scipy.sparse.csgraph import connected_components

from .batch_asw import silhouette_batch
from .batch_graph_connectivity import graph_connectivity
from .batch_kbet import kBET
from .batch_ilisi import ilisi_graph
from .batch_pcr import pcr_comparison