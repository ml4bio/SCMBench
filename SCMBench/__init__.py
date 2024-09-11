r"""
GLUE (Graph-Linked Unified Embedding)
"""


from . import genomics, data, graph, models, num, plot
from .utils import config, log

try:
    from importlib.metadata import version
except ModuleNotFoundError:
    from pkg_resources import get_distribution
    version = lambda name: get_distribution(name).version


name = "SCMBench"
# __version__ = version(name)
