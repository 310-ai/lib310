""" lib310 python package """

import logging

# from ._settings import settings
from . import data as db
from . import machinelearning as ml
from . import tools as tools
from . import visualization as plot

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata
package_name = "lib310"
__version__ = importlib_metadata.version(package_name)

# settings.verbosity = logging.INFO

b10_logger = logging.getLogger("lib310")
b10_logger.propagate = False

__all__ = [
    'db', 'ml', 'tools', 'plot',
]
