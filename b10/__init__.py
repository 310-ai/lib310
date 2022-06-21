""" b10 python package """

import logging

from ._settings import settings
from . import data as db
from . import machinelearning as MLAPI
from . import tools as tl
from . import visualization as pl

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata
package_name = "b10"
__version__ = importlib_metadata.version(package_name)

settings.verbosity = logging.INFO

b10_logger = logging.getLogger("b10")
b10_logger.propagate = False

__all__ = [
    'db', 'MLAPI', 'tl', 'pl',
]
