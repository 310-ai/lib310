""" b10 python package """

import logging

from . import data as dataAPI
from . import ml as MLAPI

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata
package_name = "b10"
__version__ = importlib_metadata.version(package_name)

test_var = "test"

b10_logger = logging.getLogger("b10")
b10_logger.propagate = False

__all__ = [
    'dataAPI', 'MLAPI'
]
