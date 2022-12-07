""" lib310 python package """

import logging

from .lazy_loader import LazyLoader

# from ._settings import settings
# from . import database as db
# from . import data
# from . import machinelearning as ml
# from . import tools as tools
# from . import visualization as plot
db = LazyLoader('db', globals(), 'lib310.database')
data = LazyLoader('data', globals(), 'lib310.data')
ml = LazyLoader('ml', globals(), 'lib310.machinelearning')
tools = LazyLoader('tools', globals(), 'lib310.tools')
plot = LazyLoader('plot', globals(), 'lib310.visualization')

try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata
try:
    package_name = "lib310"
    __version__ = importlib_metadata.version(package_name)
except:
    __version__ = 'dev'

# settings.verbosity = logging.INFO

b10_logger = logging.getLogger("lib310")
b10_logger.propagate = False

__all__ = [
    'db', 'data', 'ml', 'tools', 'plot',
]
__author__ = '310.ai team ("Saman Fekri <saman@310.ai>",Mohsen Naghipourfar <naghipourfar@berkeley.edu>", "Ismail Naderi <inaderi@310.ai>)'
