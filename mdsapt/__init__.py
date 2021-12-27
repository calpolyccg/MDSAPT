"""An MDA-kit for calcuating quantum interactions in psi4."""

# Add imports here
from . import log

# Import core classes
from .reader import InputReader
from .optimizer import Optimizer
from .sapt import TrajectorySAPT
from .viewer import Viewer

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions


def create_logger(logfile='mdsapt.log'):
    logger = log.create('MDSAPT', logfile)
    return logger

def log_banner():
    """Log program name and licence at INFO level"""
    logger.info(f"MDSAPT {__version__} starting")
    logger.info("Copyright (c) 2021 Alia Lescoulie")
    logger.info("Released under MIT Licence")

logger = create_logger()
log_banner()