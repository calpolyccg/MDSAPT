"""An MDA-kit for calcuating quantum interactions in psi4."""

# Add imports here
from . import log

# Import core items
from .config import Config, load_from_yaml_file
from .sapt import TrajectorySAPT, DockingSAPT
from .viewer import Viewer

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions


def create_logger(logfile='mdsapt.log'):
    logger = log.create('mdsapt', logfile)
    return logger


def log_banner():
    """Log program name and licence at INFO level"""
    logger.info(f"MD-SAPT {__version__} starting")
    logger.info("Copyright (c) 2021-2022 Alia Lescoulie, Astrid Yu, and Ashley Ringer McDonald")
    logger.info("Released under GPLv3 License")


logger = create_logger()
log_banner()
