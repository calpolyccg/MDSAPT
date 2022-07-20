# pylint: skip-file
"""An MDA-kit for calcuating quantum interactions in psi4."""
import logging

from . import log

# Import core classes and methods
from .config import Config, load_from_yaml_file
from .sapt import TrajectorySAPT, DockingSAPT

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions


def create_logger(logfile='mdsapt.log') -> logging.Logger:
    """Starts logger"""
    md_sapt_log = log.create('mdsapt', logfile)
    return md_sapt_log


def log_banner() -> None:
    """Log program name and licence at INFO level"""
    logger.info("MD-SAPT %s starting" % __version__)
    logger.info("Copyright (c) 2021-2022 Alia Lescoulie, Astrid Yu, and Ashley Ringer McDonald")
    logger.info("Released under GPLv3 License")


logger = create_logger()
log_banner()
