# Copied from Gromacs Wrapper by Oliver Beckstein
# https://github.com/Becksteinlab/GromacsWrapper/blob/master/gromacs/utilities.py

import os

import errno

from contextlib import contextmanager

import logging

logger = logging.getLogger('mdsapt.utils')


@contextmanager
def in_dir(directory, create=True):
    """Context manager to execute a code block in a directory.
    * The directory is created if it does not exist (unless
      create=False is set)
    * At the end or after an exception code always returns to
      the directory that was the current directory before entering
      the block.
    """
    startdir = os.getcwd()
    try:
        try:
            os.chdir(directory)
            logger.debug("Working in %r...", directory)
        except OSError as err:
            if create and err.errno == errno.ENOENT:
                os.makedirs(directory)
                os.chdir(directory)
                logger.info("Working in %r (newly created)...", directory)
            else:
                logger.exception("Failed to start working in %r.", directory)
                raise
        yield os.getcwd()
    finally:
        os.chdir(startdir)
