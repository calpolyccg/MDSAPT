# pylint: skip-file
"""
Based on log.py https://github.com/Becksteinlab/MDPOW/blob/develop/mdpow/log.py

Logs issues throughout MD-SAPT
"""
import logging


def create(log_name: str, logfile: str) -> logging.Logger:
    """
    Create the logger for MD-SAPT
    """
    logger = logging.getLogger(log_name)

    logger.setLevel(logging.DEBUG)

    logfile = logging.FileHandler(logfile)
    logger_formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    logfile.setFormatter(logger_formatter)
    logger.addHandler(logfile)

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logger.addHandler(console)

    return logger
