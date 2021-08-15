# Based on log.py https://github.com/Becksteinlab/MDPOW/blob/develop/mdpow/log.py
import logging


def create(logname: str, logfile: str):
    logger = logging.getLogger(logname)

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
