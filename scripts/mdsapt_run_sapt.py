#!/usr/bin/env python

from argparse import ArgumentParser

import logging

logger = logging.getLogger('mdsapt')

import mdsapt

if __name__ == '__main__':
    parser = ArgumentParser(usage=__doc__)
    parser.add_argument('filename', metavar='I', type=str, nargs='?',
                        default='input.yaml',
                        help='Filename for input yaml file')
    parser.add_argument('outfile', metavar='O', type=str, nargs='?',
                        default='out.csv',
                        help='Name for output csv file')

    filename = vars(parser.parse_args())['filename']

    settings = mdsapt.InputReader(filename)

    optimizer = mdsapt.Optimizer(settings)

    if optimizer.num_failed_residue != 0:
        logger.error('optimization failed see log for list of failed residues')

    sapt = mdsapt.TrajectorySAPT(settings, optimizer).run(settings.start, settings.stop, settings.step)

    logger.info('saving results to CSV')
    sapt.results.to_csv(filename)
