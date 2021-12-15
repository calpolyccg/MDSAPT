#!/usr/bin/env python

import logging
logger = logging.getLogger('mdsapt')

from argparse import ArgumentParser


if __name__ == '__main__':
    parser = ArgumentParser(usage=__doc__)
    parser.add_argument('filename', metavar='F', type='str', nargs='?',
                        default='input.yaml',
                        help='Filename for generated input file input.yaml')

    filename = vars(parser.parse_args())['filename']

    with open('./mdsapt/data/template_input.yaml', 'r') as template:
        template_data = template

    with open(f'{filename}.yaml', 'w+') as new_file:
        new_file.write(template_data)
    logger.info(f'Generated template input file {filename}')
