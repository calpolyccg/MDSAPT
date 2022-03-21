import os
import logging

import click


logger = logging.getLogger(__name__)

_dir_path = os.path.dirname(os.path.realpath(__file__))
"""Location of the mdsapt package."""


@click.command()
def cli():
    """
    MDSAPT - Molecular Dynamics Symmetry-Adapted Perturbation Theory

    By Alia Lescoulie, Astrid Yu, and Ashley Ringer McDonald.

    This command-line interface lets you easily do common MDSAPT-related tasks.
    """


@cli.command()
@click.argument(
    'filename',
    default='input.yaml',
    help="Path of the input file to generate."
)
@click.option(
    '-t', '--template',
    dest='template',
    help="Template to generate from.",
    type=click.Choice(['trajectory', 'docking'], case_sensitive=False),
    default='trajectory',
)
@click.option(
    '-f', '--force',
    help="If provided, overwrites existing files.",
    default=False,
)
def generate(filename: str, template: str, force: bool):
    """
    Generate a template input file.
    """
    ensure_safe_to_overwrite(filename, force)

    # TODO: make a wizard for these templates
    template_path = os.path.join(_dir_path, 'mdsapt', 'data', f'{template}_template.yaml')

    with open(template_path, 'r') as template:
        template_data = template.read()
    with open(filename, 'w') as new_file:
        new_file.write(template_data)

    logger.info(f'Generated template input file {filename}')


@cli.command()
@click.argument(
    'input',
    dest='in_file',
    default='input.yaml',
    help="Input file. Use `mdsapt generate` to generate a template.",
)
@click.argument(
    'output',
    dest='out_file',
    default='out.csv',
    help="Output CSV file to save results to."
)
@click.option('-f', '--force', help="If provided, overwrites existing files.", default=False)
def run(in_file: str, out_file: str, force: bool):
    """
    Run a SAPT calculation.
    """
    ensure_safe_to_overwrite(out_file, force)

    settings = mdsapt.InputReader(in_file)

    optimizer = mdsapt.Optimizer(settings)

    if optimizer.num_failed_residue != 0:
        logger.error('optimization failed see log for list of failed residues')

    sapt = mdsapt.TrajectorySAPT(settings, optimizer).run(settings.start, settings.stop, settings.step)

    logger.info('saving results to CSV')
    sapt.results.to_csv(out_file)


def ensure_safe_to_overwrite(path: str, force: bool):
    """
    Helper function to ensure that it's safe to overwrite the given file, and
    halts the program if not.
    """
    if not os.path.exists(out_file):
        return

    if force:
        logger.warning("will overwrite existing CSV %s", out_file)
        return

    logger.error("Halting, file already exists: %s", out_file)
    logger.error("If you want to overwrite that file, add the -f flag")
    os.exit(-1)

