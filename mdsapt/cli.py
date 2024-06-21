# pylint: disable=missing-module-docstring
import logging
import os
import sys

import click

# Note that we do not import MDSAPT. This is a speed optimization; it is imported later.


logger = logging.getLogger(__name__)

_dir_path = os.path.dirname(os.path.realpath(__file__))
"""Location of the mdsapt package to be used in resolving templates."""


@click.group()
def cli():
    # pylint: disable=line-too-long
    """
    MDSAPT - Molecular Dynamics Symmetry-Adapted Perturbation Theory, by Alia Lescoulie, Astrid Yu, and Ashley Ringer McDonald.

    This command-line interface lets you easily do common MDSAPT-related tasks.
    """
    # pylint: disable=line-too-long


@cli.command()
@click.argument(
    'filename',
    default='input.yaml',
)
@click.option(
    '-t', '--template', 'template',
    help="Template to generate from. By default, trajectory.",
    type=click.Choice(['trajectory', 'docking'], case_sensitive=False),
    default='trajectory',
)
@click.option(
    '-f', '--force',
    help="If provided, overwrites existing files.",
    is_flag=True,
)
def generate(filename: str, template: str, force: bool):
    """
    Generate a template input file at filename.
    """
    # pylint: disable=fixme
    ensure_safe_to_overwrite(filename, force)

    # TODO: make a wizard for these templates
    template_path = os.path.join(_dir_path, 'data', f'{template}_template.yaml')

    with open(template_path, 'r', encoding='utf-8') as template:
        template_data = template.read()
    with open(filename, 'w', encoding='utf-8') as new_file:
        new_file.write(template_data)

    logger.info('Generated template input file %s', filename)
    # pylint: enable=fixme

@cli.command()
@click.argument(
    'in_file',
    default='input.yaml',
)
@click.argument(
    'out_file',
    default='out.csv',
)
@click.option(
    '-f', '--force',
    help="If provided, overwrites existing files.",
    is_flag=True,
)
def run(in_file: str, out_file: str, force: bool):
    """
    Run a SAPT calculation using the configuration in in_file. Outputs will be written to
    out_file.
    """
    import mdsapt  # pylint: disable=import-outside-toplevel

    ensure_safe_to_overwrite(out_file, force)

    config = mdsapt.load_from_yaml_file(in_file)
    if isinstance(config.analysis, mdsapt.config.TrajectoryAnalysisConfig):
        sapt = mdsapt.TrajectorySAPT(config)
        frames = config.analysis.frames
        sapt.run(frames.start, frames.stop, frames.step)
    elif isinstance(config.analysis, mdsapt.config.DockingAnalysisConfig):
        sapt = mdsapt.DockingSAPT(config)
        sapt.run()

    logger.info('saving results to CSV')
    sapt.results.to_csv(out_file)


def ensure_safe_to_overwrite(path: str, force: bool):
    """
    Helper function to ensure that it's safe to overwrite the given file, and
    halts the program if not.
    """
    if not os.path.exists(path):
        return

    if force:
        logger.warning("will overwrite existing CSV %s", path)
        return

    logger.error("Halting, file already exists: %s", path)
    logger.error("If you want to overwrite that file, add the -f flag")
    sys.exit(-1)
