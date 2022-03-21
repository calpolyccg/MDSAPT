"""
Main entrypoint for the CLI.
"""

# Note that we import directly from the CLI instead of the full package. This is an optimization
# since some CLI tasks do not need the entire package.
from .cli import cli


if __name__ == '__main__':
    cli()

