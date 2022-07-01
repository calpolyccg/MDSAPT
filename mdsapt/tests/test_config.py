from pathlib import Path

import pytest
from pydantic import ValidationError

from mdsapt.config import load_from_yaml_file


resources_dir = Path(__file__).parent / 'testing_resources'

def test_template_can_be_loaded():
    load_from_yaml_file(resources_dir / "test_input.yaml")
