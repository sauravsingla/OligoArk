from importlib.metadata import version

import oligoark


def test_package_metadata_version_matches_public_version() -> None:
    assert version("oligoark") == oligoark.__version__ == "0.4.0"
