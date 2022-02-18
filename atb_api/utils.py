import importlib
from .exceptions import InvalidTokenError


def require_package(package, installation=""):
    installation_messages = {
        "rdkit": "conda install -c conda-forge rdkit",
        "psi4": "conda install -c psi4 psi4",
        "openff.toolkit": "conda install -c conda-forge openff-toolkit",
        "mdanalysis": "conda install -c conda-forge mdanalysis",
    }
    try:
        importlib.import_module(package)
    except ImportError:
        err = (
            f"This function requires the {package} package, "
            "but it could not be imported. "
            "Please make sure it is installed"
        )
        if not installation and package in installation_messages:
            installation = installation_messages[package]
        if installation:
            err += f", or install it with `{installation}`"
        raise ImportError(err) from None


def read_from_string_or_file(file: str) -> str:
    file = str(file)
    try:
        with open(file, "r") as f:
            file = f.read()
    except:
        pass
    return file
