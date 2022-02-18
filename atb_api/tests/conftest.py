import pytest
from atb_api.api import ATBApi
from atb_api.molecule import Molecule
from .datafiles import MOL_903922_JSON


@pytest.fixture(scope="session")
def api():
    # if "ATB_API_TOKEN" not in os.environ:
    #     pytest.skip("api token not available")
    api_ = ATBApi()  # api_token=os.environ["ATB_API_TOKEN"])
    return api_


@pytest.fixture(scope="session")
def mol903922():
    return Molecule.parse_file(MOL_903922_JSON)
