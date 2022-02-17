import os
import pytest
from atb_api.api import ATBApi


@pytest.fixture(scope="session")
def api():
    if "ATB_API_TOKEN" not in os.environ:
        pytest.skip("api token not available")
    api_ = ATBApi(api_token=os.environ["ATB_API_TOKEN"])
    return api_
