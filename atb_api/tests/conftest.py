import os
import pytest
from atb_api.api import ATBApi


# @pytest.mark.skipif("ATB_API_TOKEN" not in os.environ, reason="api token not available")
@pytest.fixture(scope="session")
def api():
    api_ = ATBApi()  # api_token=os.environ["ATB_API_TOKEN"])
    return api_
