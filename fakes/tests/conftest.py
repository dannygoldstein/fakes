import pytest
import requests
import tempfile
import pathlib

CACHE_DIR = tempfile.mkdtemp()
URL = 'https://portal.nersc.gov/cfs/astro250/fakes/ztf_20181120510683_000718_zg_c01_o_q1_sciimg.fits'


@pytest.fixture
def science_image():
    r = requests.get(URL)
    outname = pathlib.Path(CACHE_DIR) / URL.split('/')[-1]
    with open(outname, 'wb') as f:
        f.write(r.content)

    return outname
