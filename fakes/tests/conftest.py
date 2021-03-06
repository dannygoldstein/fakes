import pytest
import requests
import tempfile
import pathlib

CACHE_DIR = tempfile.mkdtemp()
URL = 'https://portal.nersc.gov/cfs/astro250/fakes/ztf_20200512312280_' \
      '000718_zg_c01_o_q1_sciimg_ra213.5613_dec38.2261_asec500.fits'

PSF_URL = 'https://portal.nersc.gov/cfs/astro250/fakes/ztf_2020051231228' \
          '0_000718_zg_c01_o_q1_sciimgdaopsfcent.fits'


@pytest.fixture
def science_image():
    r = requests.get(URL)
    outname = pathlib.Path(CACHE_DIR) / URL.split('/')[-1]
    with open(outname, 'wb') as f:
        f.write(r.content)

    return f'{outname.absolute()}'


@pytest.fixture
def ipac_psf_cutout():
    r = requests.get(PSF_URL)
    outname = pathlib.Path(CACHE_DIR) / PSF_URL.split('/')[-1]
    with open(outname, 'wb') as f:
        f.write(r.content)

    return f'{outname.absolute()}'
