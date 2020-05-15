import fakes
import numpy as np
from astropy.coordinates import SkyCoord

rng = np.random.RandomState(1234)

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u

N_FAKE = 4


def test_fakes(science_image):
    with fits.open(science_image) as hdul:
        header = hdul[0].header
        wcs = WCS(header)
        footprint = wcs.calc_footprint()

    # realize fakes randomly throughout the footprint
    rx = rng.uniform(size=N_FAKE)
    ry = rng.uniform(size=N_FAKE)

    # keep things bright
    mag = rng.uniform(low=15, high=18, size=N_FAKE)

    minra, mindec = footprint.min(axis=0)
    maxra, maxdec = footprint.max(axis=0)

    coord = SkyCoord(minra + (maxra - minra) * rx,
                     mindec + (maxdec - mindec) * ry,
                     unit='deg')

    fakes.inject_psf(science_image, mag, coord)

    with fits.open(science_image) as hdul:
        catalog = hdul[-1].data

    # call sextractor and set persist=True to produce a catalog
    fakes.measure_psf(science_image, persist=True)
    catname = science_image.replace('.fits', '.cat')

    with fits.open(catname) as hdul:
        meascat = hdul[-1].data

    # match fakes to measured stuff
    insky = SkyCoord(catalog['fake_ra'], catalog['fake_dec'], unit='deg')
    outsky = SkyCoord(meascat['X_WORLD'], meascat['Y_WORLD'], unit='deg')
    i1, i2, sep, _ = insky.search_around_sky(outsky, seplimit=0.5 * u.arcsec)

    # assert that all fakes are detected and that measured fake positions
    # agree with the input positions to within 0.5"
    assert len(i1) == N_FAKE
