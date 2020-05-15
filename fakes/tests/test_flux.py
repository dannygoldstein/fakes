import fakes
import numpy as np

rng = np.random.RandomState(1234)

from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import photutils

from astropy.stats import SigmaClip
from photutils import MedianBackground


N_FAKE = 3


def test_recovered_flux_aperture(science_image):

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
        data = hdul[-2].data
        header = hdul[0].header
        table = hdul[-1].data

    # subtract off the background
    sigma_clip = SigmaClip(sigma=3.0)
    bkg = MedianBackground(sigma_clip)
    bkg_val = bkg.calc_background(data)
    bkgsub = data - bkg_val

    APERTURE_RADIUS = 3 * u.pixel

    for row in table:
        ra, dec, mag = row['fake_ra'], row['fake_dec'], row['fake_mag']
        coord = SkyCoord(ra, dec, unit='deg')
        apertures = photutils.SkyCircularAperture(coord, r=APERTURE_RADIUS)
        phot_table = photutils.aperture_photometry(bkgsub, apertures, wcs=wcs)
        flux = phot_table['aperture_sum'][0]

        assert abs(-2.5 * np.log10(flux) + header['MAGZP'] + header['APCOR4'] - mag) < 0.05
