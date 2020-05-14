import os
import pathlib
import tempfile
import subprocess
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


__all__ = ['measure_psf', 'inject_psf']


# this needs to be odd
NPIX = 31


def measure_psf(image, persist=False):
    """Measure the PSF on ZTF science image with path `image`, returning
    a galsim DES_PSFEx object representing the PSF.  Sextractor and PSFex
    are used to model the PSF.

    Options
    -------
    persist: boolean,
        persist intermediate products produced by sextractor and psfex
        instead of deleting them immediately after this function is called

    """

    from galsim.des import DES_PSFEx

    # set up configuration parameters
    confdir = pathlib.Path(__file__).parent / 'astromatic'
    sexconf = confdir / 'sex.conf'
    sexparm = confdir / 'sex.param'
    nnwname = confdir / 'default.nnw'
    convnam = confdir / 'default.conv'
    psfconf = confdir / 'psfex.conf'

    # write the output to a temporary directory

    if persist:
        outdir = os.getcwd()
    else:
        outdir = tempfile.mkdtemp()

    catname = pathlib.Path(outdir) / pathlib.Path(image).name.replace(
        'fits', 'cat'
    )

    psfname = f'{catname.absolute()}'.replace('cat', 'psf')

    cmd = f'sex -c {sexconf} -CATALOG_NAME {catname} -PARAMETERS_NAME {sexparm} ' \
          f'-STARNNW_NAME {nnwname} -FILTER_NAME {convnam} {image}'

    # run sextractor
    subprocess.check_call(cmd.split())

    # now do psfex
    cmd = f'psfex -c {psfconf} {catname}'
    subprocess.check_call(cmd.split())

    # return the psf object
    psf = DES_PSFEx(psfname, image)
    return psf


def inject_psf(image, mag, coord, psf=None, seed=None):
    """Realize the DES_PSFEx PSF model `psf` at location `coord` on
    ZTF science image `image` (image path) with magnitude `mag` in
    the AB system, fluctuated by Poisson noise.

    If `image` is passed as an HDUlist, must be opened in update mode.
    """

    # initialize the random number generator
    import galsim
    rng = galsim.BaseDeviate(seed)

    # handle both scalar and vector inputs
    mag = np.atleast_1d(mag)
    if coord.isscalar:
        coord = coord.reshape([1])

    with fits.open(image, mode='update') as hdul:

        # read in the WCS
        header = hdul[0].header
        wcs = WCS(header=header)

        # measure the PSF using PSFEx if not already specified
        if psf is None:
            psf = measure_psf(image)

        # load the image into galsim
        gimage = galsim.fits.read(hdu_list=hdul)

        # convert the world coordinates to pixel coordinates
        ipos = wcs.all_world2pix(
            [[pos.ra.deg, pos.dec.deg] for pos in coord], 1
        )

        for mag, pos in zip(mag, ipos):

            # calculate the flux of the object
            flux = 10 ** (-0.4 * (mag - header['MAGZP']))

            image_pos = galsim.PositionD(*pos)

            # store the center of the nearest integer pixel
            iimage_pos = galsim.PositionI(*tuple(map(round, pos)))

            # get the noise
            noise = galsim.PoissonNoise(rng)

            # realize the psf at the coordinates
            realization = psf.getPSF(image_pos).withFlux(flux)

            # get the local wcs
            lwcs = psf.getLocalWCS(iimage_pos)

            # calculate the offset between the stamp center
            # and the profile center
            offset = image_pos - iimage_pos

            # draw the image
            imout = realization.drawImage(wcs=lwcs, offset=offset,
                                          nx=NPIX, ny=NPIX)

            # add the noise
            imout.addNoise(noise)

            # shift the image to the right spot
            imout.setCenter(iimage_pos)

            # get the bounds for inpainting
            bounds = imout.bounds

            # add the photons
            gimage[bounds] = gimage[bounds] + imout

        # save it as a new hdu
        galsim.fits.write(gimage, hdu_list=hdul)
