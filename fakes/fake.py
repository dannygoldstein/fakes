import os
import pathlib
import tempfile
import subprocess
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table


__all__ = ['measure_psf', 'inject_psf']


# this needs to be even
NPIX = 16


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
        outdir = os.path.dirname(image)
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


def inject_psf(image, mags, coord, psf=None, seed=None):
    """Realize the DES_PSFEx PSF model `psf` at location `coord` on
    ZTF science image `image` (image path) with magnitude `mag` in
    the AB system, fluctuated by Poisson noise.

    If `image` is passed as an HDUlist, must be opened in update mode.
    """

    # initialize the random number generator
    import galsim
    rng = galsim.BaseDeviate(seed)

    # handle both scalar and vector inputs
    mags = np.atleast_1d(mags)
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

        for mag, pos in zip(mags, ipos):

            # calculate the measured flux of the object
            flux = 10 ** (-0.4 * (mag - header['MAGZP']))

            image_pos = galsim.PositionD(*pos)

            # store the center of the nearest integer pixel
            iimage_pos = galsim.PositionI(*tuple(map(round, pos)))

            ix, iy = iimage_pos.x, iimage_pos.y

            # calculate the offset between the stamp center
            # and the profile center
            offset = image_pos - iimage_pos

            # create an output stamp onto which to draw the fluctuated
            # psf

            bounds = galsim.BoundsI(ix - NPIX,
                                    ix + NPIX,
                                    iy - NPIX,
                                    iy + NPIX)

            # check that there is at least some overlap between the
            bounds = bounds & gimage.bounds
            if not bounds.isDefined():
                raise RuntimeError('No overlap between PSF stamp and image. '
                                   'Is the object coordinate contained by the '
                                   'image?')

            # get the noise
            noise = galsim.PoissonNoise(rng)

            # realize the psf at the coordinates
            realization = psf.getPSF(image_pos).withFlux(flux)

            # get the local wcs
            lwcs = psf.getLocalWCS(image_pos)

            # draw the image
            imout = realization.drawImage(wcs=lwcs, offset=offset,
                                          nx=NPIX * 2 + 1, ny=NPIX * 2 + 1)

            # add the noise
            imout.addNoise(noise)

            # shift the image to the right spot
            imout.setCenter(iimage_pos)

            # add the photons
            gimage[bounds] = gimage[bounds] + imout[bounds]

        # save it as a new hdu
        galsim.fits.write(gimage, hdu_list=hdul)

        # propagate original WCS to output extension
        wcskeys = wcs.to_header(relax=True)
        hdul[-1].header.update(wcskeys)

        # add a record of the fakes as a bintable
        record = {'fake_mag': mags,
                  'fake_ra': coord.ra.deg,
                  'fake_dec': coord.dec.deg,
                  'fake_x': ipos[:, 0],
                  'fake_y': ipos[:, 1]}

        # save the bintable as an extension
        table = Table(record)
        nhdu = fits.BinTableHDU(table.as_array())
        hdul.append(nhdu)
