import os
import pathlib
import tempfile
import subprocess
from astropy.io import fits
from astropy.wcs import WCS


__all__ = ['measure_psf', 'inject_psf']


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


def inject_psf(image, mag, coord, psf=None):
    """Realize the DES_PSFEx PSF model `psf` at location `coord` on
    ZTF science image `image` (either an HDUList or image path) with
    magnitude `mag` in the AB system, fluctuated by Poisson noise.

    If `image` is passed as an HDUlist, must be opened in update mode.
    """

    import galsim
    with fits.open(image, mode='update') as hdul:

        header = hdul[0].header
        wcs = WCS(header)

        if psf is None:
            psf = measure_psf(image)

        # calculate the flux of the object
        flux = 10**(-0.4 * (mag - header['MAGZP']))

        # get the image coordinates of the psf
        ix, iy = wcs.all_world2pix([[coord.ra.deg, coord.dec.deg]], 1)[0]
        image_pos = galsim.PositionD(ix, iy)

        # realize the psf at the coordinates
        realization = psf.getPSF(image_pos).withFlux(flux)

        # get the local wcs
        lwcs = psf.getLocalWCS(image_pos)

        #noise = galsim.PoissonNoise()

        # now load the image into galsim
        gimage = galsim.fits.read(hdu_list=hdul)
        realization.drawImage(image=gimage, wcs=lwcs, add_to_image=True)

        # save it as a new hdu
        galsim.fits.write(gimage, hdu_list=hdul)
